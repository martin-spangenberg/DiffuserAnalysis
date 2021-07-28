from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import itertools
from utilities.readscan import readscan
from utilities.waveform import Waveform, DEFAULT_FMAX

class XVarOptions:
    ANGLE = "angle"
    COSANGLE = "cos"
    HEIGHT = "height"
    TIME = "time"
    ALL_CHOICES = [ANGLE, COSANGLE, HEIGHT, TIME]

def plot(data, labels, output=None, xlim=None, ylim=None, anglelim=None, heightlim=None, yvar="area", xvar=XVarOptions.ANGLE, fmax=DEFAULT_FMAX, dohemisphere=False):
    fig = plt.figure(figsize=(9,18))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    Ncol = len(data)
    colmap = itertools.cycle(range(Ncol))
    for i, (d, l) in enumerate(zip(data, labels)):
        color = matplotlib.cm.hot(float(next(colmap))/float(Ncol))
        _plot_series(ax1, ax2, d.scanpoints, color, l, yvar=yvar, xlim=xlim, anglelim=anglelim, heightlim=heightlim, xvar=xvar, fmax=fmax, dohemisphere=dohemisphere)
    #finalise the plot
    ax1.legend()
    ax2.legend()
    if xlim and not xvar == XVarOptions.TIME:
        for ax in [ax1, ax2]:
            ax.set_xlim(xlim)
    if ylim:
        ax2.set_ylim(ylim)
    fig.tight_layout()
    if output is None:
        plt.show()
    else:
        _writefig(fig, output)
    return

def _plot_series(ax1, ax2, data, color, label, yvar="area", xlim=None, anglelim=None, heightlim=None, xvar=XVarOptions.ANGLE, fmax=DEFAULT_FMAX, dohemisphere=False):
    temperature = 0.
    humidity = 0.
    if xvar == XVarOptions.COSANGLE:
        xtransform = lambda x: math.cos(x.coord_angle*math.pi/180.)
    elif xvar == XVarOptions.ANGLE:
        xtransform = lambda x: x.coord_angle
    elif xvar == XVarOptions.HEIGHT:
        xtransform = lambda x: x.coord_y
    elif xvar == XVarOptions.TIME:
        xtransform = lambda x: x.num_entry
    get_yvar = lambda x: getattr(Waveform(x.axis_time, x.samples_PMT, fmax), yvar)
    X_list, Y_list = [], []
    for d in data:
        if rangecheck(xtransform(d), xlim) and rangecheck(d.coord_angle, anglelim) and rangecheck(d.coord_y, heightlim):
            X_list.append(xtransform(d))
            Y_list.append(get_yvar(d) * hemisphere_correction(d.coord_angle, yvar, dohemisphere))
            temperature += d.lab_temp
            humidity += d.lab_humid

    temperature /= len(X_list)
    humidity /= len(X_list)
    print("Average temperature: {:.2f}C".format(temperature))
    print("Average humidity: {:.2f}%".format(humidity))

    X = np.array(X_list, dtype=float)
    Y = np.array(Y_list, dtype=float)
    #plot data points
    ax1.plot(X, Y, "o", color=color, alpha=0.5)
    #plot mean and sigma
    meanX, meanY, errY = _get_angle_scan_moments(X_list, Y_list)
    ax1.plot(meanX, meanY, "-", color=color, label=label)
    ax1.fill_between(meanX, meanY-errY, meanY+errY, alpha=0.5, color=color)
    ax1.set_ylabel(_get_label(yvar))
    if xvar == XVarOptions.COSANGLE:
        ax1.set_xlabel("cos(theta)")
    elif xvar == XVarOptions.ANGLE:
        ax1.set_xlabel("angle [degrees]")
    elif xvar == XVarOptions.TIME:
        ax1.set_xlabel("measurement number")
    #plot fractional change from sigma
    if xvar == XVarOptions.COSANGLE:
        norm_x_point = np.max(X)
        if norm_x_point != 1.0:
            print("WARNING: data contains no point at cos(theta) == 1.0, normalising to {}".format(norm_x_point))
        offset = np.mean([y for x, y in zip(X, Y) if x == 1.0])
    else:
        norm_x_point = X[np.argmin(np.abs(X))]
        if norm_x_point != 0.0:
            print("WARNING: data contains no point at theta == 0.0, normalising to {}".format(norm_x_point))
        offset = np.mean([y for x, y in zip(X, Y) if x == norm_x_point])
    meanY = ((meanY / offset) - 1.0) * 100.0
    errY = ((errY / offset) * 100.0)
    ax2.plot(meanX, meanY, "-", color=color, label=label)
    ax2.fill_between(meanX, meanY-errY, meanY+errY, alpha=0.5, color=color)
    ax2.set_ylabel(r"relative difference [%]")
    ax2.set_xlabel("angle [degrees]")
    return

def rangecheck(var, lim):
    return lim == None or lim[0] <= var <= lim[1]

def _writefig(fig, output):
    try:
        os.makedirs(os.path.dirname(output))
    except:
        pass
    fig.savefig(output)
    return

def _get_label(var):
    try:
        return {"power":"power [pW]",
                "angle": "angle [degrees]",
                "height" : "PMT pulse height [V]",
                "fwhm" : "PMT FWHM [pulse area [Vs]",
            }[var]
    except KeyError:
        return str(var)

def _get_angle_scan_moments(X, Y):
    Z = {}
    for x, y in zip(X, Y):
        if not x in Z:
            Z[float(x)] = []
        Z[float(x)].append(float(y))
    moments = [(x, np.mean(z), np.std(z)) for (x, z) in sorted(Z.items())]
    newX, meanY, errY = zip(*moments)
    return np.array(newX), np.array(meanY), np.array(errY)

def hemisphere_correction(angle, var, dohemisphere):
    if dohemisphere and var in ["area", "amplitude"]:
        theta = float(angle)*np.pi/180. # Convert to radians
        return 1.0 / (1.0 - (abs(theta) / np.pi))
    else:
        return 1.0

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("-o", "--output", help="Output filename", default=None)
    parser.add_argument("--labels", help="Comma separated list of names", default=None)
    parser.add_argument("filename", nargs="+", help="Input filename.")
    parser.add_argument("--xlim", type=str, help="Set x-axis limits.", default=None)
    parser.add_argument("--ylim", type=str, help="Set y-axis limits.", default=None)
    parser.add_argument("--anglelim", type=str, help="Set limits on angle.", default=None)
    parser.add_argument("--heightlim", type=str, help="Set limits on PMT height.", default=None)
    parser.add_argument("--yvar", type=str, default="area", help="Choose y-axis variable")
    parser.add_argument("--xvar", type=str, default=XVarOptions.ANGLE, choices=XVarOptions.ALL_CHOICES)
    parser.add_argument("--hemisphere-correction", "-c", action="store_true", help="Apply correction (1-theta/pi) to intensity measurments to account for hemisphere geometry.")
    parser.add_argument("--filter", action="store_true", help="Switch on low-pass filter.")
    parser.add_argument("--fmax", help="Set frequency cut for low pass filter.", default=DEFAULT_FMAX, type=float)
    return parser.parse_args()

def main():
    args = parsecml()
    matplotlib.rc('legend', fontsize="x-small")
    data = [readscan(f, "diffuser") for f in args.filename]
    if args.labels:
        labels = args.labels.split(",")
    else:
        labels = ["series %s" % n for n in range(len(data))]
    xlim, ylim, anglelim, heightlim = None, None, None, None
    if args.xlim:
        low, high = args.xlim.split(",")
        xlim = (float(low), float(high))
    if args.ylim:
        low, high = args.ylim.split(",")
        ylim = (float(low), float(high))
    if args.anglelim:
        low, high = args.anglelim.split(",")
        anglelim = (float(low), float(high))
    if args.heightlim:
        low, high = args.heightlim.split(",")
        heightlim = (float(low), float(high))
    fmax=args.fmax if args.filter else None
    plot(data, labels, args.output, xlim=xlim, ylim=ylim, anglelim=anglelim, heightlim=heightlim, yvar=args.yvar, xvar=args.xvar, fmax=fmax, dohemisphere=args.hemisphere_correction)
    return

if __name__ == "__main__":
    main()