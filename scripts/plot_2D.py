from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import itertools
from utilities.readscan import readscan
from utilities.waveform import Waveform, DEFAULT_FMAX

class XVarOptions:
    ANGLE = "angle"
    COSANGLE = "cos"
    ALL_CHOICES = [ANGLE, COSANGLE]

def plot(data, labels, output=None, xlim=None, ylim=None, zlim=None, xvar=XVarOptions.ANGLE, zvar="area", fmax=DEFAULT_FMAX, dohemisphere=False):
    fig = plt.figure(figsize=(9,18))
    ax1 = fig.add_subplot(2, 1, 1, projection="3d")
    ax2 = fig.add_subplot(2, 1, 2, projection="3d")
    Ncol = len(data)
    colmap = itertools.cycle(range(Ncol))
    for i, (d, l) in enumerate(zip(data, labels)):
        color = matplotlib.cm.hot(float(next(colmap))/float(Ncol))
        _plot_series(ax1, ax2, d.scanpoints, color, l, xlim=xlim, ylim=ylim, zlim=zlim, xvar=xvar, zvar=zvar, fmax=fmax, dohemisphere=dohemisphere)
    #finalise the plot
    ax1.legend()
    ax2.legend()
    if xlim:
        for ax in [ax1, ax2]:
            ax.set_xlim(xlim)
    if ylim:
        for ax in [ax1, ax2]:
            ax.set_ylim(ylim)
    # if zlim:
    #     ax2.set_ylim(ylim)
    fig.tight_layout()
    if output is None:
        plt.show()
    else:
        _writefig(fig, output)
    return

def _plot_series(ax1, ax2, data, color, label, xlim=None, ylim=None, zlim=None, xvar=XVarOptions.ANGLE, zvar="area", fmax=DEFAULT_FMAX, dohemisphere=False):
    if xvar == XVarOptions.COSANGLE:
        xtransform = lambda x: math.cos(x.coord_angle*math.pi/180.)
    elif xvar == XVarOptions.ANGLE:
        xtransform = lambda x: x.coord_angle
    get_zvar = lambda x: getattr(Waveform(x.axis_time, x.samples_PMT, fmax), zvar)

    X_list, Y_list, Z_list = [], [], []
    for d in data:
        if rangecheck(xlim, ylim, d.coord_angle, d.coord_y):
            X_list.append(xtransform(d))
            Y_list.append(d.coord_y)
            Z_list.append(get_zvar(d) * (hemisphere_correction(d.coord_angle, zvar) if dohemisphere else 1))

    X = np.array(X_list, dtype=float)
    Y = np.array(Y_list, dtype=float)    
    Z = np.array(Z_list, dtype=float)

    #ax1.scatter(X, Y, Z, color=color, alpha=0.5)

    newX, newY, meanZ, errZ = _get_angle_scan_moments(X_list, Y_list, Z_list)
    #ax1.scatter(newX, newY, meanZ, "-", color=color, label=label)
    #ax1.fill_between()

def _get_angle_scan_moments(X, Y, Z):
    M = {}
    for x, y, z in zip(X, Y, Z):
        if not (x, y) in M:
            M[(float(x), float(y))] = []
        M[(float(x), float(y))].append(float(z))
    moments = [(coord[0], coord[1], np.mean(m), np.std(m)) for (coord, m) in sorted(M.items())]
    newX, newY, meanZ, errZ = zip(*moments)
    return np.array(newX), np.array(newY), np.array(meanZ), np.array(errZ)

def rangecheck(xlim, ylim, xcoord, ycoord):
    x_in_range = xlim is None or xlim[0] <= xcoord <= xlim[1]
    y_in_range = ylim is None or ylim[0] <= ycoord <= ylim[1]
    return x_in_range and y_in_range

def hemisphere_correction(angle, var):
    if var in ["area", "amplitude"]:
        theta = float(angle)*np.pi/180. # Convert to radians
        return 1.0 / (1.0 - (abs(theta) / np.pi))


def parsecml():
    parser = ArgumentParser()
    parser.add_argument("-o", "--output", help="Output filename", default=None)
    parser.add_argument("--labels", help="Comma separated list of names", default=None)
    parser.add_argument("filename", nargs="+", help="Input filename.")
    parser.add_argument("--xlim", type=str, help="Set x-axis limits.", default=None)
    parser.add_argument("--ylim", type=str, help="Set y-axis limits.", default=None)
    parser.add_argument("--zlim", type=str, help="Set z-axis limits.", default=None)
    parser.add_argument("--xvar", type=str, default=XVarOptions.ANGLE, choices=XVarOptions.ALL_CHOICES)
    parser.add_argument("--zvar", type=str, default="area", help="Choose z-axis variable")
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
    xlim, ylim, zlim = None, None, None
    if args.xlim:
        low, high = args.xlim.split(",")
        xlim = (float(low), float(high))
    if args.ylim:
        low, high = args.ylim.split(",")
        ylim = (float(low), float(high))
    if args.zlim:
        low, high = args.zlim.split(",")
        zlim = (float(low), float(high))
    fmax=args.fmax if args.filter else None
    plot(data, labels, args.output, xlim=xlim, ylim=ylim, zlim=zlim, xvar=args.xvar, zvar=args.zvar, fmax=fmax, dohemisphere=args.hemisphere_correction)
    return

if __name__ == "__main__":
    main()