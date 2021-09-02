from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import itertools
import math
from utilities.readscan import readscan
from utilities.waveform import Waveform, DEFAULT_FMAX

class XVarOptions:
    ANGLE = "angle"
    COSANGLE = "cos"
    ALL_CHOICES = [ANGLE, COSANGLE]

class YVarOptions:
    HEIGHT = "height"
    VERTANGLE = "vertangle"
    ALL_CHOICES = [HEIGHT, VERTANGLE]

def plot(data, labels, output=None, xlim=None, ylim=None, zlim=None, centre_height=None, pmt_dist=None, xvar=XVarOptions.ANGLE, yvar=YVarOptions.HEIGHT, zvar="area", fmax=DEFAULT_FMAX, dohemisphere=False):
    fig = plt.figure(figsize=(9,18))
    ax1 = fig.add_subplot(2, 1, 1, projection="3d")
    ax2 = fig.add_subplot(2, 1, 2, projection="3d")
    Ncol = len(data)
    colmap = itertools.cycle(range(Ncol))
    for i, (d, l) in enumerate(zip(data, labels)):
        color = matplotlib.cm.hot(float(next(colmap))/float(Ncol))
        _plot_series(ax1, ax2, d.scanpoints, color, l, xlim=xlim, ylim=ylim, zlim=zlim, centre_height=centre_height, pmt_dist=pmt_dist, xvar=xvar, yvar=yvar, zvar=zvar, fmax=fmax, dohemisphere=dohemisphere)
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

def _plot_series(ax1, ax2, data, color, label, xlim=None, ylim=None, zlim=None, centre_height=None, pmt_dist=None, xvar=XVarOptions.ANGLE, yvar=YVarOptions.HEIGHT, zvar="area", fmax=DEFAULT_FMAX, dohemisphere=False):
    if xvar == XVarOptions.COSANGLE:
        xtransform = lambda x: math.cos(x.coord_angle*math.pi/180.)
    elif xvar == XVarOptions.ANGLE:
        xtransform = lambda x: x.coord_angle

    if yvar == YVarOptions.HEIGHT:
        ytransform = lambda x: x.coord_y
    elif yvar == YVarOptions.VERTANGLE:
        ytransform = lambda x: np.sign(x.coord_y-centre_height) * math.atan(abs(x.coord_y-centre_height) / pmt_dist) * 360 / (2*math.pi)

    if yvar == YVarOptions.VERTANGLE:
        ztransform = lambda x: getattr(Waveform(x.axis_time, x.samples_PMT, fmax), zvar) * fluxtransform(pmt_dist, x.coord_y-centre_height)
    else:
        ztransform = lambda x: getattr(Waveform(x.axis_time, x.samples_PMT, fmax), zvar)

    X_list, Y_list, Z_list = [], [], []
    for d in data:
        if rangecheck(xlim, ylim, d.coord_angle, d.coord_y):
            X_list.append(xtransform(d))
            Y_list.append(ytransform(d))
            Z_list.append(ztransform(d) * (hemisphere_correction(d.coord_angle, zvar) if dohemisphere else 1))

    X = np.array(X_list, dtype=float)
    Y = np.array(Y_list, dtype=float)    
    Z = np.array(Z_list, dtype=float)

    ax1.scatter(X, Y, Z, color=color, alpha=0.5)
    #ax1.plot_trisurf(X, Y, Z, color=color, alpha=0.5)


    ax1.set_xlabel("horizontal angle [degrees]")
    ax1.set_ylabel("vertical angle [degrees]")
    ax1.set_zlabel("area")


    newX, newY, meanZ, errZ = _get_angle_scan_moments(X_list, Y_list, Z_list)

    #ax1.errorbar(newX, newY, meanZ, zerr=0.2, linestyle='none')



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

# Transform incoming signal (flux) to match angle in vertical plane at constant distance pmt_dist
# This is needed in order to plot signal vs. angle in vertical plane, since experimental setup only moves PMT straight up/down
def fluxtransform(d, h):
    return math.pow(math.pow(d, 2) + math.pow(h, 2), 3/2) / math.pow(d, 3)

def hemisphere_correction(angle, var):
    if var in ["area", "amplitude"]:
        theta = float(angle)*np.pi/180. # Convert to radians
        return 1.0 / (1.0 - (abs(theta) / np.pi))


def parsecml():
    parser = ArgumentParser()
    parser.add_argument("-o", "--output", help="Output filename", default=None)
    parser.add_argument("--labels", help="Comma separated list of names", default=None)
    parser.add_argument("filename", nargs="+", help="Input filename.")
    parser.add_argument("--combinefiles", action="store_true", help="Combine input files and plot them as a single series.")
    parser.add_argument("--xlim", type=str, help="Set x-axis limits.", default=None)
    parser.add_argument("--ylim", type=str, help="Set y-axis limits.", default=None)
    parser.add_argument("--zlim", type=str, help="Set z-axis limits.", default=None)
    parser.add_argument("--xvar", type=str, default=XVarOptions.ANGLE, choices=XVarOptions.ALL_CHOICES)
    parser.add_argument("--yvar", type=str, default=YVarOptions.HEIGHT, choices=YVarOptions.ALL_CHOICES)
    parser.add_argument("--zvar", type=str, default="area", help="Choose z-axis variable")
    parser.add_argument("--centre_height", type=str, default="0.0", help="Set height where diffuser is centred vs. PMT [mm]. Only used when yvar=vertangle")
    parser.add_argument("--pmt_dist", type=str, help="Set PMT distance from diffuser [mm]. Only used when yvar=vertangle")
    parser.add_argument("--hemisphere-correction", "-c", action="store_true", help="Apply correction (1-theta/pi) to intensity measurments to account for hemisphere geometry.")
    parser.add_argument("--filter", action="store_true", help="Switch on low-pass filter.")
    parser.add_argument("--fmax", help="Set frequency cut for low pass filter.", default=DEFAULT_FMAX, type=float)
    return parser.parse_args()

def main():
    args = parsecml()
    matplotlib.rc('legend', fontsize="x-small")
    data = [readscan(f, "diffuser") for f in args.filename]
    if args.combinefiles:
        data = [combinefiles(data)]
    if args.labels:
        labels = args.labels.split(",")
    else:
        labels = ["series %s" % n for n in range(len(data))]
    xlim, ylim, zlim, centre_height, pmt_dist = None, None, None, None, None
    if args.xlim:
        low, high = args.xlim.split(",")
        xlim = (float(low), float(high))
    if args.ylim:
        low, high = args.ylim.split(",")
        ylim = (float(low), float(high))
    if args.zlim:
        low, high = args.zlim.split(",")
        zlim = (float(low), float(high))
    if args.centre_height:
        centre_height = float(args.centre_height)
    if args.pmt_dist:
        pmt_dist = float(args.pmt_dist)
    fmax=args.fmax if args.filter else None
    plot(data, labels, args.output, xlim=xlim, ylim=ylim, zlim=zlim, centre_height=centre_height, pmt_dist=pmt_dist, xvar=args.xvar, yvar=args.yvar, zvar=args.zvar, fmax=fmax, dohemisphere=args.hemisphere_correction)
    return

if __name__ == "__main__":
    main()
