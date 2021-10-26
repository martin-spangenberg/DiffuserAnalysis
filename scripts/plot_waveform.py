from argparse import ArgumentParser
import itertools
import numpy as np
from collections import namedtuple
import matplotlib
import matplotlib.pyplot as plt
from scripts.utilities.readscan import readscan, ScanPoint
from scripts.utilities.waveform import Waveform, DEFAULT_FMAX

def plot(data, labels, channelname, output=None, xlim=None, ylim=None, height=0.0, angles=[0.0], peaknorm=False, areanorm=False, statistics=False, showfilter=True, fmax=DEFAULT_FMAX):
    fig = plt.figure(figsize=(9,9))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)    
    Ncol = len(data) * len(angles)
    colmap = itertools.cycle(range(Ncol))
    for d, l in itertools.zip_longest(data, labels):
        for angle in angles:
            color = matplotlib.cm.hot(float(next(colmap))/float(Ncol))
            waveform = _pick_waveform(d.scanpoints, channelname, height, angle, peaknorm, areanorm)
            lwitha = l + str(r"$\theta=%.0f$"%angle)
            ax1.plot(waveform.X, waveform.Y, label=lwitha, color=color)
            if statistics:
                _plot_stats(waveform, ax1, lwitha, color)
            if showfilter:
                waveform_filter = _pick_waveform(d.scanpoints, channelname, height, angle, peaknorm, areanorm, fmax)
                ax1.plot(waveform_filter.X, waveform_filter.Y, label="filtered")
            fX, fY = waveform.fft
            ax2.plot(fX, np.abs(fY), label=lwitha, color=color)
    #finalise the plot
    ax1.legend()
    ax2.legend()
    if xlim:
        ax1.set_xlim(xlim)
    if ylim:
        ax1.set_ylim(ylim)
    ax1.set_xlabel("time [ns]")
    ax1.set_ylabel("voltage [V]")
    ax2.set_xlabel("frequency [GHz]")
    ax2.set_xlim((0.0, 1.0))
    ax2.set_ylabel("1/voltage [V^-1]")
    fig.tight_layout()
    if output is None:
        plt.show()
    else:
        _writefig(fig, output)
    return

def _pick_waveform(scanpoints, channelname, height, angle, peaknorm=False, areanorm=False, fmax=None):
    for scanpoint in scanpoints:
        if scanpoint.coord_y == height:
            if scanpoint.coord_angle == angle:
                waveform = Waveform(scanpoint.axis_time, getattr(scanpoint, channelname), fmax, peaknorm, areanorm)
                return waveform
    raise Exception("no matching record for angle", angle, "and height", height)

def _plot_stats(wave, ax1, label, color):
    _print_waveform(wave, label)
    ax1.plot([min(wave.X), max(wave.X)], [wave.pedestal, wave.pedestal], "--", color=color)
    ax1.plot([wave.peak_time, wave.peak_time], [wave.pedestal, wave.pedestal - wave.amplitude], "--", color=color)
    ax1.plot([wave.peak_time - wave.rise, wave.peak_time], [-wave.amplitude/2.0 + wave.pedestal, -wave.amplitude/2.0 + wave.pedestal], "--", color=color)
    ax1.plot([wave.peak_time, wave.peak_time+wave.fall], [-wave.amplitude/2.0 + wave.pedestal, -wave.amplitude/2.0  + wave.pedestal], "--", color=color)
    ax1.fill_between([wave.X[wave.integralindex_low], wave.X[wave.integralindex_high]], [-wave.amplitude, -wave.amplitude], [wave.pedestal, wave.pedestal], alpha=0.25)
    return

def _print_waveform(wave, label):
    print("\"%s\"" % label)
    for var in ["pedestal", "amplitude", "peak_time", "fwhm", "rise", "fall", "area"]:
        print("%s : %.3e" % (var, getattr(wave, var)))
    return

def _writefig(fig, output):
    try:
        os.makedirs(os.path.dirname(output))
    except:
        pass
    fig.savefig(output)
    return

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("-o", "--output", help="Output filename", default=None)
    parser.add_argument("--channel", type=int, default=1, choices=[1, 2], help="Select which digitizer input channel to plot.")
    parser.add_argument("--height", help="Height to plot.", default="0", type=float)
    parser.add_argument("--angle", help="Comma separated list of angles to plot.", default="0", type=str)
    parser.add_argument("--labels", help="Comma separated list of names", default=None)
    parser.add_argument("--xlim", type=str, help="Set x-axis limits.", default=None)
    parser.add_argument("--ylim", type=str, help="Set y-axis limits.", default=None)
    parser.add_argument("--peak-norm", action="store_true", help="Normalise to peak.")
    parser.add_argument("--area-norm", action="store_true", help="Normalise to unit area.")
    parser.add_argument("--stats", action="store_true", help="Show waveform statistics")
    parser.add_argument("--filter", action="store_true", help="Switch on low-pass filter.")
    parser.add_argument("--fmax", help="Set frequency cut for low pass filter.", default=DEFAULT_FMAX, type=float)
    parser.add_argument("filename", nargs="+", help="Input filename.")
    return parser.parse_args()

def main():
    args = parsecml()
    plt.rc('legend', fontsize="x-small")
    data = [readscan(f, "diffuser") for f in args.filename]
    channelvars = ["samples_PMT", "samples_PD"]
    channelname = channelvars[args.channel-1]
    if args.labels:
        labels = args.labels.split(",")
    else:
        labels = ["series %s" % n for n in range(len(data))]
    xlim, ylim = None, None
    if args.xlim:
        low, high = args.xlim.split(",")
        xlim = (float(low), float(high))
    if args.ylim:
        low, high = args.ylim.split(",")
        ylim = (float(low), float(high))
    height = args.height
    angles = [float(x) for x in args.angle.split(",")]
    plot(data, labels, channelname, args.output,
         xlim=xlim, ylim=ylim, height=height, angles=angles,
         peaknorm=args.peak_norm,
         areanorm=args.area_norm,
         statistics=args.stats,
         showfilter=args.filter,
         fmax=args.fmax,
    )
    return

if __name__ == "__main__":
    main()
