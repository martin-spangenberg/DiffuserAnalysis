import ROOT
import numpy as np
import datetime
from collections import namedtuple
from collections import defaultdict
from .waveform import Waveform

class ScanPoint:
    def __init__(self, time_epoch, coord_angle, coord_y, lab_temp, lab_humid, dt, samples_PMT, samples_PD, num_entry):
        self.time_epoch = time_epoch
        self.coord_angle = coord_angle
        self.coord_y = coord_y
        self.lab_temp = lab_temp
        self.lab_humid = lab_humid
        self.samples_PMT = np.array(list(samples_PMT), dtype=float)
        self.samples_PD = np.array(list(samples_PD), dtype=float)
        self.axis_time = np.arange(len(samples_PMT), dtype=float) * dt * 1.e9
        self.num_entry = num_entry

    def getWaveformPMT(self, fmax=None, peaknorm=False, areanorm=False):
      return Waveform(self.axis_time, self.samples_PMT, fmax, peaknorm, areanorm)
    
    def getWaveformPD(self, fmax=None, peaknorm=False, areanorm=False):
      return Waveform(self.axis_time, self.samples_PD, fmax, peaknorm, areanorm)

    def getDatetime(self):
      return datetime.datetime.fromtimestamp(self.time_epoch)

# Tuple containing metadata and scanpoints for a full scan
DiffuserScan = namedtuple("DiffuserScan", ["ID_diffuser",
                                           "ID_PMT",
                                           "ID_PD",
                                           "ID_lightsource",
                                           "ID_experimentalist",
                                           "notes",
                                           "pulse_rate",
                                           "pulse_N",
                                           "scanpoints"])

def readscan(filename, treename):
    file = ROOT.TFile.Open(filename, "READ")
    tree = getattr(file, treename)

    tree.GetEntry(0)

    formatversion = (
        int(tree.version_major),
        int(tree.version_minor),
        int(tree.version_patch),
    )

    # Check that input file format is implemented
    if formatversion[:2] == (1, 0):
        diffuserscan = DiffuserScan(
            str(tree.ID_diffuser),
            str(tree.ID_PMT),
            str(tree.ID_PD),
            str(tree.ID_lightsource),
            str(tree.ID_experimentalist),
            str(tree.notes),
            float(tree.pulse_rate),
            int(tree.pulse_N),
            [],
        )

        for i, entry in enumerate(tree):
            scanpoint = ScanPoint(entry.time_epoch,
                                  entry.coord_angle,
                                  entry.coord_y,
                                  entry.lab_temp,
                                  entry.lab_humid,
                                  entry.dt,
                                  entry.waveform_PMT,
                                  entry.waveform_PD,
                                  i)

            diffuserscan.scanpoints.append(scanpoint)

    else:
        raise NotImplementedError("Input format version not implemented")

    return diffuserscan