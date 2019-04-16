# -*- coding: utf-8 -*-
"""Test program using Python wrapper for timsdata.dll, tested with Python 3.6"""

from bdal.io.timsdata import TimsData
import os, sys
import numpy as np

if len(sys.argv) < 2:
    raise RuntimeError("need arguments: tdf_directory")

analysis_dir = sys.argv[1]

tdf = TimsData(analysis_dir, True)

frame_id=1
ids_query = tdf.conn.execute("SELECT Id, NumScans FROM Frames where Id={}".format(frame_id)).fetchall()

for id, num_scans in ids_query:
    # read msms spectra for MS1 frame with specific Frame.Id
    msms_spectra = tdf.readPasefMsMsForFrame(id)
    for precursor_id, msms_spectrum in msms_spectra.items():
        print ("precursor {}, number msms peaks {}".format(precursor_id, len(msms_spectrum[0])))
        print (msms_spectrum[0])
        print (msms_spectrum[1])