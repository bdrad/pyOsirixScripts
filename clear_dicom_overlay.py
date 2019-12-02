# Author: Aarash Heydari with help from https://groups.google.com/forum/#!topic/horos-project/pwz9aSWhYJ8

import osirix
import dicom
import numpy as np
import shutil
import os

vc = osirix.frontmostViewer()
pixs = vc.pixList()

wait = vc.startWaitProgressWindow('Deleting all DICOM Overlays…', len(pixs))
for pix in pixs:
    ds = dicom.read_file(pix.sourceFile)

    for j in range(16):
        try:
            del ds[0x60000010 + j*0x20000]
            del ds[0x60000011 + j*0x20000]
            del ds[0x60000022 + j*0x20000] 
            del ds[0x60000040 + j*0x20000] 
            del ds[0x60000050 + j*0x20000]
            del ds[0x60000100 + j*0x20000]
            del ds[0x60000102 + j*0x20000]
            del ds[0x60003000 + j*0x20000]
            print("Deleted ROI at index ", j, " for dicom at slice ", ds.ImagePositionPatient[2])
                  
        except Exception:
	        pass

    ds.save_as(pix.sourceFile)
    wait.incrementBy(1.0)

vc.endWaitWindow(wait)
