###############################################
# Author: Aarash Heydari (aheyd@berkeley.edu) #
###############################################
import osirix
import dicom
import numpy as np
import shutil
import os
import subprocess
import csv

csv_file_path = "./Users/admin/Downloads/annotations.csv"

vc = osirix.frontmostViewer()
pixs = vc.pixList()
instanceToPixPosition = {pix.imageObj().instanceNumber:i for i, pix in enumerate(pixs)}

def import_from_csvs():
    vc = osirix.frontmostViewer()
    for scans, dcmseriesuid in [load_scan(x.paths()) for x in osirix.currentBrowser().databaseSelection()]:
        spacing, origin = scans.values()[0].PixelSpacing, scans.values()[0].ImagePositionPatient
        with open(csv_file_path) as csvfile:
            csvreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
            columns = next(csvreader) # Read the first row, which is names of columns rather than data.

            # This dict records number of overlays already in a slice, where key = slice Z position
            NextOverlayIndex = {} 
            for i, row in enumerate(csvreader):
                try:
                     s = row[0].split(",")
                     csvseriesuid, coordX, coordY, coordZ, diameter_mm = s[0], s[1], s[2], s[3].strip("\'"), float(s[4])
                except:
                     print("Failed with len row: ", len(row[0].split(",")), "at i = %d" % i)
                     print(row[0])
                     continue
                if dcmseriesuid == csvseriesuid:
                     print("equality at csvrow = %d" % i)
                     vX, vY = world_to_voxel(float(coordX), float(coordY), spacing, origin)
                     print("ROI at voxel coordinates", coordX, coordY, vX, vY, coordZ)
                     Zlist = [z for z in scans.keys()]
                     scanZ = sorted(Zlist, key=lambda x: np.abs(float(x) - np.round(float(coordZ.rstrip('0')))))[0]
                     print(scanZ)
                     ds = scans[scanZ]
                     overlayIndex = NextOverlayIndex.get(scanZ, 0)
                     mask = create_circular_mask(512, 512, vX, vY, radius=diameter_mm/2*(1/ds.PixelSpacing[0]))
                     embed(mask, ds, overlayIndex, OVERLAY_CSV_NAME_SUFFIX)
                     roiNew = osirix.ROI(itype='tPlain',buffer=mask.T.astype(np.bool),name="ROI read from CSV row %d" % i,ipixelSpacing=spacing, iimageOrigin=origin[:2])
                     vc.setROI(roiNew,position=instanceToPixPosition[ds.InstanceNumber])
                     NextOverlayIndex[scanZ] = overlayIndex + 1
                     ds.save_as(pixs[instanceToPixPosition[ds.InstanceNumber]].sourceFile)
                     print("saved file in %s for instance %f with overlayIndex = %d, radius=%f\n\n\n" % (pixs[instanceToPixPosition[ds.InstanceNumber]].sourceFile, ds.InstanceNumber, overlayIndex, diameter_mm))

#############
### Utils ###
#############
OVERLAY_CSV_NAME_SUFFIX = "Annotation read from CSV"
# Embeds a mask into a DICOM file at the given overlay index
def embed(mask, dcm_file, overlay_index, name):
    if overlay_index > 16:
        return # cannot support more than 16 overlays.
    # packbits converts array of integer 1s and 0s to array of numbers.
    # It interprets a run of 8 bits as a number.
    # i.e. [1, 0, 0, 0, 0, 0, 0, 0] -> [128]
    # The reshaping and flattening is to accomodate DICOM Header's expected format.
    w = np.packbits(mask.reshape(-1,8)[:,::-1].flatten('C'))
    dcm_file.add_new(0x60000010 + overlay_index*0x20000 , 'US', 512)
    dcm_file.add_new(0x60000011 + overlay_index*0x20000 , 'US', 512)
    dcm_file.add_new(0x60000022 + overlay_index*0x20000 , 'LO', name)
    dcm_file.add_new(0x60000040 + overlay_index*0x20000 , 'CS', 'R')
    dcm_file.add_new(0x60000050 + overlay_index*0x20000 , 'SS', [1,1])
    dcm_file.add_new(0x60000100 + overlay_index*0x20000 , 'US', 1)
    dcm_file.add_new(0x60000102 + overlay_index*0x20000 , 'US', 0)
    dcm_file.add_new(0x60003000 + overlay_index*0x20000 , 'OW', w)

# Source: https://stackoverflow.com/questions/44865023/circular-masking-an-image-in-python-using-numpy-arrays
def create_circular_mask(h, w, centerX, centerY, radius):
    # Grid of indices
    Y, X = np.ogrid[:h, :w]
    # Grid of distances
    dist_from_center = np.sqrt((X - centerX)**2 + (Y-centerY)**2)
    # Grid of 1s and 0s for distances <= radius
    mask = np.array([[1 if dist <= radius else 0 for dist in dist_row] for dist_row in dist_from_center])
    return mask

# Returns the list of DICOM files for the patient sorted by Z position
def load_scan(paths):
   paths2 = [s for s in paths if '.dcm' in s and hasattr(dicom.read_file(s, force=True), "InstanceNumber") and hasattr(dicom.read_file(s, force=True), "ImagePositionPatient")]
   slices = [dicom.read_file(s, force=True) for s in paths2]
   slices2 = [x for x in slices if hasattr(x, "InstanceNumber") and hasattr(x, "ImagePositionPatient")]
   slices=slices2
   indices = sorted(list(range(len(paths2))), key=lambda x: slices[x].InstanceNumber)
   sortedPaths = np.array(paths)[indices]
   slices.sort(key=lambda x: int(x.InstanceNumber))
   slices2 = {float(str(x.ImagePositionPatient[2]).rstrip('0')): x for x in slices}
   pathDict = {float(str(x.ImagePositionPatient[2]).rstrip('0')): sortedPaths[count] for count, x in enumerate(slices)}
   return slices2, slices[0].SeriesInstanceUID

def add_row_to_csv_file(csvpath, seriesuid, x, y, z, diam_mm):
    columns = ["seriesuid", "coordX", "coordY", "coordZ", "diameter_mm"]
    with open(csvpath, 'ab') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([seriesuid, x, y, z, diam_mm])

def area_to_diameter(area):
    return 2 * np.sqrt(area / np.pi)

def world_to_voxel(x, y, spacing, origin):
    vY = abs((-origin[1] + y)*(1/spacing[1]))
    vX = abs((-origin[0] + x)*(1/spacing[0]))
    return int(np.round(vX)), int(np.round(vY))

if __name__ == "__main__":
    import_from_csvs()
    vc.needsDisplayUpdate()