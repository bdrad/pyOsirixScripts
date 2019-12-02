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

csv_file_path  = "./Users/admin/Downloads/annotations.csv"
DEFAULT_DIAM   = 6 # in milimeters
DEFAULT_COLOR = (252, 15, 41)
DEFAULT_OPACITY = 10.691605567932129
DEFAULT_THICKNESS = 0.5

def embed_rois():
    vc   = osirix.frontmostViewer()
    pixs = vc.pixList()
    rois = vc.roiList(movieIdx=0)
    for i in range(len(rois)):
        if len(rois[i]):
            print(pixs[i].sourceFile)
            ds = dicom.read_file(pixs[i].sourceFile)
# Supported ROIs are types that support the getMapFromROI() function. This pixelwise bitmask can be directly put into the DICOM header.
# Unsupported ROIs include ovals, points, and rectangles. For ovals we embed a circle based on the average of the major and minor axis. For points we infer a 6mm diameter. For rectangles, we embed the rectangle into the 
            supported_rois = ()
            unsupported_rois = ()
            for j in range(len(rois[i])):
                print(i, j)
                print("roi type: ", rois[i][j].type)
                print("Centroid: ", rois[i][j].centroid())
                if (rois[i][j].type == 'tPlain') or (rois[i][j].type == 'tPencil') or (rois[i][j].type == 'tCPolygon') or (rois[i][j].type == 'tOPolygon'):
                    supported_rois += rois[i][j],
                else:
                    unsupported_rois += rois[i][j],
            print("\n\nEmbedding rois…")
            for j in range(len(supported_rois)):
                roi = supported_rois[j]
                # roi.Area() is units of cm^2. We want units of mm^2 because we want mm diameter.
                # Thus we mutiply area by the square of mm per cm, which is 10^2 = 100.
                roiDiam = 2*area_to_radius(roi.roiArea()*100) 
                if roiDiam == 0:
                    roiDiam = DEFAULT_DIAM
                mask = get_bitmask_from_roi(pixs[i], roi)
                # Embed pixel mask into DICOM Header
                embed(mask, ds, j, roi.name+'/'+OVERLAY_PYOSIRIX_NAME_SUFFIX)
                # Add new row to CSV
                wX, wY = voxel_to_world(roi.centroid()[0], roi.centroid()[1], ds.PixelSpacing, ds.ImagePositionPatient)
                add_row_to_csv_file(csv_file_path, ds.SeriesInstanceUID, wX, wY, str(ds.ImagePositionPatient[2]).strip("\'"), roiDiam)
                beautify_roi(roi)
                print("\n\n")
    
            for j in range(len(supported_rois), len(supported_rois) + len(unsupported_rois)):
                roi = unsupported_rois[j - len(supported_rois)]
                roiDiam = 2*area_to_radius(roi.roiArea()*100)
                if roiDiam == 0:
                    roiDiam = DEFAULT_DIAM
                centX, centY = roi.centroid()[0], roi.centroid()[1]
                # For points and ovals, draw a circle using the radius.
                if len(roi.points) == 1 or roi.type == 't2DPoint' or roi.type == 'tOval':
                    mask = create_circular_mask(512, 512, roi.centroid()[0], roi.centroid()[1], radius=(roiDiam/2)*(1/ds.PixelSpacing[0]))
                elif roi.type == 'tROI': # Rectangle
                    centX, centY = np.mean(roi.points[:,0]), np.mean(roi.points[:, 1])
                    print("tROI:",centX, centY)
                    mask = np.zeros((512, 512)).astype(int)
                    points = roi.points
                    xs = points[:, 0]
                    ys = points[:, 1]
                    for x_i in range(int(np.round(min(xs))), int(np.round(max(xs)))):
                        for y_i in range(int(np.round(min(ys))), int(np.round(max(ys)))):
                            mask[x_i, y_i] = 1
                    mask = np.swapaxes(mask, 0, 1)
                else:
                    print("Ignoring an Unsupported ROI type: %s\n\n" % roi.type)
                    continue
                    
                embed(mask, ds, j, roi.name+'/'+OVERLAY_PYOSIRIX_NAME_SUFFIX)
                beautify_roi(roi)
                # Add new row to CSV
                wX, wY = voxel_to_world(centX, centY, ds.PixelSpacing, ds.ImagePositionPatient)
                add_row_to_csv_file(csv_file_path, ds.SeriesInstanceUID, wX, wY, str(ds.ImagePositionPatient[2]).strip("\'"), roiDiam)
                print("\n\n")
            ds.save_as(pixs[i].sourceFile)
            print("saved file to sourceFile %s \n\n\n" % pixs[i].sourceFile)

##########################
#         UTILS          #
##########################
OVERLAY_PYOSIRIX_NAME_SUFFIX = "ROI embedded using pyOsirix"
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

def add_row_to_csv_file(csvpath, seriesuid, x, y, z, diam_mm):
    columns = ["seriesuid", "coordX", "coordY", "coordZ", "diameter_mm"]
    with open(csvpath, 'ab') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([seriesuid, x, y, z, diam_mm])
        print("wrote to csv", x, y, z, diam_mm)

def area_to_radius(area):
    return (np.sqrt(area / np.pi))

def get_bitmask_from_roi(pix, roi):
	mask = np.swapaxes(pix.getMapFromROI(roi).astype(np.uint8), 0, 1)
	return mask
	
def voxel_to_world(vX, vY, spacing, origin):
    wX = (vX * spacing[0]) + origin[0]
    wY = (vY * spacing[1]) + origin[1]
    return wX, wY

def beautify_roi(roi):
    roi.color = DEFAULT_COLOR
    roi.thickness = DEFAULT_THICKNESS
    roi.opacity = DEFAULT_OPACITY

if __name__ == "__main__":
    embed_rois()
    vc.needsDisplayUpdate()