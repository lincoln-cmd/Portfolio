'''
<CT intensities and Hounsfield units>

 In the 'Hounsfield scale, the intensity of air is -1000 and one of water is 0. This is an absolute scale.
 - lung area: around -500 HU
 - fats: -200 ~ -100 Hu
 - soft tissue: 30 ~ 45 HU
 - bone area: excess 500 HU
'''



'''
<CT data visualization: level and window>

 The medical image convention to clip the Housenfield range is by choosing a central intensity, called level and a window.

max = level + window / 2
min = level - window / 2 

'''



import matplotlib.pyplot as plt
import numpy as np

def show_slice_window(slice, level, window):
    '''
    Function to display an image slice
    "Input is a numpy 2d array
    '''
    max = level + window / 2
    min = level - window / 2
    slice = slice.clip(min, max)
    plt.figure()
    plt.imshow(slice.T, cmap = 'gray', origin = 'lower')
    plt.savefig('L' + str(level) + 'W' + str(window))
    
    
# lung segmentation based on intensity values
import nibabel as nib
ct_img = nib.load('/workspace/free/CT_Images/slice001.nii.gz')
print(ct_img.header)


# 1. find pixel dimensions to calculate the area in mm^2

def find_pix_dim(ct_img):
    '''
    Get the pixdim of the CT image.
    A general solution that gets the pixdim indicated from the image dimensions. From the last 2 image dimensions, we get their pixel dimension.
    Args: ct_img: nib image
    
    Returns: List of the 2 pixel dimensions
    '''
    pix_dim = ct_img.header['pixdim'] # example [1, 2, 1.5, 1, 1] 
    dim = ct_img.header['dim'] # example [1, 512, 512, 1, 1]
    max_indx = np.argmax(dim)
    pixdimX = pix_dim[max_indx]
    dim = np.delete(dim, max_indx)
    pix_dim = np.delete(pix_dim, max_indx)
    max_indy = np.argmax(dim)
    pixdimY = pix_dim[max_indy]
    return [pixdimX, pixdimY] #example [2, 1.5]

# 2. Binarize image using intensity thresholding

# 3. Contour finding
def intensity_seg(ct_numpy, min = 1000, max = -300):
    slipped = clip_ct(ct_numpy, min, max)
    return measure.find_contours(clipped, 0.95)

# 4. Find the lung area from a set of possible contours
'''
 extract a convex polygon from the contour using scipy.
 2 constraints:
  - the contour of the lungs must be a close dset(always true)
  - the contour must have a minimum volume of 2000pixels to represent the lung.
 '''

def find_lugs(contours):
    '''
    Chooses the contours that correspond to the lungs and the body. First, we exclude non-closed sets-contours. Then we assume some min area and volume to exclude small contours. Then the body is excluded as the highest volume closed set. The remaining areas correspond to the lungs.
    Args:
        contours: all the detected contours
    
    Returns: contours that correspond to the lung area.
    '''
    body_and_lung_contours = []
    vol_contours = []
    
    for contour in contours:
        hull = ConvexHull(contour)
        
        # set some constraints for the volume
        if hull.volume > 2000 and set_is_closed(contour):
            body_and_lung_contours.append(contour)
            vol_contours.append(hull.colume)
            
    # Discard body contour
    if len(body_and_lung_contours) == 2:
        return body_and_lung_contours
    elif len(body_and_lung_contours) > 2:
        vol_contours, body_and_lung_contours = (list(t) for t in zip(*sorted(zip(vol_contours, body_and_lung_contours))))
        body_and_lung_contours.pop(-1) # discard body
        return body_and_lung_contours # only lungs left
    
# 5. Contour to binary mask

'''
 Utilize pillow python library that draws a polygon and creates a binary image mask to save nifty file.
'''

import numpy as np
from PIL import Image, ImageDraw

def create_mask_from_polygon(image, contours):
    '''
    Creates a binary mask with the dimensions of the image and converts the list of polygon-contours to binary masks and merges them together.
    Args:
        image: the image that the contours refer to
        contours: list of contours
        
    Returns:
    '''
    
    lung_mask = np.array(Image.new('L', image.shape, 0))
    for contour in contours:
        x = contour[:, 0]
        y = contour[:, 1]
        polygon_tuple = list(zip(x, y))
        img = Image.new('L', image.shape, 0)
        ImageDraw.Draw(img).polygon(polygon_tuple, outline = 0, fill = 1)
        mask = np.array(img)
        lung_mask += mask
        
    lung_mask[lung_mask > 1] = 1 # sanity check to make 100% sure that the mask is binary
    return lung_mask.T # transpose it to be aligned with the image dims


def save_nifty(img_np, name, affine):
    '''
    Binary masks should be converted to 255 so it can be displayed in a nii viewer. We pass the affine of the inital image to make sure it exits in the same image coorinate space.
    Args:
        img_np: the binary mask
        name: output name
        affine: 4x4np array
    
    Returns:
    '''
    
    img_np[img_np == 1] = 255
    ni_img = nib.Nifti1Image(img_np, affine)
    nib.save(ni_img, name + '.nii.gz')
    
    
'''
 Download free medical images: https://www.aliza-dicom-viewer.com/download
'''

# Segment the main vessels and compute the vessels over lung area ratio

