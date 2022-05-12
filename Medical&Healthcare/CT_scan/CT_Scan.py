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
