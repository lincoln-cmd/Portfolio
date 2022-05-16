# Introduction to 3D medical imaging for machine learning: preprocessing and augmentations

# 3D medical imaging for machine learning
# - preprocessing
# - augmentations


import matplotlib as plt


# The images will be shown in 3 planes: sagittal, coronal, axial looking from left to right.


# Two-dimensional planes visualization

def show_mid_slice(img_numpy, title = 'img'):
    '''
     Accepts an 3D numpy array and shows median slices in all three planes.
    '''
    assert img_numpy.ndim == 3
    n_i, n_j, n_k = img_numpy.shape
    
    # sagittal (left image)
    center_i1 = int((n_i - 1) / 2)
    # coronal (center image)
    center_j1 = int((n_j - 1) / 2)
    # axial slice (right image)
    center_k1 = int((n_k - 1) / 2)
    
    show_slices([img_numpy[center_i1, :, :], img_numpy[:, center_j1, :], img_numpy[:, : center_k1]])
    
    plt.suptitle(title)
    

def show_slices(slices):
    '''
     Function to display a row of image slices. 
     Input is a list of numpy 2D image slices.
    '''
    fig, axes = plt.subplots(1, len(slices))
    for i, slice in enumerate(slices):
        axes[i].imshow(slice.T, cmap = 'gray', origin = 'lower')
        

# show_mid_slice(epi_img_numpy, 'first image')
# show_mid_slice(anatomy_img_numpy, 'second image')

