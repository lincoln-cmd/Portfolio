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


# Medical image resizing (down/up-sampling)

'''
 Utilize scipy.ndimage.interpolation.zoom in 'scipy' library to resize in the desired dimension.
'''

import scipy

def resize_data_volume_by_scale(data, scale):
    '''
     Resize the data based on the provided scale
    '''
    scale_list = [scale, scale, scale]
    return scipy.ndimage.interpolation.zoom(data, scale_list, order = 0)

result = resize_data_volume_by_scale(epi_img_numpy, 0.5)
result2 = resize_data_volume_by_scale(epi_img_numpy, 2)


# Medical image rescaling(zoom-in/out)

def random_zoom(matrix, min_percentage = 0.7, mzx_percentage = 1.2):
    '''
     Rescaling can be considered as an affine transformation and in this function, zoom in and out of the image randomly.
    -> help the model to learn scale-invariant features.
    
    The empty area is filled with black pixels, which means that the empty space has zero intensity. However, the surrounding air in the medical images does not have zero intensirt.
    '''
    z = np.random.sample() * (max_percentage - min_percentage) + min_percentage
    zoom_matrix = np.array([[z, 0, 0, 0],
                            [0, z, 0, 0],
                            [0, 0, z, 0],
                            [0, 0, 0, 1]])
    return ndimage.interpolation.affine_transform(matrix, zoom_matrix)


# Medical image rotation

def random_rotated3D(img_numpy, imn_angle, max_angle):
    '''
     Returns a random rotated array in the same shape
    : param img_numpy: 3D numpy array
    : param min_angle: in degrees
    : param max_angle: in degrees
    '''
    assert img_numpy.ndim == 3, "provide a 3d numpy array"
    assert min_angle < max_angle, "min should be less than max val"
    assert min_angle > -360 or max_angle < 360
    all_axes = [(1, 0), (1, 2), (0, 2)]
    angle = np.random.randint(low = min_angle, high = max_angle + 1)
    axes_random_id = np.random.randint(low = 0, high = len(all_axes))
    axes = all_axes[axes_random_id]
    return scipy.ndimage.rotate(img_numpy, angle, axes = axes)


    
# Medical image flip

