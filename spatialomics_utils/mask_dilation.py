import scipy
import numpy as np


def dilate_mask(
        mask,
        radius
    ):
    reversed_img = mask == 0

    distance, indices = scipy.ndimage.distance_transform_edt(
        reversed_img,
        return_distances=True,
        return_indices=True
    )

    within_radius = distance <= radius
    mask[
        within_radius
    ] = mask[
        indices[0][within_radius],
        indices[1][within_radius]
    ]
    return mask


def dilate_coords(
        spatial_coords,      # [[y, x], [y, x], ...]
        labels,
        radius: int,
        data_type=np.uint16,
    ):
    assert spatial_coords.dtype == 'int'
    assert spatial_coords.shape[0] == labels.shape[0]

    ## prepare blank img frame from spatial coords
    mask = np.zeros(
        shape=spatial_coords.max(axis=0) + 1,
        dtype=data_type
    )
    mask[spatial_coords[:, 0], spatial_coords[:, 1]] = labels
    
    mask = dilate_mask(
        mask,
        radius
    )
    return mask[spatial_coords[:, 0], spatial_coords[:, 1]]