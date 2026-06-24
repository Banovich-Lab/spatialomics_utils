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
    ):
    assert spatial_coords.shape[0] == labels.shape[0]

    unknown_idx = np.where(labels == 0)[0]
    unknown_coords = spatial_coords[unknown_idx]

    known_idx = np.where(labels != 0)[0]
    known_coords = spatial_coords[known_idx]

    tree = scipy.spatial.cKDTree(known_coords)
    distances, nearest_pos = tree.query(unknown_coords, k=1)

    nearest_known_idx = known_idx[nearest_pos]
    nearest_known_labels = labels[nearest_known_idx]
    nearest_known_labels[
        distances > radius
    ] = 0

    dilated_labels = labels.copy()
    dilated_labels[unknown_idx] = nearest_known_labels
    return dilated_labels