import skimage


def mask_he_img(img):
    reversed_img = 255 - img
    reversed_img = reversed_img.max(axis=2)
    thresh = skimage.filters.threshold_triangle(
        reversed_img
    )
    in_tissue_mask = reversed_img > thresh
    return in_tissue_mask