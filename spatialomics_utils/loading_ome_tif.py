import zarr
import tifffile         
import numpy as np


def lazy_load(
        tif_path: str, 
        level: str, 
        ystart: int, 
        yend: int, 
        xstart: int, 
        xend: int
    ) -> np.ndarray:
    """
    Lazily load an OME-TIFF file using zarr.

    Parameters
    ----------
    tif_path : str
        Path to the OME-TIFF file.
    level : str
        The level of the pyramid to load.
    ystart : int
        The starting y-coordinate.
    yend : int
        The ending y-coordinate.
    xstart : int
        The starting x-coordinate.
    xend : int
        The ending x-coordinate.

    Returns
    -------
    np.ndarray
        An array representing the OME-TIFF data.
    """
    try:
        store = tifffile.imread(
            tif_path,
            return_as='zarr'
        )
        lazy_img = zarr.open(store, mode='r')

    except Exception as e:
        print(f"Error loading OME-TIFF file: {e}")
        print("Suggested version: tifffile==2026.6.1 zarr==3.2.1")
        return None

    # Check the groups
    print("Ome TIF Group", lazy_img.tree()) 

    square_crop = lazy_img[level][..., ystart:yend, xstart:xend]
    return square_crop
