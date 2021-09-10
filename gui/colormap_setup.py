import numpy as np
from silx.gui.colors import Colormap
def setup_colormap(data):
    try:
        data = data.astype(np.float)
        vmin = np.nanmedian(data) - np.nanstd(data)*3
        if vmin < np.nanmin(data):
            vmin = np.nanmin(data)
        vmax = np.nanmedian(data) + np.nanstd(data)*3
        if vmax > np.nanmax(data):
            vmax = np.nanmax(data)
        colormap = Colormap(name  = 'gray',
                    normalization = 'linear',
                    vmin          = vmin,
                    vmax          = vmax
                    )
    except Exception as e:
        print(e)
        colormap = Colormap(name = 'gray',
                    normalization= 'linear',
                    vmin         = 0,
                    vmax         = 1,
                    )
    return colormap

