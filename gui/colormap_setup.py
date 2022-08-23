import numpy as np
from silx.gui.colors import Colormap
def setup_colormap(data,vmin=None,vmax=None,normalization='linear'):
    try:
        data = data.astype(np.float)
        if isinstance(vmin,type(None)):
            vmin = np.nanmedian(data) - np.nanstd(data)*3
            if vmin < np.nanmin(data):
                vmin = np.nanmin(data)
        else:
            vmin = vmin
        if isinstance(vmax,type(None)):
            vmax = np.nanmedian(data) + np.nanstd(data)*3
            if vmax > np.nanmax(data):
                vmax = np.nanmax(data)
        else:
            vmax = vmax
        if np.isnan(vmin):
            vmin = 0
        if np.isnan(vmax):
            vmax = 1
        colormap = Colormap(name  = 'inferno',
                    normalization = normalization,
                    vmin          = vmin,
                    vmax          = vmax
                    )
    except Exception as e:
        print(e)
        colormap = Colormap(name = 'inferno',
                    normalization= normalization,
                    vmin         = 0,
                    vmax         = 1,
                    )
    return colormap

