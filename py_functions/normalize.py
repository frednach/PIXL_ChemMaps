import numpy as np

# define function to normalize oxide values to 0...1 and flatten the extreme 2 percentiles on each end
def normalize_percentile(series, low=2, high=98):
    #map extreme values below "low"[%] threshold to 0
    vmin = np.nanpercentile(series, low)
    #map extreme values above "high"[%] threshold to 1
    vmax = np.nanpercentile(series, high)

    normalized = (series - vmin) / (vmax - vmin)

    return np.clip(normalized, 0, 1)