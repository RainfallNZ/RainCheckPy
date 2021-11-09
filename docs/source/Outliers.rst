Outliers
========

``rain_outliers()`` function determines, for each time step, the ratio of the rain amount to the 99th percentile of the non-zero observations.
A large positive number indicates the observation is very much greater than the 99th percentile, and so is likely to be an outlier. This number may be considered an outlier index.
