Single Series checks
====================

Each of these functions accepts a pandas timeseries dataframe as input, and provides a dataframe as the output.

There are five single-series checks:

Outliers
--------

``rain_outliers()`` function determines, for each time step, the ratio of the rain amount to the 99th percentile of the non-zero observations.
A large positive number indicates the observation is very much greater than the 99th percentile, and so is likely to be an outlier. This number may be considered an outlier index.

Impossibles
-----------

The ``impossibles()`` function checks for negative numbers and non-numbers in the rainfall amount values.

It returns a boolean timeseries where TRUE indicates an impossible value.

Date-Time Issues
----------------

The ``DateTimeIssues()`` function checks for any duplicate date-times.

It returns a boolean timeseries where TRUE indicates a duplicate.

Dry Spells
----------

The ``DrySpells()`` function calculates the length of time (in days) of the dry spell that each zero rain observation is within. This function should be applied to daily data.

Repeated Values
---------------

The ``RepeatedValues()`` function returns the number of items of a repeated sequence that each rain observation is within.

Homogeneous Values
------------------

The ``Homogeneity()`` function tests for unusual changes in the rainfall timeseries. It uses the pyHomogeneity package and applies the Pettitt non-parametric change-point detection algorithm at the annual level.

It returns a boolean time series, with TRUE for the latest homogeneous section of the rain data time series.

Sub-freezing Rain
-----------------

The ``SubFreezingRain()`` function identifies rain observations when the temperature is below 0 degrees C.

This function requires a second input time series with temperautre (in degrees C).

It returns a boolean time series where TRUE is when the temperature is below 0 degrees C.
