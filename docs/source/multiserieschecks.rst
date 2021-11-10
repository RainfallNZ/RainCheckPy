Multi-series Checks
===================

All the multi-series check functions take two pandas timeseries dataframes as their input.

affinity
--------

The ``affinity()`` function determines how similar two series are. This is used to determine if a site is suitable for comparisons, but it can also be used to see if a data series is an exact duplicate.
The ``affinity()`` function is based on the implementation of Lewis et al. (2018). It is the likelihood of both sites having the same indication of rainfall ocurrence (i.e. both show rain or both show not-rain).
It returns a number from 0 (no affinity) to 1 (identical).
An affinity greater than 0.9 is expected for sites in similar areas

Spearman Rank Correlation
-------------------------

The ``spearman()`` function calculates the Spearman rank correlation. This is a nonparametric test for rank correlation. This may be used to determine if a site is suitable for comparisons, but it can also be used to see if a data series is an exact duplicate.
It returns a number from -1 (oposite ranking) to 0 (the ranks do not correlate at all) to 1 (identical).
A rank correlation over about 0.65 is high.

Neighborhood Divergence
-----------------------

The ``neighborhoodDivergence()`` function calculates how different an observation at the test site is to the reference site.
For each time step the difference in rainfall is found between the test and reference site. The 95th percentile of these differences is determined. The NeigborhoodDivergence index is the ratio of the positive difference to the 95th percentile or the negative differences to the 5th percentile. The numbers range from 0 to inf. A value of 1, indicates the difference between the observations at the test and reference sites is equal to the outer 5th percentile of all test-to-reference differences. A value of 4 indicates the difference is four times the outer 5th percentile The function generates a dataframe with two columns, "HighCFOutliers" and "LowCFOutliers".

Dry Spell Divergence
--------------------

The ``DrySpellDivergence()`` function calculates how dry the test site is compared to the reference site.
For both sites the proportion-of-dry-spells-over-the-previous-15-days is calculated. The difference of these proportions between the test and reference site is calculated, and the 95th perentile of that series is found.
For each time step, the ratio of the difference to the 95th percentile is determined.
This is the DrySpellDivergence index.
The possible values range from 0 to inf. A 0 indicates no dry spells in the 15 days. A value of 1 indicates the dry-spell-proportion difference equals the 95th percentile. A very large number indicates the dry-spell-proportion difference is much larger than the 95th percentile. The index is only calculated if there are observations in both the test and refernce site throughout the 15 days. The function generates a dataframe of the DrySpellDivergence index.
