
## Results
> For the analyses on the ProTherm&lowast; dataset see [O2567_analyses](notes/O2567_analyses.md)

> For the initial analyses on the ProTherm&lowast; dataset see [ProTherm-star_initial_analyses](notes/ProTherm-star_initial_analyses.md)

> For the second analyses on the ProTherm&lowast; dataset see [ProTherm-star_second_analyses](notes/ProTherm-star_analyses.md)


The results were consistent between datasets, therefore herein are reported the final analyses.
Also the alteration of settings resulted in only mild difference in scores:

![errors](../images/O_distro_errors.png)

One issue is that there are severe outliers. The fences under most settings are less than ±4 kcal/mol,
but the outliers can go upto 126 kcal/mol.
This is driven by the inability to un-distort a clash requiring large scale backbone corrections.
Venus frontend shows ∆∆G greater than 10 kcal/mol as "> 10 kcal/mol".

![outliers](../images/O_outliers.png)

For analytical purposes, median is used herein instead of mean because:

* it is not sensitive to outliers —especially given that the amount of values beyond 5 kcal/mol are unimportant
* some samples show large deviation, non-zero skew and non-three kurtosis

Standard deviation is based on the mean, the median equivalent is the
[median absolute deviation (MAD)](https://en.wikipedia.org/wiki/Median_absolute_deviation).
MAD in a normal distribution does not match the deviation, as it is ~0.67 the size of the latter 
—median and mean differ depending on skewness (positive skew => lower median than mean), which likewise affect MAD.

Pearson &rho; is the ratio of covariance over the product of deviations and
I cannot manage to rigorously/justifiably [convert it to a median form](median_Pearson.md). 

![O_conditions](../images/O_conditions.png)

The percentage of samples correctly allocated to the >2 kcal/mol bin is roughly 75% for all settings used,
The median absolute error is around 0.9 kcal/mol and the mean absolute error is 1.3 kcal/mol.
The major discriminant is time.

![O_time](../images/O_time.png)

This utterly eliminates cartesian scorefunction and the best scorefunction appears to the ref2015 series.
Based on the split of different conditions, ref2015 12 Å x2 appears ideal.