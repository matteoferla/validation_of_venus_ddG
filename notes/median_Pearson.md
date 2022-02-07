## Test: median Pearson

Ignoring outliers (i.e. everything beyond the fence) or setting an arbitrary cutoff,
would skew the data for the calculation of the Pearson &rho;.

The Pearson &rho; is the ratio of the covariance over the product of the standard deviations.
A standard deviation is the root of the variance, which is 
the mean of the square of the difference of the residual of each value from the mean.
The covariance is calculated like a regular variance,
but with the term averaged is the product of the residuals from the respective means.
Here I am talking about vector vs. vector, so the covariance is a number or
pedantically a 2x2 matrix with tha same value repeat twice in the antidiagonal and the two different variances in the diagonal
—I mention this because covariance matrix is generally what is discussed in places.

Now, stuff gets funky with the median. The median of the squares of the residuals gives a rubbish value,
hence the median absolute deviation.

However, I tried naïvely extending this to covariance and that proved tricky.
There is a paper online about
median convariant matrices (MCM) but it's rather intense and intended for matrices.

Using the same formula blindly with the product results in a mess:
the problem stems from using the absolute of the product.
This leaves two options, not using the absolute as if it were a vanilla (co)variance,
or using the square root of the product and squaring the median of that array, 
which assumes the residuals are as they should be for a vanilla variance.

The covariance of random values has got to be zero. 

```python
x, y = np.random.standard_normal(1000), np.random.standard_normal(1000)
correct:float = np.cov(x, y)[0][1]
naive:float = np.median(  np.abs( (x - np.median(x)) * (y - np.median(y)) )  )
hacked:float = np.median(  ( (x - np.median(x)) * (y - np.median(y)) )  )
hacked2:float = np.median(  ( np.abs((x - np.median(x)) * (y - np.median(y)) )**1/2)  )**2
results = {'regular': correct, 'naive': expected, 'hacked': hacked, 'double-hacked': hacked2}
print(pd.Series(results).round(2).to_markdown())
```

|               |     0 |
|:--------------|------:|
| regular       |  0.02 |
| naive         |  0.33 |
| hacked        | -0    |
| double-hacked |  0.03 |

The naive approach is screamingly wrong. But the other two are reasonable...
So diving it by the product of the deviations (ie. getting the 0,1 scaled version aka. Karl Pearson's &rho;):

```python
def median_rho(x, y, calc_cov: Callable) -> float:
    """
    Screw you Karl Pearson, ignore my outliers please.
    """
    x = np.array(x)
    y = np.array(y)
    calc_mad : Callable[[np.array], float] = lambda a: np.median( np.abs(a - np.median(a)) )
    cov = calc_cov(x, y)  # this is likely the wrong step
    return cov / (calc_mad(x) * calc_mad(y))

results = []

calc_naive = lambda x, y: np.median(  np.abs( (x - np.median(x)) * (y - np.median(y)) )  )
calc_hacked = lambda x, y: np.median(  ( (x - np.median(x)) * (y - np.median(y)) )  )
calc_hacked2  = lambda x, y: np.median(  ( np.abs((x - np.median(x)) * (y - np.median(y)) )**1/2)  )**2


for gen in (np.random.standard_normal, np.arange):
    x, y = gen(1000), gen(1000)
    results.append({'name': gen.__name__,
                    'regular': pearsonr(x,y)[0],
                    'naive': median_rho(x,y, calc_naive),
                    'hacked': median_rho(x,y, calc_hacked ),
                    'double-hacked': median_rho(x,y, calc_hacked2 )
                   })
print(pd.DataFrame(results).round(2).to_markdown())
```
|    | name            |   regular |   naive |   hacked |   double-hacked |
|---:|:----------------|----------:|--------:|---------:|----------------:|
|  0 | standard_normal |     -0.02 |    0.82 |       -0 |            0.08 |
|  1 | arange          |      1    |    1    |        1 |        15625.1  |

* The real Pearson predicts no correlation for a random sample, but does for an perfect one.
* The naive approach obviously fails with the uncorrelated random sample as seen for the covariance.
* The double hack (=square of the MAD of the root of the element-wise products) goes off the rails.
* The no-absolute hack instead predicts both correctly.

Testing that last function on the actual data, gives results consistent with expectations!

```python
def median_rho(x, y) -> float:
    """
    Screw you Karl Pearson, ignore my outliers please.
    """
    x = np.array(x)
    y = np.array(y)
    calc_mad : Callable[[np.array], float] = lambda a: np.median( np.abs(a - np.median(a)) )
    cov = np.median(  ( (x - np.median(x)) * (y - np.median(y)) )  )  # this is likely the wrong step
    return cov / (calc_mad(x) * calc_mad(y))

c = pd.pivot_table(cleanscores, values=['ddG', 'flipped_experimental_ddG'], index=['condition'], aggfunc=list)
c['rho'] = c.apply(lambda row: pearsonr(row.ddG, row.flipped_experimental_ddG)[0], 1)
c['median_rho'] = c.apply(lambda row: median_rho(row.ddG, row.flipped_experimental_ddG), 1)
print(c[['rho', 'median_rho']].sort_values('rho', ascending=False).round(2).to_markdown())
```
| condition        |   rho |   median_rho |
|:-----------------|------:|-------------:|
| cart_extra       |  0.36 |         0.41 |
| beta_12x2        |  0.34 |         0.4  |
| talaris2014_9x3  |  0.34 |         0.46 |
| radius9x3        |  0.34 |         0.36 |
| beta_9x3         |  0.32 |         0.4  |
| radius9x1        |  0.32 |         0.35 |
| talaris2014_12x1 |  0.32 |         0.37 |
| full_12x3        |  0.31 |         0.4  |
| full_12x2        |  0.31 |         0.4  |
| radius12x2       |  0.31 |         0.38 |
| radius12x3       |  0.3  |         0.35 |
| talaris2014_9x1  |  0.3  |         0.42 |
| default          |  0.29 |         0.34 |
| radius12x1       |  0.29 |         0.34 |
| incrementor      |  0.29 |         0.34 |
| full_score       |  0.28 |         0.34 |
| full_10x1        |  0.28 |         0.35 |
| con_12x3         |  0.28 |         0.39 |
| full_con_12x3    |  0.27 |         0.39 |
| radius15x1       |  0.27 |         0.31 |
| outer_con        |  0.26 |         0.4  |
| talaris2014_12x3 |  0.23 |         0.41 |
| incr_12x2        |  0.22 |         0.41 |
| full_12x4        |  0.22 |         0.39 |
| incr_12x3        |  0.22 |         0.39 |
| full_11x1        |  0.21 |         0.34 |
| radius15x3       |  0.2  |         0.34 |
| radius20x3       |  0.2  |         0.34 |
| full_20x5        |  0.19 |         0.36 |
| cartesian_12x3   |  0.12 |         0.37 |
| cartesian        |  0.12 |         0.31 |

That is the two columns are similar (the rho of the two rho is 0.35...) and the biggest difference is for the outlier-heavy
samples.

However, this is utter speculatative maths, so not usable.

