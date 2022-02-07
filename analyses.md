# Analysis

For previous iterations of analysis see [notes](notes).
These will contain specific tests _etc._, whereas this is the code for a summary table.

## Function definitions

## Import

```python3
import json
import pandas as pd

pd.options.mode.chained_assignment = None
import numpy as np
from IPython.display import display, HTML


def entitle(*args: str, header: str = 'h1'):
    content = ' &mdash; '.join(map(str, args))
    display(HTML(f'<{header}>{content}</{header}>'))


def read_db(db_name: str, scaled: bool = False) -> pd.DataFrame:
    from sqlitedict import SqliteDict
    db_scores = SqliteDict(db_name, encode=json.dumps, decode=json.loads, autocommit=True)
    dict_scores = {k: v for k, v in db_scores.items() if isinstance(v, dict)}
    scores = pd.DataFrame.from_dict(dict_scores, orient='index')
    scores['has_pro'] = scores.mutation.str.contains('P')
    scores['time'] = scores.time.astype(float)
    scores['ddG'] = scores.ddG.astype(float)
    if not scaled:
        # scores['ddG'] is not 0.239 scaled kJ <=> kcal
        # 0.239 scaled kJ <=> kcal
        scores['ddG'] = scores.ddG * 0.239
    return scores
```

### Outlier elimination

```python3
from typing import Tuple
import pandas as pd


def wo_outliers(scores: pd.DataFrame, k: float = 1.5) -> Tuple[pd.DataFrame, Tuple[float, float]]:
    """
    Tukey fence.
    """
    q25 = scores.ddG.quantile(0.25)
    q75 = scores.ddG.quantile(0.75)
    iqr = q75 - q25
    fences = (q25 - k * iqr, q75 + k * iqr)
    return scores.loc[(fences[0] < scores.ddG) & (scores.ddG < fences[1])], fences
```

## Stats

```python3
import numpy as np
import pandas as pd


def get_confusion_allocation(row: pd.DataFrame, cutoff: float = 0, precision: int = 1) -> str:
    """
    https://academic.oup.com/bib/article/22/6/bbab184/6289890#312128384
    uses inverted potentials, but appears to be a typical confusion matrix
    
    This could have been done without a function.
    The perfect zero is impossible for empirical, but not for
    calculated, which has generally 1 digit precision.
    One option is adding `np.finfo(float).eps` but that is a bit weird.
    Anyway there are only 5 .2f values, and 54 .1f values in `experimental_ddG`
    So it should not affect the matrix.
    """
    calc = np.sign(np.round(row.ddG, precision) - cutoff)
    emp = np.sign(np.round(row.experimental_ddG, precision) - cutoff)
    # np.sign gives -1, 0, +1
    if np.isnan(row.ddG):
        return 'NA'
    elif calc + emp == 2:
        return 'TP'
    elif calc + emp == -2:
        return 'TN'
    elif calc == 0 and emp == 0:
        return 'TN'
    elif calc == 1:
        return 'FP'
    else:
        return 'FN'


def get_confused(experimental_ddG: pd.Series, calculated_ddG: pd.Series, cutoff) -> dict:
    from scipy.stats import contingency
    from functools import partial
    confusion = (pd.DataFrame({'experimental_ddG': experimental_ddG,
                               'ddG':              calculated_ddG})
                 .apply(partial(get_confusion_allocation, cutoff=cutoff),
                        axis='columns')
                 .value_counts()
                 .to_dict())
    confusion_matrix = np.array([[confusion[truth + sign] for sign in 'PN'] for truth in 'TF'])
    return dict(
        sensitivity=confusion['TP'] / (confusion['TP'] + confusion['FN']),
        specificity=confusion['TN'] / (confusion['TN'] + confusion['FP']),
        false_discovery_rate=confusion['FP'] / (confusion['FP'] + confusion['TP']),
        fallout=confusion['FP'] / (confusion['FP'] + confusion['TN']),
        accuracy=(confusion_matrix[0, :].sum()) / (confusion_matrix[:].sum()),
        phi_coefficient=contingency.association(confusion_matrix, method='pearson')
    )


def get_linregression(experimental_ddG: pd.Series, calculated_ddG: pd.Series) -> dict:
    """
    Note that the intercept is not set to zero.
    """
    from scipy.stats import linregress
    lin = linregress(calculated_ddG, experimental_ddG)
    return dict(slope=lin.slope,
                slope_se=lin.stderr,
                intercept=lin.intercept,
                intercept_se=lin.intercept_stderr,
                )


def get_general(experimental_ddG: pd.Series, calculated_ddG: pd.Series) -> dict:
    from scipy.stats import pearsonr
    return dict(pcc=pearsonr(experimental_ddG, calculated_ddG)[0],
                rmse=sum((experimental_ddG - calculated_ddG) ** 2 / len(experimental_ddG)) ** 0.5,
                mae=(experimental_ddG - calculated_ddG).abs().mean(),
                mdae=(experimental_ddG - calculated_ddG).abs().median(),
                )


def get_ddG_stats(experimental_ddG: pd.Series, calculated_ddG: pd.Series) -> dict:
    add_suffix = lambda d, s: {k + s: v for k, v in d.items()}
    confused_zero = get_confused(experimental_ddG, calculated_ddG, cutoff=0)
    confused_two = get_confused(experimental_ddG, calculated_ddG, cutoff=1.9)
    return {
        **get_general(experimental_ddG, calculated_ddG),
        **get_linregression(experimental_ddG, calculated_ddG),
        **add_suffix(confused_zero, '_@0'),
        **add_suffix(confused_two, '_@+2')
    }
```

## Run

```python3

# -----------------------------------------------------------------------------
scores = read_db('O2567-redux-scores.db')  # noqa
condition = '12x2_reg'
subscores = scores.loc[scores.condition == condition]
inscores, fences = wo_outliers(subscores)  # noqa
n_outliers = len(subscores) - len(inscores)
print(f'fences @ k=1.5, {fences}, {n_outliers} outliers out of {len(subscores)}')

details = {'dataset':     'O2567',
           'condition':   condition,
           'n_outliers':  n_outliers,
           'lower_fence': fences[0],
           'upper_fence': fences[1],
           **get_ddG_stats(inscores.experimental_ddG, inscores.ddG),  # noqa
           **{'time_' + k: v for k, v in subscores.time.describe().to_dict().items()}
           }
```