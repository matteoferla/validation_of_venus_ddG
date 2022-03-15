## Settings files

```python
import json
# code fixes the defaults, but I want a table
default = dict(radius=12, 
               cycles=2, 
               scorefxn_name='ref2015',
               prevent_acceptance_of_incrementor=True,
               neighbour_only_score=False,
               outer_constrained=False,
               single_chain=True, 
               remove_ligand=True, 
               use_pymol_for_neighbours=False)

settings = {'8x1_reg':  {**default, 'radius': 8, 'cycles': 1},
            '12x1_reg': {**default, 'radius': 12, 'cycles': 1},
            '12x2_reg': {**default, 'radius': 12, 'cycles': 2},
            '12x3_reg': {**default, 'radius': 12, 'cycles': 3},
            '15x3_reg': {**default, 'radius': 15, 'cycles': 3},
            '12x2_con':  {**default, 'radius': 12, 'cycles': 2, 'outer_constrained': True},
            '12x2_neigh':  {**default, 'radius': 12, 'cycles': 2, 'neighbour_only_score': True},
            '12x2_con':  {**default, 'radius': 12, 'cycles': 2, 'prevent_acceptance_of_incrementor': False},
            '12x2_cart': {**default, 'radius': 12, 'cycles': 2, 'scorefxn_name': 'ref2015_cart'},
            '12x2_beta': {**default, 'radius': 12, 'cycles': 2, 'scorefxn_name': 'beta_nov16'},
            '12x2_betacart': {**default, 'radius': 12, 'cycles': 2, 'scorefxn_name': 'beta_nov16_cart'},
           }

with open('settings.json', 'w') as fh:
    json.dump(obj=settings, fp=fh)
```

For :octocat: output

```python
import os, json


def get_data(filename):
    with open('settings.json') as fh:
        return json.load(fh)

import pandas as pd

print(pd.DataFrame({k: v for filename in os.listdir() if 'settings.json' in filename for k, v in get_data(filename).items()}).transpose().to_markdown())
```