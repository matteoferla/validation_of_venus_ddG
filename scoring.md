# Scoring

For previous iterations of scoring see [notes](notes).
In this version the major differences are:

* the reference tables are made consistent as opposed to the analyses changing the names for them
* More distributed runs



## Rectification

The standaridised datasets end in `.standarised.csv` (NB :uk: not :us:) and are in [standardised_inputs folder](standardised_inputs).


### Init michelanglo data

Adapted from `prepare.ipynb`.

```python3
from michelanglo_protein import ProteinCore, global_settings, Structure, ProteinAnalyser, Mutation
from michelanglo_protein.generate import ProteinGatherer
global_settings.startup('/users/brc/matteo/michelanglo/protein-data')

```
### Download safeguard function

The safeguard decorator-class from [this gist](https://gist.github.com/matteoferla/24d9a319d05773ae219dd678a3aa11be)
was used initially.

```
def retrieve_from_gist(gist_sha: str, gist_filename: str, wanted_variable: str):
    import requests, importlib
    codeblock = requests.get(f'https://api.github.com/gists/{gist_sha}').json()['files'][gist_filename]['content']
    faux_global = {**globals(),  'Any': importlib.import_module('typing').Any,   'defaultdict': importlib.import_module('collections').defaultdict   }
    exec(codeblock, faux_global)
    return faux_global[wanted_variable]

Safeguard = retrieve_from_gist(gist_sha = '24d9a319d05773ae219dd678a3aa11be', gist_filename = 'safeguard.py', wanted_variable='Safeguard')
```

