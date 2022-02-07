# Scoring

For previous iterations of scoring see [notes](notes).
The major difference is that in this version the reference tables are made consistent
as opposed to the analyses changing the names for them...

## Rectification

### Init michelanglo data

```python3
from michelanglo_protein import ProteinCore, global_settings, Structure, ProteinAnalyser, Mutation
from michelanglo_protein.generate import ProteinGatherer
global_settings.startup('/users/brc/matteo/michelanglo/protein-data')

```
### Download safeguard function
```

def retrieve_from_gist(gist_sha = '24d9a319d05773ae219dd678a3aa11be', gist_filename = 'safeguard.py', wanted_variable='Safeguard'):
    import requests, importlib
    codeblock = requests.get(f'https://api.github.com/gists/{gist_sha}').json()['files'][gist_filename]['content']
    faux_global = {**globals(),  'Any': importlib.import_module('typing').Any,   'defaultdict': importlib.import_module('collections').defaultdict   }
    exec(codeblock, faux_global)
    return faux_global[wanted_variable]

# https://gist.github.com/matteoferla/24d9a319d05773ae219dd678a3aa11be
Safeguard = retrieve_from_gist()
```