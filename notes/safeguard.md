## Safeguard

Some files _may_ have a decorator class called `Safeguard`.
This is from [this gist](https://gist.github.com/matteoferla/24d9a319d05773ae219dd678a3aa11be)

```
def retrieve_from_gist(gist_sha: str, gist_filename: str, wanted_variable: str):
    import requests, importlib
    codeblock = requests.get(f'https://api.github.com/gists/{gist_sha}').json()['files'][gist_filename]['content']
    faux_global = {**globals(),  'Any': importlib.import_module('typing').Any,   'defaultdict': importlib.import_module('collections').defaultdict   }
    exec(codeblock, faux_global)
    return faux_global[wanted_variable]

Safeguard = retrieve_from_gist(gist_sha = '24d9a319d05773ae219dd678a3aa11be', gist_filename = 'safeguard.py', wanted_variable='Safeguard')
```

