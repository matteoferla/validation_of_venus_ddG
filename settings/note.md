For :octocat:

```python
import os, json


def get_data(filename):
    with open('settings.json') as fh:
        return json.load(fh)

import pandas as pd

print(pd.DataFrame({k: v for filename in os.listdir() if 'settings.json' in filename for k, v in get_data(filename).items()}).transpose().to_markdown())
```