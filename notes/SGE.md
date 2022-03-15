## SGE

common functions:

```python
# qsub in chunks
from mako.template import Template
import re, os
import pandas as pd
from typing import (Optional)

with open('sge.mako', 'r') as fh:
    template = fh.read()
    # print(', '.join(re.findall(r'\$\{(.*?)}', template)))

def split_up(dataset_name, n_chunks):
    reference = pd.read_csv(f'{dataset_name}.standarised.csv', index_col=0)
    import numpy as np
    for i, split in enumerate(np.array_split(reference, n_chunks)):
        split.to_csv(f'{dataset_name}.standarised_{i}.csv')
        
def qsub_chunks(dataset_name,
                n_chunks,
                dataset_filename: Optional[str]=None,
                protein_filename: Optional[str]=None,
                settings_filename='settings.json'):
    if not protein_filename:
        protein_filename=f'{dataset_name}_protein.p'
    for i in range(n_chunks):
        if not dataset_filename:
            subdataset_filename=f'{dataset_name}.standarised_{i}.csv'
        else:
            subdataset_filename=dataset_filename.format(i)
        with open(f'{dataset_name}_{i}.sh', 'w') as fh:
            makoed = Template(template).render(jobname=f'{dataset_name}_{i}',
                                               cores=5,
                                               queue='short.qc',
                                               py_filename='analyse.py',
                                               py_args=f'{dataset_name} {subdataset_filename} {protein_filename} {settings_filename} 5'
                                               )
            fh.write(makoed)
        os.system(f'qsub {dataset_name}_{i}.sh')
```

These got run as:

```python
%%script false --no-raise-error
dataset_name = 'O2567-AF2'
n_chunks = 10
split_up(dataset_name, n_chunks)
qsub_chunks(dataset_name, n_chunks, settings_filename='minimal_settings.json')
```

## Magic

These may appear in places.

The cell magic called `rescomp` is to run the contents of the cell as a submitted job.

```python
%%rescomp jobname=test queue=test.qc cores=1

with open('greeting.txt', 'w') as fh:
    fh.write('hello world')
import sys
import time
time.sleep(10)
sys.exit(0)
```

```python
from IPython.core.magic import (register_cell_magic)
import re
from mako.template import Template
@register_cell_magic
def rescomp(line, cell) -> int:
    """
    Returns job id
    https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/cluster-usage
    
    %%rescomp jobname=test queue=test.qc cores=1
    """
    jobname = re.search(r'jobname=(\S*)', line).group(1).strip()
    queue = re.search(r'queue=(\S*)', line).group(1).strip()
    cores = re.search(r'cores=(\d+)', line).group(1).strip()
    py_filename = f'temp_{jobname}.py'
    with open(py_filename, 'w') as fh:
        fh.write(cell)
    with open('sge.mako', 'r') as fh:
        template = fh.read()
    makoed = Template(template).render(jobname=jobname,
                                        cores=cores,
                                        queue=queue,
                                       py_filename=py_filename,
                                       py_args=''
                                       )
    with open(f'temp_{jobname}.sh', 'w') as fh:
            fh.write(makoed)
    import subprocess
    p = subprocess.Popen(["qsub", f'temp_{jobname}.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = p.stdout.read().decode()
    err = p.stderr.read().decode()
    if err:
        raise RuntimeError(err)
    job_id = int(re.search('job (\d+)', out).group(1))
    if 'JobStatusChecker' in globals():
        globals()['JobStatusChecker'].jobs.append(job_id)
    return out

    # print(', '.join(re.findall(r'\$\{(.*?)}', template)))

del rescomp
```

Other code that _may_ be called is the following (do not copy for inspiration as `drmaa` module is annoying as it uses v1).

```python

def job_ids():
    import subprocess
    p = subprocess.Popen(["qstat"], stdout=subprocess.PIPE)
    out = p.stdout.read()
    return [int(line.strip().split()[0]) for line in out.decode().strip().split('\n')[2:] ]

import os
os.environ['DRMAA_LIBRARY_PATH'] = '/mgmt/uge/8.6.8/lib/lx-amd64/libdrmaa.so.1.0'
import drmaa

class JobStatusChecker:
    session = None
    jobs = []
    
    def __init__(self):
        if not self.session:
            self.session = drmaa.Session()
            self.session.initialize()
            self.__class__.session = self.session
            
    def __del__(self):
        self.session.exit()
        
    def __getitem__(self, job_id):
        try:
            return self.session.jobStatus(str(job_id))
        except drmaa.InvalidJobException:
            return 'disappeared'
    
    def last(self):
        assert len(self.jobs), 'No known jobs'
        job_id = self.jobs[-1]
        return job_id, self[job_id]
```



## No interent
There is no internet in the cores so this was used to test the calls, which disables the use of sockets
Note that a subprocess is a socket, ae!

```python
import socket
class guard:
    """
    This blocks the internet, but also blocks subprocesses unfortunately.
    """
    def __init__(*args, **kwargs):
        raise Exception("I told you not to use the Internet!")
socket.socket = guard
```