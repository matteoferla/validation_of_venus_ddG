### O2567

> TODO move code over for the gen of O2567.corrected.csv

```python
import os, json
import pandas as pd
from sqlitedict import SqliteDict
from michelanglo_protein import Structure
dataset_name = 'O2567-redux'

reference = pd.read_csv('O2567.corrected.csv', index_col=0)
reference['dataset'] = 'O2567'
reference['Experimental ddG'] = - reference['Experimental ddG']
# the columns labels differ between every dataset
columns = reference.columns.values.tolist()
reference.columns = reference.columns.map({**dict(zip(columns, columns)),
                                           'Residue number':   'resi',
                                           'PDB code':         'PDB_code',
                                           'Chain':            'chain',
                                           'Experimental ddG': 'experimental_ddG',
                                           'Wild':             'from_resn',
                                           'Mutated':          'to_resn',})
reference['protein_id'] = reference['PDB_code'] + '_' + reference['chain']
reference.to_csv('O2567.standarised.csv')

import pickle
models = {}
for protein_id in reference.protein_id.unique():
    pdb_code, chain = protein_id.split('_')
    models[protein_id] = Structure(id=chain,
                                  description='',
                                  x=0, y=99999,
                                  code=pdb_code,
                                  type='rcsb',
                                  chain=chain,
                                  offset=0)
    models[protein_id].get_coordinates()
    
with open('O2567_protein.p', 'wb') as fh:
    pickle.dump(obj=models, file=fh)
```

### Star

```python
import os, json
import pandas as pd
from sqlitedict import SqliteDict
from michelanglo_protein import Structure
dataset_name = 'star'

reference = pd.read_csv('star.csv', index_col=0)
reference['dataset'] = 'Protherm*'
reference['Experimental DDG'] = reference['Experimental DDG']
# the columns labels differ between every dataset
columns = reference.columns.values.tolist()
key_fixes = {'Residue Number':   'resi',
             'PDB ID':         'PDB_code',
             'Chain':            'chain',
             'Experimental DDG': 'experimental_ddG',
             'Wild Type':             'from_resn',
             'Mutation':          'to_resn',}
assert all([k in columns for k in key_fixes]), ([k for k in key_fixes if k not in columns], columns)
reference.columns = reference.columns.map({**dict(zip(columns, columns)),
                                           **key_fixes,
                                          })
reference['protein_id'] = reference['PDB_code'] + '_' + reference['chain']
reference.to_csv('star.standarised.csv')

import pickle
models = {}
for protein_id in reference.protein_id.unique():
    pdb_code, chain = protein_id.split('_')
    models[protein_id] = Structure(id=chain,
                                  description='',
                                  x=0, y=99999,
                                  code=pdb_code,
                                  type='rcsb',
                                  chain=chain,
                                  offset=0)
    models[protein_id].get_coordinates()
    
with open('star_protein.p', 'wb') as fh:
    pickle.dump(obj=models, file=fh)
```

### AlphaFold2

> ToDo Copy code over for Uniprot id allocation

```python
from michelanglo_protein import ProteinAnalyser
reference = pd.read_csv('O2567.standarised.csv', index_col=0)

import pickle
models = {}
missing = []
for protein_id in reference.protein_id.unique():
    uniprot = reference.loc[reference.protein_id == protein_id].uniprot.values[0]
    protein = ProteinAnalyser(uniprot=uniprot)
    protein.add_alphafold2()
    if not len(protein.alphafold2):
        missing.append(protein_id) # impossible with v2
        continue
    structure = protein.alphafold2[0]
    structure.id = protein_id
    try:
        structure.get_coordinates()
        models[protein_id] = structure
    except ConnectionError as error:
        missing.append(protein_id)
    
with open('O2567_AF2_protein.p', 'wb') as fh:
    pickle.dump(obj=models, file=fh)
```
alt.
```python
fixed_reference = pd.read_csv('../O2567.uniprot.csv', index_col=0)
uniprot_resis = dict(zip((fixed_reference['PDB_chain']+ '_'+fixed_reference['mutation']).values, fixed_reference['uniprot_residue_number'].values.astype(int)))

reference = pd.read_csv('O2567.standarised.csv', index_col=0)

def convert_resi(row: pd.DataFrame) -> int:
    p_id = row['PDB_chain']+ '_'+row['mutation']
    if p_id not in uniprot_resis:
        return 0
    elif uniprot_resis[p_id] <1:
        return 0
    else:
        return uniprot_resis[p_id]
reference['pdb_resi'] = reference.resi
reference['pdb_mutation'] = reference.mutation
reference['resi'] = reference.apply(convert_resi, axis=1)
reference['mutation'] = reference.from_resn+reference.resi.astype(str)+reference.to_resn
reference['dataset'] = 'O2567-AF2'
reference = reference.loc[(reference.resi > 0) & (~reference.protein_id.isin(missing))]
reference.to_csv('O2567-AF2.standarised.csv')
```

## Silent

```python
reference = pd.read_csv('O2567.standarised.csv', index_col=0)
reference['to_resn'] = reference['from_resn']
reference['experimental_ddG'] = 0
reference['note'] = 'Silent mutation dataset'
reference.to_csv('O2567-silent.standarised.csv')
```

## S1342

```python
import os
from michelanglo_protein import global_settings

global_settings.startup(os.environ['MICHELANGLO_PROTEIN_DATA'])


def get_pdb2uni():
    """
    reference_data['PDB-CODE + _ + CHAIN'].map(pdb2uni).fillna('P404')
    """
    import os, json
    with open(os.path.join(global_settings.dictionary_folder, 'uniprot2pdb.json')) as fh:
        uni2pdb = json.load(fh)
    # format: uni2pdb['P24565'] = ['1PNB_A']
    pdb2uni = {}  # format: pdb2uni['1PNB_A'] = 'P24565'
    for uni, pdbs in uni2pdb.items():
        for pdb in pdbs:
            pdb2uni[pdb] = uni
    return pdb2uni


reference = pd.read_csv('S1342_130-Proteins.csv')
reference['dataset'] = 'S1342'
reference['Experimental'] = - reference['Experimental']
# the columns labels differ between every dataset
columns = reference.columns.values.tolist()
reference.columns = reference.columns.map({**dict(zip(columns, columns)),
                                           'Position':     'resi',
                                           'PDB':          'PDB_code',
                                           'Chain':        'chain',
                                           'Experimental': 'experimental_ddG',
                                           'Wild.1':       'from_resn',
                                           'Mutant.1':     'to_resn', })
reference['protein_id'] = reference['PDB_code'] + '_' + reference['chain']
reference['mutation'] = reference['from_resn'] + reference['resi'].astype(str) + reference['to_resn']
pdb2uni = get_pdb2uni()
reference['uniprot'] = reference['protein_id'].map(pdb2uni).fillna('P404')
reference = reference.sort_values('PDB_code')
reference.to_csv('S1342.standarised.csv')
reference
```
and
```python
import pickle
models = {}
for protein_id in reference.protein_id.unique():
    pdb_code, chain = protein_id.split('_')
    models[protein_id] = Structure(id=chain,
                                  description='',
                                  x=0, y=99999,
                                  code=pdb_code,
                                  type='rcsb',
                                  chain=chain,
                                  offset=0)
    models[protein_id].get_coordinates()
    
with open('S1342_protein.p', 'wb') as fh:
    pickle.dump(obj=models, file=fh)
```

## SwissModel

```python
from michelanglo_protein import ProteinAnalyser
reference = pd.read_csv('O2567-AF2.standarised.csv', index_col=0)
import json
import pickle
models = {}
proteins = {}
missing = []
bad_uniprots = []
# protein_id is a totally disconnected pdb-code _ chain
for protein_id in reference.protein_id.unique():
    subset = reference.loc[reference.protein_id == protein_id]
    uniprot = subset.uniprot.values[0]
    if uniprot in bad_uniprots:
        continue
    elif uniprot in proteins:
        protein = proteins[uniprot]
    else:
        protein = ProteinAnalyser(uniprot=uniprot)
        try:
            protein.retrieve_structures_from_swissmodel()
            protein.pdbs = []
            protein.alphafold2 = []
            proteins[uniprot] = protein
        except json.JSONDecodeError as error:
            bad_uniprots.append(uniprot)
            break
    try:
        for model in protein.swissmodel:
            for uniprot_residue_number in subset.resi.unique():
                # just in case there are no matches for one.
                if model.x < uniprot_residue_number < model.y:
                    model.get_coordinates()
                    models[protein_id] = structure
                    raise ZeroDivisionError
    except ConnectionError as error:
        missing.append(protein_id)
    except ZeroDivisionError as error:  # double break.
        pass
    
with open('O2567-swissmodel_protein.p', 'wb') as fh:
    pickle.dump(obj=models, file=fh)
```