## Dataset generation

## System
The dataset of scores was done on a cluster node (x32 Intel Xeon (Haswell, x86_64) CPU E5-2640 v3 @ 2.60GHz)
with PyRosetta4 release 2020.49 (Dec2020) for python37 on Linux.

The data was stored with SqliteDict in case of interuption, but was not actually needed
(SqliteDict is basically a dict with permanence).

## Prep
Start the michelanglo_protein files —not needed...
```python
from michelanglo_protein import ProteinCore, global_settings, Structure, ProteinAnalyser, Mutation
from michelanglo_protein.generate import ProteinGatherer
global_settings.startup('/users/brc/matteo/michelanglo/protein-data')
```

Correct logging

```python
import logging, sys

log = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
handler.set_name('stdout')
handler.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
log.addHandler(handler)
```

```python
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7579412/#TS1
import pandas as pd
import os, json

reference_data = pd.read_csv('Data_Sheet_2.CSV')

with open(os.path.join(global_settings.dictionary_folder, 'uniprot2pdb.json')) as fh:
    uni2pdb = json.load(fh)
# format: uni2pdb['P24565'] = ['1PNB_A']

pdb2uni = {}   # format: pdb2uni['1PNB_A'] = 'P24565'
for uni, pdbs in uni2pdb.items():
    for pdb in pdbs:
        pdb2uni[pdb] = uni

reference_data['PDB_chain'] = reference_data['PDB ID'] + '_' + reference_data.Chain
reference_data['mutation'] = reference_data['Wild Type'] + reference_data['Residue Number'] + reference_data['Mutation']
reference_data['uniprot'] = reference_data['PDB_chain'].map(pdb2uni).fillna('P404')
```

Define methods:

```python
from michelanglo_protein import Structure
import time, pymol2
from sqlitedict import SqliteDict
results = SqliteDict('neoscores.db', encode=json.dumps, decode=json.loads, autocommit=True)

def check_sequence(protein):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.read_pdbstr(protein.pdbblock, 'xxx')
        atom = pymol.cmd.get_model(f'resi {protein.mutation.residue_index} and name CA').atom
        if atom:
            mapping = {name3.upper(): name1 for name1, name3, fullname in protein.mutation.names}
            assert mapping[atom[0].resn] == protein.mutation.from_residue, 'mismatch'
        else:
            raise ValueError(f'resi {protein.mutation.residue_index} and name CA does not exist')
        

def get_protein(row, silent=False):
    protein = ProteinAnalyser(uniprot=row['uniprot'])
    if silent:
        mut = row['Wild Type'] + str(row['Residue Number']) + row['Wild Type']
    else:
        mut = row['mutation']
    protein.mutation = Mutation(mut)
    
    model = Structure(id=row['Protherm ID'],
                      description='',
                      x=0, y=99999,
                      code=row['PDB ID'],
                      type='rcsb',
                      chain=row['Chain'],
                      offset = 0)

    protein.pdbs.append(model)
    protein.analyse_structure(no_conservation=True)
#     if not protein.check_mutation():
#         raise ValueError(f"Discrepancy for {row['PDB ID']}")
    check_sequence(protein)
    return protein

def analyse(protein,
            acc,
            radius, cycles, scorefxn_name,
            use_pymol_for_neighbours = False,
            neighbour_only_score = True,
            outer_constrained = True,
            debug=False
           ):
    protein.use_pymol_for_neighbours = use_pymol_for_neighbours
    protein.radius = radius
    protein.scorefxn_name = scorefxn_name
    protein.cycles = cycles
    tick = time.time()
    protein.energetics = None
    data = protein.analyse_FF(spit_process=not debug,
                              neighbour_only_score=neighbour_only_score,
                              outer_constrained=outer_constrained)
    tock = time.time()
    with open(f'structures/{acc}.native.pdb', 'w') as fh:
        fh.write(data['native'])
    with open(f'structures/{acc}.mutant.pdb', 'w') as fh:
        fh.write(data['mutant'])
    reply = dict(time=tock - tick)
    for k in ['ddG', 'scores', 'rmsd', 'dsol', 'score_fxn', 'neighbours',
              'cycles', 'radius', 'neighbouring_ligand', 'n_constraints']:
        reply[k] = data[k]
    if 'ddG' not in data:
        print(data.keys())
        reply['error'] = data['error']
    else:
        for k in ('ddG', 'neighbours'):
            reply[k] = data[k]
    return reply

def analyse_row(row,
                neighbour_only_score=True,
                outer_constrained=True,
                silent=False, 
                cycleses=(3,), # (1, 2, 3)
                radii=(12,), # (8, 10, 12)
                scorefxn_names=('ref2015', ), # ('ref2015', 'beta_nov16', 'ref2015_cart', 'beta_nov16_cart')
                debug=False):
    if debug:
        Exceptionality = ()
    else:
        Exceptionality = Exception
    try:
        protein = get_protein(row, silent=silent)
        for cycles in cycleses:
            for radius in (12,):
                for scorefxn_name in scorefxn_names:
                    #acc = f"{row['Protherm ID']}_{scorefxn_name}_{radius}_{cycles}"
                    silent_mark = 's' if silent else 'S'
                    neigh_mark = 'n' if neighbour_only_score else 'N'
                    con_mark = 'c' if outer_constrained else 'C'
                    acc = f"{row['Protherm ID']}_{scorefxn_name}_{radius}_{cycles}_{neigh_mark}{silent_mark}{con_mark}"
                    results = SqliteDict('neoscores.db', encode=json.dumps, decode=json.loads, autocommit=True)
                    if acc in results:
                        if debug:
                            print(f'{acc}_present')
                        continue
                    collected = {**row.to_dict(),
                                    'silent': silent,
                                    'constraints': outer_constrained,
                                    'neighbor_score': neighbour_only_score,
                                    **analyse(protein=protein,
                                              acc=acc,
                                              radius = radius, cycles=cycles,
                                              scorefxn_name = scorefxn_name,
                                              outer_constrained=outer_constrained,
                                              neighbour_only_score=neighbour_only_score
                                             )}
                    # weird numpy bug
                    collected['Protherm ID'] = int(collected['Protherm ID'])
                    collected['SASA'] = float(collected['SASA'])
                    collected['Experimental DDG'] = float(collected['Experimental DDG'])
                    
                    if debug:
                        for k in collected:
                            print(k, type(collected[k]), collected[k])
                    results[acc] = collected
    except Exceptionality as error:
        print(error.__class__.__name__, str(error))
```

## Calculate
Run test as appropriate:

```python
partial_analyses = partial(analyse_row, silent=False,
                             neighbour_only_score=True, outer_constrained=False,
                            cycleses=(3,),
                             radii=(12,),
                            )
                            
with ThreadPoolExecutor(max_workers = 10) as pool:
      pool.map(partial_analyses, [row for i, row in reference_data.iterrows()])
```
