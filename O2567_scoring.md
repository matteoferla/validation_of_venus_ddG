## Read dataset
The `ci0c00591_si_002` Excel spreadsheet from the SI of the O2567 dataset is an Excel spreadsheet that
I could not convert, so saved it as a csv in Excel, which screwed up some fields —see next cells.
```python
import pandas as pd

# This has a weird line.
#reference = pd.read_excel('ci0c00591_si_002.xlsx', engine='openpyxl')
reference = pd.read_csv('../O2567.csv')
# fix the excel BS.
reference.at[(reference['PDB code'] == '1.00E+65'), 'PDB code'] = '1E65'
reference['PDB code'] = reference['PDB code'].astype(str)
unique = reference['PDB code'].unique()
print(f'There are {len(unique)} unique structures')
unique
```

    There are 106 unique structures
    array(['1A2J', '1A2P', '1A6M', '1AEP', '1AIE', '1AY7', '1B26', '1C2R',
           '1C9O', '1CQM', '1CYC', '1CYO', '1DIV', '1E65', '1EGP', '1EY0',
           '1FEP', '1FNF', '1FTG', '1FXA', '1G4I', '1GVP', '1H3T', '1H7M',
           '1HG7', '1HUW', '1IDS', '1IFC', '1IG5', '1IHB', '1IRK', '1ISP',
           '1JBE', '1LMB', '1LNI', '1MJC', '1ARR', '1OIA', '1OMP', '1QGV',
           '1QHJ', '1QLP', '1RCF', '1RG8', '1RHG', '1RN1', '1SCE', '1U5P',
           '1V7Y', '1WAA', '1XPB', '1YCC', '1YTC', '2BLS', '2CBE', '2CI2',
           '2CRK', '2H61', '2HFT', '2IMM', '2LZM', '2PPN', '2Q98', '2RJX',
           '2TRX', '2TS1', '3B0K', '3BLM', '3GP6', '3I40', '3MVT', '3PF4',
           '3PGK', '3SIL', '3SNF', '3UA6', '3UTO', '3VSY', '4BLM', '4HIL',
           '4IPY', '4NX6', '4P2P', '4POC', '5C0Z', '5CE4', '5CGQ', '5LVE',
           '5T4L', '9ILB', '1HL4', '1PGA', '1OTR', '1FNA', '1BI7', '2RN2',
           '4AXT', '1UZC', '1UBQ', '2ZCC', '1PTA', '451C', '1BPI', '1HNG',
           '1TEN', '1APS'], dtype=object)

## Init Mike
```python
from michelanglo_protein import ProteinCore, global_settings, Structure, ProteinAnalyser, Mutation
from michelanglo_protein.generate import ProteinGatherer
global_settings.startup('/users/brc/matteo/michelanglo/protein-data')
```

## Polish dataset
The dataset is 
```python
import os, json

with open(os.path.join(global_settings.dictionary_folder, 'uniprot2pdb.json')) as fh:
    uni2pdb = json.load(fh)
# format: uni2pdb['P24565'] = ['1PNB_A']

pdb2uni = {}   # format: pdb2uni['1PNB_A'] = 'P24565'
for uni, pdbs in uni2pdb.items():
    for pdb in pdbs:
        pdb2uni[pdb] = uni

reference['Chain'] = reference.Chain.astype(str).map({'-': 'A', 'A': 'A', '_': 'A', 'nan': 'A',
                                                    'B': 'B', '3': '3', 'I': 'I'})
reference['PDB_chain'] = reference['PDB code'] + '_' + reference.Chain
reference['uniprot'] = reference['PDB_chain'].map(pdb2uni).fillna('P404')
reference['mutation'] = reference['Wild'] + reference['Residue number'].astype(str) + reference['Mutated']
reference.to_csv('O2567.corrected.csv')
reference
```
|      | PDB code   | Wild   |   Residue number | Mutated   | Chain   |   Experimental ddG |    RSA | Temperature   |     pH | Method    | PDB_chain   | uniprot   | mutation   |
|-----:|:-----------|:-------|-----------------:|:----------|:--------|-------------------:|-------:|:--------------|-------:|:----------|:------------|:----------|:-----------|
|    0 | 1A2J       | C      |               33 | S         | A       |         -1.1       |   0.3  | 25            |   7.5  | Urea      | 1A2J_A      | P0AEG4    | C33S       |
|    1 | 1A2J       | C      |               30 | S         | A       |         -1.8       |   2.5  | 25            |   7.5  | Urea      | 1A2J_A      | P0AEG4    | C30S       |
|    2 | 1A2P       | H      |               18 | Q         | A       |         -1.44      |  39.1  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | H18Q       |
|    3 | 1A2P       | T      |                6 | P         | A       |         -3.3       |  51.4  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | T6P        |
|    4 | 1A2P       | D      |                8 | G         | A       |         -1.3       |  62.4  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | D8G        |
|    5 | 1A2P       | D      |                8 | S         | A       |         -0.7       |  62.4  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | D8S        |
|    6 | 1A2P       | D      |               12 | G         | A       |         -0.2       |  49.6  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | D12G       |
|    7 | 1A2P       | D      |               12 | S         | A       |         -0.9       |  49.6  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | D12S       |
|    8 | 1A2P       | T      |               16 | G         | A       |         -1.3       |  56.5  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | T16G       |
|    9 | 1A2P       | T      |               16 | A         | A       |         -0.5       |  56.5  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | T16A       |
|   10 | 1A2P       | Y      |               17 | G         | A       |         -4         |  60.8  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | Y17G       |
|   11 | 1A2P       | Y      |               17 | S         | A       |         -2         |  60.8  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | Y17S       |
|   12 | 1A2P       | H      |               18 | S         | A       |         -2.1       |  39.1  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | H18S       |
|   13 | 1A2P       | H      |               18 | N         | A       |         -1.8       |  39.1  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | H18N       |
|   14 | 1A2P       | H      |               18 | D         | A       |         -2.2       |  39.1  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | H18D       |
|   15 | 1A2P       | H      |               18 | R         | A       |         -1         |  39.1  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | H18R       |
|   16 | 1A2P       | T      |               26 | D         | A       |         -0.08      |  24    | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | T26D       |
|   17 | 1A2P       | S      |               28 | G         | A       |         -0.5       |  60.8  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | S28G       |
|   18 | 1A2P       | S      |               28 | A         | A       |          0.6       |  60.8  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | S28A       |
|   19 | 1A2P       | E      |               29 | S         | A       |         -1.2       |  49.5  | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | E29S       |
|   20 | 1A2P       | Q      |               31 | G         | A       |         -1.1       |  53    | 25            |   6.3  | Urea      | 1A2P_A      | P00648    | Q31G       |

## Code for analysis

```python
from michelanglo_protein import ProteinCore, global_settings, Structure, ProteinAnalyser, Mutation

import time, pymol2, json
from sqlitedict import SqliteDict
from concurrent.futures import ThreadPoolExecutor
from functools import partial

results = SqliteDict('O2567-scores.db', encode=json.dumps, decode=json.loads, autocommit=True)

def analyse_table(reference, settings, max_workers):
    w_settings = partial(analyse_row, settings=settings)
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        pool.map(w_settings, [row for i, row in reference.iterrows()])

def analyse_row(row, settings, silent=False, debug=False):
    if debug:
        Exceptionality = ()
    else:
        Exceptionality = Exception
    try:
        protein = get_protein(row, silent=silent)
        for setting_name, setting in settings.items():
            acc = f'{protein.pdbs[0].id}_{setting_name}'
            analysed = analyse(pdbblock=protein.pdbblock,
                              mutation=protein.mutation,
                              acc=acc,
                              **setting)
            collected = {**row.to_dict(),
                         'silent': silent,
                          'condition': setting_name,
                         'apriori': protein.mutation.apriori_effect,
                         **analysed}
            # weird numpy issue (Object of type int64 is not JSON serializable)
            # collected['Protherm ID'] = int(collected['Protherm ID'])
            for k in ['Experimental ddG',
                      'RSA', 'Temperature', 'pH']:
                collected[k] = float(collected[k])
            for k in ['Residue number']:  # 'Protherm ID'
                collected[k] = int(collected[k])
            if debug:
                for k in collected:
                    print(k, type(collected[k]), collected[k])
            results = SqliteDict('O2567-scores.db', encode=json.dumps, decode=json.loads, autocommit=True)
            results[acc] = collected
    except Exceptionality as error:
        print(error.__class__.__name__, str(error))

def check_sequence(protein):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.read_pdbstr(protein.pdbblock, 'xxx')
        atom = pymol.cmd.get_model(f'resi {protein.mutation.residue_index} and name CA').atom
        if atom:
            mapping = {name3.upper(): name1 for name1, name3, fullname in protein.mutation.names}
            assert mapping[atom[0].resn] == protein.mutation.from_residue, 'mismatch'
        else:
            raise ValueError(f'resi {protein.mutation.residue_index} and name CA does not exist')


def get_protein(row, silent=False,
                code_col='PDB code',
                wild_type_col='Wild',
                res_num_col='Residue number',
                mutated_col='Mutated',
                chain_col='Chain'):
    protein = ProteinAnalyser(uniprot=row['uniprot'])
    if silent:
        mut = row[wild_type_col] + str(row[res_num_col]) + row[wild_type_col]
    else:
        mut = row[wild_type_col] + str(row[res_num_col]) + row[mutated_col]
    identifier = row[code_col] + '_' + mut
    protein.mutation = Mutation(mut)
    model = Structure(id=identifier,
                      description='',
                      x=0, y=99999,
                      code=row[code_col],
                      type='rcsb',
                      chain=row[chain_col],
                      offset=0)
    protein.pdbs.append(model)
    protein.analyse_structure(no_conservation=True)
    #     if not protein.check_mutation():
    #         raise ValueError(f"Discrepancy for {row['PDB ID']}")
    check_sequence(protein)
    return protein

def analyse(pdbblock,
            mutation,
            acc,
            radius=12,
            cycles=1,
            scorefxn_name='ref2015',
            use_pymol_for_neighbours=False,
            neighbour_only_score=True,
            outer_constrained=True,
            remove_ligands=True,
            single_chain=True,
            prevent_acceptance_of_incrementor=False,
            debug=False
            ):
    protein = ProteinAnalyser(scorefxn_name=scorefxn_name, cycles=cycles, radius=radius, use_pymol_for_neighbours=use_pymol_for_neighbours)
    protein.mutation = mutation
    protein.use_pymol_for_neighbours = use_pymol_for_neighbours
    tick = time.time()
    protein.energetics = {'native': pdbblock}
    data = protein.analyse_FF(spit_process=not debug,
                              neighbour_only_score=neighbour_only_score,
                              outer_constrained=outer_constrained,
                              remove_ligands=remove_ligands,
                              prevent_acceptance_of_incrementor=prevent_acceptance_of_incrementor,
                              single_chain=single_chain)
    tock = time.time()
    if 'native' in data:
        with open(f'O2567-structures/{acc}.native.pdb', 'w') as fh:
            fh.write(data['native'])
    if 'mutant' in data:
        with open(f'O2567-structures/{acc}.mutant.pdb', 'w') as fh:
            fh.write(data['mutant'])
    reply = dict(time=tock - tick)
    if 'ddG' not in data:
        if len(data) > 1:
            print('keys present', data.keys())
        elif 'error' in data:
            reply['error'] = data['error']
            print(data['error'])
    for k in ['ddG', 'scores', 'rmsd', 'dsol', 'score_fxn', 'neighbours',
              'cycles', 'radius', 'neighbouring_ligand', 'n_constraints']:
        if k in data:
            reply[k] = data[k]
    return reply
```

## Test analysis setup

```python
# import logging
# log = logging.getLogger()
# handler = logging.handlers.
# log.addHandler(handler)

default = {'default':
                      {'radius': 12,
                       'cycles': 1,
                       'scorefxn_name': 'talaris2014',
                       'neighbour_only_score': True,
                       'outer_constrained': False,
                       'remove_ligands': True,
                       'single_chain': True,
                       'prevent_acceptance_of_incrementor': False,
                       'bogus':True}
            }
analyse_row(reference.iloc[0], default, silent=False, debug=True)
```
    definitionless structure: 1A2J
    dict_keys(['error'])
    PDB code <class 'str'> 1A2J
    Wild <class 'str'> C
    Residue number <class 'int'> 33
    Mutated <class 'str'> S
    Chain <class 'str'> A
    Experimental ddG <class 'float'> -1.1
    RSA <class 'float'> 0.3
    Temperature <class 'float'> 25.0
    pH <class 'float'> 7.5
    Method <class 'str'> Urea
    PDB_chain <class 'str'> 1A2J_A
    uniprot <class 'str'> P0AEG4
    mutation <class 'str'> C33S
    silent <class 'bool'> False
    condition <class 'str'> default
    apriori <class 'str'> The mutation changes one amino acid to another that is equally sized, more flexible.
    time <class 'float'> 15.677705526351929
    error <class 'str'> AttributeError:'Mutator' object has no attribute 'bogus'
    
## Run analysis

```python
import pandas as pd

settings = pd.read_csv('settings.csv').set_index('condition').transpose().to_dict()
analyse_table(reference, settings, max_workers=30)
```