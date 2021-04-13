## Step 1. Scoring

> This is a modified Markdown-exported version of the notebook (see [data](data) for actual notebook)

```python
from michelanglo_protein import ProteinCore, global_settings, Structure, ProteinAnalyser, Mutation
from michelanglo_protein.generate import ProteinGatherer
global_settings.startup('/users/brc/matteo/michelanglo/protein-data')
```

    Folder path set to /users/brc/matteo/michelanglo/protein-data





    <michelanglo_protein.settings_handler.GlobalSettings at 0x7f45f821d410>




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
reference_data
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>PDB ID</th>
      <th>Protherm ID</th>
      <th>Residue Number</th>
      <th>Chain</th>
      <th>Wild Type</th>
      <th>Mutation</th>
      <th>SASA</th>
      <th>Experimental DDG</th>
      <th>Classifiers</th>
      <th>PDB_chain</th>
      <th>mutation</th>
      <th>uniprot</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.224900</td>
      <td>0.800000</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.500000</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.600000</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1TUP</td>
      <td>2256</td>
      <td>273</td>
      <td>A</td>
      <td>R</td>
      <td>H</td>
      <td>0.468525</td>
      <td>0.350000</td>
      <td>'positive to non-charged polar' 'buried'</td>
      <td>1TUP_A</td>
      <td>R273H</td>
      <td>P04637</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1LZ1</td>
      <td>1076</td>
      <td>110</td>
      <td>A</td>
      <td>V</td>
      <td>L</td>
      <td>0.429977</td>
      <td>-0.071702</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1LZ1_A</td>
      <td>V110L</td>
      <td>P61626</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>763</th>
      <td>1STN</td>
      <td>2014</td>
      <td>78</td>
      <td>A</td>
      <td>K</td>
      <td>C</td>
      <td>0.792768</td>
      <td>0.400000</td>
      <td>'large to small' 'involves cysteine' 'surface'</td>
      <td>1STN_A</td>
      <td>K78C</td>
      <td>P00644</td>
    </tr>
    <tr>
      <th>764</th>
      <td>1STN</td>
      <td>2046</td>
      <td>82</td>
      <td>A</td>
      <td>T</td>
      <td>C</td>
      <td>0.524742</td>
      <td>0.100000</td>
      <td>'involves cysteine' 'surface'</td>
      <td>1STN_A</td>
      <td>T82C</td>
      <td>P00644</td>
    </tr>
    <tr>
      <th>765</th>
      <td>1SUP</td>
      <td>2207</td>
      <td>206</td>
      <td>A</td>
      <td>Q</td>
      <td>C</td>
      <td>0.328425</td>
      <td>-1.250000</td>
      <td>'large to small' 'involves cysteine' 'buried'</td>
      <td>1SUP_A</td>
      <td>Q206C</td>
      <td>P00782</td>
    </tr>
    <tr>
      <th>766</th>
      <td>1TUP</td>
      <td>2253</td>
      <td>242</td>
      <td>A</td>
      <td>C</td>
      <td>S</td>
      <td>0.119379</td>
      <td>2.940000</td>
      <td>'involves cysteine' 'buried'</td>
      <td>1TUP_A</td>
      <td>C242S</td>
      <td>P04637</td>
    </tr>
    <tr>
      <th>767</th>
      <td>1VQB</td>
      <td>2338</td>
      <td>19</td>
      <td>A</td>
      <td>V</td>
      <td>C</td>
      <td>0.570471</td>
      <td>0.300000</td>
      <td>'involves cysteine' 'surface'</td>
      <td>1VQB_A</td>
      <td>V19C</td>
      <td>P69543</td>
    </tr>
  </tbody>
</table>
<p>768 rows × 12 columns</p>
</div>




```python
from michelanglo_protein import Structure
import time, pymol2

def check_sequence(protein):
    # The residue numbers are PDB residue numbers, but one can never be too careful.
    with pymol2.PyMOL() as pymol:
        pymol.cmd.read_pdbstr(protein.pdbblock, 'xxx')
        atom = pymol.cmd.get_model(f'resi {protein.mutation.residue_index} and name CA').atom
        if atom:
            mapping = {name3.upper(): name1 for name1, name3, fullname in protein.mutation.names}
            assert mapping[atom[0].resn] == protein.mutation.from_residue, 'mismatch'
        else:
            raise ValueError(f'resi {protein.mutation.residue_index} and name CA does not exist')
        

def get_protein(row):
    protein = ProteinAnalyser(uniprot=row['uniprot'])
    protein.mutation = Mutation(row['mutation'])

    model = Structure(id=row['Protherm ID'],
                      description='',
                      x=0, y=99999,
                      code=row['PDB ID'],
                      type='rcsb',
                      chain=row['Chain'],
                      offset = 0)

    protein.pdbs.append(model)
    protein.analyse_structure()
#     if not protein.check_mutation():
#         raise ValueError(f"Discrepancy for {row['PDB ID']}")
    check_sequence(protein)
    return protein

def analyse(protein, radius, scorefxn_name, use_pymol_for_neighbours = False):
    protein.use_pymol_for_neighbours = use_pymol_for_neighbours
    protein.radius = radius
    protein.scorefxn_name = scorefxn_name
    tick = time.time()
    protein.energetics = None
    data = protein.analyse_FF()
    tock = time.time()
    if 'ddG' not in data:
        print(data.keys())
    return dict(radius=protein.radius,
                 use_pymol_for_neighbours=protein.use_pymol_for_neighbours,
                 ddG=data['ddG'],
                 neighs = len(data['neighbours']),
                 scorefxn_name=protein.scorefxn_name,
                 time=tock - tick)

def analyse_row(row):
    try:
        protein = get_protein(row)
        for radius in (8, 10, 12):
            for scorefxn_name in ('ref2015', 'beta_nov16', 'ref2015_cart', 'beta_nov16_cart'):
                acc = f"{row['Protherm ID']}_{scorefxn_name}_{radius}"
                results = SqliteDict('scores.db', encode=json.dumps, decode=json.loads, autocommit=True)
                if acc in results:
                    continue
                results[acc] = {**row.to_dict(),
                                **analyse(protein, radius = radius, scorefxn_name = scorefxn_name)}
    except Exception as error:
        print(error.__class__.__name__, str(error))
```


```python
from sqlitedict import SqliteDict
results = SqliteDict('data/scores.db', encode=json.dumps, decode=json.loads, autocommit=True)
```


```python
# 'Protherm ID' is indeed unique.
reference_data['Protherm ID'].nunique() == len(reference_data['Protherm ID'])
```




    True




```python
# from multiprocessing import Pool
# # AssertionError daemonic processes are not allowed to have children

# with Pool(10) as p:
#     p.map(analyse_row, [row for i, row in reference_data.iloc[120:150].iterrows()])
```


```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(max_workers = 15) as pool:
      pool.map(analyse_row, [row for i, row in reference_data.iterrows()])
```



```python
pd.DataFrame(dict(results.items())).transpose().head(20)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>PDB ID</th>
      <th>Protherm ID</th>
      <th>Residue Number</th>
      <th>Chain</th>
      <th>Wild Type</th>
      <th>Mutation</th>
      <th>SASA</th>
      <th>Experimental DDG</th>
      <th>Classifiers</th>
      <th>PDB_chain</th>
      <th>mutation</th>
      <th>uniprot</th>
      <th>radius</th>
      <th>use_pymol_for_neighbours</th>
      <th>ddG</th>
      <th>neighs</th>
      <th>scorefxn_name</th>
      <th>time</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2199_ref2015_8</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.2249</td>
      <td>0.8</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
      <td>8</td>
      <td>False</td>
      <td>2.5621</td>
      <td>11</td>
      <td>ref2015</td>
      <td>9.55795</td>
    </tr>
    <tr>
      <th>2199_beta_nov16_8</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.2249</td>
      <td>0.8</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
      <td>8</td>
      <td>False</td>
      <td>-0.0440791</td>
      <td>11</td>
      <td>beta_nov16</td>
      <td>9.33401</td>
    </tr>
    <tr>
      <th>2199_ref2015_10</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.2249</td>
      <td>0.8</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
      <td>10</td>
      <td>False</td>
      <td>1.79129</td>
      <td>17</td>
      <td>ref2015</td>
      <td>11.8185</td>
    </tr>
    <tr>
      <th>2199_beta_nov16_10</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.2249</td>
      <td>0.8</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
      <td>10</td>
      <td>False</td>
      <td>1.34752</td>
      <td>17</td>
      <td>beta_nov16</td>
      <td>11.038</td>
    </tr>
    <tr>
      <th>2199_ref2015_12</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.2249</td>
      <td>0.8</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
      <td>12</td>
      <td>False</td>
      <td>1.89625</td>
      <td>28</td>
      <td>ref2015</td>
      <td>14.157</td>
    </tr>
    <tr>
      <th>2199_beta_nov16_12</th>
      <td>1STN</td>
      <td>2199</td>
      <td>122</td>
      <td>A</td>
      <td>E</td>
      <td>F</td>
      <td>0.2249</td>
      <td>0.8</td>
      <td>'polar to hydrophobic' 'negative to hydrophobi...</td>
      <td>1STN_A</td>
      <td>E122F</td>
      <td>P00644</td>
      <td>12</td>
      <td>False</td>
      <td>-0.696345</td>
      <td>28</td>
      <td>beta_nov16</td>
      <td>15.73</td>
    </tr>
    <tr>
      <th>1711_ref2015_8</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.5</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
      <td>8</td>
      <td>False</td>
      <td>0.765898</td>
      <td>12</td>
      <td>ref2015</td>
      <td>9.74753</td>
    </tr>
    <tr>
      <th>1711_beta_nov16_8</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.5</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
      <td>8</td>
      <td>False</td>
      <td>3.29626</td>
      <td>12</td>
      <td>beta_nov16</td>
      <td>9.09404</td>
    </tr>
    <tr>
      <th>1711_ref2015_10</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.5</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
      <td>10</td>
      <td>False</td>
      <td>0.514132</td>
      <td>23</td>
      <td>ref2015</td>
      <td>11.7462</td>
    </tr>
    <tr>
      <th>1711_beta_nov16_10</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.5</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
      <td>10</td>
      <td>False</td>
      <td>1.60971</td>
      <td>23</td>
      <td>beta_nov16</td>
      <td>11.5977</td>
    </tr>
    <tr>
      <th>1711_ref2015_12</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.5</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
      <td>12</td>
      <td>False</td>
      <td>-0.882054</td>
      <td>30</td>
      <td>ref2015</td>
      <td>13.0197</td>
    </tr>
    <tr>
      <th>1711_beta_nov16_12</th>
      <td>1STN</td>
      <td>1711</td>
      <td>139</td>
      <td>A</td>
      <td>I</td>
      <td>A</td>
      <td>0.326741</td>
      <td>3.5</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1STN_A</td>
      <td>I139A</td>
      <td>P00644</td>
      <td>12</td>
      <td>False</td>
      <td>2.99445</td>
      <td>30</td>
      <td>beta_nov16</td>
      <td>13.7258</td>
    </tr>
    <tr>
      <th>3_ref2015_8</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.6</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
      <td>8</td>
      <td>False</td>
      <td>1.65583</td>
      <td>12</td>
      <td>ref2015</td>
      <td>9.00631</td>
    </tr>
    <tr>
      <th>3_beta_nov16_8</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.6</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
      <td>8</td>
      <td>False</td>
      <td>1.27046</td>
      <td>12</td>
      <td>beta_nov16</td>
      <td>9.3355</td>
    </tr>
    <tr>
      <th>3_ref2015_10</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.6</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
      <td>10</td>
      <td>False</td>
      <td>-0.544857</td>
      <td>19</td>
      <td>ref2015</td>
      <td>12.2666</td>
    </tr>
    <tr>
      <th>3_beta_nov16_10</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.6</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
      <td>10</td>
      <td>False</td>
      <td>-3.05145</td>
      <td>19</td>
      <td>beta_nov16</td>
      <td>13.251</td>
    </tr>
    <tr>
      <th>3_ref2015_12</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.6</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
      <td>12</td>
      <td>False</td>
      <td>-0.615064</td>
      <td>34</td>
      <td>ref2015</td>
      <td>17.3526</td>
    </tr>
    <tr>
      <th>3_beta_nov16_12</th>
      <td>1AAR</td>
      <td>3</td>
      <td>45</td>
      <td>A</td>
      <td>F</td>
      <td>W</td>
      <td>0.150026</td>
      <td>-0.6</td>
      <td>'hydrophobic to hydrophobic' 'buried'</td>
      <td>1AAR_A</td>
      <td>F45W</td>
      <td>P0CH28</td>
      <td>12</td>
      <td>False</td>
      <td>0.35734</td>
      <td>34</td>
      <td>beta_nov16</td>
      <td>17.6672</td>
    </tr>
    <tr>
      <th>2256_ref2015_8</th>
      <td>1TUP</td>
      <td>2256</td>
      <td>273</td>
      <td>A</td>
      <td>R</td>
      <td>H</td>
      <td>0.468525</td>
      <td>0.35</td>
      <td>'positive to non-charged polar' 'buried'</td>
      <td>1TUP_A</td>
      <td>R273H</td>
      <td>P04637</td>
      <td>8</td>
      <td>False</td>
      <td>0.465832</td>
      <td>13</td>
      <td>ref2015</td>
      <td>20.152</td>
    </tr>
    <tr>
      <th>2256_beta_nov16_8</th>
      <td>1TUP</td>
      <td>2256</td>
      <td>273</td>
      <td>A</td>
      <td>R</td>
      <td>H</td>
      <td>0.468525</td>
      <td>0.35</td>
      <td>'positive to non-charged polar' 'buried'</td>
      <td>1TUP_A</td>
      <td>R273H</td>
      <td>P04637</td>
      <td>8</td>
      <td>False</td>
      <td>2.05686</td>
      <td>13</td>
      <td>beta_nov16</td>
      <td>21.0362</td>
    </tr>
  </tbody>
</table>
</div>




```python

```
