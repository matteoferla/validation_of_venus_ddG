## GoF and LoF dataset

An additional dataset was considered: [the dataset](
https://www.liebertpub.com/doi/suppl/10.1089/gtmb.2010.0036/suppl_file/Supp_Data.pdf)
from [the following paper](https://www.liebertpub.com/doi/10.1089/gtmb.2010.0036):

    Sarah E. Flanagan, Ann-Marie Patch, and Sian Ellard.Genetic Testing and Molecular Biomarkers.Aug 2010.533-537.

This dataset contains GoF and LoF variants for three protein.
However, the GoF could be explained nicely by visual inspection —exactly what Venus is for!—,
but the data needed too much fiddlying to prove non-visually and 
also MF did not have time to read the literature discussing this.

As a result the code here is just dumped.

```python
import pandas as pd

data = pd.read_csv('rough.csv')

def is_yes(v:str):
    if not isinstance(v, str):
        return None
    return 'Yes' in v

data['Functional studies'] = data['Functional studies'].apply(is_yes)
data['Response to sulphonylureas'] = data['Response to sulphonylureas'].apply(is_yes)
data['SIFT prediction'].unique()
data['SIFT_tolerance'] = data['SIFT prediction'].map({'tolerated': True, 'not_tolerated': False})
data['PolyPhen_damaging'] = data['PolyPhen prediction'].map({'benign': False, 'possibly_damaging': True, 'probably_damaging': True})
data['protein'] = data['group'].apply(lambda v: v.split('-')[0])
data['lof'] = data['group'].str.contains('Loss')
data['gof'] = data['group'].str.contains('Gain')
data['uniprot'] = data['protein'].map({'ABCC8': 'Q09428', 'KCNJ11': 'Q14654', 'GCK': 'P35557'})
data['mutation'] = data['Mutation name']
del data['Mutation name']
del data['SIFT prediction']
del data['PolyPhen prediction']
data.to_csv('polished.csv')
```

...

```python
import os
from michelanglo_protein import ProteinCore, Structure, ProteinAnalyser, Mutation, global_settings
global_settings.startup(os.environ['MICHELANGLO_PROTEIN_DATA'])

import pickle

with open('ddGs.p', 'rb') as fh:
    ddG_map = pickle.load(fh)

# direct copypaste
variants_folder_name = 'variants'

from michelanglo_protein import Variant, ProteinAnalyser
from typing import (Union)

def get_variant(protein, grouping: str) -> Union[Variant, None]:
    for variant in getattr(protein, grouping):
        if variant.mutation == protein.mutation.mutation:
            return variant
    else:
        None
        
from functools import partial

get_clinvar = partial(get_variant, grouping='clinvar')
get_gnomAD = partial(get_variant, grouping='gnomAD')
        
def get_clinvar_details(protein) -> dict:
    clinvar = get_clinvar(protein)
    if clinvar:
        return dict(clinvar_id=clinvar.id,
                    clinvar_description=clinvar.description,
                    clinvar_homozygous=clinvar.homozygous,
                    clinvar_consequence=clinvar.consequence,
                    clinvar_impact = clinvar.impact,
                    clinvar_N=clinvar.N
                   )
    return {}

def get_gnomAD_details(protein) -> dict:
    gnomAD = get_gnomAD(protein)
    if gnomAD:
        return dict(gnomAD_id=gnomAD.id,
                    gnomAD_description=gnomAD.description,
                    gnomAD_homozygous=gnomAD.homozygous,
                    gnomAD_consequence=gnomAD.consequence,
                    # gnomAD_impact = gnomAD.impact, pointless.
                    gnomAD_N=gnomAD.N,
                    gnomAD_frequency=gnomAD.frequency,
                    gnomAD_frequencies=gnomAD.frequencies
                   )
    return {}

def get_summary(protein) -> dict:
    general = dict(gene_name=protein.gene_name,
                uniprot=protein.uniprot,
                mutation=str(protein.mutation),
                ddG=protein.energetics['ddG'],
                dsol=protein.energetics['dsol'],
                rsa=protein.structural.RSA,
                features=protein.get_features_at_position(),
                distance_to_closest_ligand=protein.structural.distance_to_closest_ligand,
                closest_ligand=protein.structural.closest_ligand,
                neighbor_ptm_details=[n for n in protein.structural.neighbours if len(n['ptms'])],
                neighbor_other_chain_details=[n for n in protein.structural.neighbours if n['other_chain'] and len(n['resn']) == 1], # filter for peptide
                neighbor_gnomAD_details=[n for n in protein.structural.neighbours for g in n['gnomads']],
                neighbor_clinvar_details=[n for n in protein.structural.neighbours for g in n['clinvars']],
                neighbor_ligand_details=[n for n in protein.structural.neighbours if len(n['resn']) > 1],
               )
    return {**general, **get_gnomAD_details(protein), **get_clinvar_details(protein)}

import copy

import os
# ddG_map :Dict[Tuple[str, str], float] imported.

def get_path_name(row:pd.Series):
    subfolder_name = row.uniprot[:3]
    file_path = os.path.join(variants_folder_name, subfolder_name, row.uniprot + '_'+row.mutation+'.p')
    if not os.path.exists(file_path):
        return None
    return file_path

def get_protein(file_path) -> Union[ProteinAnalyser, None]:
    if not file_path:
        return {}
    with open(file_path, 'rb') as fh:
        return pickle.load(fh)

def get_datum(file_path):
    """
    Loads the protein and returns the get_summary dict
    """
    protein :Union[ProteinAnalyser, None] = get_protein(file_path)
    if protein:
        return get_summary(protein)
    return {}

        
def analyse_mutation(protein, mutation):
    try:
        mutant = copy.deepcopy(protein)
        mutant.mutation = Mutation(mutation)
        if not mutant.check_mutation():
            return 'discrepant'
        mutant.predict_effect()
        # protein.get_best_model() --> none if no structure
        mutant.analyse_structure(no_conservation=True) # protein.structural.structure
        # protein.structural = StructureAnalyser
        mutant.analyse_FF(spit_process=False) # protein.energetics is filled, e.g. protein.energetics['ddG']
        if not os.path.exists(f'variants/{protein.uniprot[:3]}'):
            os.mkdir(f'variants/{protein.uniprot[:3]}')
        with open(f'variants/{protein.uniprot[:3]}/{protein.uniprot}_{mutation}.p', 'wb') as fh:
            pickle.dump(obj=mutant, file=fh)
        return 'successful'
    except Exception as error:
        print(f'{type(error)}: {error}')
        return 'error'

def analyse_this(protein, mutation):
    protein.mutation = Mutation(mutation)
    if not protein.check_mutation():
        print(protein.mutation_discrepancy())
        return
        #raise ValueError()
    analyse_mutation(protein, mutation)
```

Unfortunately, one protein is running off the non-canonical isoform

```python
uniprot = 'Q09428'
taxid = 9606
protein = ProteinAnalyser(uniprot=uniprot, taxid=taxid).load()

iso2 = '''MPLAFCGSENHSAAYRVDQGVLNNGCFVDALNVVPHVFLLFITFPILFIGWGSQSSKVHI
HHSTWLHFPGHNLRWILTFMLLFVLVCEIAEGILSDGVTESHHLHLYMPAGMAFMAAVTS
VVYYHNIETSNFPKLLIALLVYWTLAFITKTIKFVKFLDHAIGFSQLRFCLTGLLVILYG
MLLLVEVNVIRVRRYIFFKTPREVKPPEDLQDLGVRFLQPFVNLLSKGTYWWMNAFIKTA
HKKPIDLRAIGKLPIAMRALTNYQRLCEAFDAQVRKDIQGTQGARAIWQALSHAFGRRLV
LSSTFRILADLLGFAGPLCIFGIVDHLGKENDVFQPKTQFLGVYFVSSQEFLANAYVLAV
LLFLALLLQRTFLQASYYVAIETGINLRGAIQTKIYNKIMHLSTSNLSMGEMTAGQICNL
VAIDTNQLMWFFFLCPNLWAMPVQIIVGVILLYYILGVSALIGAAVIILLAPVQYFVATK
LSQAQRSTLEYSNERLKQTNEMLRGIKLLKLYAWENIFRTRVETTRRKEMTSLRAFAIYT
SISIFMNTAIPIAAVLITFVGHVSFFKEADFSPSVAFASLSLFHILVTPLFLLSSVVRST
VKALVSVQKLSEFLSSAEIREEQCAPHEPTPQGPASKYQAVPLRVVNRKRPAREDCRGLT
GPLQSLVPSADGDADNCCVQIMGGYFTWTPDGIPTLSNITIRIPRGQLTMIVGQVGCGKS
SLLLAALGEMQKVSGAVFWSSSLPDSEIGEDPSPERETATDLDIRKRGPVAYASQKPWLL
NATVEENIIFESPFNKQRYKMVIEACSLQPDIDILPHGDQTQIGERGINLSGGQRQRISV
ARALYQHANVVFLDDPFSALDIHLSDHLMQAGILELLRDDKRTVVLVTHKLQYLPHADWI
IAMKDGTIQREGTLKDFQRSECQLFEHWKTLMNRQDQELEKETVTERKATEPPQGLSRAM
SSRDGLLQDEEEEEEEAAESEEDDNLSSMLHQRAEIPWRACAKYLSSAGILLLSLLVFSQ
LLKHMVLVAIDYWLAKWTDSALTLTPAARNCSLSQECTLDQTVYAMVFTVLCSLGIVLCL
VTSVTVEWTGLKVAKRLHRSLLNRIILAPMRFFETTPLGSILNRFSSDCNTIDQHIPSTL
ECLSRSTLLCVSALAVISYVTPVFLVALLPLAIVCYFIQKYFRVASRDLQQLDDTTQLPL
LSHFAETVEGLTTIRAFRYEARFQQKLLEYTDSNNIASLFLTAANRWLEVRMEYIGACVV
LIAAVTSISNSLHRELSAGLVGLGLTYALMVSNYLNWMVRNLADMELQLGAVKRIHGLLK
TEAESYEGLLAPSLIPKNWPDQGKIQIQNLSVRYDSSLKPVLKHVNALIAPGQKIGICGR
TGSGKSSFSLAFFRMVDTFEGHIIIDGIDIAKLPLHTLRSRLSIILQDPVLFSGTIRFNL
DPERKCSDSTLWEALEIAQLKLVVKALPGGLDAIITEGGENFSQGQRQLFCLARAFVRKT
SIFIMDEATASIDMATENILQKVVMTAFADRTVVTIAHRVHTILSADLVIVLKRGAILEF
DKPEKLLSRKDSVFASFVRADK'''.replace('\n', '')


import pyrosetta_help as ph

alignments = ph.get_alignment(target=iso2, template=protein.sequence)

def make_map(al, pose_offset=0):
    # ungapped to w/ gapped
    # pose index (fortran-style) to MSA index (C++-style)
    gap_map = [i for i, r in enumerate(al) if r != '-']
    return {pose_offset + ungap_i: gap_i for ungap_i, gap_i in enumerate(gap_map)}

isotwo2al = make_map(alignments['target'])  #iso2
canon2al = make_map(alignments['template']) # canon
al2canon = dict(zip(canon2al.values(), canon2al.keys()))
isotwo2canon = {k: al2canon[g] for k, g in isotwo2al.items() if g in al2canon}

data['original_position'] = data.mutation.str.extract(r'\w(\d+)\w').astype(int)

import re

def convert(row):
    if row.uniprot == 'Q09428':
        return re.sub(r'(\w)(\d+)(\w)', lambda match: f'{match.group(1)}{isotwo2canon[int(match.group(2))]}{match.group(3)}', row.mutation)
    return row.mutation

data['original_mutation'] = data['mutation']
data['mutation'] = data.apply(convert, axis=1)

def is_missing(row:pd.Series):
    subfolder_name = row.uniprot[:3]
    file_path = os.path.join(variants_folder_name, subfolder_name, row.uniprot + '_'+row.mutation+'.p')
    if not os.path.exists(file_path):
        return True
    return False

len(data.apply(is_missing, axis=1)) #.protein.value_counts()
```

```python
%%script false --no-raise-error

## Add data!
# Waters and chain needs fixing though!
# does `get_coordinates` do that?

taxid = 9606
uniprot = 'Q14654'
protein = ProteinAnalyser(uniprot=uniprot, taxid=taxid)
if not protein.exists():
    raise ValueError
protein.load()
protein.add_alphafold2()
protein.retrieve_structures_from_swissmodel()
for model in protein.pdbs + protein.swissmodel + protein.alphafold2:
    model.get_coordinates()
protein.dump()
```

...
```python
%%script false --no-raise-error

taxid = 9606

for uniprot in {'ABCC8': 'Q09428', 'KCNJ11': 'Q14654', 'GCK': 'P35557'}.values():
    protein = ProteinAnalyser(uniprot=uniprot, taxid=taxid).load()
    data.loc[data.apply(is_missing, axis=1) & (data.uniprot == uniprot)].mutation.apply(lambda mutation: analyse_this(protein, mutation))
```
Analysis
```python
data['data_path'] = data.apply(get_path_name, axis=1)
summaries = pd.DataFrame(data.data_path.apply(get_datum).tolist(), index=data.index)
del summaries['uniprot']
del summaries['mutation']
analysis = pd.concat([data, summaries], axis='columns', verify_integrity=True)
analysis
```
Adding key values (see analyses.md)
```python
def map_into_ddG(row: pd.Series) -> None:
    get_ddG = lambda variant_name: ddG_map[row.uniprot, variant_name] if (row.uniprot, variant_name) in ddG_map else float('nan')
    if str(row.neighbor_gnomAD_details) == 'nan':
        return
    for neigh in row.neighbor_gnomAD_details: #type: dict
        assert isinstance(neigh, dict), f'{type(neigh)} is not dict'
        neigh['gnomads_ddG'] = [get_ddG(variant_name)  for variant_name in neigh['gnomads']]
    for neigh in row.neighbor_clinvar_details: #type: dict
        neigh['clinvars_ddG'] = [get_ddG(variant_name)  for variant_name in neigh['clinvars']]

# trinary
def is_benign(impact) -> Union[None, bool]:
    if str(impact) == 'nan':
        return None
    elif 'pathogenic' in impact.lower():
        return False
    elif 'benign' in impact.lower():
        return True
    else:
        return None
    
def get_distance_to_closest_gnomAD(neighs):
    distances = [neigh['distance'] for neigh in neighs]
    if len(distances):
        return sorted(distances)[0]
    return float('nan')
        
from functools import partial

def get_distance_to_closest_ddG_gnomAD(neighs, cutoff=1):
    # gnomads_ddG is list, gnomads is list
    distances = [neigh['distance'] for neigh in neighs if max(neigh['gnomads_ddG']) > cutoff]
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

import json

with open('hmz.json', 'r') as fh:
    hmzs = list(map(tuple, json.load(fh)))
    
def get_distance_to_hmz_gnomAD(row: pd.Series):
    if not isinstance(row.neighbor_gnomAD_details, list):
        return float('nan')
    distances = [neigh['distance'] for neigh in row.neighbor_gnomAD_details for gnomad in neigh['gnomads'] if (row.uniprot, gnomad) in hmzs]
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

def get_distance_to_closest_other_chain(neighs):
    distances = [neigh['distance'] for neigh in neighs]
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

def get_distance_to_closest_ptm(neighs):
    distances = [neigh['distance'] for neigh in neighs]
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

def get_distance_to_closest_phospho(neighs):
    distances = [neigh['distance'] for neigh in neighs if 'phosphorylated' in neigh['ptms']]
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

from functools import partial

def get_distance_to_closest_hmz_ddG_gnomAD(row: pd.Series, cutoff=1):
    # gnomads_ddG is list, gnomads is list
    distances = []
    if not isinstance(row.neighbor_gnomAD_details, list):
        return float('nan')
    for neigh in row.neighbor_gnomAD_details:
        for gnomad, ddG in zip(neigh['gnomads'], neigh['gnomads_ddG']):
            if ddG > cutoff and (row.uniprot, gnomad) in hmzs:
                distances.append(neigh['distance'])
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

def get_stability_cat(ddG):
    if ddG >= 2:
        return 'destable'
    elif ddG >= 1:
        return 'mildy_destable'  # ie. error ambiguous.
    else:
        return 'stable'
    
al_cutoff = 5e-4
print(al_cutoff**2*5e7) # UK population is 6.7e7
al = analysis.loc[analysis.gnomAD_frequency >= al_cutoff]
als = pd.Series(zip(al.uniprot.values, al.mutation.values)).to_list()

def get_distance_to_closest_common_ddG_gnomAD(row: pd.Series, cutoff=1):
    distances = []
    if not isinstance(row.neighbor_gnomAD_details, list):
        return float('nan')
    for neigh in row.neighbor_gnomAD_details:
        for gnomad, ddG in zip(neigh['gnomads'], neigh['gnomads_ddG']):
            if ddG > cutoff and (row.uniprot, gnomad) in als:
                distances.append(neigh['distance'])
    if len(distances):
        return sorted(distances)[0]
    return float('nan')

distance_cutoff = 10  # chosen arbitrarily, less than 12 Å.
is_within_distance = lambda d: d <= distance_cutoff

def get_stability_cat(ddG):
    if ddG >= 2:
        return 'destable'
    elif ddG >= 1:
        return 'mildy_destable'  # ie. error ambiguous.
    else:
        return 'stable'

nan2list = lambda d: d if isinstance(d, list) else []

def tweak_analysis(analysis: pd.DataFrame):
    analysis.apply(map_into_ddG, axis=1)
    analysis['N_neighbor_gnomAD'] = analysis.neighbor_gnomAD_details.apply(nan2list).apply(len)
    analysis['N_neighbor_clinvar'] = analysis.neighbor_clinvar_details.apply(nan2list).apply(len)
    analysis['N_neighbor_other_chain'] = analysis.neighbor_other_chain_details.apply(nan2list).apply(len)
    analysis['N_neighbor_ptm_details'] = analysis.neighbor_ptm_details.apply(nan2list).apply(len)
    # analysis['is_benign'] = analysis.clinvar_impact.apply(is_benign)
    # analysis['is_benign_cat'] = analysis.is_benign.map({True: 'benign', False: 'pathogenic', None: 'indeterminate'})
    analysis['distance_to_closest_gnomAD'] = analysis.neighbor_gnomAD_details.apply(nan2list).apply(get_distance_to_closest_gnomAD)
    analysis['distance_to_closest_one_destabilising_gnomAD'] = analysis.neighbor_gnomAD_details.apply(nan2list).apply(partial(get_distance_to_closest_ddG_gnomAD, cutoff=1))
    analysis['distance_to_closest_two_destabilising_gnomAD'] = analysis.neighbor_gnomAD_details.apply(nan2list).apply(partial(get_distance_to_closest_ddG_gnomAD, cutoff=2))
    analysis['distance_to_closest_phospho'] = analysis.neighbor_ptm_details.apply(nan2list).apply(get_distance_to_closest_phospho)
    analysis['distance_to_hmz_gnomAD'] = analysis.apply(get_distance_to_hmz_gnomAD, axis=1)
    analysis['distance_to_closest_other_chain'] = analysis.neighbor_other_chain_details.apply(nan2list).apply(get_distance_to_closest_other_chain)
    analysis['distance_to_closest_ptm'] = analysis.neighbor_ptm_details.apply(nan2list).apply(get_distance_to_closest_ptm)
    analysis['distance_to_closest_hmz_one_destabilising_gnomAD'] = analysis.apply(partial(get_distance_to_closest_hmz_ddG_gnomAD, cutoff=1), axis=1)
    analysis['distance_to_closest_hmz_two_destabilising_gnomAD'] = analysis.apply(partial(get_distance_to_closest_hmz_ddG_gnomAD, cutoff=2), axis=1)
    analysis['distance_to_closest_common_one_destabilising_gnomAD'] =analysis.apply(partial(get_distance_to_closest_common_ddG_gnomAD, cutoff=1), axis=1)
    analysis['distance_to_closest_common_two_destabilising_gnomAD'] = analysis.apply(partial(get_distance_to_closest_common_ddG_gnomAD, cutoff=2), axis=1)
    analysis['is_mildly_unstable'] = analysis.ddG.apply(lambda v: v >= 1)
    analysis['is_unstable'] = analysis.ddG.apply(lambda v: v >= 2)
    analysis['has_mildly_unstable_gnomAD'] = analysis.distance_to_closest_one_destabilising_gnomAD.apply(is_within_distance)
    analysis['has_unstable_gnomAD'] = analysis.distance_to_closest_two_destabilising_gnomAD.apply(is_within_distance)
    analysis['has_mildly_unstable_hmz_gnomAD'] = analysis.distance_to_closest_hmz_one_destabilising_gnomAD.apply(is_within_distance)
    analysis['has_unstable_hmz_gnomAD'] = analysis.distance_to_closest_hmz_two_destabilising_gnomAD.apply(is_within_distance)
    analysis['has_mildly_unstable_common_gnomAD'] = analysis.distance_to_closest_common_one_destabilising_gnomAD.apply(is_within_distance)
    analysis['has_unstable_common_gnomAD'] = analysis.distance_to_closest_common_two_destabilising_gnomAD.apply(is_within_distance)
    analysis['has_gnomAD'] = analysis.distance_to_closest_gnomAD.apply(is_within_distance)
    analysis['has_other_chain'] = analysis.distance_to_closest_other_chain.apply(is_within_distance)
    analysis['has_ptm'] = analysis.distance_to_closest_ptm.apply(is_within_distance)
    #clean['has_close_ptm'] = clean.distance_to_closest_ptm.apply(lambda d: d <= 5)
    analysis['has_ligand'] = analysis.distance_to_closest_ligand.apply(is_within_distance)
    analysis['has_phospho'] = analysis.distance_to_closest_phospho.apply(is_within_distance)
    def has_features(name, features) -> bool:
        if not isinstance(features, list): # NaN!
            return False
        return any([True if entry['type'] == name else False for entry in features])
    
    analysis['is_transmembrane'] = analysis['features'].apply(lambda features: has_features('transmembrane region', features))
    analysis['is_other_mod'] = analysis['features'].apply(lambda features: has_features('modified residue', features))
    analysis['is_NTP'] = analysis['features'].apply(lambda features: has_features('nucleotide phosphate-binding region', features))
    analysis['is_ROI'] = analysis['features'].apply(lambda features: has_features('region of interest', features))
    # unique() returns np.array hence the tolist not to_list
    # and the closest ligand is a selctor `[ATP]4.H11:_`
    get_ligands_for_uniprot = lambda uniprot: analysis.loc[analysis.uniprot == uniprot].closest_ligand.str.extract(r'\[(.*)\]')[0].unique().tolist()
    ligands = {uniprot: get_ligands_for_uniprot(uniprot)  for uniprot in analysis.uniprot.unique()}
    max_n_ligands = max(map(len, ligands.values()))
    ligand_name = analysis.closest_ligand.str.extract(r'\[(.*)\]')[0]
    for n in range(max_n_ligands):
        # this is a simplification and will fail if ligand X is close but Y is closer...
        analysis[f'has_ligand_{n}'] = ligand_name == analysis.uniprot.apply(lambda u: ligands[u][n] if len(ligands[u]) > n else float('nan'))
        # overly distant: remove
        analysis.at[~analysis['has_ligand'], f'has_ligand_{n}'] = False
```

...
```python
# This for human viewing, therefore the variables are purposefully in British English.

import numpy as np
from IPython.display import display, HTML
from typing import Tuple

def summarise(df,
              column2nice, 
              indices=['uniprot', 'effect'], 
              aggfunc=sum,
              ratio:Union[None, Tuple[str, str]]=None):
    
    columns = ['counts'] + list(column2nice.keys())
    xtabbed = pd.pivot_table(df,
                             values=list(column2nice.keys()),
                             index=indices,
                             aggfunc= aggfunc,
                            )
    display(HTML(f"<h3>aggfunc='{aggfunc.__name__}'</h3>"))
    xtabbed['counts'] = pd.pivot_table(df,values='ddG',index=indices, aggfunc=len)['ddG']
    xtabbed = xtabbed[columns]
    xtabbed.index.name = None
    xtabbed.columns = xtabbed.columns.to_series().map({'counts': 'Counts', **column2nice})
    txtabbed = xtabbed.transpose()
    if ratio:
        nominator = ratio[0]
        denominator = ratio[1]
        txtabbed[f'{nominator} (fraction of total)'] = txtabbed[nominator]/txtabbed.loc['Counts', nominator]
        txtabbed[f'{denominator} (fraction of total)'] = txtabbed[denominator]/txtabbed.loc['Counts', denominator]
        txtabbed['Ratio'] = (txtabbed[f'{nominator} (fraction of total)'] / txtabbed[f'{denominator} (fraction of total)'])
        txtabbed.columns = txtabbed.columns.to_series().str.capitalize()
    return txtabbed.round(2)


get_ligands_for_uniprot = lambda uniprot: analysis.loc[analysis.uniprot == uniprot].closest_ligand.str.extract(r'\[(.*)\]')[0].unique().tolist()
ligands = {uniprot: get_ligands_for_uniprot(uniprot)  for uniprot in analysis.uniprot.unique()}
max_n_ligands = max(map(len, ligands.values()))

def distance_to_specific_ligand(row, n):
    if not isinstance(row.neighbor_ligand_details, list):
        return float('nan')
    ds = [neigh['distance'] for neigh in row.neighbor_ligand_details if neigh['resn'] == ligands[row.uniprot][n]]
    if ds:
        return min(ds)
    return float('nan')

from functools import partial

for n in range(max_n_ligands):
    dsln = partial(distance_to_specific_ligand, n=n)
    analysis[f'distance_to_ligand_{n}'] = analysis.apply(dsln, axis=1)
    analysis[f'has_ligand_{n}'] = analysis[f'distance_to_ligand_{n}'].apply(lambda d: d <= 10)

column2nice = {'is_unstable': 'Destabilising (> 2 kcal/mol)',
           'is_mildly_unstable': 'Destabilising (> 1 kcal/mol)',
                         'has_ligand': 'Close to ligand',
                         'has_ligand_0': 'Close to ligand_0',
                         'has_ligand_1': 'Close to ligand_1',
                         'has_ligand_2': 'Close to ligand_2',
               'is_transmembrane': 'Transmembrane',
               'is_other_mod': 'modified residue',
               'is_NTP': 'in nucleotide phosphate-binding region',
               'is_ROI': 'in "Region of Interest"',
         'has_other_chain': 'Close to interface',
         'has_ptm': 'Close to post-translational-modification',
         #'has_close_ptm': 'Very close to post-translational-modification',
        'has_phospho': 'Close to phosphorylation',
         'has_unstable_hmz_gnomAD': 'Close to destabilising homozygous gnomADs (> 2 kcal/mol)',
         'has_mildly_unstable_hmz_gnomAD': 'Close to destabilising homozygous gnomADs (> 1 kcal/mol)',
         'has_unstable_common_gnomAD': f'Close to destabilising common gnomADs ({al_cutoff}, > 2 kcal/mol)',
         'has_mildly_unstable_common_gnomAD': f'Close to destabilising common gnomADs ({al_cutoff}, > 1 kcal/mol)',
         'has_mildly_unstable_gnomAD': 'Close to destabilising gnomADs (> 1 kcal/mol)',
         'has_unstable_gnomAD': 'Close to destabilising gnomADs  (> 2 kcal/mol)',}

clean = analysis.loc[~analysis.ddG.isna()]
summarise(clean, column2nice)
```
...
```python
def task(args):
    uniprot, mutation = args
    taxid = 9606
    try:
        if os.path.exists(os.path.join(variants_folder_name, uniprot[:3], uniprot + '_'+mutation+'.p')):
            print(uniprot, mutation, 'previously done')
            return
        protein = ProteinAnalyser(uniprot=uniprot, taxid=taxid).load()
        protein.mutation = Mutation(mutation)
        if protein.check_mutation():
            analyse_mutation(protein, mutation)
        else:
            print(protein.mutation_discrepancy())
    except Exception as error:
        print(error)
            
from multiprocessing import Pool

with Pool(20) as p:
    p.map(task, neighbors)
```
...
```python
ddG_map = {}

def warn_about(uniprot, mutation):
    print()
    protein = ProteinAnalyser(uniprot=uniprot, taxid=taxid).load()
    try:
        protein.mutation = Mutation(mutation)
        if protein.check_mutation():
            print('No', uniprot, mutation, '--- Should have worked')
        else:
            print('No', uniprot, mutation, protein.mutation_discrepancy())
    except Exception as error:
        print('No', uniprot, mutation, error)

for uniprot, mutation in neighbors:
    subfolder_name = uniprot[:3]
    file_path = os.path.join(variants_folder_name, subfolder_name, uniprot + '_'+mutation+'.p')
    if not os.path.exists(file_path):
        warn_about(uniprot, mutation)
        continue
    protein :Union[ProteinAnalyser, None] = get_protein(file_path)
    if not protein:
        warn_about(uniprot, mutation)
        continue
    datum = get_summary(protein)
    ddG_map[datum["uniprot"], datum["mutation"]] = datum["ddG"]
```