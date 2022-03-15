"""
Score dataset.

    $ python analyse.py test test.csv O2567_protein.p settings.json 10
    # dataset_name dataset_filename protein_filename settings_filename max_workers
"""
import argparse
import json
import os
import pickle
import pymol2
import time
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from functools import partial
from typing import Any
import pandas as pd
from michelanglo_protein import ProteinCore, ProteinAnalyser, Mutation  # noqa   It is there.
from michelanglo_protein.analyse import StructureAnalyser  # noqa   It is there.
from sqlitedict import SqliteDict

def check_sequence(protein: ProteinCore):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.read_pdbstr(protein.pdbblock, 'xxx')
        atom = pymol.cmd.get_model(f'resi {protein.mutation.residue_index} and name CA').atom
        if atom:
            mapping = {name3.upper(): name1 for name1, name3, fullname in protein.mutation.names}
            assert mapping[atom[0].resn] == protein.mutation.from_residue, 'mismatch'
        else:
            raise ValueError(f'resi {protein.mutation.residue_index} and name CA does not exist')

class TableScorer:
    wanted_keys = ['ddG', 'scores', 'rmsd', 'dsol', 'score_fxn', 'neighbours',
                   'cycles', 'radius', 'neighbouring_ligand', 'n_constraints']

    errors = defaultdict(list)
    error_to_catch = Exception  # note that thread pool catches error already: setting to () will be couterproductive
    default_setting = dict(debug=False,
                           radius=8,
                           cycles=1,
                           scorefxn_name='ref2015',
                           prevent_acceptance_of_incrementor=True,
                           neighbour_only_score=False,
                           outer_constrained=False,
                           single_chain=True,
                           remove_ligands=True,
                           use_pymol_for_neighbours=False)

    def __init__(self,
                 dataset_name: str,
                 reference: pd.DataFrame,
                 settings: dict,
                 results: SqliteDict,
                 outfolder: str,
                 models: dict,
                 max_workers: int = 1):
        """
        This does not bar define the attributes.
        call again to run.

        :param dataset_name:
        :param reference:
        :param settings:
        :param results:
        :param outfolder:
        :param models:
        :param max_workers:
        """
        self.reference = reference
        self.settings = settings
        self.max_workers = max_workers
        self.outfolder = outfolder
        self.dataset_name = dataset_name
        self.models = models
        self.results = results
        if not os.path.exists(self.outfolder):
            os.mkdir(self.outfolder)

    # --------------------------------------------------------------------------------------------------------
    def __call__(self):
        """
        Run in subthreads for saner globals.
        """
        with ThreadPoolExecutor(max_workers=self.max_workers) as pool:
            pool.map(self.safely_analyse_row, [row for i, row in self.reference.iterrows()])
        return self

    # --------------------------------------------------------------------------------------------------------
    def analyse_row(self, row, silent=False) -> None:
        # silent is different as it affects the reference values.
        if isinstance(row, pd.Series):
            row = row.to_dict()  # make into dictionary
        if silent:
            row['to_resn'] = row['from_resn']
            row['empirical_ddG'] = 0.
        for setting_name in self.settings.keys():
            self.results[f'{row["protein_id"]}_{row["mutation"]}_{setting_name}'] = {**row,
                                                                                     'dataset_name': self.dataset_name,
                                                                                     'condition':    setting_name,
                                                                                     }
        protein = self.get_protein(row)
        for setting_name, setting in self.settings.items():
            acc = f'{row["protein_id"]}_{row["mutation"]}_{setting_name}'
            analysed = self.analyse(pdbblock=protein.pdbblock,
                                    mutation=protein.mutation,
                                    acc=acc,
                                    **setting)
            collected = {**row,
                         'dataset_name': self.dataset_name,
                         'silent':       silent,
                         'condition':    setting_name,
                         'apriori':      protein.mutation.apriori_effect
                         }
            self.results[acc] = {**collected, **analysed}

    def safely_analyse_row(self, row, silent=False):
        """
        Silent is silent mutation A201A
        Disable via ``ts.error_to_catch = ()``

        :param row:
        :param silent:
        :return:
        """
        try:
            self.analyse_row(row, silent=silent)
        except self.error_to_catch as error:
            msg = f'{error} processing {row.protein_id}'
            self.errors[error.__class__.__name__].append(msg)

    def get_protein(self, row: pd.Series):
        protein = ProteinAnalyser(uniprot=row['uniprot'])  # do not retrieve.
        mut = row['from_resn'] + str(row['resi']) + row['to_resn']
        protein.mutation = Mutation(mut)
        model = self.models[row['protein_id']]
        protein.pdbs.append(model)
        protein.analyse_structure(no_conservation=True)
        if protein.structural is None:  # force it by mockery
            # protein.pdbblock is a dynamic attribute that reads protein.structural.coordinates (or the minimised...)
            protein.structural = StructureAnalyser.__new__(StructureAnalyser)
            protein.structural.coordinates = model.coordinates
        check_sequence(protein)  # operates on protein.mutation
        return protein

    def analyse(self,
                pdbblock,
                mutation,
                acc,
                **setting
                ):
        # the coordinates are passed via a backdoor n pdbbock
        self.fix_setting(setting)
        protein = ProteinAnalyser(scorefxn_name=setting['scorefxn_name'],
                                  cycles=setting['cycles'],
                                  radius=setting['radius'])
        protein.mutation = mutation
        tick = time.time()
        # `protein.pdbblock` checks first `protein.energetics['native']`
        # before `protein.structural.coordinates`
        protein.energetics = {'native': pdbblock}
        data = protein.analyse_FF(spit_process=not setting['debug'],
                                  neighbour_only_score=setting['neighbour_only_score'],
                                  outer_constrained=setting['outer_constrained'],
                                  remove_ligands=setting['remove_ligands'],
                                  prevent_acceptance_of_incrementor=setting['prevent_acceptance_of_incrementor'],
                                  use_pymol_for_neighbours=setting['use_pymol_for_neighbours'],
                                  single_chain=setting['single_chain'])
        tock = time.time()
        # ## save:
        for form in ('native', 'mutant'):
            if form in data:
                with open(f'{self.outfolder}/{acc}.{form}.pdb', 'w') as fh:
                    fh.write(data[form])
        # ## ready output
        reply = dict(time=tock - tick)
        if 'error' in data:
            reply['error'] = data['error']
        for k in self.wanted_keys:
            if k in data:
                reply[k] = data[k]
        return reply

    @classmethod
    def fix_setting(cls, setting: dict):
        for key, value in cls.default_setting.items():
            if key in setting:
                setting[key] = type(value)(setting[key])
            else:
                setting[key] = value

    # ------------------------------ utilities ---------------------------------------

    def add_errors_to_pickle(self, filename: str):
        """
        Lazily coded makes it a prepend, but shmeh
        """
        if os.path.exists(filename):
            with open(filename, 'rb') as fh:
                errors = pickle.load(file=fh)
                for error_name, msgs in errors.items():
                    self.errors[error_name].extend(msgs)
        with open(filename, 'wb') as fh:
            pickle.dump(obj=self.errors, file=fh)

    @staticmethod
    def read_protein_file(filename):
        with open(filename, 'rb') as fh:
            return pickle.load(file=fh)

    @staticmethod
    def read_settings_file(filename):
        with open(filename, 'r') as fh:
            return json.load(fp=fh)

    @staticmethod
    def read_db_file(filename):
        return SqliteDict(filename,
                          encode=json.dumps,
                          decode=json.loads,
                          autocommit=True)

    def to_df(self):
        dict_scores = {k: v for k, v in self.results.items() if isinstance(v, dict)}
        return pd.DataFrame.from_dict(dict_scores, orient='index')

    @classmethod
    def test_sample(cls,
                    dataset_name,
                    dataset_filename,
                    database_filename,
                    protein_filename,
                    settings_filename='minimal_settings.json'):
        """
        Test a subset of the dataset_filename

        :param dataset_name: experiment name
        :param dataset_filename: csv file with the data to score
        :param database_filename: sqlite3 filename for the results
        :param protein_filename: pickle of the protein structures
        :param settings_filename: json of the settings to test
        :return:
        """
        from IPython.display import display
        reference = pd.read_csv(dataset_filename, index_col=0)
        results = cls.read_db_file(database_filename)
        self = TableScorer(f'{dataset_name}-test',
                         reference=reference.sample(5),
                         settings=cls.read_settings_file(settings_filename),
                         results=results,
                         outfolder=f'{dataset_name}-test',
                         models=cls.read_protein_file(protein_filename),
                         max_workers=1)()
        print(ts.errors)
        scores = ts.to_df()
        display(scores)

    @classmethod
    def test_one(cls,
                    dataset_name,
                    dataset_filename,
                    database_filename,
                    protein_filename,
                    settings_filename='minimal_settings.json'):
        """
        Test a single row at random with error raised off.

        :param dataset_name: experiment name
        :param dataset_filename: csv file with the data to score
        :param database_filename: sqlite3 filename for the results
        :param protein_filename: pickle of the protein structures
        :param settings_filename: json of the settings to test
        :return:
        """
        from IPython.display import display
        reference = pd.read_csv(dataset_filename, index_col=0)
        results = cls.read_db_file(database_filename)
        self = TableScorer(f'{dataset_name}-test',
                         reference=reference.sample(1),
                         settings=cls.read_settings_file(settings_filename),
                         results=results,
                         outfolder=f'{dataset_name}-test',
                         models=cls.read_protein_file(protein_filename),
                         max_workers=1)()
        # analyse row is behind the safety net and thread pool (which secretly supresses errors...)
        self.analyse_row(self.reference.iloc[0])
        return self.to_df()


# ---------- MAIN ----------------------
if __name__ == '__main__':
    # ## prep
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('dataset_name', type=str, help='dataset name')
    parser.add_argument('dataset_filename', type=str, help='dataset filename')
    parser.add_argument('protein_filename', type=str, help='protein filename')
    parser.add_argument('settings_filename', type=str, help='settings filename')
    parser.add_argument('max_workers', type=int, help='max workers')

    args = parser.parse_args()
    for as_expected in ['.csv' in args.dataset_filename,
                        '.p' in args.protein_filename,
                        '.json' in args.settings_filename]:
        assert as_expected, 'wrong argument order.'

    dataset_name = args.dataset_name
    # ## run
    ts = TableScorer(dataset_name,
                     reference=pd.read_csv(args.dataset_filename, index_col=0),
                     settings=TableScorer.read_settings_file(args.settings_filename),
                     results=TableScorer.read_db_file(f'{dataset_name}-scores.db'),
                     outfolder=f'{dataset_name}-structures',
                     models=TableScorer.read_protein_file(args.protein_filename),
                     max_workers=args.max_workers)
    ts()
    # ## parse errors
    error_filename = f'{dataset_name}-errors.p'
    ts.add_errors_to_pickle(error_filename)
    print('complete')
