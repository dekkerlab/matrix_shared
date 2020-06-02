import os
import sys
import glob
import numpy as np
import pandas as pd
from cooler import Cooler
import matplotlib
import matplotlib.pyplot as plt
import h5py
import seaborn as sns
import shelve
from collections import defaultdict, Iterable

# cooler_path = '/net/levsha/share/lab/dekkerU54/new_files/'
# cooler_paths = ['/net/levsha/share/lab/U54/2019_mapping_hg38/U54_matrix/cooler_library/',
#                 '/net/levsha/share/lab/U54/2019_mapping_hg38/U54_deep/cooler_library_group/']

# dot_paths = ['/net/levsha/share/lab/U54/2019_mapping_hg38/U54_matrix/snakedots/',
#              '/net/levsha/share/lab/U54/2019_mapping_hg38/U54_deep/snakedots/']

# analysis_path = '/net/levsha/share/sameer/U54/hic_matrix/'
# db_path = '/net/levsha/share/sameer/U54/hic_matrix/metadata/U54_matrix_info'
hela_chroms = ['chr4', 'chr14', 'chr17', 'chr18', 'chr20', 'chr21']


class Database:
    
    def __init__(self, db_path):
        self.db_path = db_path
        
        if os.path.exists(f'{db_path}.dat'):
            with shelve.open(db_path, flag='r') as db:
                self.metadata = db['metadata']
                self.cooler_paths = db['cooler_paths']
                self.analysis_path = db['analysis_path']
                self.dot_paths = db['dot_paths']
                keys = list(db.keys())
                keys.remove('metadata')
                keys.remove('cooler_paths')
                keys.remove('analysis_path')
                keys.remove('dot_paths')
                self.keys  = keys
        else:
            self.metadata = None
            self.keys = []
            self.cooler_paths = ''
            self.analysis_path = ''
            self.dot_paths = ''
        
        
    def create_dataset(self, table, cool_paths, analysis_path, dot_paths):
        
        if os.path.exists(f'{self.db_path}.dat'):
            print(f'Database already exists at {self.db_path}')
            raise FileExistsError
        else:
            assert np.all([s in table.columns for s in ['lib_name', 'celltype', 'xlink', 'enzyme', 'cycle', 'seq']]), 'table does not contain the required columns'
            with shelve.open(self.db_path, flag='n') as db:
                db['cooler_paths'] = cool_paths
                self.cooler_paths = cool_paths
                
                db['analysis_path'] = analysis_path
                self.analysis_path = analysis_path
                
                db['dot_paths'] = dot_paths
                self.dot_paths = dot_paths
                
                db['metadata'] = table
                self.metadata = table
            
            
    def get_tables(self, keys=None):
        
        if self.metadata is None:
            print('Database not initialized')
            raise
        else:
            result = self.metadata
            
        if keys is None:
            return result
        elif isinstance(keys, str):
            keys = [keys]
            
        with shelve.open(self.db_path, flag='r') as db:
            for key in keys:
                assert key in self.keys, "Key not found in database"

                df = db[key]
                
                result = result.merge(df, on='lib_name', how='outer')
        
        return result

    
    def add_table(self, key, table):
        assert 'lib_name' in table.columns, "Please pass table with lib_names columns in it"
        table_lib_names = table['lib_name'].values
        
        with shelve.open(self.db_path, flag='w') as db:
            assert key not in self.keys, "Key already exists. If you wish to modify this, please use modify_table() method"
            meta_lib_names = db['metadata']['lib_name'].values
            assert np.all(meta_lib_names == table_lib_names), 'List of libraries does not match those in metadata'
            
            db[key] = table
            self.keys.append(key)
     
    
    def remove_table(self, key):
        bad_keys = ['metadata', 'cooler_paths', 'analysis_path', 'dot_paths']
        assert key not in bad_keys, f"The following keys should not be deleted: {bad_keys}"
            
        with shelve.open(self.db_path, flag='w') as db:
            assert key in self.keys, "Key not found in database"
            
            del db[key]
            self.keys.remove(key)
            
            
    def modify_table(self, key, new_table):
        assert 'lib_name' in new_table.columns, "Please pass table with lib_names columns in it"
        
        table_lib_names = new_table['lib_name'].values
        meta_lib_names = self.metadata['lib_name'].values
        assert np.all(meta_lib_names == table_lib_names), 'List of libraries does not match those in metadata'
        
        with shelve.open(self.db_path, flag='w') as db:
            assert key in self.keys, "Key not found in database. If you want to add a table, please use add_table() method"

            del db[key]
            db[key] = new_table
                
    
    def get_coolers(self, table, res=1000000):

        names = table['lib_name'].values
        cool_dict = defaultdict(list)

        for name in names:
            if name not in self.metadata['lib_name'].values:
                print(f'Name: {name} not found in metadata. Skipping')
                continue
                
            cool_dict['lib_name'].append(name)
            flag = True
            for cpath in self.cooler_paths:
                if f'{name}.hg38.mapq_30.1000.mcool' in os.listdir(cpath):
                    flag = False
                    cool = Cooler(cpath+f'{name}.hg38.mapq_30.1000.mcool::/resolutions/{res}')
                    cool_dict[f'cooler_{res}'].append(cool)

            if flag:
                print(f'Cooler not found matching {name}. Appending np.nan to appropriate row')
                cool_dict[f'cooler_{res}'].append(np.nan)

        df = pd.DataFrame(cool_dict)
        df = table.copy(deep=True).merge(df, on='lib_name', how='outer')
        return df
    
    
    def get_eigendecomps(self, table, res=1000000, subdir='eigdecomp/'):

        comp_path = self.analysis_path+f'{subdir}/{res}/'
        if not os.path.isdir(comp_path):
            print(f'{comp_path} is not a valid directory')
            raise
        
        names = table['lib_name'].values
        keys = ['lams', 'vectors']
        comp_dict = defaultdict(list)

        for name in names:
            if name not in self.metadata['lib_name'].values:
                print(f'Name: {name} not found in metadata. Skipping')
                continue
                
            comp_dict['lib_name'].append(name)

            for k in keys:
                if f'{name}.hdf5' in os.listdir(comp_path):
                    comp_dict[f'{k}_{res}'].append(
                        pd.read_hdf(comp_path+f'{name}.hdf5', key=k))
                else:
                    comp_dict[f'{k}_{res}'].append(np.nan)


        df = pd.DataFrame(comp_dict)
        df = table.copy(deep=True).merge(df, on='lib_name', how='outer')

        return df

    
    def get_scalings(self, table, subdir='scalings/global/', trans=False):
        scale_path = self.analysis_path+subdir
        if not os.path.isdir(scale_path):
            print(f'{scale_path} is not a valid directory')
            raise
            
        names = data['lib_name'].values
        if trans:
            keys = ['scaling','trans_lvl']
        else:
            keys = ['scaling']
        scale_dict = defaultdict(list)
        for name in names:
            if name not in self.metadata['lib_name'].values:
                print(f'Name: {name} not found in metadata. Skipping')
                continue
                
            scale_dict['lib_name'].append(name)

            if f'{name}.hdf5' in os.listdir(scale_path):
                scale_dict['scaling'].append(Scaling(scale_path+f'{name}.hdf5', keys))
            else:
                scale_dict['scaling'].append(np.nan)

        df = pd.DataFrame(scale_dict)
        df = data.copy(deep=True).merge(df, on='lib_name', how='outer')
    
        return df


    def get_pileups(self, table, subdir='pileup/dots/5000/', col_name='dots'):
        
        pileup_path = self.analysis_path+subdir
        if not os.path.isdir(pileup_path):
            print(f'{pileup_path} is not a valid directory')
            raise
            
        names = table['lib_name'].values
        pileup_dict = defaultdict(list)
        for name in names:
            if name not in self.metadata['lib_name'].values:
                print(f'Name: {name} not found in metadata. Skipping')
                continue
                
            pileup_dict['lib_name'].append(name)
            pileup_dict[col_name].append(Pileup(pileup_path, f'{name}.npy'))

        df = pd.DataFrame(pileup_dict)
        df = table.copy(deep=True).merge(df, on='lib_name', how='outer')

        return df


    def get_insulations(self, table, subdir='insulation/100kb_win_20px/', col_name='insulation'):
    
        ins_path = self.analysis_path+subdir
        if not os.path.isdir(ins_path):
            print(f'{ins_path} is not a valid directory')
            raise
            
        names = table['lib_name'].values
        ins_dict = defaultdict(list)
        
        for name in names:
            if name not in self.metadata['lib_name'].values:
                print(f'Name: {name} not found in metadata. Skipping')
                continue
                
            ins_dict['lib_name'].append(name)

            if f'{name}.txt' in os.listdir(ins_path):
                ins_dict[col_name].append(pd.read_csv(ins_path+f'{name}.txt', sep='\t'))
            else:
                ins_dict[col_name].append(np.nan)

        df = pd.DataFrame(ins_dict)
        df = table.copy(deep=True).merge(df, on='lib_name', how='outer')

        return df
    
    
class Scaling():
    # I created this for the same reason I created a Pileup class.
    def __init__(self, path, keys):
        self.path = path
        self.keys = keys
                
    def load(self):
        if not os.path.isfile(self.path):
            return np.nan
        if len(self.keys) == 1:
            return pd.read_hdf(self.path, key=self.keys[0])
            
        return tuple(pd.read_hdf(self.path, key=k) for k in self.keys)

    
class Pileup():
    # I save pileups as .npy files. I created this pileup class because loading all ~70 pileups for the matrix is very slow and memory intensive. So instead, when I ask the Dataset class to get me a pileup, it returns this Pileup object. Each pileup other can be loaded as needed.
    
    def __init__(self, path, file_name):
        if path[-1] != '/':
            path = f'{path}/'
        self.path = path
        self.name = file_name
                
    def load(self):
        if self.name not in os.listdir(self.path):
            return np.nan
        else:
            return np.load(f'{self.path}{self.name}')
        
        
def filter_data(df, filter_dict):
    for key in filter_dict.keys():
        assert key in df.columns, f'Column named {key} not found in DataFrame'
        
    out_df = df.copy()
    for dict_item in filter_dict.items():
        out_df = _find_matches(out_df, dict_item)
    
    out_df = out_df.sort_values(['celltype','xlink','enzyme','seq','cycle'])
    return out_df

def _find_matches(in_df, dict_item):
    col, val = dict_item
    
    if isinstance(val, str):
        out_df =  in_df[in_df[col] == val]
    elif isinstance(val, Iterable):
        out_df = []
        for item in val:
            out_df.append(_find_matches(in_df, (col, item)))
        out_df = pd.concat(out_df)
    else:
        out_df =  in_df[in_df[col] == val]
    
    return out_df