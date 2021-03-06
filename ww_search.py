from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
import os
import shutil
import multiprocessing as mp
from tqdm import tqdm
import time
import sys
import gzip
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT, GetUSR

from functools import partial



#debug = 'OFF'

#min_sim = 0.3
headers = "smiles\tidnumber\tsimilarity_usr_cat\tREF_smiles\tREF_idnumber\n"
#base_folder = '/media/dudenko/STORAGE/Projects/CS-Services/Shape-Based-Similarity/DEV/Base_TEST_MMFF94s_v4_default'
#result_file = f'results_xray-active-ligands.MMFF94s.v4.default.REF.txt'

#SDFile = 'RIV_STI_AQD_08K.sdf'

def set_var(min_sims, base_folders, SDFiles, result_files, debugs = 'OFF'):

    global min_sim
    #global base_folder
    #global SDFile
    #global result_file
    #global debug

    min_sim=min_sims
    base_folder=base_folders
    SDFile=SDFiles
    result_file=result_files
    debug = 'OFF'
    print(min_sim)
    print(result_files)

    return min_sim #, base_folder, SDFile, result_file, debug



def initializerSD(SDFile):
    global qshape2s
    global qsmiles
    global qids
    qshape2s = list()
    qsmiles = list()
    qids = list()
    suppl = Chem.SDMolSupplier(SDFile, removeHs=True)
    for mol_search in suppl:
        try:
            qshape2s.append(GetUSRCAT(mol_search))
            qsmiles.append(Chem.MolToSmiles(mol_search))
            qids.append(mol_search.GetProp("_Name"))
            print("reading ", Chem.MolToSmiles(mol_search), mol_search.GetProp("_Name"))
        except:
            print("Error with:", Chem.MolToSmiles(mol_search), mol_search.GetProp("_Name"))



def list_to_line(_list, sep: str = '\t') -> str:
    return f'{sep.join(_list)}\n'


def file_list(folder: str):
    fl = list()
    for root, dirs, files in os.walk(folder):
        for file in files:
            name, ext = os.path.splitext(os.path.basename(file))
            if ext.lower() in ['.gz', '.bz2']:
                fl.append(os.path.join(root, file))
    return fl


def human_format(num: int) -> str:
    num = float(f'{num:.3g}')
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    num = f'{num:.1f}'.rstrip('0').rstrip('.')
    let = ['', 'K', 'M', 'B', 'T'][magnitude]
    return f'{num}{let}'


def add_count_to_file(file: str, count: int) -> str:
    name, ext = os.path.splitext(os.path.basename(file))
    if ext.lower() in ['.gz', '.bz2']:
        ext = '.' + name.split('.')[-1] + ext
        name = '.'.join(name.split('.')[:-1])
    folder = os.path.dirname(file)
    new_file = os.path.join(folder, f'{name}_{human_format(count)}{ext}')
    os.replace(file, new_file)
    return new_file


def worker(file, min_sim,debug):
    results = list()
    usrcat_list = list()
    with gzip.open(file, 'rb') as f:
        while True:
            try:
                line, usrcat_list = pickle.load(f)
                for qshape2, qsmi, qid in zip(qshape2s, qsmiles, qids):
                    for usrcat_i in usrcat_list:
                        similarity_usrcat = GetUSRScore(qshape2, usrcat_i)
                        if (similarity_usrcat >= min_sim):
                            r = line.rstrip('\n').split('\t')
                            print("FOUND")
                            print(list_to_line([r[0], r[1], f'{similarity_usrcat:.4f}', qsmi, qid]))
                            results.append(([r[0], r[1], f'{similarity_usrcat:.4f}', qsmi, qid]))
            except EOFError:
                print("EOF Message after ", line)
                break
    results = sorted(results, key=lambda x:x[2], reverse = True)
    res_list = []
    res_mols = []
    for res_line in results:
        if res_line[1] in res_mols:
            if debug == 'ON':
                res_list.append(list_to_line(res_line))
        else:
            res_mols.append(res_line[1])
            res_list.append(list_to_line(res_line))
    results = res_list

    return results





def do_search(min_sims, base_folders, SDFiles, result_files, debugs = 'OFF'):

    set_var(min_sims, base_folders, SDFiles, result_files, debugs = 'OFF')
    print(min_sim)
    result_file=result_files
    base_folder=base_folders
    SDFile=SDFiles

    c = 0
    start = time.time()
    with open(result_file, 'wt') as f, mp.Pool(processes=1, initializer=initializerSD, initargs=(SDFile,)) as pool:
        f.write(headers)
        for res in tqdm(pool.imap(partial(worker, min_sim=min_sim, debug='OFF'), sorted(file_list(base_folder))), total=len(file_list(base_folder))):
            if res:
                c += len(res)
                f.writelines(res)
    add_count_to_file(result_file, c)
    end = time.time()
    print(f'Done, {c:_} results, {end - start:.1f} seconds')

if __name__ == '__main__':
    c = 0
    start = time.time()
    with open(result_file, 'wt') as f, mp.Pool(processes=1, initializer=initializerSD, initargs=(SDFile,)) as pool:
        f.write(headers)
        for res in tqdm(pool.imap(worker, sorted(file_list(base_folder))), total=len(file_list(base_folder))):
            if res:
                c += len(res)
                f.writelines(res)
    add_count_to_file(result_file, c)
    end = time.time()
    print(f'Done, {c:_} results, {end - start:.1f} seconds')
