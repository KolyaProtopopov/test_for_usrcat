from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
from rdkit.Chem import DataStructs
import os
import sys
import time
import shutil
import multiprocessing as mp
from tqdm import tqdm
import gzip
from rdkit import RDLogger
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT, GetUSR
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from urllib.request import urlopen
from tempfile import NamedTemporaryFile
from shutil import unpack_archive

param_gen = AllChem.ETKDGv3()
param_gen.pruneRmsThresh=0.5
param_gen.numThreads=0
param_gen.maxIters=8192
param_gen.onlyHeavyAtomsForRMS=True
param_gen.useSmallRingTorsions = True
param_gen.useRandomCoords = True

remover = SaltRemover()
uncharger = rdMolStandardize.Uncharger()



db_index = sys.argv[1]
numconf = 300
etrunc = 500

base_index = 'CS_stock_'
file_in = 'xray-active-ligands.smi'
folder_out = ''
file_failed = '3D_failed.smi'

logfile = open(base_index + ".log", 'at')
print("Start LOGGING:", file = logfile)
logfile.close()


chunk_size = 10*1024

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)


def gen_for_worker(file1, file2):
    n = 0
    with open(file1, 'rb', buffering=chunk_size) as f:
#        next(f)
        while True:
            lines = f.readlines(chunk_size)
            if not lines:
                break
            n += 1
            yield lines, n, file2


def worker(line_n):
    lines, n, folder_out = line_n
    logging = list()
    failed = list()

    file_out = os.path.join(folder_out, base_index + "_" + f'{n:08d}.bin.gz')
    print(file_out)
    print("FILE_OUT", file_out)
    size = 0
    with gzip.open(file_out, 'wb') as f:
        for line in lines:
            try:
                size += len(line)
                line = line.decode().rstrip('\n')
                smiles = line.split('\t', 1)[0]
                id = line.split('\t', 1)[1]
                moll = Chem.MolFromSmiles(smiles)
                moll.SetProp('_Name', str(id))
                moll = remover.StripMol(moll)
                moll = uncharger.uncharge(moll)

                nr = int(AllChem.CalcNumRotatableBonds(moll))
                if nr <= 3:
                    nc = 50
                elif nr > 6:
                    nc = 300
                else:
                    nc = nr**3

                if numconf:
                    nc = numconf

                moll = Chem.AddHs(moll)
                cids_moll = AllChem.EmbedMultipleConfs(moll, numConfs = 300, params = param_gen)
                mp_moll = AllChem.MMFFGetMoleculeProperties(moll, mmffVariant = 'MMFF94s')
                AllChem.MMFFOptimizeMoleculeConfs(moll, maxIters = 8192, numThreads = 0, mmffVariant = 'MMFF94s')
                logging.append('Mol ID: ' + str(id) + ' and SMILES: ' + str(smiles) + '\n')
                logging.append('Initial nc: ' + str(nc) + ' and the pruned nc: ' + str(len(cids_moll)) + '\n')

                res = []
                res_trunc = []
                ground_e = None
                for cid in cids_moll:
                    ff = AllChem.MMFFGetMoleculeForceField(moll, mp_moll, confId=cid)
                    e = ff.CalcEnergy()
                    res.append((cid, e))
                sorted_res = sorted(res, key=lambda x:x[1])
                for num, energy in sorted_res:
                    if ground_e == None:
                        ground_e = energy
                        logging.append('Ground state E:' + str(ground_e) + '\n')
                    if (energy - ground_e) <= etrunc:
                        res_trunc.append((num, energy))
                sorted_res = res_trunc
                logging.append('NUMBER of truncated conformers ' + str(len(sorted_res)) + '\n')
                sorted_res = sorted_res[0:min(len(sorted_res), nc)]
                logging.append('NUMBER of stored conformers ' + str(len(sorted_res)) + '\n')
                logging.append('TOP10 lowest conformers ' + str(sorted_res[0:9]) + '\n')
                logging.append('TOP10 highest conformers ' + str(sorted_res[-10:]) + '\n')

                moll = Chem.RemoveHs(moll)
                pickle.dump((line, tuple(GetUSRCAT(moll, confId=num) for num, energy in sorted_res)), f)
            except:
                failed.append(line + '\n')
#               failed.append(smiles + '\t' + id + '\n')
                print("Failed at: ", line)
    return size, logging, failed


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

def download_sets(web_path, path_to_base):
    if web_path == 2:
        zipurl =''
    elif web_path == 3:
        zipurl =''
    elif web_path == 4:
        zipurl =''
    elif web_path == 5:
        zipurl =''
    elif web_path == 6:
        zipurl =''
    elif web_path == 7:
        zipurl =''
    with urlopen(zipurl) as zipresp, NamedTemporaryFile() as tfile:
        tfile.write(zipresp.read())
        tfile.seek(0)
        unpack_archive(tfile.name, path_to_base, format = 'zip')


def do_base(files_in, base_folder, fail_file):
    global file_in
    global folder_out
    global file_failed

    file_in = files_in
    folder_out = base_folder
    file_failed = fail_file



    numconf = 300
    etrunc = 500

    start = time.time()
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    err_counter = 0
    res_counter = 0
    with mp.Pool(processes=1) as pool, open(file_failed, 'at') as f, open(base_index + ".log", 'at') as l, \
        tqdm(unit_scale=True, mininterval=1.0, total=os.path.getsize(file_in), unit='B', unit_divisor=1024) as t:
        for res in pool.imap(worker, gen_for_worker(file_in, folder_out)):
            t.update(res[0])
            l.writelines(res[1])
            res_counter += res[0]
            err_counter += len(res[2])
            f.writelines(res[2])
    add_count_to_file(file_failed, err_counter)
    end = time.time()
    total_time = (end - start) / 3600

    logfile = open(base_index + ".log", 'at')
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours', file = logfile)
    logfile.close()
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours')

if __name__ == '__main__':

    start = time.time()
    if not os.path.exists(folder_out):
        os.makedirs(folder_out)
    err_counter = 0
    res_counter = 0
    with mp.Pool(processes=1) as pool, open(file_failed, 'at') as f, open(base_index + ".log", 'at') as l, \
        tqdm(unit_scale=True, mininterval=1.0, total=os.path.getsize(file_in), unit='B', unit_divisor=1024) as t:
        for res in pool.imap(worker, gen_for_worker(file_in)):
            t.update(res[0])
            l.writelines(res[1])
            res_counter += res[0]
            err_counter += len(res[2])
            f.writelines(res[2])
    add_count_to_file(file_failed, err_counter)
    end = time.time()
    total_time = (end - start) / 3600

    logfile = open(base_index + ".log", 'at')
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours', file = logfile)
    logfile.close()
    print(f'Done. {res_counter:_} molecules processed and {err_counter:_} molecules failed. Finished in {total_time:.1f} hours')
