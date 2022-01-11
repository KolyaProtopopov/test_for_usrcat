import click
import os
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

#@click.command()
#@click.option('--input', '-i', help='input file SMI or SDF', required=True)
#@click.option('--output', '-o', help='output file SDF', default='gen_confs.sdf')
#@click.option('--numconf', required=False, type=int)
#@click.option('--etrunc', default=500, required=False, type=float)
#@click.option('--add_ref', '-r', default=False, type=bool)
def ww_confgen(input, output, numconf, etrunc, add_ref):
    param = AllChem.ETKDGv3()
    param.pruneRmsThresh = 0.5
    param.numThreads=0
    param.maxIters=8192
    param.onlyHeavyAtomsForRMS=True
    param.useSmallRingTorsions = True
    param.useRandomCoords = True

    remover = SaltRemover()
    uncharger = rdMolStandardize.Uncharger()



    name, ext = os.path.splitext(os.path.basename(input))
    mols = []
    if ext.lower() in ['.smi']:
        print("SMILES are being processed")
#        mols = Chem.SmilesMolSupplier(input, " \t", 0, 1, titleLine=True, sanitize=True)
        with open(input, 'rb') as f:
            lines = f.readlines()
        for line in lines:
            line = line.decode().rstrip('\n')
            id = line.split('\t', 1)[1]
            smi = line.split('\t', 1)[0]
            try:
                mol = Chem.MolFromSmiles(smi)
                mol.SetProp('_Name', str(id))
                mols.append(mol)
            except:
                print("Error while reading", smi, "with id:", id)
    elif ext.lower() in ['.sdf', 'sd']:
        print("SDF are being processed")
        mols = Chem.SDMolSupplier(input, removeHs=True)
    else:
        print("Uknown extension")
        exit()

    for mol in mols:
        mol = remover.StripMol(mol)
        mol = uncharger.uncharge(mol)

        nr = int(AllChem.CalcNumRotatableBonds(mol))
        if nr <= 3:
            nc = 50
        elif nr > 6:
            nc = 300
        else:
            nc = nr**3

        if numconf:
            nc = numconf

        label = mol.GetProp("_Name")
        print("applying NC", nc, "for ", label, "conformers")

        mol = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(mol, numConfs = 300, params = param)
        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, maxIters=8192, mmffVariant='MMFF94s')

        sdfile = open(output, 'at')
        w = Chem.SDWriter(sdfile)
        res = []
        res_trunc = []
        ground_e = None
        for cid in cids:
            ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
            e = ff.CalcEnergy()
            res.append((cid, e))
        sorted_res = sorted(res, key=lambda x:x[1])
        print("NUMBER of sorted_res conformers", len(sorted_res))
        for num, energy in sorted_res:
            if ground_e == None:
                ground_e = energy
                print("Ground state E:", ground_e)
            if (energy - ground_e) <= etrunc:
                res_trunc.append((num, energy))
        sorted_res = res_trunc
#        print("NUMBER of truncated conformers", len(sorted_res))
        sorted_res = sorted_res[0:min(len(sorted_res), nc)]
#        print("NUMBER of stored conformers", len(sorted_res))
#        print("TOP10 lowest conformers", sorted_res[0:9])
#        print("TOP10 highest conformers", sorted_res[-10:])
#        print("NUMBER of requested conformers (NC) ", nc)
        rdMolAlign.AlignMolConformers(mol)
        for cid, e in sorted_res:
            mol.SetProp('_Name', str(label))
            mol.SetProp('idnumber', str(label))
            mol.SetProp('CID', str(cid))
            mol.SetProp('Energy', str(e))
            w.write(mol, confId=cid)
        w.close()

if __name__=='__main__':
    ww_confgen()
