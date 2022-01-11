def make_import():
    #OS libs
    import click
    import os

    #RDKIT Libs
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdMolAlign
    from rdkit.Chem.SaltRemover import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize

    #CS_scripts
    from ww_confgen import ww_confgen
    from ww_Base import gen_for_worker, human_format, add_count_to_file
    from ww_Base import worker as wb
    from ww_Base import do_base
    from ww_search import initializerSD, list_to_line, file_list
    from ww_search import worker as ws
    from ww_interface import ref_smi, refSDF_base_path, energy_cut_off, num_conf, confgen_menu
    from ww_interface import minimum_similarity, refSDF_base_path, usrcat_base_path, results_file_path, serch_menu
    from ww_interface import Failed_smi, Base_menu, Base_smi
    from ww_interface import ww_menu, tabs_box, tab_menu

    #widgets
    import ipywidgets as widgets

    return click, os, Chem, rdDistGeom,  AllChem, rdMolAlign, SaltRemover, rdMolStandardize, ww_confgen, gen_for_worker, human_format, add_count_to_file, wb, do_base, initializerSD, list_to_line, file_list, ws, ref_smi, refSDF_base_path, energy_cut_off, num_conf, confgen_menu, minimum_similarity, refSDF_base_path, usrcat_base_path, results_file_path, serch_menu, Failed_smi, Base_menu, Base_smi, ww_menu, tabs_box, tab_menu, widgets

#click, os, Chem, rdDistGeom,  AllChem, rdMolAlign, SaltRemover, rdMolStandardize, ww_confgen, gen_for_worker, human_format, add_count_to_file, wb, do_base, initializerSD, list_to_line, file_list, worker, ref_smi, refSDF_base_path, energy_cut_off, num_conf, confgen_menu, minimum_similarity, refSDF_base_path, usrcat_base_path, results_file_path, serch_menu, Failed_smi, Base_menu, Base_smi, ww_menu, tabs_box, tab_menu, widgets
