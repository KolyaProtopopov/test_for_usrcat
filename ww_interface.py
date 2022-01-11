import ipywidgets as widgets


# CONFGEN MENU

ref_smi = widgets.Text(
    value='',
    placeholder='Type path to your molecules .smi',
    description='Path:',
    disabled=False
)

refSDF_base_path = widgets.Text(
    value='',
    placeholder='Type path to output file .sdf',
    description='Path:',
    disabled=False
)

energy_cut_off = widgets.FloatSlider(
    value=7.5,
    min=0,
    max=500.0,
    step=0.1,
    description= 'Energy cut of:',
    disabled=False,
    continuous_update=False,
    orientation= 'horizontal',
    readout=True,
    readout_format='.1f',
)

num_conf = widgets.IntSlider(
    value=42,
    min=0,
    max=300,
    step=1,
    description='Number of generated conformers:',
    disabled=False,
    continuous_update=False,
    orientation='horizontal',
    readout=True,
    readout_format='d'
)

confgen_menu = widgets.VBox([ref_smi, refSDF_base_path, energy_cut_off, num_conf ])

# SEARCH MENU

minimum_similarity = widgets.FloatSlider(
    value= 0.30,
    min=0,
    max=1.00,
    step=0.01,
    description= 'Minimum similariry value:',
    disabled=False,
    continuous_update=False,
    orientation= 'horizontal',
    readout=True,
    readout_format='.1f',
)

usrcat_base_path = widgets.Text(
    value='./TEST',
    placeholder='Type path to USRCAT fingerprints base folder',
    description='Path to USRCAT FP:',
    disabled=False
)

results_file_path = widgets.Text(
    value='',
    placeholder='Type path to results file',
    description='Path to USRCAT FP:',
    disabled=False
)

serch_menu = widgets.VBox([minimum_similarity, refSDF_base_path, usrcat_base_path, results_file_path ])

ww_menu = widgets.HBox([confgen_menu, serch_menu])

# BASE PREPARATION MENU

Base_smi = widgets.Text(
    value='./test.smi',
    placeholder='Type path to your base molecules .smi',
    description='Path:',
    disabled=False
)

Failed_smi = widgets.Text(
    value='./FAILED.smi',
    placeholder='Type path to your failed molecules .smi',
    description='Path:',
    disabled=False
)

select_CS_base = widgets.Dropdown(
    options=[('None', 1), ('Screening Compound Catalog', 2), ('GPCRs Targeted Library', 3), ('Protein Kinases Targeted Libraries', 4),
    ('Ion Channels Targeted Library', 5), ('ChEMBL analogs', 6), ('Lead-Like Diversity Library', 7)],
    value=2,
    description='CS Library:',
    disabled=False,
)

#Select your library
#GPCRs Library

#Screening Compound Catalog
#In-Stock Lead-Like compounds
#ChEMBL analogs


lable_path = widgets.Label(
    value='Select path to USRCAT fingerprints base:'
)

label_select_base = widgets.Label(
    value='Select library from Chemspace:'
)

#label_chose_base = widgets.Label(
#    value='Or chose your own library:'
#)

lable_smi_to_create = widgets.Label(
    value='Select *.SMI file for USRCAT fingerprints creation:'
)

Base_menu = widgets.VBox([lable_path, usrcat_base_path, label_select_base, select_CS_base, lable_smi_to_create, Base_smi, Failed_smi])

# CONFIG MENU

tabs_box = [confgen_menu, serch_menu, Base_menu]
tab_menu = widgets.Tab()
tab_menu.children = tabs_box

tab_menu.set_title(0, 'Confgen')
tab_menu.set_title(1, 'Searching properties')
tab_menu.set_title(2, 'Base preporation')
