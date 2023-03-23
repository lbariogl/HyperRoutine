
import os
import yaml

run_period = "LHC_22r"
hyperloop_path = '/alice/cern.ch/user/a/alihyperloop/outputs/70921/7453/Stage_1/'
files_to_download = ['AO2D.root', 'AnalysisResults.root']

target_dir = '/home/fmazzasc/train_out'

download = True
merge = True


if download:

    if not os.path.exists(f'{target_dir}/{run_period}'):
            os.makedirs(f'{target_dir}/{run_period}')
    else:
        os.system(f'rm -rf {target_dir}/{run_period}/*')

    os.system(f'alien.py ls {hyperloop_path} > {target_dir}/{run_period}/out.txt')

    with open(f'{target_dir}/{run_period}/out.txt') as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    for index,line in enumerate(lines):
        input_path = hyperloop_path + line
        if not os.path.exists(f'{target_dir}/{run_period}/{line}'):
                os.makedirs(f'{target_dir}/{run_period}/{line}')

        for file in files_to_download:
            os.system(f'alien.py cp {input_path}/{file} file:{target_dir}/{run_period}/{line}/.')

if merge:
    print(f'Merging files of {run_period}')
    for file in files_to_download:
        ## add run period before .root extension
        file_out = file.split('.')[0] + '_' + run_period + '.root'
        if(file == 'AO2D.root'):
            file_out = file.split('.')[0] + '_' + run_period + '_temp.root'
        os.system(f'hadd -f {target_dir}/{run_period}/{file_out} {target_dir}/{run_period}/*/{file}')
        ## remove directory structure
        if(file == 'AO2D.root'):
            file_out_merged = file.split('.')[0] + '_' + run_period + '.root'
            os.system(f"""root -l -b -q MergeTrees.cc'("{target_dir}/{run_period}/{file_out}", "{target_dir}/{run_period}/{file_out_merged}")'""")
            os.system(f'rm {target_dir}/{run_period}/{file_out}')

    # clean directories, list all of them and check if they are a directory
    for dir in os.listdir(f'{target_dir}/{run_period}'):
        if os.path.isdir(f'{target_dir}/{run_period}/{dir}'):
            os.system(f'rm -rf {target_dir}/{run_period}/{dir}')
