
import os
import yaml

import argparse

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('download_or_merge',
                    help="Choose between `download` and `merge`.")
parser.set_defaults(download_or_merge='')
parser.add_argument('--hyperloop', dest='hyperloop_path',
                    help="Path to the files on the hyperloop.")
parser.set_defaults(hyperloop_path='')
parser.add_argument('--local', dest='local_path',
                    help="Path to the files in local.")
parser.set_defaults(local_path='')
parser.add_argument('--period', dest='run_period',
                    help="Run period.")
parser.set_defaults(run_period='')
parser.add_argument('--config-file', dest='config_file',
                    help="path to the YAML file with configuration.")
parser.set_defaults(config_file='')
args = parser.parse_args()

download_or_merge = args.download_or_merge
if download_or_merge != 'download' and download_or_merge != 'merge':
    raise ValueError('Choose between `download` and `merge`. Tertium non datur.')

hyperloop_path = args.hyperloop_path
local_path = args.local_path
run_period = args.run_period

if args.config_file != "":
    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    if download_or_merge == 'download':
        hyperloop_path = config['hyperloop_path']
    local_path = config['local_path']
    run_period = config['run_period']

files_to_download = ['AO2D.root', 'AnalysisResults.root']

if download_or_merge == 'download':

    if(local_path == ''):
        raise ValueError('Parameter `local_path` not set.')
    if(run_period == ''):
        raise ValueError('Parameter `run_period` not set.')
    if(hyperloop_path == ''):
        raise ValueError('Parameter `hyperloop_path` not set.')

    if not os.path.exists(f'{local_path}/{run_period}'):
        os.makedirs(f'{local_path}/{run_period}')
    else:
        os.system(f'rm -rf {local_path}/{run_period}/*')

    os.system(
        f'alien.py ls {hyperloop_path} > {local_path}/{run_period}/out.txt')

    with open(f'{local_path}/{run_period}/out.txt') as file:
        lines = file.readlines()
        lines = [line.rstrip() for line in lines]

    for index, line in enumerate(lines):
        input_path = hyperloop_path + line
        if not os.path.exists(f'{local_path}/{run_period}/{line}'):
            os.makedirs(f'{local_path}/{run_period}/{line}')

        for file in files_to_download:
            os.system(
                f'alien.py cp {input_path}/{file} file:{local_path}/{run_period}/{line}/.')

else:

    if(local_path == ''):
        raise ValueError('Parameter `local_path` not set.')
    if(run_period == ''):
        raise ValueError('Parameter `run_period` not set.')

    print(f'Merging files of {run_period}')

    for file in files_to_download:
        # add run period before .root extension
        file_out = file.split('.')[0] + '_' + run_period + '.root'
        if (file == 'AO2D.root'):
            file_out = file.split('.')[0] + '_' + run_period + '_temp.root'
        os.system(
            f'hadd -f {local_path}/{run_period}/{file_out} {local_path}/{run_period}/*/{file}')
        # remove directory structure
        if (file == 'AO2D.root'):
            file_out_merged = file.split('.')[0] + '_merged_' + run_period + '.root'
            os.system(
                f"""root -l -b -q merge/MergeTrees.cc++'("{local_path}/{run_period}/{file_out}", "{local_path}/{run_period}/{file_out_merged}")'""")
    # clean directories, list all of them and check if they are a directory
    for dir in os.listdir(f'{local_path}/{run_period}'):
        if os.path.isdir(f'{local_path}/{run_period}/{dir}'):
            os.system(f'rm -rf {local_path}/{run_period}/{dir}')
