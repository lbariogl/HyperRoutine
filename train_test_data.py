import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler
import hipe4ml.analysis_utils as au
import hipe4ml.plot_utils as pu

import mplhep as mpl
import xgboost as xgb

import yaml

matplotlib.use('pdf')
plt.style.use(mpl.style.ALICE)


###############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.", default='')
args = parser.parse_args()

config_file = args.config_file

with open(config_file, 'r') as stream:
    try:
        params = yaml.full_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


do_training = params['do_training']
do_application = params['do_application']

input_data_path = params['input_data_path']
input_mc_path = params['input_mc_path']
output_dir = params['output_dir']

training_preselections = params['training_preselections']
training_variables = params['training_variables']
test_set_size = params['test_set_size']
bkg_fraction = params["background_over_signal"]
random_state = params["random_state"]
hyperparams = params["hyperparams"]


### create output directory if it does not exist
if not os.path.exists(output_dir):
        os.makedirs(output_dir)

##  ml_plots_dir
figures_ml_path = output_dir + "/figures_ML"
if not os.path.exists(figures_ml_path):
        os.makedirs(figures_ml_path)



print('**********************************')
print('    Running train_test_data.py    ')
print('**********************************')

if do_training:

        signalH = TreeHandler(input_mc_path, "O2mchypcands")
        bkgH = TreeHandler(input_data_path, "O2datahypcands")
        
        ## select background by taking the sidebands of the mass distribution
        if training_preselections != '':
                signalH.apply_preselections(training_preselections)
                bkgH.apply_preselections(f"(fMassH3L<2.95 or fMassH3L>3.02) and {training_preselections}")
        else:
                bkgH.apply_preselections(f"(fMassH3L<2.95 or fMassH3L>3.02)")

        if bkg_fraction!=None:
                bkgH.shuffle_data_frame(size=bkg_fraction*len(signalH), inplace=True, random_state=random_state)
        
        print("Signal events: ", len(signalH))
        print("Background events: ", len(bkgH))

        train_test_data = au.train_test_generator([signalH, bkgH], [1,0], test_size=test_set_size, random_state=random_state)

        ### create ML output directory if it does not exist
        figures_ml_path = output_dir + "/figures_ML"
        if not os.path.exists(figures_ml_path):
                os.makedirs(figures_ml_path)

        distr = pu.plot_distr([bkgH, signalH], training_variables, bins=63, labels=['Signal',"Background"],colors=["blue","red"], log=True, density=True, figsize=(18, 13), alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(figures_ml_path + "/features_distributions.png", bbox_inches='tight')
        corr = pu.plot_corr([signalH,bkgH], training_variables + ["fMassH3L"], ['Signal',"Background"])
        corr[0].savefig(figures_ml_path + "/correlations.png",bbox_inches='tight')

        print("---------------------------------------------")
        print("Data loaded. Training and testing ....")

        model_hdl = ModelHandler(xgb.XGBClassifier(), training_variables)
        model_hdl.set_model_params(hyperparams)
        y_pred_test = model_hdl.train_test_model(train_test_data, True, True)

        print("Model trained and tested. Saving results ...")

        bdt_out_plot = pu.plot_output_train_test(model_hdl, train_test_data, 100, True, ["Signal", "Background"], True, density=True)
        bdt_out_plot.savefig(figures_ml_path + "/bdt_output.png")
        feature_importance_plot = pu.plot_feature_imp(train_test_data[2], train_test_data[3], model_hdl)
        feature_importance_plot[0].savefig(figures_ml_path + "/feature_importance_1.png")
        feature_importance_plot[1].savefig(figures_ml_path + "/feature_importance_2.png")


        ## dump model handler and efficiencies vs score
        model_hdl.dump_model_handler(output_dir + "/model_hndl.pkl")
        eff_arr = np.round(np.arange(0.5,0.99,0.01),2)
        score_eff_arr = au.score_from_efficiency_array(train_test_data[3], y_pred_test, eff_arr)
        np.save(output_dir + "/efficiency_arr.npy", eff_arr)
        np.save(output_dir + "/score_efficiency_arr.npy",score_eff_arr)

        print("Training done")
        print("---------------------------------------------")
        del signalH, bkgH


if do_application:

        print("---------------------------------------------")
        print("Starting application: ..")

        dataH = TreeHandler(input_data_path, "O2datahypcands")
        if training_preselections != '':
                dataH.apply_preselections(training_preselections)

        bdt_eff_arr = np.load(output_dir + "/efficiency_arr.npy")
        score_eff_arr = np.load(output_dir + "/score_efficiency_arr.npy")

        model_hdl = ModelHandler()
        model_hdl.load_model_handler(output_dir + "/model_hndl.pkl")

        dataH.apply_model_handler(model_hdl, column_name="model_output")

        print("Application done. Saving results ...")
        dataH.write_df_to_parquet_files("dataH.parquet", output_dir)






