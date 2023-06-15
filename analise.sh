python3 analyse_tree.py --config-file "configs/config_antimatter.yaml"
python3 analyse_tree.py --config-file "configs/config_matter.yaml"
rm ../results/HypertritonResults_tot_pass4.root
hadd ../results/HypertritonResults_tot_pass4.root ../results/HypertritonResults_matter_pass4.root ../results/HypertritonResults_antimatter_pass4.root
