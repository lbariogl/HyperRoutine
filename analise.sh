python3 analyse_tree.py --config-file "config_antimatter.yaml"
python3 analyse_tree.py --config-file "config_matter.yaml"
rm ../results/HypertritonResults_tot.root
hadd ../results/HypertritonResults_tot.root ../results/HypertritonResults_matter.root ../results/HypertritonResults_antimatter.root
