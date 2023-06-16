python3 analyse_tree.py --config-file "configs/config_antimatter.yaml"
python3 QAplots.py --input-file ../results/HypertritonResults_antimatter_pass4.root --output-dir ../results/plots/matter
python3 analyse_tree.py --config-file "configs/config_matter.yaml"
python3 QAplots.py --input-file ../results/HypertritonResults_matter_pass4.root
rm ../results/HypertritonResults_tot_pass4.root
hadd ../results/HypertritonResults_tot_pass4.root ../results/HypertritonResults_matter_pass4.root ../results/HypertritonResults_antimatter_pass4.root
python3 QAplots.py --input-file ../results/HypertritonResults_tot_pass4.root --output-dir ../results/plots/tot
