python3 analyse_tree.py --config-file "configs/analysis/config_antimatter.yaml" --dump-out-tree
python3 QAplots.py --input-file ../results/HypertritonResults_antimatter.root --output-dir ../results/plots/antimatter
python3 analyse_tree.py --config-file "configs/analysis/config_matter.yaml" --dump-out-tree
python3 QAplots.py --input-file ../results/HypertritonResults_matter.root --output-dir ../results/plots/matter
rm ../results/HypertritonResults_tot.root
hadd ../results/HypertritonResults_tot.root ../results/HypertritonResults_matter.root ../results/HypertritonResults_antimatter.root
python3 QAplots.py --input-file ../results/HypertritonResults_tot.root --output-dir ../results/plots/tot
