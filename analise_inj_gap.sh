python3 analyse_tree.py --config-file "configs/config_antimatter_inj_gap.yaml"
python3 analyse_tree.py --config-file "configs/config_matter_inj_gap.yaml"
rm ../results/HypertritonResults_inj_gap_tot.root
hadd ../results/HypertritonResults_inj_gap_tot.root ../results/HypertritonResults_inj_gap_matter.root ../results/HypertritonResults_inj_gap_antimatter.root
