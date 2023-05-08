GENERAL_FILE="../results/AnalysisResults.root"
if [ -f "$GENERAL_FILE" ]; then
  echo "Removing $GENERAL_FILE"
  rm $GENERAL_FILE
fi
hadd $GENERAL_FILE /data/lbariogl/hyp_mc/round_*/AnalysisResults.root
python3 analyse_tree.py --config-file "configs/config_antimatter_inj_gap.yaml"
python3 analyse_tree.py --config-file "configs/config_matter_inj_gap.yaml"
NEW_TOT_FILE="../results/HypertritonResults_inj_gap_tot.root"
if [ -f "$NEW_TOT_FILE" ]; then
  echo "Removing $NEW_TOT_FILE"
  rm $NEW_TOT_FILE
fi
hadd $NEW_TOT_FILE ../results/HypertritonResults_inj_gap_matter.root ../results/HypertritonResults_inj_gap_antimatter.root
