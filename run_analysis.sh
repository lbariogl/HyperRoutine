
## integrated analysis (for performance studies)
python3 analyse_tree.py --config-file configs/analysis/config_antimatter.yaml
python3 analyse_tree.py --config-file configs/analysis/config_mc.yaml
python3 signal_extraction.py --config-file configs/signal_extraction/config_signal_extraction_antimat.yaml

## pt analysis
python3 pt_analysis.py --config-file configs/analysis/config_pt_bins.yaml

## ct analysis
python3 ct_analysis.py --config-file configs/analysis/config_ct_bins.yaml

