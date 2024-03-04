# HyperRoutine
Analysis routine for light hypernuclei

### Dependencies
- ROOT > 6.26
- hipe4ml, possibly in dev mode. To install the in-development package:

```bash
git clone https://github.com/hipe4ml/hipe4ml.git
```

Then, from the repository base directory
```bash
pip install -e .[dev]
```



### How to run the analyses
- Inspect the output tree and produce basic histograms:
```bash
python3 analyse_tree.py --config-file config/analyse_tree/your_config.yaml
```
- Extract the raw and corrected pt spectrum:
```bash
python3 pt_analysis.py --config-file config/pt_analysis/your_config.yaml
```

- Extract the raw and corrected ct spectrum:
```bash
python3 ct_analysis.py --config-file config/ct_analysis/your_config.yaml
```
