import os

pdg0 = 211
pdg1 = 1000020030
mother_pdg = 1010010030
debug = "false"
path = "/data/fmazzasc/its_data/sim/hyp_gap_trig_2/"
 
suffix = "def"

os.system(f"""root -l -b -q DauTreeBuilder.C+'({pdg0}, {pdg1}, {mother_pdg}, {debug}, "{path}", "{suffix}")' """)
os.system(f""" python3 afterburner_eff_study.py --dau_pdg {pdg0} --suffix "{suffix}" """)
os.system(f""" python3 afterburner_eff_study.py --dau_pdg {pdg1} --suffix "{suffix}" """)