import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None

files_train = sorted(glob.glob("results/train/achrom/*.csv"))
files_test = sorted(glob.glob("results/test/achrom/*.csv"))

for i, data in enumerate(files_train):
  print(files_train[i])
  
  df = pd.read_csv(files_train[i], sep = "\t")
  pred = pd.read_csv(files_test[i], sep = "\t")
  
  # Calibration
  # RCs = achrom.get_RCs_vary_lcp(sequences = df['Mascot_seq'].tolist(), RTs=df['RT'].tolist(), lcp_range=(-1.0,1.0), term_aa=False, rcond=None)
  RCs = achrom.get_RCs(sequences = df['Mascot_seq'].tolist(), RTs=df['RT'].tolist(), term_aa=False)
  
  # Prediction
  pep = pred['Mascot_seq'].tolist()

  RTpred = numpy.zeros(len(pep))
  for j in range(0,len(pep)):
      RTpred[j] = achrom.calculate_RT(pep[j], RCs, raise_no_mod=False)
  
  out = pd.DataFrame({'Mascot_seq': pep,'RT': RTpred})    
  out.to_csv(files_train[i].replace("train", "predict", 1))
