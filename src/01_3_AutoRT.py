import glob
import pandas as pd
from pyteomics import achrom
import numpy
rcond = None

files_train = sorted(glob.glob("../results/train/AutoRT/*.tsv"))[0:2]
files_test = sorted(glob.glob("../results/test/AutoRT/*.tsv"))[0:2]

for i, data in enumerate(files_train):
  print(files_train[i])
  
  ## Training
  python autort.py train -i files_train[i] -o results/train/AutoRT/ -e 40 -b 64 -u s -m ../models/base_models_PXD006109/model.json --add_ReduceLROnPlateau --early_stop_patience 10
  
  ## Prediction
  #python autort.py predict -t files_test[i] -s tf_model/model.json -o tf_prediction/ -p test
