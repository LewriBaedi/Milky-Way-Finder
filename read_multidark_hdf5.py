import os
import pandas as pd
import numpy as np
import h5py as h5

os.chdir(os.path.dirname(__file__))
sim = input('Please enter the path to the hdf dataset:')

desired_columns = ['X','Y','Z','X_Velocity','Y_Velocity','Z_Velocity','Mvir']
with h5.File(sim,'r') as f:
    d = {}
    for i in list(f.keys()):
        if i in desired_columns:
            d[i] = np.asarray(f[i])
    df = pd.DataFrame(data=d)

df.to_pickle('data_mdpl_select.pkl')