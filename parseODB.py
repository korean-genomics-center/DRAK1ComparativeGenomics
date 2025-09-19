# %%
import glob
import os

import orthodb

# %%
api = orthodb.OdbAPI()
print('orthologs...')
t_orthologs =  api.orthologs("",species="9606_0")



# %%
