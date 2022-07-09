#!/usr/bin/env python
# coding: utf-8

# # Run notebooks

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import os
from glob import glob
import papermill as pm
import util


# In[2]:


cluster, client = util.get_ClusterClient()
cluster.scale(5)
client


# In[ ]:


output_dir = '_computed-notebooks'
get_ipython().system(' mkdir _computed-notebooks')
from tqdm import tqdm 
for key in tqdm(['SSH_2','SST']):
    for nb in tqdm(['tracer_snapshot']):
        input_path = f'{nb}.ipynb'
        output_path = (f'{output_dir}/{nb}-{key}.ipynb')
        nb_api = pm.inspect_notebook(input_path)
        print(f'executing {input_path}')
        o = pm.execute_notebook(
            input_path=input_path,
            output_path=output_path,
            kernel_name='TPAC_O2',
            parameters=dict(variable_id=key)
        )


# In[18]:


cluster.close()

