# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import numpy as np

# %%
df = pd.read_csv('./output/model_vectors.txt',header=None)

# %%
df

# %%
df.shape

# %%

# %%
sum_col = np.sum(df,axis=1)

# %%
sum_col.min()

# %%
np.argmin(sum_col)

# %%
df.iloc[1264125,:].values

# %%
pd.DataFrame(sum_col).value_counts()

# %% tags=[]
# 27    412159
# 26    353774
# 25    236589
# 24    120231
# 23     44626
# 22     11328
# 21      1740
# 20       120

# %%
np.where(sum_col == 27)[0][1]

# %%
df.iloc[1535,:].values

# %%
df22 = df.loc[sum_col ==22].copy()

# %%
df22.to_csv('./output/df22.csv')

# %%
multi_histo = pd.read_csv('./output/results_histo_rand_prior.txt',index_col=0,header=None)

# %%
multi_histo

# %%
np.sum(multi_histo[1] > 0.025)

# %%
vec = multi_histo[1] > 0.025

# %%
multi_histo.loc[vec,:]

# %%
vec.values

# %%
sum_col[vec.values].min()

# %%
vec_best = [1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,1,1,1,0]

# %%
sum(vec_best)

# %%
(df == vec_best).all(1).any()

# %%
df24.iloc[0,:].values

# %%
