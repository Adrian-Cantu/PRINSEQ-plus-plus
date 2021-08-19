# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import pandas as pd

# %%
df=pd.read_csv('timing_run2',sep='\t',header=None, names=['reads','software','time'])
ftr = [3600,60,1]
df['time']=['0'+x if x[1] else x for x in df['time']]
df['time']=['00:'+x if x.count(':')==1 else x for x in df['time']]
df['time']=[sum([a*b for a,b in zip(ftr, map(float,TT.split(':')))]) for TT in df['time']]

# %%

#col3 =[]

# %%
import seaborn as sns
import matplotlib.pyplot as plt
unique = df["software"].append(df["software"]).unique()
palette = dict(zip(unique, sns.color_palette()))
palette.update({"02_NE_fastp":palette['01_NE_fastp']})
palette.update({"04_NE_fastp":palette['01_NE_fastp']})
palette.update({"02_NE_prinseq++":palette['01_NE_prinseq++']})
palette.update({"04_NE_prinseq++":palette['01_NE_prinseq++']})
palette.update({"02_trimmomatic":palette['01_trimmomatic']})
palette.update({"04_trimmomatic":palette['01_trimmomatic']})

# %%
import seaborn as sns
import matplotlib.pyplot as plt
unique = df["software"].append(df["software"]).unique()
palette = dict(zip(unique, sns.color_palette()))
palette.update({"01_NE_fastp":(1,0,0)})
palette.update({"02_NE_fastp":(1,0,0)})
palette.update({"04_NE_fastp":(1,0,0)})

palette.update({"01_NE_prinseq++":(0,1,0)})
palette.update({"02_NE_prinseq++":(0,1,0)})
palette.update({"04_NE_prinseq++":(0,1,0)})

palette.update({"01_trimmomatic":(0,0,1)})
palette.update({"02_trimmomatic":(0,0,1)})
palette.update({"04_trimmomatic":(0,0,1)})


# %%

sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.lineplot(x="reads", y="time", hue="software",data=df)
#ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",
#                  hue_order = ['01_prinseq++','prinseq_lite','01_fastp','01_NE_fastp','fastqc','01_NE_prinseq++',
#                              '01_trimmomatic'],ci=99, data=df)
ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",ci=99, data=df)



ax.set_ylabel('Time (s)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/60 for x in ax.get_ybound()])
ax2.set_ylabel('time (minutes)')
ax2.grid(None)
fig.tight_layout()
fig.savefig('MM_time_comp_all')

# %%
sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.lineplot(x="reads", y="time", hue="software",data=df)
#ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",
#                  hue_order = ['01_prinseq++','prinseq_lite','01_fastp','01_NE_fastp','fastqc','01_NE_prinseq++',
#                              '01_trimmomatic'],ci=99, data=df)
ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",ci=99, data=df,palette=palette)



ax.set_ylabel('Time (s)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/60 for x in ax.get_ybound()])
ax2.set_ylabel('time (minutes)')
ax2.grid(None)
fig.tight_layout()
plt.show()
fig.savefig('MM_time_comp_all_cc')

# %%
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.pointplot(x="reads", y="time", hue="software",data=df)
ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",
                  hue_order = ['prinseq_lite','01_NE_fastp','fastqc','01_NE_prinseq++',
                              '01_trimmomatic'],ci=99, data=df, palette=palette)
ax.set_ylabel('Time (s)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/60 for x in ax.get_ybound()])
ax2.set_ylabel('time (minutes)')
ax2.grid(None)
fig.tight_layout()
fig.savefig('time_comp_all')

# %%
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.pointplot(x="reads", y="time", hue="software",data=df)
ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",
                  hue_order = ['01_fastp','01_NE_fastp','01_NE_prinseq++',
                              '01_trimmomatic'],ci=99, data=df)
ax.set_ylabel('Time (s)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/60 for x in ax.get_ybound()])
ax2.set_ylabel('time (minutes)')
ax2.grid(None)
fig.tight_layout()
