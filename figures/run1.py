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

# %%
df=pd.read_csv('timing_run1',sep='\t',header=None, names=['reads','software','time'])
ftr = [3600,60,1]
df['time']=['0'+x if x[1] else x for x in df['time']]
df['time']=['00:'+x if x.count(':')==1 else x for x in df['time']]
df['time']=[sum([a*b for a,b in zip(ftr, map(float,TT.split(':')))]) for TT in df['time']]

# %%
df['reads'].apply(pd.to_numeric)

# %%
#some_array=['fastqc','01_fastp']
#df.loc[df['software'].isin(some_array)]

# %%
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.pointplot(x="reads", y="time", hue="software",data=df)
ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",
                  hue_order = ['01_prinseq++','prinseq_lite','01_fastp','01_NE_fastp','fastqc','01_NE_prinseq++',
                              '01_trimmomatic'],ci=99, data=df)
ax.set_ylabel('Time (s)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/60 for x in ax.get_ybound()])
ax2.set_ylabel('time (minutes)')
ax2.grid(None)
fig.tight_layout()

plt.show()
fig.savefig('time_comp_all')

# %%
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.pointplot(x="reads", y="time", hue="software",data=df)
ax = sns.lineplot(x="reads", y="time", hue="software",err_style="bars",
                  hue_order = ['01_fastp','01_NE_fastp','fastqc','01_NE_prinseq++',
                              '01_trimmomatic'],ci=99, data=df)
ax.set_ylabel('Time (s)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/60 for x in ax.get_ybound()])
ax2.set_ylabel('time (minutes)')
ax2.grid(None)
fig.tight_layout()
plt.show()
fig.savefig('time_comp_few')

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
plt.show()

# %%
df_ram=pd.read_csv('ram_run1',sep='\t',header=None, names=['reads','software','ram'])

# %%
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="darkgrid")
fig, ax = plt.subplots()
fig.set_size_inches(8, 6)
#ax = sns.pointplot(x="reads", y="time", hue="software",data=df)
ax = sns.lineplot(x="reads", y="ram", hue="software",err_style="bars",
                  hue_order = ['01_prinseq++','prinseq_lite','01_fastp','01_NE_fastp','fastqc','01_NE_prinseq++',
                              '01_trimmomatic'],ci=99, data=df_ram)
ax.set_ylabel('Memory (Kbytes)')    
ax.set_xlabel('Pair-end Sequences (millions)')
ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax2.set(ylim=[x/1024 for x in ax.get_ybound()])
ax2.set_ylabel('Memory (Mbytes)')
ax2.grid(None)
fig.tight_layout()
plt.show()
fig.savefig('RAM_comp_all')

# %%
#df_ram.loc{df_ram['software']=='01_prinseq++'])
#df_ram['software']==''
