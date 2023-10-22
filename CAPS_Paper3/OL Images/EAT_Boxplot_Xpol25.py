# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%%

# Import libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
 
#%%
path="\\coe-fs.engr.tamu.edu\Grads\sohom070994\Documents\GitHub\\"
Dyn_8_4=np.array([4.3900923408717,3.77488001533228,3.65683655951233,3.560663810963,3.93497745724651,3.78715171367794,3.71638013974579,3.65372693982899,3.52223458104721,3.72646739497319,3.51914314339282,4.00882569538859,3.41755270563801,3.54478357256661,3.69950712960179,3.83018352035592,4.23478155733677,3.67909747006514,3.74248339674853,3.21112951537752,3.7995996194742,3.84928842512908,3.5253108877206,4.13694345836506,3.79583224106091,4.39071038243269,3.73860166730726,3.12237008304566,4.04502890722104,3.30624261906793,4.22221098124165,3.95306589263621])
Uni_8_4=np.array([6.4875435216886,5.55963653364009,5.5378669573472,5.18787773827957,5.92522477765109,5.15471471909679,4.9269437837455,4.7676154910001,5.30130605898911,5.04381370739986,5.38239469599147,5.76615498703939,5.38133233216579,4.09790256645617,4.7189511666607,5.61862442069422,6.22082946610881,5.15970783036801,4.70214339723045,4.45847256045855,4.79253597494149,6.3519689150777,5.42062187961622,5.53536139501712,5.15812447813332,4.58188825010073,5.28018807585611,5.08168084332822,5.20952923252659,4.54009423017166,6.12067819983476,5.61950288807162])
#%%
data=[Uni_8_4,Dyn_8_4]

fig = plt.figure(figsize=(6,4))

#change fatness
xl=0  
xu=0.75

#change distance from ends
s=3.5

# Creating axes instance
ax = fig.add_axes([0.1, 0.1, .85, .85])#if stating two args 1 then axis will not show

# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
# Creating plot
bp = plt.boxplot(data,patch_artist=True,positions=[(xu-xl)/s,xu-(xu-xl)/s],showmeans=True,showfliers=False,
                 meanprops={"marker":"x","markeredgecolor":"black"},
                 medianprops=dict(linestyle='dashed', linewidth=1, color='black')
                 )


colors = ['#808080', '#D3D3D3']
 
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)


plt.xlim(xl,xu)
plt.grid(axis = 'y',alpha=0.3)

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) 

# plt.xlabel('categories')
plt.ylabel('Expected access time (days)')

ax.legend(bp["boxes"], ["25% division policy", "Optimal"], loc='upper right')

fig1 = plt.gcf()

# show plot
# plt.plot()
#%%
fig1.savefig('\\coe-fs.engr.tamu.edu\Grads\sohom070994\Documents\GitHub\\EAT-Xpol25.pdf',dpi=100)  
fig1.show()
# %%
