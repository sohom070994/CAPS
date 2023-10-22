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
path="C:\\Users\\Sohom\\Desktop\\CAPS_Paper3\\OL Images\\Matplotlib_Images\\"
Dyn_8_4=np.array([4.3900923408717,3.77488001533228,3.65683655951233,3.560663810963,3.93497745724651,3.78715171367794,3.71638013974579,3.65372693982899,3.52223458104721,3.72646739497319,3.51914314339282,4.00882569538859,3.41755270563801,3.54478357256661,3.69950712960179,3.83018352035592,4.23478155733677,3.67909747006514,3.74248339674853,3.21112951537752,3.7995996194742,3.84928842512908,3.5253108877206,4.13694345836506,3.79583224106091,4.39071038243269,3.73860166730726,3.12237008304566,4.04502890722104,3.30624261906793,4.22221098124165,3.95306589263621])
Uni_8_4=np.array([10.4084843955387,9.34468010582045,9.58522839415061,9.86592668056372,8.12378833380826,8.00270252709804,7.25914497577021,7.22583256182129,8.38989287773111,8.95226766801557,8.54625508948945,7.09279700398968,7.57587299354791,7.08236937912519,7.59951761018208,10.0571159587179,8.56575206623095,7.95940841406109,9.10289959793726,8.05234666948941,9.73823536430406,10.3210706130871,7.937269045141,7.72885190858107,7.23685411652824,7.69867331243855,7.56599775257374,7.47079743808427,8.23294804174074,7.22915273832231,8.27247789265235,8.69905751352684])
Smart_21 = np.array([6.4875435216886,5.55963653364009,5.5378669573472,5.18787773827957,5.92522477765109,5.15471471909679,4.9269437837455,4.7676154910001,5.30130605898911,5.04381370739986,5.38239469599147,5.76615498703939,5.38133233216579,4.09790256645617,4.7189511666607,5.61862442069422,6.22082946610881,5.15970783036801,4.70214339723045,4.45847256045855,4.79253597494149,6.3519689150777,5.42062187961622,5.53536139501712,5.15812447813332,4.58188825010073,5.28018807585611,5.08168084332822,5.20952923252659,4.54009423017166,6.12067819983476,5.61950288807162])
Smart_31 = np.array([4.96649495285826,5.26375123319136,5.97258392607759,5.15225385830926,5.4945439051705,5.3380206828457,5.00169536107784,5.15414934777285,5.09990954747907,5.42101882242714,4.74864051620516,4.55767436486415,4.53172356348204,4.45493463809609,4.236716074228,6.21400963723209,5.83076558608209,5.14271140452254,4.80120688913196,4.92787266434161,4.32263054990519,5.41866037592087,4.98959151438791,4.77377625113532,5.10828898888937,5.01069488539622,4.46625915024491,5.34391507231251,4.41128183658691,4.69219621099882,6.02339489731034,5.38107303714497])
#%%
data=[Uni_8_4,Smart_21,Smart_31,Dyn_8_4]

fig = plt.figure(figsize=(6,4))

# #change fatness
xl=0  
xu=1

#change distance from ends
# s=3
# positions=[(xu-xl)/s,(((xu-xl)/s)+(xu-(xu-xl)/s))/2,xu-(xu-xl)/s]
# Creating axes instance
ax = fig.add_axes([0.1, 0.1, .85, .85])#if stating two args 1 then axis will not show

# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))
# Creating plot
bp = plt.boxplot(data,patch_artist=True,positions=[0.2,0.4,0.6,0.8],showmeans=True,showfliers=False,
                 meanprops={"marker":"x","markeredgecolor":"black"},
                 medianprops=dict(linestyle='dashed', linewidth=1, color='black')
                 )


colors = ['#d3d3d3','#8b8b8b', '#626262','#36454f']
 
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)

# bp['boxes'][1].set(hatch = '--')

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

ax.legend(bp["boxes"], ["Uniform","Mid-semester (2:1)", "Mid-semester (3:1)","Dynamic"], loc='upper right')

fig1 = plt.gcf()

# show plot
plt.plot()
#%%
import os
directory_of_python_script = os.path.dirname(os.path.abspath(__file__))
#%%
fig1.savefig(os.path.join(directory_of_python_script, "EAT-8pct_BP.pdf"),dpi=100)  
fig1.show()
# %%
