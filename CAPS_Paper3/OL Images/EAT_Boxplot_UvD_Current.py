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


colors = ['red', 'green']
 
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

ax.legend(bp["boxes"], ["Uniform", "Dynamic"], loc='upper right')

fig1 = plt.gcf()

# show plot
# plt.plot()
#%%
fig1.savefig(path+'EAT-8pct_BP.pdf',dpi=100)  
fig1.show()