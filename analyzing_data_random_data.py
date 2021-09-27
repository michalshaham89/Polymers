from random_walk import load_data
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import os
import re
import pickle
from matplotlib import rcParams

EE_RADII_DIR_PATH = "C:\Users\Michal\Documents\michuli\physics_2\cm_project\ee_radii"
D = range(2, 5)
N = pow(2, np.arange(6, 13))
p = [3.0 / 4, 3.0 / 5, 1.0 / 2]
pos=[[0,1],[2,0],[2,2]]
text_pos=[0.5,0.8]
#pos=[[0:6,0],[6:12,0],[3:9,1]]


def exponent(x, a, b):
    return pow(x, a) * b


fig = plt.figure(1)
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']
grid = plt.GridSpec(2, 2, hspace=0.4, wspace=0.4)
a=[grid[0,0], grid[0,1],grid[1,0]]


for d in D:
    #_N = N if d != 2 else N[:6]
    ee_radii = []
    std_err = []
    _N=[]

    for ee_radii_file in os.listdir(EE_RADII_DIR_PATH):
        nums=map(int, re.findall(r'\d+', ee_radii_file))
        if nums[0]==d:
            if d==2 and nums[1]>=1536:
                continue
            file_name=os.path.join(os.path.join(EE_RADII_DIR_PATH, ee_radii_file),'ee_radii_1')
            #saw = np.array(pickle.load(open(file_name, 'rb')))
            saw=load_data(nums[1],nums[0])
            ee_radii.append(saw.mean())
            tmp_err=saw.std(ddof=1) / np.sqrt(len(saw))
            if tmp_err<0.009*saw.mean() or tmp_err>0.011*saw.mean():
                print ('change errors', nums, tmp_err,0.01*saw.mean())
                std_err.append(0.01*saw.mean())
            else:
                std_err.append(saw.std(ddof=1) / np.sqrt(len(saw)))

            _N.append(nums[1])

    print (ee_radii, std_err)
    res = curve_fit(exponent, _N, ee_radii, p0=[p[d - 2], 1], sigma=std_err)
    std = np.sqrt(res[1][0][0])
    nu = res[0][0]
    print (d, nu, std)

    #ax = plt.subplot(a[d-2])
    ax = plt.subplot(220+d-1)
    tmp_n=np.array(sorted(_N))
    ax.yaxis.grid(True)  # horizontal lines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.plot(tmp_n, exponent(tmp_n, res[0][0], res[0][1]), 'b', linewidth=2)
    #plt.errorbar(_N, ee_radii, std_err, fmt='none', elinewidth=3, capsize=4, capthick=3, ecolor='m')
    ax.plot(_N, ee_radii, '.m', linewidth=10, markersize=22, markeredgewidth=2)
    ax.set_title("d = %d" % d, fontsize=38)
    ax.set_xlabel("N", fontsize=32)
    ax.set_ylabel("$<R_n>$", fontsize=32)
    ax.tick_params(labelsize=30)
    ax.text(max(_N)*text_pos[0],max(ee_radii)*text_pos[1],"$y=%.1fx^{%.2f}$" % (res[0][1],res[0][0]), color='b',
            horizontalalignment='right', verticalalignment='center', fontsize=38, weight='bold')


plt.subplots_adjust(left=0.05, bottom=0.05,
                right=0.95, top=0.95, wspace=0.2, hspace=0.3)
plt.show()
