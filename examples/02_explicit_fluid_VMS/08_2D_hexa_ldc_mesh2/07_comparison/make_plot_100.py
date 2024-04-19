import numpy as np
import matplotlib.pyplot as plt

fs=12
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('text', usetex=True)

true_hor = np.loadtxt('re_100_horizontal.txt',skiprows=1)
true_ver = np.loadtxt('re_100_vertical.txt',skiprows=1)

sim_hor05_100 = np.loadtxt('hor05_100.csv',skiprows=1,delimiter=',')

sim_vert05_100 = np.loadtxt('vert05_100.csv',skiprows=1,delimiter=',')

# x=8,y=9
sv_x = 8
# vx=3,vy=4
sv_y = 4

f, ax = plt.subplots(1,2, figsize=(6, 3))
ax[0].set_title(r'$v_y$, y=0.5',fontsize=fs)
ax[0].plot(true_ver[:,0],true_ver[:,1],'ks',label='true',markersize=3)
ax[0].plot(sim_vert05_100[:,sv_x],sim_vert05_100[:,sv_y],'bs',label='simulated',markersize=3)
ax[0].set_xlabel(r'x',fontsize=fs)
ax[0].set_ylabel(r'$v_y$',fontsize=fs)

# x=8,y=9
sv_x = 9
# vx=3,vy=4
sv_y = 3

ax[1].set_title(r'$v_x$, x=0.5',fontsize=fs)
ax[1].plot(true_hor[:,0],true_hor[:,1],'ks',markersize=3)
ax[1].plot(sim_hor05_100[:,sv_x],sim_hor05_100[:,sv_y],'bs',markersize=3)
ax[1].set_xlabel(r'y',fontsize=fs)
ax[1].set_ylabel(r'$v_x$',fontsize=fs)

ax[0].legend()
plt.tight_layout()
# plt.show()
plt.savefig('re100.pdf')



