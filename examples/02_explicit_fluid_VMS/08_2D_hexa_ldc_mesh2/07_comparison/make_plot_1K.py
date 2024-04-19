import numpy as np
import matplotlib.pyplot as plt

fs=12
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')
plt.rc('text', usetex=True)

true_hor = np.loadtxt('re_1000_horizontal.txt',skiprows=1)
true_ver = np.loadtxt('re_1000_vertical.txt',skiprows=1)

sim_hor005_1K = np.loadtxt('hor005_1K.csv',skiprows=1,delimiter=',')
sim_hor01_1K = np.loadtxt('hor01_1K.csv',skiprows=1,delimiter=',')
sim_hor05_1K = np.loadtxt('hor05_1K.csv',skiprows=1,delimiter=',')
sim_hor09_1K = np.loadtxt('hor09_1K.csv',skiprows=1,delimiter=',')
sim_hor095_1K = np.loadtxt('hor095_1K.csv',skiprows=1,delimiter=',')

sim_vert005_1K = np.loadtxt('vert005_1K.csv',skiprows=1,delimiter=',')
sim_vert01_1K = np.loadtxt('vert01_1K.csv',skiprows=1,delimiter=',')
sim_vert05_1K = np.loadtxt('vert05_1K.csv',skiprows=1,delimiter=',')
sim_vert09_1K = np.loadtxt('vert09_1K.csv',skiprows=1,delimiter=',')
sim_vert095_1K = np.loadtxt('vert095_1K.csv',skiprows=1,delimiter=',')

# x=8,y=9
sv_x = 8
# vx=3,vy=4
sv_y = 4

f, ax = plt.subplots(2,5, figsize=(13, 5))
ax[0,0].set_title(r'$v_y$, y=0.05',fontsize=fs)
ax[0,0].plot(true_ver[:,0],true_ver[:,1],'ks',label='true',markersize=2)
ax[0,0].plot(sim_hor005_1K[:,sv_x],sim_hor005_1K[:,sv_y],'bs',label='simulated',markersize=2)
ax[0,0].set_xlabel(r'x',fontsize=fs)
ax[0,0].set_ylabel(r'$v_y$',fontsize=fs)

ax[0,1].set_title(r'$v_y$, y=0.1',fontsize=fs)
ax[0,1].plot(true_ver[:,0],true_ver[:,2],'ks',markersize=2)
ax[0,1].plot(sim_hor01_1K[:,sv_x],sim_hor01_1K[:,sv_y],'bs',markersize=2)
ax[0,1].set_xlabel(r'x',fontsize=fs)
ax[0,1].set_ylabel(r'$v_y$',fontsize=fs)

ax[0,2].set_title(r'$v_y$, y=0.5',fontsize=fs)
ax[0,2].plot(true_ver[:,0],true_ver[:,3],'ks',markersize=2)
ax[0,2].plot(sim_hor05_1K[:,sv_x],sim_hor05_1K[:,sv_y],'bs',markersize=2)
ax[0,2].set_xlabel(r'x',fontsize=fs)
ax[0,2].set_ylabel(r'$v_y$',fontsize=fs)

ax[0,3].set_title(r'$v_y$, y=0.9',fontsize=fs)
ax[0,3].plot(true_ver[:,0],true_ver[:,4],'ks',markersize=2)
ax[0,3].plot(sim_hor09_1K[:,sv_x],sim_hor09_1K[:,sv_y],'bs',markersize=2)
ax[0,3].set_xlabel(r'x',fontsize=fs)
ax[0,3].set_ylabel(r'$v_y$',fontsize=fs)

ax[0,4].set_title(r'$v_y$, y=0.95',fontsize=fs)
ax[0,4].plot(true_ver[:,0],true_ver[:,5],'ks',markersize=2)
ax[0,4].plot(sim_hor095_1K[:,sv_x],sim_hor095_1K[:,sv_y],'bs',markersize=2)
ax[0,4].set_xlabel(r'x',fontsize=fs)
ax[0,4].set_ylabel(r'$v_y$',fontsize=fs)

# x=8,y=9
sv_x = 9
# vx=3,vy=4
sv_y = 3

ax[1,0].set_title(r'$v_x$, x=0.05',fontsize=fs)
ax[1,0].plot(true_hor[:,0],true_hor[:,1],'ks',markersize=2)
ax[1,0].plot(sim_vert005_1K[:,sv_x],sim_vert005_1K[:,sv_y],'bs',markersize=2)
ax[1,0].set_xlabel(r'y',fontsize=fs)
ax[1,0].set_ylabel(r'$v_x$',fontsize=fs)

ax[1,1].set_title(r'$v_x$, x=0.1',fontsize=fs)
ax[1,1].plot(true_hor[:,0],true_hor[:,2],'ks',markersize=2)
ax[1,1].plot(sim_vert01_1K[:,sv_x],sim_vert01_1K[:,sv_y],'bs',markersize=2)
ax[1,1].set_xlabel(r'y',fontsize=fs)
ax[1,1].set_ylabel(r'$v_x$',fontsize=fs)

ax[1,2].set_title(r'$v_x$, x=0.5',fontsize=fs)
ax[1,2].plot(true_hor[:,0],true_hor[:,3],'ks',markersize=2)
ax[1,2].plot(sim_vert05_1K[:,sv_x],sim_vert05_1K[:,sv_y],'bs',markersize=2)
ax[1,2].set_xlabel(r'y',fontsize=fs)
ax[1,2].set_ylabel(r'$v_x$',fontsize=fs)

ax[1,3].set_title(r'$v_x$, x=0.09',fontsize=fs)
ax[1,3].plot(true_hor[:,0],true_hor[:,4],'ks',markersize=2)
ax[1,3].plot(sim_vert09_1K[:,sv_x],sim_vert09_1K[:,sv_y],'bs',markersize=2)
ax[1,3].set_xlabel(r'y',fontsize=fs)
ax[1,3].set_ylabel(r'$v_x$',fontsize=fs)

ax[1,4].set_title(r'$v_x$, x=0.95',fontsize=fs)
ax[1,4].plot(true_hor[:,0],true_hor[:,5],'ks',markersize=2)
ax[1,4].plot(sim_vert095_1K[:,sv_x],sim_vert095_1K[:,sv_y],'bs',markersize=2)
ax[1,4].set_xlabel(r'y',fontsize=fs)
ax[1,4].set_ylabel(r'$v_x$',fontsize=fs)

ax[0,0].legend()
plt.tight_layout()
# plt.show()
plt.savefig('re1K.pdf')



