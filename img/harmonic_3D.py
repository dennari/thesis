from mpl_toolkits.mplot3d import Axes3D
from mpltools import style
from scipy.io import loadmat
import mympltools as util
import numpy as np
import matplotlib.pyplot as plt

style.use('dippa3D')
d=loadmat('harmonic_3D.mat')

keys = ('lh_em_n','est_em_n','lh_bfgs_n','est_bfgs_n')
colors = ('#348ABD','#E24A33')
labels = ('$\sigma_\omega$','$\sigma_x$','$\ell$')
w = 426.79134/72.27
#h = w
#fig = plt.figure(figsize=(w,2*w),facecolor='w')
zllim = 1.5e4
for j in range(2):
	fig,ax = util.getpadfigure(figw=w/1.2,figh=w/1.2,is3D=True)
	#fig = plt.figure(figsize=(w,w),facecolor='w')
	#ax = fig.add_subplot(111,projection='3d')
	lh = d[keys[j*2]]
	est = d[keys[j*2+1]]
	for k in range(100): 
		zi = lh[:,k]>zllim
		ax.plot(est[0,zi,k],est[1,zi,k],lh[zi,k],
		alpha=0.5,
		color=colors[j],
		lw=0.9,
		marker=None,
		markersize=3)
	ax.plot(est[0,-1,:],est[1,-1,:],lh[-1,:],
		alpha=0.7,
		color=colors[j],
		markerfacecolor='black',#colors[j],
		markeredgewidth=0.,
		lw=0.5,
		marker='*',
		markersize=4)

	#ax.set_zlim([1.6e4, 2.0e4])
	ax.set_xlim([-1.5, 0.5])
	ax.view_init(55,-40)
	#plt.axis('tight')
	ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
	ax.set_zlabel(labels[2])
	fig.patch.set_alpha(0.)
	fig.savefig('harmonic_3D%d.pdf'%j)
	# ax.xaxis.set_ticklabels([])
	# ax.yaxis.set_ticklabels([])
	# ax.zaxis.set_ticklabels([])


