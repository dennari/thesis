from mpl_toolkits.mplot3d import Axes3D
from mpltools import style
from scipy.io import loadmat

style.use('dippa3D')
d=loadmat('data/Harmonic_qx_qw_python.mat')
ax = (figure()).add_subplot(111,projection='3d')
lh = d['lh_em_n']
est = d['est_em_n']
for k in range(10): 
	ax.plot(est[0,:,k],est[1,:,k],lh[:,k],alpha=1.0,color='#348ABD',lw=3)
	#ax.plot(est[0,:,k],est[1,:,k],lh[:,k],'.',mfc='black',ms=4,alpha=0.4)


ax = (figure()).add_subplot(111,projection='3d')
lh = d['lh_bfgs_n']
est = d['est_bfgs_n']
for k in range(10): 
	ax.plot(est[0,:,k],est[1,:,k],lh[:,k],alpha=1.0,color='#E24A33',lw=3)
	#ax.plot(est[0,:,k],est[1,:,k],lh[:,k],'.',mfc='black',ms=4,alpha=0.4)
