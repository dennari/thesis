#!/usr/local/bin/python
import numpy as np
import mympltools as util
import matpyplot

if __name__ == "__main__":
	ax = matpyplot.draw('ar1_ex_em.mat')
	lines = ax.get_lines();
	
	# put better data limits
	lh = lines[0]
	x,y = lh.get_xdata(),lh.get_ydata()
	alpha = np.diag([0.05, 0.05, 0.1, 0.1])
	T = np.array([[1,-1,0,0],[-1,1,0,0],[0,0,1,-1],[0,0,-1,1]]).dot(alpha)
	ax.axis((T+np.eye(4)).dot([x.min(),x.max(),y.min(),y.max()]))

	# add vertical lines
	lbmax = lines[-1].get_xdata()
	lines[-1].set_visible(False)
	texts = (r'$\theta_{k}$',r'$\theta_{k+1}$',r'$\theta_{k+2}$')
	for k in range(len(lbmax)):
		ax.axvline(lbmax[k],lw=0.4,color='black',alpha=0.4)
		ax.annotate(texts[k],xy=(lbmax[k],1),xytext=(-5,3),
			xycoords=('data','axes fraction'),
			textcoords='offset points',size='small')
	
	ax.get_figure().savefig('ar1_ex_em.pdf')	
