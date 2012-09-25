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
	ax.axvline(lbmax[0],lw=0.8,color='#E24A33',alpha=0.5)
	ax.axvline(lbmax[1],lw=0.8,color='#E24A33',alpha=0.5)
	ax.get_figure().savefig('ar1_ex_em.pdf')	
