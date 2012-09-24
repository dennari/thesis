from mpltools import style
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from util import *

def draw():
	style.use('dippa')
	beta = (np.sqrt(5)+1)/2 # golden ratio
	cycle = plt.rcParams["axes.color_cycle"]
	defc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

	inp = loadmat('ar1_ex.mat',squeeze_me=True,struct_as_record=False)
	K,x = (inp["data"][0][0],inp["data"][0][1])
	K,y = (inp["data"][1][0],inp["data"][1][1])
	a,lhs = (inp["data"][2][0],inp["data"][2][1])

	alpha = 0.05;
	beta = (np.sqrt(5)+1)/2 # golden ratio
	w = 426.79134/72.27; # textwidth in in inches
	h = w/(2*beta);

	fig = plt.figure(figsize=(w,h),facecolor='w')

	ax = fig.add_subplot(121)
	l1,l2 = ax.plot(K,x,K,y,'kx',ms=3,mec='black')
	l2.set_alpha(0.9);
	axstretch(ax,alpha)
	xlabel = ax.set_xlabel(r'$k$')
	lg = ax.legend((r'$x$',r'$y$'))
	#legendtolabelcolor(lg,xlabel)
	ax.set_title(r'AR(1) simulation',family='serif')

	ax = fig.add_subplot(122);
	mi = lhs.argmax();
	#ax.plot(a,lhs,a[mi],lhs[mi],'o',mec='none',mfc='black')
	ax.plot(a,lhs)
	axstretch(ax,alpha)
	xlabel = ax.set_xlabel(r'$a$')
	ax.set_title(r'Likelihood of $a$',family='serif')
	#legendtolabelcolor(lg,xlabel)

	return(fig)

if __name__ == "__main__":
	(draw()).savefig('ar1_ex.pdf')