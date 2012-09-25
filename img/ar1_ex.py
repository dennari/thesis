from mpltools import style
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import util

def draw():
	style.use('dippa')
	beta = (np.sqrt(5)+1)/2 # golden ratio
	cycle = plt.rcParams["axes.color_cycle"]
	defc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

	inp = loadmat('ar1_ex.mat',squeeze_me=True,struct_as_record=False)
	K,x = (inp["data"][0][0],inp["data"][0][1])
	K,y = (inp["data"][1][0],inp["data"][1][1])
	a,lhs = (inp["data"][2][0],inp["data"][2][1])

	# size computations
	alpha = 0.05;
	beta = (np.sqrt(5)+1)/2 # golden ratio
	tw = 426.79134/72.27; # latex textwidth in in inches
	figw = tw/2.0 # no scaling in latex
	leftmargin = 0.35 # inches
	axwidth = figw-leftmargin # inches
	axheight = axwidth/beta #inches
	bottommargin = 0.35 # inches
	topmargin = 0.22 # inches
	figh = topmargin+axheight+bottommargin # inches

	fig1 = plt.figure(figsize=(figw,figh),facecolor='w')
	ax = fig1.add_axes((leftmargin/figw,bottommargin/figh,axwidth/figw,axheight/figh))
	l1,l2 = ax.plot(K,x,K,y,'kx',ms=3,mec='black')
	l2.set_alpha(0.8);
	util.axstretch(ax,alpha)
	xlabel = ax.set_xlabel(r'$k$')
	lg = ax.legend((r'$x$',r'$y$'))
	#legendtolabelcolor(lg,xlabel)
	ax.set_title(r'AR(1) simulation',family='serif')
	#fig.savefig('ar1_ex_a.pdf')

	fig2 = plt.figure(figsize=(figw,figh),facecolor='w')
	ax = fig2.add_axes((leftmargin/figw,bottommargin/figh,axwidth/figw,axheight/figh))
	mi = lhs.argmax();
	#ax2.plot(a,lhs,a[mi],lhs[mi],'o',mec='none',mfc='black')
	ax.plot(a,lhs)
	util.axstretch(ax,alpha)
	xlabel = ax.set_xlabel(r'$a$')
	ax.set_title(r'Likelihood of $a$',family='serif')
	#fig.savefig('ar1_ex_b.pdf')
	#legendtolabelcolor(lg,xlabel)

	return((fig1,fig2))

if __name__ == "__main__":
	fig1,fig2 = draw()
	#inchs = fig.dpi_scale_trans.inverted()
	#extent = ax1.get_window_extent().transformed(inchs).expanded(1.2,1.4)
	fig1.savefig('ar1_ex_a.pdf')	
	#extent = ax2.get_window_extent().transformed(inchs).expanded(1.2,1.4)
	fig2.savefig('ar1_ex_b.pdf')
	#fig.savefig('ar1_ex.pdf')