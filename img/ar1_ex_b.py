#!/usr/local/bin/python
from mpltools import style
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import mympltools as util
import matpyplot

if __name__ == "__main__":
    ax = matpyplot.draw('ar1_ex_b.mat')
    x,y = ax.get_lines()[0].get_xdata(),ax.get_lines()[0].get_ydata()
    mi = y.argmax()
    ax.hold(True)
    #ax.plot(x[mi], y[mi],'k*')
    #xt = ax.get_xticks()
    #xt.append(x[mi])
    ax.axvline(x[mi],lw=0.4,color='black',alpha=1.0)
    ax.annotate(r'MAP',xy=(x[mi],1),xytext=(-8,2),
        xycoords=('data','axes fraction'),
        textcoords='offset points',size='x-small')
	# ax.annotate('ML/MAP',
 #                    xy=(x[mi], y[mi]), xycoords='data',
 #                    xytext=(-45,-25), textcoords='axes points',
 #                    size='small', alpha=1.0, color='0.2',
 #                    arrowprops=dict(arrowstyle="->",
 #                                    color="0.0",
 #                                    relpos=(0,0.5),
 #                                    shrinkA=5, shrinkB=6,
 #                                    patchA=None,
 #                                    patchB=None,
 #                                    connectionstyle="angle3,angleA=110,angleB=20",
 #                                    )
 #                    )
	#ax.annotate('x',xy=(x[mi], y[mi]),color="black") 
    ax.get_figure().savefig('ar1_ex_b.pdf')	
