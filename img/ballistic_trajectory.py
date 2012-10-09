#!/usr/local/bin/python
from mpltools import style
from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import mympltools as util
import matpyplot

if __name__ == "__main__":
    ax = matpyplot.draw('ballistic_trajectory.mat')
    lim = ax.axis("image")
    ax.set_xticks(np.arange(0,lim[1]+1,10))
    ax.set_yticks(np.arange(0,lim[3]+1,10))
    util.padaxis(ax,0.03,l=lim)
    #ax.grid()

    ax.get_figure().savefig('ballistic_trajectory.pdf')	
