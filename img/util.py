import scipy.io
import numpy as np


def mattodict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, scipy.io.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict
def axstretch(ax,alpha):
	# stretch both dimensions by alpha
	T = np.array([[1,-1,0,0],[-1,1,0,0],[0,0,1,-1],[0,0,-1,1]])*alpha
	# all dimensions as 1x4 vector
	l = ax.axis("tight")
	# make the transformation (add one to make it relative to current)
	ax.axis((T+np.eye(4)).dot(l))

def legendtolabelcolor(lg,lb):
	for text in lg.get_texts():
		text.set_color(lb.get_color())