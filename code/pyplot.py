#!/usr/local/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpltools import style
from mpltools.special import errorfill
import scipy.io
from scipy.io import loadmat


def plot(argv=None):
	A = np.array # shortcut

	#fm = font_manager.FontManager()
	#p = font_manager.FontProperties("Latin Modern Roman",size=10)
	#print(fm.findfont(p))

	style.use('dippa')
	argv = sys.argv[1:] if argv is None else argv
	dataFileName = argv[0]
	try:
		outputFileName = argv[1]
	except IndexError as e:
		outputFileName = None

	beta = (np.sqrt(5)+1)/2 # golden ratio
	cycle = plt.rcParams["axes.color_cycle"]
	defc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
	
	inp = loadmat(dataFileName,squeeze_me=True,struct_as_record=False)
			
	try:
		xlabel = inp["xlabel"]
	except KeyError as e:
		xlabel = None
	try:
		ylabel = inp["ylabel"]
	except KeyError as e:
		ylabel = None
	try:
		title = inp["title"]
	except KeyError as e:
		title = None

	try:
		legend = inp["legend"]
	except KeyError as e:
		legend = None

	try:
		legendkw = _todict(inp["legendkw"])
	except KeyError as e:
		legendkw = {}

	try:
		w = inp["w"]
	except KeyError as e:
		w = 5

	try:
		h = inp["h"]
	except KeyError as e:
		h = w/beta
	
	try:
		alpha = inp["alpha"]
	except KeyError as e:
		alpha = 0.03


	fig = plt.figure(figsize=(w,h),facecolor='w')
	ax = fig.add_subplot(111)
	# plot
	for k,triplet in enumerate(inp["data"]):
		if len(triplet) < 2:
			continue
		if len(triplet) > 2 and type(triplet[2]) is not unicode:
			triplet[2] = ""
		arg = triplet[0:3]
		#print(arg)
		kw = _setcolor(triplet,cycle,defc,k)
		#print(kw)
		if kw.has_key("yerr"):
			kw["ax"] = ax
			errorfill(arg[0],arg[1],**kw)
		else:
			ax.plot(*arg,**kw)
		ax.hold(True)
	
	# stretch both dimensions by alpha
	T = A([[1,-1,0,0],[-1,1,0,0],[0,0,1,-1],[0,0,-1,1]])*alpha
	# all dimensions as 1x4 vector
	l = ax.axis("tight")
	# make the transformation (add one to make it relative to current)
	ax.axis((T+np.eye(4)).dot(l))

	if title is not None:
		ax.set_title(title)
	if xlabel is not None:
		xlabel = ax.set_xlabel(xlabel)
	if ylabel is not None:
		ylabel = ax.set_ylabel(ylabel)
	if legend is not None:
		lg = ax.legend(legend,**legendkw)
		if xlabel is not None:
			for text in lg.get_texts():
				text.set_color(xlabel.get_color())

	if outputFileName is not None:
		fig.savefig(outputFileName)



def _setcolor(triplet,cycle,defc,k):
	''' triplet[2] is a matlab '.-b' type string '''
	c = None
	try:
		 cc = [a for a in defc if triplet[2].find(a) > -1]
		 if len(cc):
		 	c = cc[0]
	except IndexError:
		pass

	try:
		kw = _todict(triplet[3])
	except IndexError:
		kw = {}
	
	try:
		c = kw["color"]
	except KeyError:
		pass

	if c is None:
		c = cycle[np.mod(k,len(cycle))]
	
	kw["color"] = c
	return kw

def _todict(matobj):
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

if __name__ == '__main__':
	plot()
