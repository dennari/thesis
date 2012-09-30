#!/usr/local/bin/python2.7
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpltools import style
#from mpltools.special import errorfill
from scipy.io import loadmat
from scipy.io.matlab.mio5_params import mat_struct
import mympltools as util


def draw(dataFileName):
	A = np.array # shortcut

	style.use('dippa')


	
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
		margins = {'margins': inp["margins"]}
	except KeyError as e:
		margins = {}
	
	try:
		alpha = inp["alpha"]
	except KeyError as e:
		alpha = 0.05

	try:
		axlim = inp["axis"]
	except KeyError as e:
		axlim = None

	try:
		ticklabels = inp["ticklabels"]
	except KeyError as e:
		ticklabels = None

	fig,ax = util.getpadfigure(w,**margins)

	# plot
	lines = []
	for k,triplet in enumerate(inp["data"]):
		if len(triplet) < 2:
			continue
		if len(triplet) > 2 and type(triplet[2]) is not unicode:
			triplet[2] = ""
		arg = triplet[0:3]
		#print(arg)
		kw = _setcolor(triplet,cycle,defc,k)
		lines.append(ax.plot(*arg,**kw)[0])

		#print(kw)
		#if kw.has_key("yerr"):
		#	kw["ax"] = ax
		#	errorfill(arg[0],arg[1],**kw)
		#else:
		#	ax.plot(*arg,**kw)
		ax.hold(True)
	
	util.padaxis(ax,alpha,l=axlim)

	if title is not None:
		ax.set_title(title,family='serif')
	if xlabel is not None:
		xlabel = ax.set_xlabel(xlabel,family='serif')
	if ylabel is not None:
		ylabel = ax.set_ylabel(ylabel,family='serif')
	if legend is not None:
		#print(legend)
		#print(lines)
		lg = ax.legend(lines,legend.flat,**legendkw)
		# if xlabel is not None:
		# 	for text in lg.get_texts():
		# 		text.set_color(xlabel.get_color())
		for text in lg.get_texts():
			text.set_family('serif')
	if ticklabels is not None:
		if not ticklabels[0]:
			ax.set_xticklabels([])
		if not ticklabels[1]:
			ax.set_yticklabels([])
			
	return(ax)



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
        if isinstance(elem, mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

if __name__ == '__main__':
	argv = sys.argv[1:]
	dataFileName = argv[0]
	outputFileName = argv[1]
	ax = draw(dataFileName)
	ax.get_figure().savefig(outputFileName)
