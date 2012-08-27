#!/usr/local/bin/python
import numpy as np
import sys
import matplotlib.pyplot as plt
from mpltools import style
from scipy.io import loadmat

def plot(argv=None):
	A = np.array # shortcut

	#fm = font_manager.FontManager()
	#p = font_manager.FontProperties("Latin Modern Roman",size=10)
	#print(fm.findfont(p))

	style.use('dippa')
	argv = sys.argv[1:] if argv is None else argv
	dataFileName = argv[0]
	outputFileName = argv[1]

	beta = (np.sqrt(5)+1)/2 # golden ratio
	inp = loadmat(dataFileName,squeeze_me=True)
			
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
	ax.hold()
	# plot
	ddim = A([np.nan]*4) 
	cycle = plt.rcParams["axes.color_cycle"]
	for k,triplet in enumerate(inp["data"]):
		try:
			ax.plot(*triplet[0:3],**triplet[3])
		except IndexError as e:
			try:
				ax.plot(*triplet[0:3])
			except IndexError as e:

		x = triplet[0]
		xmin = x.min()
		xmax = x.max()
		if np.isnan(ddim[0]) or xmin < ddim[0]:
			ddim[0] = xmin
		if np.isnan(ddim[1]) or xmax < ddim[1]:
			ddim[1] = xmax
			
		y = triplet[1]
		ymin = y.min()
		ymax = y.max()
		if np.isnan(ddim[2]) or ymin < ddim[2]:
			ddim[2] = ymin
		if np.isnan(ddim[3]) or ymax < ddim[3]:
			ddim[3] = ymax

		if len(x.shape) == 2:
			if(x.shape[0] < x.shape[1]):
				x = x.T
		if len(y.shape) == 2:
			if(y.shape[0] < y.shape[1]):
				y = y.T
		if len(triplet) > 2:
			s = triplet[2]
			if len(triplet) > 3 and type(triplet[3]) is dict:
				kw = triplet[3]
				ax.plot(x,y,s,**kw)
			else:
				ax.plot(x,y,s)
		else:
			ax.plot(x,y)
		ax.hold(True)

	if title is not None:
		ax.set_title(title)
	if xlabel is not None:
		xlabel = ax.set_xlabel(xlabel)
	if ylabel is not None:
		ylabel = ax.set_ylabel(ylabel)
	if legend is not None:
		lg = ax.legend(legend)
		if xlabel is not None:
			for text in lg.get_texts():
				text.set_color(xlabel.get_color())

	
	# stretch both dimensions by alpha
	T = A([[1,-1,0,0],[-1,1,0,0],[0,0,1,-1],[0,0,-1,1]])*alpha
	# all dimensions as 1x4 vector
	l = ax.axis("tight")
	# make the transformation (add one to make it relative to current)
	ax.axis((T+np.eye(4)).dot(l))
	#l = (T+np.eye(4)).dot(l)
	#ax.set_xlim(l[:2])
	#ax.set_ylim(l[2:])
	fig.savefig(outputFileName)

if __name__ == '__main__':
	plot()
