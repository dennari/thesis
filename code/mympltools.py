import numpy as np
import matplotlib.pyplot as plt

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

def padaxis(ax,alpha):
	# stretch both dimensions by alpha
	T = np.array([[1,-1,0,0],[-1,1,0,0],[0,0,1,-1],[0,0,-1,1]]).dot(alpha)
	# shows all data, centers, returns all dimensions as 1x4 vector
	l = ax.axis("tight")
	# make the transformation (add one to make it relative to current)
	return(ax.axis((T+np.eye(4)).dot(l)))

def getpadfigure(figw,margins=(0.18,0.0,0.35,0.35)):
	beta = (np.sqrt(5)+1)/2 # golden ratio
	leftmargin = margins[3] # inches
	rightmargin = margins[1];
	axwidth = figw-leftmargin-rightmargin # inches
	axheight = axwidth/beta #inches
	bottommargin = margins[2] # inches
	topmargin = margins[0] # inches
	figh = topmargin+axheight+bottommargin # inches

	fig = plt.figure(figsize=(figw,figh),facecolor='w')
	ax = fig.add_axes((leftmargin/figw,bottommargin/figh,axwidth/figw,axheight/figh))
	return((fig,ax))