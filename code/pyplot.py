#!/usr/local/bin/python
import numpy as np
import sys
import PyQt4
print(PyQt4)
from PyQt4 import QtCore, QtGui
import matplotlib.pyplot as plt
import tables

def plot(argv=None):
	argv = sys.argv[1:] if argv is None else argv
	dataFileName = argv[0]
	outputFileName = argv[1]

	data = tables.openFile(dataFileName)
	x = np.squeeze(data.root.x[:])
	y = np.squeeze(data.root.y[:])
	try:
		s = data.root.s[:].tostring().decode('UTF-16LE')
	except tables.NoSuchNodeError as e:
		s = None
	data.close()

	fig = plt.figure(figsize=(9,6),facecolor='w')
	ax = fig.add_subplot(111,frame_on=False)
	if s is None:
		ax.plot(x,y)
	else:
		ax.plot(x,y,s)
		#print((x,y,s))
		#ax.plot(np.linspace(0,2*np.pi,200),np.sin(np.linspace(0,2*np.pi,200)),s)
	#print(plt.rcParams)
	fig.savefig(outputFileName)

if __name__ == '__main__':
	plot()
