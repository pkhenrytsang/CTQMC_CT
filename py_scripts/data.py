import numpy as np


def paverage(data,npoints=5,order=1):

	def average(dat,startit):
		ld = len(dat)
		avg = []
		std = []
		iters=[]
		for i in np.arange(npoints,ld,1,dtype =np.int):
			avg.append(np.average(dat[i-npoints:i]))
			std.append(np.std(dat[i-npoints:i]))
			iters.append(i+npoints+startit)
		avg = np.array(avg)
		std = np.array(std)
		iters = np.array(iters)
		return iters,avg,std
	
	for iord in xrange(order):
		startit = iord*npoints #higher order average will reduce amount of usable points
		iters,avg,std = average(data,startit)
		data = avg
	return iters-npoints,avg,std
	
def pdiff(dat,it0 = 0):

	ld = len(dat)
	
	diff = []
	iters = []
	
	for i in np.arange(1,ld,1):
		diff.append(dat[i]-dat[i-1])
		iters.append(i+it0)
	diff = np.array(diff)
	iters = np.array(iters)
	return iters+1,diff
	
def idiff(dat):

	ld = len(dat)
	
	diff = []
	iters = []
	
	for i in xrange(ld):
		diff.append(dat[i]-dat[0])
		iters.append(i)
	diff = np.array(diff)
	iters = np.array(iters)
	return iters+1,diff
