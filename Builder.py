import math
import numpy as np
import Point9
import matplotlib.pyplot as plt
import h5py
import random
#import Character
from netCDF4 import Dataset 
#from PIL import Image

YF = 1920
XF = 1920

YSize = [YF,YF,YF,YF,YF,YF,YF,6400,6400]
XSize = [XF,XF,XF,XF,XF,XF,XF,6400,6400]
RSize = [606,909,1212,1515,1818,1212,606]

LENS = [0,0,0,0,0,0,92162,368642,1474562,5898242]
RADS = [0,0,0,13,25,50,100,225,int(YF/2),YF]

ZMLVL = [1,2,4,8,16,32]

LEVEL = 8

Grid = Point9.getPointList()
distancescale = math.sqrt(16.0/len(Grid))
NNList = list(open('NN'+str(LEVEL)+'.txt'))

Sunshade = list(open('Sunlight3.txt'))

IGBPTypes = ('Water','EvergreenNeedleleaf','EvergreenBroadleaf','DeciduousNeedleleaf','DeciduousBroadleaf','MixedForest','ClosedShrubland','OpenShrubland','WoodedSavanna','Savanna','Grasslands','PermanentWetlands','Croplands','Urban','SparseCrops','SnowIce','Desert')
RGNTypes = ('Flat','Hill','Mountain','Peak')
CLIMTypes = ('Grasslands','Plains','Desert','Tundra','Snow')

#	This Script loads the list of points into memory and provides a variety of methods useful for calculating and visualizing data on the sphere
#	Model Atlas for Nesting on a Truncated Icosahedral Star (MANTIS)

#	Returns the dot product  between point1 and point2
def dot(XYZ1, XYZ2):
	return (XYZ1[0]*XYZ2[0]+XYZ1[1]*XYZ2[1]+XYZ1[2]*XYZ2[2])

def dot_arr1(XYZ1,XYZ2):
	return np.sum(XYZ1[:,0]*XYZ2[0]+XYZ1[:,1]*XYZ2[1]+XYZ1[:,2]*XYZ2[2])

def dot_arr2(XYZ1,XYZ2):
	return np.sum(XYZ1[:,0]*XYZ2[:,0]+XYZ1[:,1]*XYZ2[:,1]+XYZ1[:,2]*XYZ2[:,2])

#	Returns the cross product between point1 and point2 as a vector (list)
def cross(r1,r2):
	tot = 0.0
	rx = r1[1]*r2[2]-r1[2]*r2[1]
	ry = r1[2]*r2[0]-r1[0]*r2[2]
	rz = r1[0]*r2[1]-r1[1]*r2[0]
	R = (rx,ry,rz)
	return R

def cross_arr2(r1,r2):
	r3 = np.zeros(shape=r1.shape,dtype=r1.dtype)
	r3[:,0] = r1[:,1]*r2[:,2]-r1[:,2]*r2[:,1]
	r3[:,1] = r1[:,2]*r2[:,0]-r1[:,0]*r2[:,2]
	r3[:,2] = r1[:,0]*r2[:,1]-r1[:,1]*r2[:,0]
	return r3
	

#	Returns the cross product between point1 and point2 as a string
def crossStr(r1,r2):
	tot = 0.0
	rx = r1[1]*r2[2]-r1[2]*r2[1]
	ry = r1[2]*r2[0]-r1[0]*r2[2]
	rz = r1[0]*r2[1]-r1[1]*r2[0]
	return str(rx)+' '+str(ry)+' '+str(rz)

#	Returns a vector from the corresponding string of XYZ coordinates
def toVector(strpoint):
	tmp = strpoint.split(' ')
	x0 = float(tmp[0])
	y0 = float(tmp[1])
	z0 = float(tmp[2])
	return (x0,y0,z0)	

#	Returns Euclidean distance between point1 and point2
def dEuclidStr(XYZ1s, XYZ2s):
	XYZ1 = toVector(XYZ1s)
	XYZ2 = toVector(XYZ2s)
	return math.sqrt((XYZ1[0]-XYZ2[0])*(XYZ1[0]-XYZ2[0])+(XYZ1[1]-XYZ2[1])*(XYZ1[1]-XYZ2[1])+(XYZ1[2]-XYZ2[2])*(XYZ1[2]-XYZ2[2]))

#	Returns Euclidean distance between point1 and point2
def dEuclid(XYZ1, XYZ2):
	return math.sqrt((XYZ1[0]-XYZ2[0])*(XYZ1[0]-XYZ2[0])+(XYZ1[1]-XYZ2[1])*(XYZ1[1]-XYZ2[1])+(XYZ1[2]-XYZ2[2])*(XYZ1[2]-XYZ2[2]))

def dEuclid_arr2(R1,R2):
	return np.sqrt((R1[:,0]-R2[:,0])**2+(R1[:,1]-R2[:,1])**2+(R1[:,2]-R2[:,2])**2)

def dSurf_arr2(R1,R2):
	D0 = np.sqrt((R1[:,0]-R2[:,0])**2+(R1[:,1]-R2[:,1])**2+(R1[:,2]-R2[:,2])**2)
#	D0 = qb.PositivesOrZeros(
	return np.arccos(1-D0**2/2)


#	Returns the distance on unit sphere between point1 and point 2
def dSurfStr(point1, point2):
	D0 = dEuclidStr(point1,point2)
	if(D0 < 2):
		return math.acos(1-D0*D0/2.0)
	else:
		return D0*math.pi/2.0

#	Returns the distance on unit sphere between point1 and point2
def dSurf(point1, point2):
	D0 = dEuclid(point1,point2)
	if(D0 < 2):
		return math.acos(1-D0*D0/2.0)
	else:
		return D0*math.pi/2.0

#	Returns xyz coordinate associated with latitude longitude degree cooridnate (lat0, lon0)
def toXYZDeg(lat0,lon0):
	lat = lat0*math.pi/180.0
	lon = lon0*math.pi/180.0
	x = math.cos(lon)*math.cos(lat)
	y = math.sin(lon)*math.cos(lat)
	z = math.sin(lat)
	return((x,y,z))

def toXY(strpoint,xframe,yframe):
#    tmp = toVector(strpoint)
	x0 = strpoint[0]
	y0 = strpoint[1]
	z0 = strpoint[2]
	phi = math.asin(z0/math.sqrt(x0*x0+y0*y0+z0*z0))
	lam = math.pi/2.0
	the = math.acos(-1*z0)
	Y = int(the*yframe/math.pi-0.5)
	if(x0 < 0 and y0 < 0):
		lam = math.atan(y0/x0)
	if(x0 == 0 and y0 < 0):
		lam = 1.0*math.pi/2.0
	if(x0 > 0 and y0 < 0):
		lam = math.atan(y0/x0)+math.pi
	if(x0 < 0 and y0 >= 0):
		lam = math.atan(y0/x0)
	if(x0 == 0 and y0 >= 0):
		lam = 3.0*math.pi/2.0
	if(x0 > 0 and y0 >= 0):
		lam = math.atan(y0/x0)+math.pi
	if(lam < 0):
		lam = lam+2*math.pi
	if(lam >= 2*math.pi):
		lam = lam-2*math.pi
	X = int(lam*xframe/(2*math.pi))
	if(x0 < 0 and y0 < 0):
		if(X > yframe/2):
			print('Help')
	if(X >= xframe):
		X = X-xframe
	if(X < 0):
		X = X+xframe
	if(Y < 0 ):
		Y = 0
	if(Y >= yframe ):
		Y = yframe-1
	return (X,Y)


#	Returns the length of XYZ string vector
def modStr(point):
	XYZ = toVector(point)
	return math.sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2])

#	Returns the length of XYZ vector
def mod(XYZ):
	return math.sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2])


#	Returns the center of a collection of points
def getCentralPoint(points):
	xtot = 0.0
	ytot = 0.0
	ztot = 0.0
	for pts in points:
		ptXYZ = pts.getXYZPos()
		xtot += ptXYZ[0]
		ytot += ptXYZ[1]
		ztot += ptXYZ[2]
	xNew = xtot/len(points)
	yNew = ytot/len(points)
	zNew = ztot/len(points)
	size = math.sqrt(xNew*xNew+yNew*yNew+zNew*zNew)
	xNew = xNew/size
	yNew = yNew/size
	zNew = zNew/size
	return (xNew,yNew,zNew)

#Returns a list of all point IDs within radius of point
def NNID(point, radius):
	nearby = []
	for p in Grid:
		tmp = p.getXYZPos()
		x1 = tmp[0]
		y1 = tmp[1]
		z1 = tmp[2]
		if(dSurf(point,(x1,y1,z1)) <= 1.1*radius*distancescale):
			nearby.append(p.getID())
	return nearby


#	For Cases where distances > pi/2 are not desired
def NNIDW(point, radius):
	nearby = []
	for p in Grid:
		tmp = p.getXYZPos()
		x1 = tmp[0]
		y1 = tmp[1]
		z1 = tmp[2]
		if(dSurf(point,(x1,y1,z1)) <= min(1.1*radius*distancescale, math.pi/2.0) ):
			nearby.append(p.getID())
	return nearby

#Returns a list of all point IDs within radius of point
def NNID2(PID, radius):
#	print(PID)
	nearby = [PID,]
	
	for i in range(int(radius)):
		for j in list(set(nearby)):
			for k in NNList[j].split(' '):
				nearby.append(int(k))

	return nearby

# Returns an numerical average of nearest points
def NN_NUMAVG(field, point, radius):
	tmp = point.getXYZPos()
	x1 = tmp[0]
	y1 = tmp[1]
	z1 = tmp[2]
	nearby = NNID(point, radius)
	Total = 0.0
	Count = 0
	for p in nearby:
		Total += field[p]
		Count += 1
	return Total/Count

def NN_NUMAVG2(field, point, radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	Total = 0.0
	Count = 0
	for Value in nearby:
		Total += field[Value]
		Count += 1
	return Total/Count


#	Numerical Average of Nearby
def NN_NUMAVG3(field_fname, point, radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	return ER_DSCAVG(field_fname, nearby)

def ER_DSCAVG3(FDATA, point, radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	return ER_DSCAVG(FDATA, nearby)

def ER_DSCAVG4(FDATA, point, radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	return ER_DSCAVG2(FDATA, nearby)


def ER_NUMAVG3(FDATA, point, radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	return ER_NUMAVG(FDATA, nearby)


#	Descrete Average (Bayes Classifier) of points on ER Grid
def ER_DSCAVG(FDATA, points):
	IMGDATA = plt.imread('ToXY7.png')
	switch = False
	total = 0.0
	values = []
	counts = []
	uniquevalues = []
	kos = []
	for i in range(2160):
		I = (i+0.5)/12-90.0
		Values = []
		UniqueValues = list([])
		for j in range(4320):
			J = (j+0.5)/12
			COUNT = 0
			if(any(points) == C2I(IMGDATA[i,j])):
				if(switch):
					Values.append(C2I(FDATA[i,j]))
				else:
					Values.append(FDATA[i,j])
				UniqueValues = list(set(list(Values)))
		if(len(Values) != 0):
			MaxField = Values[0]
			MaxCount = 0
			for Point in UniqueValues:
				ct = 0
				for pt in Values:
					if(Point == pt):
						ct += 1
				if(ct > MaxCount):
					MaxCount = ct
					MaxField = Point
			values.append(Values)
			counts.append(MaxCount)
			kos.append(math.cos(I*math.pi/180.0))
	minkos = 10000
	for i in range(len(kos)):
		if(kos[i] < minkos):
			minkos = kos[i]
	COSCT = []
	for i in range(len(kos)):
		COSCT.append(int(kos[i]/minkos))
	Values = []
	UniqueValues = []
	for i in range(len(COSCT)):
		for j in range(COSCT[i]):
			for k in values[i]:
				Values.append(k)
				UniqueValues = list(set(Values))
	for Point in UniqueValues:
		ct = 0
		for pt in Values:
			if(Point == pt):
				ct += 1
		if(ct > MaxCount):
			MaxCount = ct
			MaxField = Point
	return Point

#	Descrete Average (Bayes Classifier) of points on ER Grid
def ER_DSCAVG2(FDATA, points):
	IMGDATA = plt.imread('ToXY7.png')
	switch = True
	total = 0.0
	values = []
	counts = []
	uniquevalues = []
	kos = []
	for i in range(2160):
		I = (i+0.5)/12-90.0
		Values = []
		UniqueValues = list([])
		for j in range(4320):
			J = (j+0.5)/12
			COUNT = 0
			if(any(points) == C2I(IMGDATA[i,j])):
				if(switch):
					Values.append(C2I(FDATA[i,j]))
				else:
					Values.append(FDATA[i,j])
				UniqueValues = list(set(list(Values)))
		if(len(Values) != 0):
			MaxField = Values[0]
			MaxCount = 0
			for Point in UniqueValues:
				ct = 0
				for pt in Values:
					if(Point == pt):
						ct += 1
				if(ct > MaxCount):
					MaxCount = ct
					MaxField = Point
			values.append(Values)
			counts.append(MaxCount)
			kos.append(math.cos(I*math.pi/180.0))
	minkos = 10000
	for i in range(len(kos)):
		if(kos[i] < minkos):
			minkos = kos[i]
	COSCT = []
	for i in range(len(kos)):
		COSCT.append(int(kos[i]/minkos))
	Values = []
	UniqueValues = []
	for i in range(len(COSCT)):
		for j in range(COSCT[i]):
			for k in values[i]:
				Values.append(k)
				UniqueValues = list(set(Values))
	MaxCount = 0
	MaxField = Values[0]
	for Point in UniqueValues:
		ct = 0
		for pt in Values:
			if(Point == pt):
				ct += 1
		if(ct > MaxCount):
			MaxCount = ct
			MaxField = Point
	return MaxField

def ER_DSCAVG5(FDATA, points):
	IMGDATA = plt.imread('ToXY'+str(LEVEL)+'.png')
	values = []
	wgt = []
	for i in points:
		tmp = i.split(' ')
		I = (int(tmp[0])+0.5)/12-90.0
		Values = []
		J = (int(tmp[1])+0.5)/12
		values.append(FDATA[int(tmp[0]),int(tmp[1])])
		wgt.append(math.cos(I*math.pi/180.0))
	minkos = 10000
	for i in range(len(wgt)):
		if(wgt[i] < minkos):
			minkos = wgt[i]
	COSCT = []
	for i in range(len(wgt)):
		COSCT.append(int(wgt[i]/minkos))
	Values = []
	UniqueValues = []
	for i in range(len(values)):
		for j in range(COSCT[i]):
			Values.append(values[i])
			UniqueValues = list(set(Values))
	MaxCount = 0
	MaxField = Values[0]
	for Point in UniqueValues:
		ct = 0
		for pt in Values:
			if(Point == pt):
				ct += 1
		if(ct > MaxCount):
			MaxCount = ct
			MaxField = Point
	return MaxField

def ER_DSCAVG6(FDATA, points):
	IMGDATA = plt.imread('ToXY'+str(LEVEL)+'.png')
	values = []
	wgt = []
	for i in points:
		tmp = i.split(' ')
		I = (int(tmp[0])+0.5)/12-90.0
		Values = []
		J = (int(tmp[1])+0.5)/12
		values.append(C2I(FDATA[int(tmp[0]),int(tmp[1])]))
		wgt.append(math.cos(I*math.pi/180.0))
	minkos = 10000
	for i in range(len(wgt)):
		if(wgt[i] < minkos):
			minkos = wgt[i]
	COSCT = []
	for i in range(len(wgt)):
		COSCT.append(int(wgt[i]/minkos))
	Values = []
	UniqueValues = []
	for i in range(len(values)):
		for j in range(COSCT[i]):
			Values.append(values[i])
			UniqueValues = list(set(Values))
	MaxCount = 0
	MaxField = Values[0]
	for Point in UniqueValues:
		ct = 0
		for pt in Values:
			if(Point == pt):
				ct += 1
		if(ct > MaxCount):
			MaxCount = ct
			MaxField = Point
	return MaxField
		
#	Numerical Average of points on ER Grid
def ER_NUMAVG(FDATA, points):
	IMGDATA = plt.imread('ToXY7.png')
	total = 0.0
	values = []
	counts = []
	uniquevalues = []
	kos = []
	for i in range(2160):
		I = (i+0.5)/12-90.0
		Values = []
		UniqueValues = list([])
		icount = 0
		for j in range(4320):
			J = (j+0.5)/12
			itot = 0.0
			if(any(points) == C2I(IMGDATA[i,j])):
				itot += FDATA[i,j]
				icount += 1
		if(icount != 0):
			wgt = math.cos(I*math.pi/180.0)
			kos.append(wgt)
			values.append(wgt*itot/icount)
	ktot = 0.0
	vtot = 0.0
	for i in range(len(kos)):
		ktot += kos[i]
		vtot += values[i]
	if(vtot != 0):
		return ktot/vtot
	else:
		return 0.0
			
	

# Returns an numerical average of nearest points
def NN_DecAVG(field, point, radius):
	nearby = NNID(point, radius)
	fieldmix = list(set(field))
	fieldlen = len(fieldmix)
	fieldcount  = [] 
	maxfield = field[0]
	maxcount = 0
	for i in fieldmix:
		ct = 0
		for j in field:
			if(i == j):
				ct += 1 
		fieldcount.append(ct)
		if(ct > maxcount):
			maxcount = ct
			maxfield = i
	return maxfield	

# Returns an numerical average of nearest points
def NN_DecAVG2(field, point, radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	
	Values = []
	UniqueValues = list([])
	for Point in nearby:
		Values.append(field[Point])
		UniqueValues = list(set(Values))
	MaxField = Values[0]
	MaxCount = 0
	for Point in UniqueValues:
		ct = 0
		for pt in Values:
			if(Point == pt):
				ct += 1
		if(ct > MaxCount):
			MaxCount = ct
			MaxField = Point
	return Point

def ER_NUMAVG5(FDATA, points):
	IMGDATA = plt.imread('ToXY'+str(LEVEL)+'.png')
	values = []
	wgt = []
	for i in points:
		tmp = i.split(' ')
		I = (int(tmp[0])+0.5)/12-90.0
		Values = []
		J = (int(tmp[1])+0.5)/12
		values.append(FDATA[int(tmp[0]),int(tmp[1])])
		wgt.append(math.cos(I*math.pi/180.0))
	minkos = 10000
	for i in range(len(wgt)):
		if(wgt[i] < minkos):
			minkos = wgt[i]
	COSCT = []
	for i in range(len(wgt)):
		COSCT.append(int(wgt[i]/minkos))
	total = 0.0
	wgtal = 0.0
	for i in range(len(Values)):
		total += Values[i]*wgt[i]
		wgtal += wgt[i]
	if(wgtal == 0):
		return 0
	else:
		return total/wgtal

# 
def ZeroOne(field,point,radius):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)

	Result = 0
	for i in nearby:
		if(field[i] == 1 or Grid[i].getIGBP() != 0 ):
			Result = 1
	return Result

# Returns nearby points whose ID is less than a cutoff value
def NNCut(point, radius, cutoff):
	nearby = NNID(point, radius)
	realNB = []
	for p in nearby:
		if(p < cutoff):
			realNB.append(p)
	return realNB
		
def NNCut2(point,radius,cutoff):
	nearby = [point]
	for i in range(radius):
		nownearby = list(nearby)
		for j in nownearby:
#			print(NNList[j])
			tmp = NNList[j].split(' ')
			nearby.append(int(tmp[1]))
			nearby.append(int(tmp[2]))
			nearby.append(int(tmp[3]))
			nearby.append(int(tmp[4]))
			nearby.append(int(tmp[5]))
			if(len(tmp) == 7):
				nearby.append(int(tmp[6]))
		nownearby =  set(nearby)
		nearby = list(nownearby)
	anynearby = []
	for i in nearby:
		if(i < cutoff and i != point):
			anynearby.append(i)
	return anynearby



#	Returns coordinates of central location for list of points
def CentralPoint(points):
	xtot = 0.0
	ytot = 0.0
	ztot = 0.0
	for pts in points:
		tmp = pts.split(' ')
		xtot = xtot+float(tmp[4])
		ytot = ytot+float(tmp[5])
		ztot = ztot+float(tmp[6])
	xNew = xtot/len(points)
	yNew = ytot/len(points)
	zNew = ztot/len(points)
	rstr = (xNew,yNew,zNew)
	xNew = xNew/mod(rstr)
	yNew = yNew/mod(rstr)
	zNew = zNew/mod(rstr)
	furthestD = 0.0
	for pts in points:
		tmp = pts.split(' ')
		tmpD = dSurf((xNew,yNew,zNew),(float(tmp[4]),float(tmp[5]),float(tmp[6])))
		if(tmpD > furthestD):
			furthestD = tmpD
	return ((str(xNew)+' '+str(yNew)+' '+str(zNew)),furthestD)

def FindNN(pointID):
	NNIDs = []
	for i in NNList[pointID].split(' '):
		NNIDs.append(int(i))
	return NNIDs
	

#	Finds and returns the shortest path (inclusive of endpoints) between point1 and point2. Path is returned as list of pointID values
def GeodesicPath(point1, point2):
	Path = []
	xyz2 = point2.getXYZPos()
	ID1 = point1.getID()
	Path.append(ID1)
	ID2 = point2.getID()
	ct = 0
	while(Path[ct] != ID2):
		dmin = 1000000.0
		closestpoint = Path[ct]
		for ptid in FindNN(Path[ct]):
			pt = Point9.getPoint(ptid)
			xyztmp = pt.getXYZPos()
			if(dSurf(xyztmp,xyz2) <= dmin):
				dmin = dSurf(xyztmp,xyz2)
				closestpoint = ptid
		Path.append(closestpoint)
		ct += 1
	return Path

def GeodesicPathID(point1ID, point2ID):
	Path = []
	point1 = Grid[point1ID]
	point2 = Grid[point2ID]
	xyz2 = point2.getXYZPos()
	ID1 = point1.getID()
	Path.append(ID1)
	ID2 = point2.getID()
	ct = 0
	while(Path[ct] != ID2):
		dmin = 1000000.0
		closestpoint = Path[ct]
#		for ptid in FindNN(Path[ct]):
		for ptid in NNID2(Path[ct],1):
			pt = Point9.getPoint(ptid)
			xyztmp = pt.getXYZPos()
			if(dSurf(xyztmp,xyz2) <= dmin):
				dmin = dSurf(xyztmp,xyz2)
				closestpoint = ptid
		Path.append(closestpoint)
		ct += 1
	return Path

#	Returns the point on grid for a specified time of year (0-365) and hour of day (0-23)
def FindSolarNoon(TimeOfYear, TimeOfDay):
	if(TimeOfYear >= 366):
		print('Too many days in year!')
		TimeOfYear = TimeOfYear%366
	if(TimeOfDay >= 24):
		print('Too many hours in a day!')
		TimeOfDay = TimeOfDay%24
	lam = (12-TimeOfDay)*math.pi/12
	if(lam < 0):
		lam = lam+2*math.pi
	phi = -1*23.45*math.pi/180.0*math.sin((TimeOfYear-79)*math.pi/365.25)
	x0 = math.cos(phi)*math.cos(lam)
	y0 = math.cos(phi)*math.sin(lam)
	z0 = math.sin(phi)
	dmin = 1000000.0
	ptmin = Grid[0]
	for pts in Grid:
		pt = pts.getXYZPos()
		if(dSurf(pt,(x0,y0,z0)) <= dmin):
			dmin = dSurf(pt,(x0,y0,z0))
			ptmin = pts
	return ptmin

def getSolarXYZ(pt0,TimeOfYear,TimeOfDay):
	pt = Grid[pt0]
	if(TimeOfYear > 366):
		print('Too many days in year!')
		TimeOfYear = TimeOfYear%366
	if(TimeOfDay >= 24):
		print('Too many hours in a day!')
		TimeOfDay = TimeOfDay%24
	lam = (TimeOfDay-12)*math.pi/12
	if(lam < 0):
		lam = lam+2*math.pi
	phi = -1*23.45*math.pi/180.0*math.sin((TimeOfYear-79)*math.pi/365.25)
	x0 = math.cos(phi)*math.cos(lam)
	y0 = math.cos(phi)*math.sin(lam)
	z0 = math.sin(phi)
	return np.array((x0,y0,z0))

#	Returns solar zenith angle (0 = overhead, 5 = sunrise/sunset, 6 = dark)
def getSolarZenithAngle(pt0,TimeOfYear,TimeOfDay):
	pt = Grid[pt0]
	if(TimeOfYear > 366):
		print('Too many days in year!')
		TimeOfYear = TimeOfYear%366
	if(TimeOfDay >= 24):
		print('Too many hours in a day!')
		TimeOfDay = TimeOfDay%24
	lam = (TimeOfDay-12)*math.pi/12
	if(lam < 0):
		lam = lam+2*math.pi
	phi = -1*23.45*math.pi/180.0*math.sin((TimeOfYear-79)*math.pi/365.25)
	x0 = math.cos(phi)*math.cos(lam)
	y0 = math.cos(phi)*math.sin(lam)
	z0 = math.sin(phi)
	r0 = pt.getXYZPos()
	d0 = dSurf(pt.getXYZPos(),(x0,y0,z0))
#	hrs = range(0,15708+2618,2618)
#	h0 = 7
	hrs = range(0,16392,683)
	h0 = 25
	for i in range(len(hrs)-1):
		if(d0*10000.0 >= hrs[i] and d0*10000.0 < hrs[i+1]):
			h0 = i
	return h0, pt.getLatDeg(), pt.getLonDeg()

#	Shows color in different shades of light
def ColorMultiply(sunlight,color):
	newcol = np.ones(shape=(4),dtype=np.float16)
	newcol[0] = sunlight[0]*color[0]
	newcol[1] = sunlight[1]*color[1]
	newcol[2] = sunlight[2]*colos[2]
	return newcol

#	Shades color col at a point pt based on time of year ty and time of day td, generates population lighting
def SolarShade(col,pt0,ty,td,cld):
	if(random.random() < (cld)):
		col = np.array((0.95,0.95,0.95,1.0))
	HourDistance = getSolarZenithAngle(pt0,ty,td)
	
	pt = Grid[pt0]
	sunshade = Sunshade[int(HourDistance)].split(' ')
	newcol = ColorMultiply( (float(sunshade[0]),float(sunshade[1]),float(sunshade[2])),col)
	if(HourDistance >= 62):
		if(pt.getPOPL() > 7):
			newcol = np.array((1.0*pt.getPOPL()/12.0,176/255.0*pt.getPOPL()/12.0,0.0,1.0))
		else:
			newcol[0] = col[0]*0.2
			newcol[1] = col[1]*0.2
			newcol[2] = col[2]*0.2
	return newcol

def SolarShade2(col, HourDistance, cld):
	if(random.random() < cld/9.0):
		col = np.array((0.95,0.95,0.95,1.0))
#	pt = Grid[pt0]
	sunshade = Sunshade[int(HourDistance*HourDistance/8.0)].split(' ')
	newcol = ColorMultiply( (float(sunshade[0]),float(sunshade[1]),float(sunshade[2])),col)
	if(HourDistance >= 24):
#		if(pt.getPOPL() > 7):
#			newcol = np.array((1.0*pt.getPOPL()/12.0,176/255.0*pt.getPOPL()/12.0,0.0,1.0))
#		else:
			newcol[0] = col[0]*0.2
			newcol[1] = col[1]*0.2
			newcol[2] = col[2]*0.2
	return newcol

#	Averages Colors
def AvgColors(colors):
	dt = np.array((0.0,0.0,0.0))
	ct = 0
	for c in colors:
		dt += c[0:3]
		ct += 1
	if(ct == 0):
		return dt
	else:
		basecol = dt/ct
		return np.array((basecol[0],basecol[1],basecol[2],1.0))

#	Averages Colors with weights wgts
def AvgColorsWGT(colors, wgts):
	dt = np.array((0.0,0.0,0.0,1.0))
	ct = 0
	for c in colors:
		dt += c
		ct += 1
	if(ct == 0):
		return dt
	else:
		return dt/ct 

#	Bayes Classifier for Color Averaging

def commoncolors(colores):
	colors = []
	for c in colores:
		colors.append(C2I(c))
	
	dt = np.array((0.0,0.0,0.0,1.0))
	vals = list(set(colors))
	nums = np.zeros(shape=(len(vals)),dtype=np.float16)
	valmax = C2I(np.array((1.0,1.0,1.0,1.0)))
	nummax = 0
	ct = 0
	for v in vals:
		if(v == C2I(dt)):
			ct += 1
		else:
			c0 = 0
			for c in colors:
				if(c == v):
					c0 += 1
			if(c0 > nummax):
				c0 = nummax
				valmax = v
	return INT2Color(valmax)

#	Fills Image Blackness with surrounding color average

def FillImage(plane,res):
	print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	OLD = plt.imread('Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float32)
	for i in range(YF):
		for j in range(XF):
			if(C2I(OLD[i,j]) == 0):
				cols = []
				if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
					cols.append(OLD[max(i-1,0),(j-1)%XF,:])
				if(C2I(OLD[max(i-1,0),j,:]) != 0):
					cols.append(OLD[max(i-1,0),j,:])
				if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
					cols.append(OLD[max(i-1,0),(j+1)%XF,:])

				if(C2I(OLD[i,(j-1)%XF,:]) != 0):
					cols.append(OLD[i,(j-1)%XF,:])
				if(C2I(OLD[i,(j+1)%XF,:]) != 0):
					cols.append(OLD[i,(j+1)%XF,:])

				if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
				if(C2I(OLD[min(i+1,YF-1),j,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
				NEW[i,j,0:3] = AvgColors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png',arr=NEW)

def AssetFill(lat,lon,res):
	
	if(res >= 2):


		print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
		XF = XSize[res]
		YF = YSize[res]
		RF = RSize[res]
		
		OLD0 = plt.imread('Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
		NEW0 = np.ones(shape=(OLD0.shape),dtype=np.float16)
#		while(CA2I(OLD).any() == 0 ):
#		for i in range(4):
		NEW0 = FastFillAVG(OLD0)
#		NEW0 = FastFillAVG(NEW0)
#		NEW0 = FastFillAVG(NEW0)
#		NEW0 = FastFillAVG(NEW0)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW0)
	else:
		print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
		XF = XSize[res]
		YF = YSize[res]
		RF = RSize[res]
		
		OLD0 = plt.imread('Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
		NEW0 = np.ones(shape=(OLD0.shape),dtype=np.float16)
#		while(CA2I(OLD).any() == 0 ):
#		for i in range(4):
		NEW0 = FastFillAVG(OLD0)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW0)

		

def ChangeColor(OLD,colold,colnew):
	MASK = SelectValue2( CA2I(OLD) , C2I(colold))
#	NEW = np.zeros(shape=(old.shape[0],old.shape[1],4),dtype=np.float16)
#	COL = np.zeros(shape=(old.shape[0],old.shape[1]),dtype=int)
	COL = CA2I(OLD)*(1-MASK)+C2I(colnew)*MASK
	return INTS2Color(COL)


def ChangeColorI(OLD,colold,colnew):
	MASK = SelectValue2( CA2I(OLD) , colold)
#	NEW = np.zeros(shape=(old.shape[0],old.shape[1],4),dtype=np.float16)
#	COL = np.zeros(shape=(old.shape[0],old.shape[1]),dtype=int)
	COL = CA2I(OLD)*(1-MASK)+colnew*MASK
	return INTS2Color(COL)

	

def IMGFill(lat,lon,res):
	
	if(res >= 2):


		print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
		XF = XSize[res]
		YF = YSize[res]
		RF = RSize[res]
		
		OLD0 = plt.imread('IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
		OLD0 = ChangeColorI(OLD0,0,1474562)
		OLD0 = ChangeColorI(OLD0,16777215,0)
		OLD0 = ChangeColorI(OLD0,1476224,0)
		NEW0 = np.ones(shape=(OLD0.shape),dtype=np.float16)
#		while(CA2I(OLD).any() == 0 ):
#		for i in range(4):
		NEW0 = FastFillBAY(OLD0,1,0)
		NEW0 = FastFillBAY(NEW0,1,2)
		NEW0 = FastFillBAY(NEW0,0,1)
		NEW0 = FastFillBAY(NEW0,2,1)
		NEW0 = ChangeColorI(NEW0,0,16777215)#1476224)
		NEW0 = ChangeColorI(NEW0,1474562,0)
		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW0)
	else:
		print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
		XF = XSize[res]
		YF = YSize[res]
		RF = RSize[res]
		
		OLD0 = plt.imread('IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
		OLD0 = ChangeColorI(OLD0,0,1474562)
		OLD0 = ChangeColorI(OLD0,16777215,0)

		NEW0 = np.ones(shape=(OLD0.shape),dtype=np.float16)
#		while(CA2I(OLD).any() == 0 ):
#		for i in range(4):
		NEW0 = FastFillBAY(OLD0,1,0)
		NEW0 = FastFillBAY(NEW0,1,2)
		NEW0 = FastFillBAY(NEW0,0,1)
		NEW0 = FastFillBAY(NEW0,2,1)
		NEW0 = ChangeColorI(NEW0,0,1476224)
		NEW0 = ChangeColorI(NEW0,1474562,0)
		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW0)


def StarFill(lat,lon):
	
	res = 2
	if(res >= 2):


		print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
#		XF = XSize[res]
#		YF = YSize[res]
#		RF = RSize[res]
		
		OLD0 = plt.imread('Assets/STAR/Plane'+str(lat)+'x'+str(lon)+'.png')
		
		NEW0 = np.ones(shape=(OLD0.shape),dtype=np.float16)
#		while(CA2I(OLD).any() == 0 ):
#		for i in range(4):
		NEW0 = FastFillAVG(OLD0)
#		NEW0 = FastFillAVG(NEW0)
#		NEW0 = FastFillAVG(NEW0)
#		NEW0 = FastFillAVG(NEW0)
		plt.imsave(fname='Assets/STAR/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW0)
	
	

def FastFillAVG(OLD):

#	MASK0 = SelectValueC(CA2I(OLD),1)
#	OLD[:,:,0:3] = OLD[:,:,0:3]*(1-MASK0[:,:,0:3])
#	return old
	OLD = ChangeColorI(OLD,16777215,16777214)

	val = 0
	CMASK = SelectValue2( CA2I(OLD) , val)
	MASK = SelectValueC( CA2I(OLD) , val)
	NEW = np.zeros(shape=(YF+2,XF+2,4),dtype=np.float16)
	TOT = np.zeros(shape=(YF+2,XF+2),dtype=int)	
	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(YF,XF,4),dtype=np.float16)
	for i in range(3):
		for j in range(3):
			NEW[i:YF+i,j:XF+j,0:3] += OLD[:,:,0:3]
			TOT[i:YF+i,j:XF+j] += (1-CMASK[:,:])
#	BMASK = SelectValue(TOT,1)
#	AMASK = SelectValue(TOT,2)
#	NEW2[:,:,0] = OLD[:,:,0]*AMASK[1:YF+1,1:XF+1]+NEW2[1:YF+1,1:XF+1,0]*BMASK[1:YF+1,1:XF+1]
#	NEW2[:,:,2] = OLD[:,:,2]*AMASK[1:YF+1,1:XF+1]+NEW2[1:YF+1,1:XF+1,2]*BMASK[1:YF+1,1:XF+1]
#
	
	NEW2[:,:,0] = (NEW[1:YF+1,1:XF+1,0]+0.001)/(TOT[1:YF+1,1:XF+1]+0.001)
	NEW2[:,:,1] = (NEW[1:YF+1,1:XF+1,1]+0.001)/(TOT[1:YF+1,1:XF+1]+0.001)
	NEW2[:,:,2] = (NEW[1:YF+1,1:XF+1,2]+0.001)/(TOT[1:YF+1,1:XF+1]+0.001)


	OLD[:,:,0:3] += NEW2[:,:,0:3]*MASK[:,:,0:3]
#	NEW[:,:,0:3] = OLD[:,:,0:3]*MASK[:,:,0:3]

	MASK = SelectValueC(CA2I(OLD),16777215)
	OLD[:,:,0:3] = OLD[:,:,0:3]*(1-MASK[:,:,0:3])
	
	return OLD

def CheckIMGFill(lat,lon,res):
	A = plt.imread('IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
	for i in range(A.shape[0]):
		for j in range(A.shape[1]):
			if(C2I(A[i,j]) > 1476224):
				print(str(i)+' '+str(j)+' '+str(C2I(A[i,j])))

def FastFillBAY(OLD,i,j):

#	MASK0 = SelectValueC(CA2I(OLD),1)
#	OLD[:,:,0:3] = OLD[:,:,0:3]*(1-MASK0[:,:,0:3])
#	return old

	val = 0
	CMASK = SelectValue2( CA2I(OLD) , val)
#	MASK = SelectValueC( CA2I(OLD) , val)
	NEW = np.zeros(shape=(YF+2,XF+2,4),dtype=np.float16)
#	TOT = np.zeros(shape=(YF+2,XF+2),dtype=int)	
	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(YF,XF,4),dtype=np.float16)
	
	

		

#	NEW[1:YF+1,1:XF+1,0:3] += OLD[:,:,0:3]
#	TOT[1:YF+1,1:XF+1] += (1-CMASK[:,:])
	NEW[i:YF+i,j:XF+j,0:3] += OLD[:,:,0:3]
#	TOT[i:YF+i,j:XF+j] += (1-CMASK[:,:])
	NEW2[:,:,0] = OLD[:,:,0]*(1-CMASK[:,:])+NEW[1:YF+1,1:XF+1,0]*CMASK[:,:]
	NEW2[:,:,1] = OLD[:,:,1]*(1-CMASK[:,:])+NEW[1:YF+1,1:XF+1,1]*CMASK[:,:]
	NEW2[:,:,2] = OLD[:,:,2]*(1-CMASK[:,:])+NEW[1:YF+1,1:XF+1,2]*CMASK[:,:]
#	for i in range(3):
#		for j in range(3):
#			NEW[i:YF+i,j:XF+j,0:3] += OLD[:,:,0:3]
#			TOT[i:YF+i,j:XF+j] += (1-CMASK[:,:])
	
#	NEW2[:,:,0] = (NEW[1:YF+1,1:XF+1,0]+0.001)/(TOT[1:YF+1,1:XF+1]+0.001)
#	NEW2[:,:,1] = (NEW[1:YF+1,1:XF+1,1]+0.001)/(TOT[1:YF+1,1:XF+1]+0.001)
#	NEW2[:,:,2] = (NEW[1:YF+1,1:XF+1,2]+0.001)/(TOT[1:YF+1,1:XF+1]+0.001)


#	OLD[:,:,0:3] += NEW2[:,:,0:3]*MASK[:,:,0:3]
#	NEW[:,:,0:3] = OLD[:,:,0:3]*MASK[:,:,0:3]

#	MASK = SelectValueC(CA2I(OLD),16777215)
#	OLD[:,:,0:3] = OLD[:,:,0:3]*(1-MASK[:,:,0:3])
	
	return NEW2

	

def FillImageLL(lat,lon,res):
	print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	RF = RSize[res]
	OLD = plt.imread('Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float16)
	for i in range(YF):
		for j in range(XF):
			if(C2I(OLD[i,j]) == 0 and ((i-540)**2+(j-960)**2) < RF**2 ):
				cols = []
				if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
					cols.append(OLD[max(i-1,0),(j-1)%XF,:])
				if(C2I(OLD[max(i-1,0),j,:]) != 0):
					cols.append(OLD[max(i-1,0),j,:])
				if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
					cols.append(OLD[max(i-1,0),(j+1)%XF,:])

				if(C2I(OLD[i,(j-1)%XF,:]) != 0):
					cols.append(OLD[i,(j-1)%XF,:])
				if(C2I(OLD[i,(j+1)%XF,:]) != 0):
					cols.append(OLD[i,(j+1)%XF,:])

				if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
				if(C2I(OLD[min(i+1,YF-1),j,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
				NEW[i,j,0:3] = AvgColors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW)

#	Now for the IMG part

	OLD = plt.imread('IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float16)
	Zero = C2I((1.0,1.0,1.0,1.0))
	for i in range(YF):
		for j in range(XF):
			if(C2I(OLD[i,j]) == Zero and ((i-540)**2+(j-960)**2) < RF**2 ):
				cols = []
				if(C2I(OLD[max(i-1,0),max(j-1,0),:]) != Zero):
					cols.append(OLD[max(i-1,0),max(j-1,0),:])
				if(C2I(OLD[max(i-1,0),j,:]) != Zero):
					cols.append(OLD[max(i-1,0),j,:])
				if(C2I(OLD[max(i-1,0),min(j+1,XF-1),:]) != Zero):
					cols.append(OLD[max(i-1,0),min(j+1,XF-1),:])

				if(C2I(OLD[i,max(j-1,0),:]) != Zero):
					cols.append(OLD[i,max(j-1,0),:])
				if(C2I(OLD[i,(j+1)%XF,:]) != Zero):
					cols.append(OLD[i,(j+1)%XF,:])

				if(C2I(OLD[min(i+1,YF-1),max(j-1,0),:]) != Zero):
					cols.append(OLD[min(i+1,YF-1),max(j-1,0),:])
				if(C2I(OLD[min(i+1,YF-1),j,:]) != Zero):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(C2I(OLD[min(i+1,YF-1),min(j+1,XF-1),:]) != Zero):
					cols.append(OLD[min(i+1,YF-1),min(j+1,XF-1),:])
				NEW[i,j,0:3] = commoncolors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW)

def FillImageLL0(lat,lon,res):
	print('Filling Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	RF = RSize[res]
	OLD = plt.imread('Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float32)
	for i in range(YF):
		for j in range(XF):
			if(C2I(OLD[i,j]) == 0 and ((i-540)**2+(j-960)**2) < RF**2 ):
				cols = []
				if(C2I(OLD[max(i-1,0),max(j-1,0),:]) != 0):
					cols.append(OLD[max(i-1,0),max(j-1,0),:])
				if(C2I(OLD[max(i-1,0),j,:]) != 0):
					cols.append(OLD[max(i-1,0),j,:])
				if(C2I(OLD[max(i-1,0),min(j+1,XF-1),:]) != 0):
					cols.append(OLD[max(i-1,0),min(j+1,XF-1),:])

				if(C2I(OLD[i,max(j-1,0),:]) != 0):
					cols.append(OLD[i,max(j-1,0),:])
				if(C2I(OLD[i,min(j+1,XF-1),:]) != 0):
					cols.append(OLD[i,min(j+1,XF-1),:])

				if(C2I(OLD[min(i+1,YF-1),max(j-1,0),:]) != 0):
					cols.append(OLD[min(i+1,YF-1),max(j-1,0),:])
				if(C2I(OLD[min(i+1,YF-1),j,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(C2I(OLD[min(i+1,YF-1),min(j+1,XF-1),:]) != 0):
					cols.append(OLD[min(i+1,YF-1),min(j+1,XF-1),:])
				NEW[i,j,0:3] = AvgColors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW)

def FillSTAR(lat,lon):
	res = 2
	print('Filling STAR '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	RF = RSize[res]
	OLD = plt.imread('Assets/STAR/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float16)
	One = 1
	for i in range(YF):
		for j in range(XF):
			if(C2I(OLD[i,j]) == 0 or C2I(OLD[i,j]) == 1):
				cols = []
				if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0 and C2I(OLD[max(i-1,0),(j-1)%XF,:]) != One ):
					cols.append(OLD[max(i-1,0),(j-1)%XF,:])
				if(C2I(OLD[max(i-1,0),j,:]) != 0 and C2I(OLD[max(i-1,0),j%XF,:]) != One ):
					cols.append(OLD[max(i-1,0),j,:])
				if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0 and C2I(OLD[max(i-1,0),(j+1)%XF,:]) != One ):
					cols.append(OLD[max(i-1,0),(j+1)%XF,:])

				if(C2I(OLD[i,(j-1)%XF,:]) != 0 and C2I(OLD[i,(j-1)%XF,:]) != One ):
					cols.append(OLD[i,(j-1)%XF,:])
				if(C2I(OLD[i,(j+1)%XF,:]) != 0 and C2I(OLD[i,(j+1)%XF,:]) != One ):
					cols.append(OLD[i,(j+1)%XF,:])

				if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0 and C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != One ):
					cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
				if(C2I(OLD[min(i+1,YF-1),j,:]) != 0 and C2I(OLD[min(i+1,YF-1),j%XF,:]) != One ):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0 and C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != One ):
					cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
				NEW[i,j,0:3] = AvgColors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='Assets/STAR/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=NEW)

def FillSTAR2(lat,lon):
	res = 2
	print('Filling STAR '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	RF = RSize[res]

	OLD = plt.imread('Assets/STAR'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float16)
	Zero = C2I((0.0,0.0,0.0,0.0))
	One = C2I((0.0,0.0,1.0,0.0))
	for i in range(YF):
		for j in range(XF):
			if(C2I(OLD[i,j]) == Zero):
				cols = []
				if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != Zero and C2I(OLD[max(i-1,0),(j-1)%XF,:]) != One ):
					cols.append(OLD[max(i-1,0),(j-1)%XF,:])
				if(C2I(OLD[max(i-1,0),j,:]) != Zero and C2I(OLD[max(i-1,0),j%XF,:]) != One ):
					cols.append(OLD[max(i-1,0),j,:])
				if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != Zero and C2I(OLD[max(i-1,0),(j+1)%XF,:]) != One ):
					cols.append(OLD[max(i-1,0),(j+1)%XF,:])

				if(C2I(OLD[i,(j-1)%XF,:]) != Zero and C2I(OLD[i,(j-1)%XF,:]) != One ):
					cols.append(OLD[i,(j-1)%XF,:])
				if(C2I(OLD[i,(j+1)%XF,:]) != Zero and C2I(OLD[i,(j+1)%XF,:]) != One ):
					cols.append(OLD[i,(j+1)%XF,:])

				if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != Zero and C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != One ):
					cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
				if(C2I(OLD[min(i+1,YF-1),j,:]) != Zero and C2I(OLD[min(i+1,YF-1),j%XF,:]) != One ):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != Zero and C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != One ):
					cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
				NEW[i,j,0:3] = commoncolors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='IMG/STAR'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=NEW)


#	Enhances Image Sharpness

def EnhanceImage(plane,res):
	print('Enhancing Plane '+str(plane)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	OLD = plt.imread('Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float32)
	for i in range(YF):
#		print(str(i))
		for j in range(XF):
			cols = []
			if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j-1)%XF,:])
			if(C2I(OLD[max(i-1,0),j,:]) != 0):
				cols.append(OLD[max(i-1,0),j,:])
			if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j+1)%XF,:])
			if(C2I(OLD[i,(j-1)%XF,:]) != 0):
				cols.append(OLD[i,(j-1)%XF,:])
			if(C2I(OLD[i,(j+1)%XF,:]) != 0):
				cols.append(OLD[i,(j+1)%XF,:])
			if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
			if(C2I(OLD[min(i+1,YF-1),j,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),j,:])
			if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
			
			NEW[i,j,0:3] = (len(cols)+8)*OLD[i,j,0:3]
			for k in range(len(cols)):
				NEW[i,j,0:3] = NEW[i,j,0:3]-cols[k][0:3]
			NEW[i,j,0:3] = NEW[i,j,0:3]/8.0
			for k in range(3):
				NEW[i,j,k] = max(min(NEW[i,j,k],1.0),0)
	plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png',arr=NEW)

def EnhanceImageLL(lat,lon,res):
	print('Enhancing Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	OLD = plt.imread('Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float32)
	for i in range(YF):
#		print(str(i))
		for j in range(XF):
			cols = []
			if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j-1)%XF,:])
			if(C2I(OLD[max(i-1,0),j,:]) != 0):
				cols.append(OLD[max(i-1,0),j,:])
			if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j+1)%XF,:])
			if(C2I(OLD[i,(j-1)%XF,:]) != 0):
				cols.append(OLD[i,(j-1)%XF,:])
			if(C2I(OLD[i,(j+1)%XF,:]) != 0):
				cols.append(OLD[i,(j+1)%XF,:])
			if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
			if(C2I(OLD[min(i+1,YF-1),j,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),j,:])
			if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
			
			NEW[i,j,0:3] = (len(cols)+8)*OLD[i,j,0:3]
			for k in range(len(cols)):
				NEW[i,j,0:3] = NEW[i,j,0:3]-cols[k][0:3]
			NEW[i,j,0:3] = NEW[i,j,0:3]/8.0
			for k in range(3):
				NEW[i,j,k] = max(min(NEW[i,j,k],1.0),0)
	plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW)

def EnhanceSTARLL(lat,lon):
	res = 2
	print('Enhancing Plane '+str(lat)+'x'+str(lon)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	OLD = plt.imread('Assets/STAR'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float16)
	for i in range(YF):
#		print(str(i))
		for j in range(XF):
			cols = []
			if(C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j-1)%XF,:])
			if(C2I(OLD[max(i-1,0),j,:]) != 0):
				cols.append(OLD[max(i-1,0),j,:])
			if(C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j+1)%XF,:])
			if(C2I(OLD[i,(j-1)%XF,:]) != 0):
				cols.append(OLD[i,(j-1)%XF,:])
			if(C2I(OLD[i,(j+1)%XF,:]) != 0):
				cols.append(OLD[i,(j+1)%XF,:])
			if(C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
			if(C2I(OLD[min(i+1,YF-1),j,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),j,:])
			if(C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
			
			NEW[i,j,0:3] = (len(cols)+8)*OLD[i,j,0:3]
			for k in range(len(cols)):
				NEW[i,j,0:3] = NEW[i,j,0:3]-cols[k][0:3]
			NEW[i,j,0:3] = NEW[i,j,0:3]/8.0
			for k in range(3):
				NEW[i,j,k] = max(min(NEW[i,j,k],1.0),0)
	plt.imsave(fname='Assets/STAR'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=NEW)



#	Bayes Classifier for color-numbers

def commonvalues(colores):
	colors = []
	for c in colores:
		colors.append(c)
	
	dt = np.array((0.0,0.0,0.0,1.0))
	vals = list(set(colors))
	nums = np.zeros(shape=(len(vals)),dtype=np.float16)
	valmax = 1000000
	nummax = 0
	ct = 0
	for v in vals:
		if(v == C2I(dt)):
			ct += 1
		else:
			c0 = 0
			for c in colors:
				if(c == v):
					c0 += 1
			if(c0 > nummax):
				c0 = nummax
				valmax = v
	return valmax

def NearestPlane(pointID):
	point = Point9.getPoint(pointID)
	xyz = point.getXYZPos()
	ptID = point.getID()
	if(ptID < 92):
		return ptID
	else:
		dmin = 200.0
		ID = 0
		for i in range(92):
			pt = Point9.getPoint(i)
			xyztmp = pt.getXYZPos()
			if(dSurf(xyztmp,xyz) <= dmin):
				dmin = dSurf(xyztmp,xyz)
				ID = i
		return ID

def NearestPlane(lat,lon):
	phi = lat*math.pi/180.0
	lam = lon*math.pi/180.0
	x0 = math.cos(phi)*math.cos(lam)
	y0 = math.cos(phi)*math.sin(lam)
	z0 = math.sin(phi)
	xyz = (x0,y0,-1*z0)
	dmin = 2.0
	ID = 0
	for i in range(60,92,1):
		pt = Point9.getPoint(i)
		xyztmp = pt.getXYZPos()
		if(dEuclid(xyztmp,xyz) <= dmin):
			dmin = dEuclid(xyztmp,xyz)
			ID = i
	return ID

#	Builds a full world from list of XYZ points
def BuildFromCRDS(fname):
	from PIL import Image
	Pointlist = list(open(fname+'.txt'))
	IGBP0 = Dataset('sdat_10004_1_20201121_222026926.nc')
	IGBP = IGBP0.variables['Band1'][:]
	REG0 = Image.open('MAPID.png')
	REG = REG0.load()
	RGN0 = Dataset('VSmallF2.nc')
	RGN = RGN0.variables['V'][:]
	HGT0 = Dataset('ZSmallF2.nc')
	HGT = HGT0.variables['Z'][:]
	CLIM0 = Image.open('SCZ6.png')
	CLIM = CLIM0.load()
	SIC0 = Image.open('SIC.png')
	SIC = SIC0.load()
#	oldList = list(open('CRDS6.pts'))
	g = open(fname+'_2.pts','w')
	ct = 0
	for i in Pointlist:
		tmp = i.split(' ')
		x0 = float(tmp[0])
		y0 = float(tmp[1])
		z0 = float(tmp[2])
		REG_XY = toXY((x0,y0,z0),4320,2160)
		REG_col = REG[REG_XY[0],REG_XY[1]]
		REG_val = 256*256*REG_col[0]+256*REG_col[1]+REG_col[2]
		IGBP_val = IGBP[2159-REG_XY[1],REG_XY[0]]
		RGN_val = RGN[REG_XY[0],REG_XY[1]]
		HGT_val = HGT[REG_XY[0],REG_XY[1]]
		CLIM_XY = toXY((x0,y0,z0),720,360)
		CLIM_val = CLIM[CLIM_XY[0],CLIM_XY[1]]
		SIC_XY = toXY((x0,y0,z0),576,361)
		SIC_val = SIC[SIC_XY[0],SIC_XY[1]][0:3]
#		print(SIC_val)
		RGNV = 0
		if(IGBP_val == 0):
			iscoast = False
			for j in NNList[ct].split(' '):
				if(Grid[int(j)].getIGBP() != 0):
					iscoast = True
			if(iscoast):
				IGBP_val = 1
		else:
			IGBP_val += 1

		if(RGN_val >= 10000):
			if( HGT_val < 2500):
				RGNV = 1
			else:
				RGNV = 3
		elif(RGN_val >= 2500):
			RGNV = 2
		elif(RGN_val >= 600):
			RGNV = 1
		else:
			RGNV = 0
		if((IGBP_val == 0 or IGBP_val == 1) and SIC_val == (255,0,0)):
			IGBP_val = 18
#			print('Ice Shelf!')
		CLIMV = 0
		if(CLIM_val == (0,0,255) or CLIM_val ==  (0,255,0)):
			CLIMV = 0
		if(CLIM_val == (127,195,31)):
			CLIMV = 1
		if(CLIM_val == (255,255,0)):
			CLIMV = 3
		if(CLIM_val == (255,0,255)):
			CLIMV = 2
		if(CLIM_val == (255,255,255)):
			CLIMV = 4
#		tmp0 = oldList[(ct-1)%LENS[LEVEL]].split(' ')
		BRDV = '0'
#		print(float(tmp0[2])-x0)
		g.write(str(ct)+' '+str(REG_val)+' '+str(IGBP_val)+' '+str(RGNV)+' '+BRDV+' '+str(x0)+' '+str(y0)+' '+str(z0)+' '+str(CLIMV)+'\n')
		ct += 1

#	Using NN and ER Mappings

def BuildFromCRDS2(fname):
	Pointlist = list(open(fname+'.pts'))
	IGBP0 = Dataset('sdat_10004_1_20201121_222026926.nc')
	IGBP1 = IGBP0.variables['Band1'][:]
	IGBP = np.zeros(shape=(2160,4320),dtype=int)
	for i in range(0,2160):
		IGBP[i,:] = IGBP1[2159-i,:]
	REG = plt.imread('MAPID.png')
	RGN0 = Dataset('VSmallF2.nc')
	RGN1 = RGN0.variables['V'][:]
	RGN = np.zeros(shape=(2160,4320),dtype=int)
	for i in range(0,2160):
		RGN[i,:] = RGN1[:,i]
	HGT0 = Dataset('ZSmallF2.nc')
	HGT1 = HGT0.variables['Z'][:]
	HGT = np.zeros(shape=(2160,4320),dtype=int)
	for i in range(0,2160):
		HGT[i,:] = HGT1[:,i]
	CLIM = plt.imread('SCZ6_16.png')
	SIC = plt.imread('SIC_16.png')
	g = open(fname+'_2.pts','w')
	ct = 0
	
	for i in Pointlist:
		tmp = i.split(' ')
		print(tmp[0])
		pt = int(tmp[0])
		CT = ''+str(pt)
		if(pt < 10000):
			CT = '0'+str(pt)
			if(pt < 1000):
				CT = '00'+str(pt)
				if(pt < 100):
					CT = '000'+str(pt)
					if(pt < 10):
						CT = '0000'+str(pt)
		x0 = float(tmp[5])
		y0 = float(tmp[6])
		z0 = float(tmp[7])
		REG_val = ER_DSCAVG6(REG,list(open('PTLists/'+CT[0]+'/'+CT[1]+'/'+CT+'.txt')))
		IGBP_val = ER_DSCAVG5(IGBP,list(open('PTLists/'+CT[0]+'/'+CT[1]+'/'+CT+'.txt')))
		RGN_val = ER_NUMAVG5(RGN,list(open('PTLists/'+CT[0]+'/'+CT[1]+'/'+CT+'.txt')))
		HGT_val = ER_NUMAVG5(HGT,list(open('PTLists/'+CT[0]+'/'+CT[1]+'/'+CT+'.txt')))
		CLIM_val = ER_DSCAVG6(CLIM,list(open('PTLists/'+CT[0]+'/'+CT[1]+'/'+CT+'.txt')))
		SIC_val = ER_DSCAVG6(SIC,list(open('PTLists/'+CT[0]+'/'+CT[1]+'/'+CT+'.txt')))
#		print(SIC_val)
		RGNV = 0
		if(IGBP_val == 0):
			iscoast = False
			for j in NNList[ct].split(' ')[0:-1]:
				if(Grid[int(j)].getIGBP() != 0):
					iscoast = True
			if(iscoast):
				IGBP_val = 1
		else:
			IGBP_val += 1

		if(RGN_val >= 10000):
			if( HGT_val < 2500):
				RGNV = 1
			else:
				RGNV = 3
		elif(RGN_val >= 2500):
			RGNV = 2
		elif(RGN_val >= 600):
			RGNV = 1
		else:
			RGNV = 0
		if((IGBP_val == 0 or IGBP_val == 1) and SIC_val == C2I((255,0,0))):
			IGBP_val = 18
#			print('Ice Shelf!')
		CLIMV = 0
		if(CLIM_val == C2I((0,0,255)) or CLIM_val ==  C2I((0,255,0,))):
			CLIMV = 0
		if(CLIM_val == C2I((128,128,0))):
			CLIMV = 1
		if(CLIM_val == C2I((255,255,0))):
			CLIMV = 3
		if(CLIM_val == C2I((255,0,255))):
			CLIMV = 2
		if(CLIM_val == C2I((255,255,255))):
			CLIMV = 4
#		tmp0 = oldList[(ct-1)%LENS[LEVEL]].split(' ')
		BRDV = '0'
#		print(float(tmp0[2])-x0)
		g.write(str(ct)+' '+str(REG_val)+' '+str(IGBP_val)+' '+str(RGNV)+' '+BRDV+' '+str(x0)+' '+str(y0)+' '+str(z0)+' '+str(CLIMV)+'\n')
		ct += 1

#BuildFromCRDS('CRDS8')
print('Starting')
def ShowNN(pointID, res):
#	rad = min(1000.0/res,int(YF/2))
	rad = min(2000.0/(2**res),int(YF/2))
	NNList = list(open('NN8.txt'))
	TPSet = set([])
	TPSet.add(pointID)
	for i in range(int(rad)):
		for j in list(TPSet):
			for k in NNList[j].split(' '):
				TPSet.add(int(k))
	tmppts = WriteLocalCRDS(Grid[pointID].getXYZPos(),rad)
	img = np.zeros(shape=(YF,XF,4),dtype='float32')
	img[:,:] = np.array((0,0,0,1))
	skelpix, ringpix, Textures = LoadTextures(res)
	for i in tmppts:
		AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
	plt.imsave(fname='test.png',arr=img)


def GenIMG(ptID, res):
	rad = min(2000.0/(2**res),int(YF/2))
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	tmppts = WriteLocalCRDS2(ptID,min(5*2**(8-res),int(YF/2)))
	img = np.zeros(shape=(YF,XF,4),dtype='float32')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 4):
		imgID = np.zeros(shape=(YF,XF,4),dtype='float32')
		imgID[:,:] = np.array((1,1,1,1))
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
	else:
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
#Z1 = list(open('Z3_8.txt'))

def GenIMG_PTS(lat, lon, res):
	Z1 = list(open('Z'+str(res)+'_8.txt'))
#	H1 = list(open('H'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	g = Planes/Zoom
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 0):
		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile4(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,lat,lon,res)
			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)
		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)


#Z1 = list(open('Z4_8.txt'))
def GenIMG2(lat, lon, res):

	Z1 = list(open('Z'+str(res)+'_8.txt'))

#	H1 = list(open('H'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDSF(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 0):
#		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTileH_3(int(i[0]),i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,lat,lon,res)
#			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)
#		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)

def GenIMG20(lat, lon, res):

	Z1 = list(open('Z'+str(res)+'_8.txt'))

	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
	rad = 500
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
#	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	tmppts = WriteLocalCRDS(ptXYZ,420)
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 0):
		for i in tmppts:
			AssignTileH_3(int(i[0]),i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,lat,lon,res)
		plt.imsave(fname='AssetsLarge/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)
	else:
		for i in tmppts:
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='AssetsLarge/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)


def GenIMG_MOON(lat, lon, res):

	Z1 = list(open('L'+str(res)+'_8.txt'))

#	H1 = list(open('H'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDSF(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 0):
#		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(int(i[0]),i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
#			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='Assets/Moon'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)
#		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)


def GenIMG_MARS(lat, lon, res):

	Z1 = list(open('M'+str(res)+'_8.txt'))

#	H1 = list(open('H'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDSF(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 0):
#		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(int(i[0]),i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
#			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='Assets/Mars'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)
#		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)


def GenIMGT(lat, lon, res):
	Z1 = list(open('T'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 4):
		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
#			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='Topo/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)
#		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='Topo/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)

#Z1 = list(open('H50_8.txt'))

def GenIMGH(lat, lon, res):
#	Z1 = list(open('H'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 4):
		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTileH_3(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,lat,lon,res)
#			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='HGT/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)
#		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTileH(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,lat,lon,res)
		plt.imsave(fname='HGT/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=img)


def GenIMG_STAR(lat, lon,res):
#	res = 2
	Z1 = list(open('Z_STAR'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures(res)
	rad = min(2000.0/(2**res),420)
#	RADS = [450,450,450,450,450,450]
#	rad = RADS[res]
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = -1*(360-lon)%360*math.pi/180.0
#	lam = -1*lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS3(ptXYZ,min(5*2**(8-res),420))
	img = np.zeros(shape=(YF,XF,4),dtype='float16')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 4):
		imgID = np.ones(shape=(YF,XF,4),dtype='float16')
#		imgID[:,:] = np.array((1,1,1,1))
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='STAR'+str(res)+''+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)
		plt.imsave(fname='STAR'+str(res)+''+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
		plt.imsave(fname='Assets/STAR'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)#img[960-175:960+175,960-175:960+175,:])



def GenIMG3(ptID, res):
	YSize = [YF,YF,YF,YF,YF,XF,3200,6400]
	XSize = [XF,XF,XF,XF,XF,XF,3200,6400]
	rad = 450 #min(2000.0/(2**res),int(YF/2))
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
#	tmppts = WriteLocalCRDS2(ptID,min(5*2**(8-res),int(YF/2)))
	tmppts = WriteLocalCRDS2(ptID,rad)
	g = open('Planes/Plane'+str(ptID)+'.txt','w')
	img = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float32')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 0):
		imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float32')
		imgID[:,:] = np.array((1,1,1,1))
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
			AssignTile3(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignID3(i[0],i[1],i[2],imgID,skelpix,res)
			tmpi = g.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
		plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
		plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
		g.close()
	else:
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
			AssignTile3(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
		plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)


def GenIMG4(lat, lon, res):
	rad = min(2000.0/(2**res),int(YF/2),420)
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),int(YF/2),420))
	
#	img = np.zeros(shape=(YF,XF,4),dtype='float32')
#	img[:,:] = np.array((0,0,0,1))

	if(res >= 0):
		imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float32')
		imgID[:,:] = np.array((1,1,1,1))
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignID2(i[0],i[1],i[2],imgID,skelpix,res)
#		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=img)
		plt.imsave(fname='IMG/Zoom'+str(res)+'/Plane'+str(int(lat*10000))+'x'+str(int(lon*10000))+'.png',arr=imgID)
	else:
#		skelpix, ringpix, Textures = LoadTextures(res)
#		for i in tmppts:
#			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
#		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat*10000)+'x'+str(lon*10000)+'.png',arr=img)
		print('wut')

def GenIMG5(ptID,res):
	Z1 = list(open('Z'+str(res)+'_8.txt'))
	skelpix, ringpix, Textures = LoadTextures3(res)
	if(ptID != 80 and ptID != 83):
		RADS = [520,520,520,520,520,520]
#		rad = min(2000.0/(2**res),520)
		rad = RADS[res]
		tmppts = WriteLocalCRDS2(ptID,rad)
		g = open('Planes/Plane'+str(ptID)+'.txt','w')
		img = np.ones(shape=(YSize[res],XSize[res],4),dtype='float16')
		img[:,:] = np.array((0,0,0,1))
		if(res >= 0):
			imgID = np.ones(shape=(YSize[res],XSize[res],4),dtype='float16')
			for i in tmppts:
				AssignTile3(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
				AssignID3(i[0],i[1],i[2],imgID,skelpix,res)
				g.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
			plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
			plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
			g.close()
		else:
			for i in tmppts:
				AssignTile3(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
	elif(ptID == 80):
		lat = -90
		LONS = [54,126,198,342,270]
		LABS = ['60','61','64','62','63']
		RADS = [450,450,450,450,450,450]
		rad = RADS[res]
		ct = 0
		for lon in LONS:
			phi = lat*math.pi/180.0
			lam = lon*math.pi/180.0
			x0 = math.cos(lam)*math.cos(phi)
			y0 = math.sin(lam)*math.cos(phi)
			z0 = -1*math.sin(phi)
			ptXYZ = (x0,y0,z0)
			tmppts = WriteLocalCRDS(ptXYZ,rad)
			g = open('Planes/Plane'+str(ptID)+LABS[ct]+'.txt','w')
			img = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
			img[:,:] = np.array((0,0,0,1))
			if(res >= 0):
				imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
				imgID[:,:] = np.array((1,1,1,1))
				for i in tmppts:
					AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
					AssignID(i[0],i[1],i[2],imgID,skelpix,res)
					g.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
				plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=imgID)
				g.close()
			else:
				for i in tmppts:
					AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
			ct += 1
	elif(ptID == 83):
		lat = 90
		LONS = [18,90,162,234,306]
		LABS = ['74','73','72','70','71']
		RADS = [450,450,450,450,450,450]
		rad = RADS[res]
		ct = 0
		for lon in LONS:
			phi = lat*math.pi/180.0
			lam = lon*math.pi/180.0
			x0 = math.cos(lam)*math.cos(phi)
			y0 = math.sin(lam)*math.cos(phi)
			z0 = -1*math.sin(phi)
			ptXYZ = (x0,y0,z0)
			tmppts = WriteLocalCRDS(ptXYZ,rad)
			g = open('Planes/Plane'+str(ptID)+LABS[ct]+'.txt','w')
			img = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
			img[:,:] = np.array((0,0,0,1))
			if(res >= 0):
				imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
				imgID[:,:] = np.array((1,1,1,1))
				for i in tmppts:
					AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
					AssignID(i[0],i[1],i[2],imgID,skelpix,res)
					g.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
				plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
			else:
				for i in tmppts:
					AssignTile2(i[0],i[1],i[2],img,Z1[int(i[0])],skelpix,ringpix,res)
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
			ct += 1

#	Fastest

def GenIMG6(ptID,res):
	skelpix, ringpix, Textures = LoadTextures3(res)
	if(ptID != 80 and ptID != 83):
		RADS = [520,520,520,520,520,520]
#		rad = min(2000.0/(2**res),520)
		rad = RADS[res]
#		tmppts = WriteLocalCRDS2(ptID,rad)
		g = open('Planes/Plane'+str(ptID)+'.txt')
		img = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
		img[:,:] = np.array((0,0,0,1))
		if(res >= 0):
			imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
			imgID[:,:] = np.array((1,1,1,1))
#			skelpix, ringpix, Textures = LoadTextures(res)
			for I in g:
				i = I.split(' ')
				AssignTile3(int(i[0]),float(i[1]),float(i[2]),img,Textures,skelpix,ringpix,res)
				AssignID3(int(i[0]),float(i[1]),float(i[2]),imgID,skelpix,res)
#				g.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
			plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
			plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
#			g.close()
		else:
			skelpix, ringpix, Textures = LoadTextures(res)
			for I in g:
				i = I.split(' ')
				AssignTile3(int(i[0]),float(i[1]),float(i[2]),img,Textures,skelpix,ringpix,res)
			plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=img)
	elif(ptID == 80):
		lat = -90
		LONS = [54,126,198,342,270]
		LABS = ['60','61','64','62','63']
		RADS = [450,450,450,450,450,450]
		rad = RADS[res]
		ct = 0
		for lon in LONS:
			phi = lat*math.pi/180.0
			lam = lon*math.pi/180.0
			x0 = math.cos(lam)*math.cos(phi)
			y0 = math.sin(lam)*math.cos(phi)
			z0 = -1*math.sin(phi)
			ptXYZ = (x0,y0,z0)
#			tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),rad))
			g = open('Planes/Plane'+str(ptID)+LABS[ct]+'.txt')
			img = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
			img[:,:] = np.array((0,0,0,1))
			if(res >= 0):
				imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
				imgID[:,:] = np.array((1,1,1,1))
#				skelpix, ringpix, Textures = LoadTextures(res)
				for I in g:
					i = I.split(' ')
					AssignTile(int(i[0]),float(i[1]),float(i[2]),img,Textures,skelpix,ringpix,res)
					AssignID(int(i[0]),float(i[1]),float(i[2]),imgID,skelpix,res)
#					g.write(str(i[0])+' '+str(i[1])+' '+str(i[2])+'\n')
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
				plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
#				g.close()
			else:
				skelpix, ringpix, Textures = LoadTextures(res)
				for I in g:
					i = I.split(' ')
					AssignTile(int(i[0]),float(i[1]),float(i[2]),img,Textures,skelpix,ringpix,res)
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
			ct += 1
	elif(ptID == 83):
		lat = 90
		LONS = [18,90,162,234,306]
		LABS = ['74','73','72','70','71']
		RADS = [450,450,450,450,450,450]
		rad = RADS[res]
		ct = 0
		for lon in LONS:
			phi = lat*math.pi/180.0
			lam = lon*math.pi/180.0
			x0 = math.cos(lam)*math.cos(phi)
			y0 = math.sin(lam)*math.cos(phi)
			z0 = -1*math.sin(phi)
			ptXYZ = (x0,y0,z0)
#			tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),rad))
			g = open('Planes/Plane'+str(ptID)+LABS[ct]+'.txt')
			img = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
			img[:,:] = np.array((0,0,0,1))
			if(res >= 0):
				imgID = np.zeros(shape=(YSize[res],XSize[res],4),dtype='float16')
				imgID[:,:] = np.array((1,1,1,1))
#				skelpix, ringpix, Textures = LoadTextures(res)
				for I in g:
					i = I.split(' ')
					AssignTile(int(i[0]),float(i[1]),float(i[2]),img,Textures,skelpix,ringpix,res)
					AssignID(int(i[0]),float(i[1]),float(i[2]),imgID,skelpix,res)
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
				plt.imsave(fname='IMG32/Zoom'+str(res)+'/Plane'+str(ptID)+'.png',arr=imgID)
			else:
				skelpix, ringpix, Textures = LoadTextures(res)
				for I in g:
					i = I.split(' ')
					AssignTile(int(i[0]),float(i[1]),float(i[2]),img,Textures,skelpix,ringpix,res)
				plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(ptID)+LABS[ct]+'.png',arr=img)
			ct += 1

#	For selected Plane & points only

def GenIMG7(ptIDs,planes,res):
	skelpix, ringpix, Textures = LoadTextures(res)
	for plane in planes:
		img = plt.imread('Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png')
		RADS = [520,520,520,520,520]
#		rad = min(2000.0/(2**res),520)
		rad = RADS[res]
		g = list(open('Planes/P2_'+str(plane)+'.txt'))
		for ptID in ptIDs:
			tmp = g[ptID].split(' ')
			AssignTile3(ptID,float(tmp[0]),float(tmp[1]),img,Textures,skelpix,ringpix,res)
		plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png',arr=img)
#	Generates New lat lon IMG from 

def GenOnEachPlane(zoom):
	YF = YSize[zoom]
	XF = XSize[zoom]
	G = list(open('XY2PT.txt'))
	ct = 0 
	for data in G:
		tmp = data.split(' ')
		plane = Builder.NearestPlane(int(tmp[0]),int(tmp[1]))
		strplane = str(plane)
		if(plane == 80):
			if(int(tmp[1]) <= 18 or int(tmp[1]) > 306):
				strplane = '8062'
			elif(int(tmp[1]) <= 90 and int(tmp[1]) > 18):
				strplane = '8060'
			elif(int(tmp[1]) <= 162 and int(tmp[1]) > 90):
				strplane = '8061'
			elif(int(tmp[1]) <= 234 and int(tmp[1]) > 162):
				strplane = '8064'
			elif(int(tmp[1]) <= 306 and int(tmp[1]) > 234):
				strplane = '8063'
		if(plane == 83):
			if(int(tmp[1]) <= 54 or int(tmp[1]) > 342):
				strplane = '8374'
			elif(int(tmp[1]) <= 126 and int(tmp[1]) > 54):
				strplane = '8373'
			elif(int(tmp[1]) <= 198 and int(tmp[1]) > 126):
				strplane = '8372'
			elif(int(tmp[1]) <= 270 and int(tmp[1]) > 198):
				strplane = '8370'
			elif(int(tmp[1]) <= 342 and int(tmp[1]) > 270):
				strplane = '8371'
#		print(strplane)
		imgpix = plt.imread('Assets32/Zoom'+str(zoom)+'/Plane'+strplane+'.png')
		imgdat = list(open('Planes/P2_'+strplane+'.txt'))
		ptid = int(tmp[2])
		center = imgdat[ptid].split(' ')
		x0 = int(float(center[0])*(2**zoom+1))+int(XF/2)
		y0 = int(float(center[1])*(2**zoom+1))+int(YF/2)
#		print(str(y0-499)+' to '+str(y0+450))
#		print(str(x0-799)+' to '+str(x0+800))
		imgnew = imgpix[y0-449:y0+450,x0-799:x0+800]
		plt.imsave(fname='AssetsF/Zoom'+str(zoom)+'/Plane'+tmp[0]+'x'+tmp[1]+'.png',arr=imgnew)
		ct += 1

def GenOnEachPlane2(lat, lon, zoom):
	YF = YSize[zoom]
	XF = XSize[zoom]
	G = list(open('XY2PT.txt'))
	ct = 0 
	ptat = 1000000000
	for data in G:
		tmp = data.split(' ')
		if(int(tmp[0]) == int(lat) and int(tmp[1]) == int(lon)):
			ptat = ct
		ct += 1
	ct = 0
	if(ptat != 1000000000):
		data = G[ptat]
		tmp = data.split(' ')
		plane = Builder.NearestPlane(int(tmp[0]),int(tmp[1]))
		strplane = str(plane)
		if(plane == 80):
			if(int(tmp[1]) <= 18 or int(tmp[1]) > 306):
				strplane = '8062'
			elif(int(tmp[1]) <= 90 and int(tmp[1]) > 18):
				strplane = '8060'
			elif(int(tmp[1]) <= 162 and int(tmp[1]) > 90):
				strplane = '8061'
			elif(int(tmp[1]) <= 234 and int(tmp[1]) > 162):
				strplane = '8064'
			elif(int(tmp[1]) <= 306 and int(tmp[1]) > 234):
				strplane = '8063'
		if(plane == 83):
			if(int(tmp[1]) <= 54 or int(tmp[1]) > 342):
				strplane = '8374'
			elif(int(tmp[1]) <= 126 and int(tmp[1]) > 54):
				strplane = '8373'
			elif(int(tmp[1]) <= 198 and int(tmp[1]) > 126):
				strplane = '8372'
			elif(int(tmp[1]) <= 270 and int(tmp[1]) > 198):
				strplane = '8370'
			elif(int(tmp[1]) <= 342 and int(tmp[1]) > 270):
				strplane = '8371'
#		print(strplane)
		imgpix = plt.imread('Assets32/Zoom'+str(zoom)+'/Plane'+strplane+'.png')
		imgdat = list(open('Planes/P2_'+strplane+'.txt'))
		ptid = int(tmp[2])
		center = imgdat[ptid].split(' ')
		x0 = int(float(center[0])*(2**zoom+1))+int(XF/2)
		y0 = int(float(center[1])*(2**zoom+1))+int(YF/2)
#		print(str(y0-499)+' to '+str(y0+450))
#		print(str(x0-799)+' to '+str(x0+800))
		imgnew = imgpix[y0-449:y0+450,x0-799:x0+800]
		plt.imsave(fname='AssetsF/Zoom'+str(zoom)+'/Plane'+tmp[0]+'x'+tmp[1]+'.png',arr=imgnew)
		ct += 1


#
def SetGrid(PointID,aff,igbp,rgn,isborder,CLIM):
	p = Point9.getPoint(PointID)
	p.setProvince(aff)
	p.setIGBP(igbp)
	p.setRGN(rgn)
	p.setBoarder(isboarder)
	p.setCLIM(CLIM)

def SetPRV(PointID,aff):
	p = Point9.getPoint(PointID)
	p.setProvince(aff)

def SetIGBP(PointID,igbp):
	p = Point9.getPoint(PointID)
	p.setIGBP(igbp)

def SetHGT(PointID,hgt):
	p = Point9.getPoint(PointID)
	p.setHGT(hgt)

def SetPOPL(PointID,popl):
	p = Point9.getPoint(PointID)
	p.setPOPL(popl)

def SetRGN(PointID, rgn):
	p = Point9.getPoint(PointID)
	p.setRGN(rgn)
	
def SetBoarder(PointID, brd):
	p = Point9.getPoint(PointID)
	p.setBoarder(brd)

def SetCLIM(PointID,CLIM):
	p = Point9.getPoint(PointID)
	p.setCLIM(CLIM)

def SaveGrid():
	Point9.SavePoints()

def QuickIMG(lat, lon, res):
	rad = min(2000.0/(2**res),100)
#	NNList = list(open('NN8.txt'))
#	TPSet = set([])
#	TPSet.add(ptID)
#	for i in range(int(rad)):
#		for j in list(TPSet):
#			for k in NNList[j].split(' '):
#				TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),rad)
	phi = lat*math.pi/180.0
	lam = lon*math.pi/180.0
	x0 = math.cos(lam)*math.cos(phi)
	y0 = math.sin(lam)*math.cos(phi)
	z0 = -1*math.sin(phi)
	ptXYZ = (x0,y0,z0)
	tmppts = WriteLocalCRDS(ptXYZ,min(5*2**(8-res),int(YF/2)))
	img = np.zeros(shape=(YF,XF,4),dtype='float32')
	img[:,:] = np.array((0,0,0,1))
	if(res >= 4):
		imgID = np.zeros(shape=(YF,XF,4),dtype='float32')
		imgID[:,:] = np.array((1,1,1,1))
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
			AssignID(i[0],i[1],i[2],imgID,skelpix,res)
		plt.imsave(fname='tmp.png',arr=img)
		plt.imsave(fname='tmpimg.png',arr=imgID)
	else:
		skelpix, ringpix, Textures = LoadTextures(res)
		for i in tmppts:
			AssignTile(i[0],i[1],i[2],img,Textures,skelpix,ringpix,res)
		plt.imsave(fname='Assets/Zoom'+str(res)+'/Plane'+str(lat)+'x'+str(lon)+'.png',arr=img)

def WriteLocalCRDS(point, radius):
	nearbyNUM = NNIDW(point,radius)
	x0 = point[0]
	y0 = point[1]
	z0 = point[2]
	if(z0 == 0):
		z0 = 0.0001
	North = (0.0,0.0,1000)
	phi = math.acos(z0)
	the = phi-math.pi/2.0
	if(math.cos(the) == 0):
		North = (0.0,1.0,0.0)
	else:
		North = (0.0,0.0,1.0/math.cos(the))
	ey0 = (0.0,0.0,1.0)
	if(math.sin(the) != 0 and math.cos(the) != 0):
		ey0 = cross(North,point)
	elif(the == 0):
		ey0 = (0.0,1.0,0.0)
	elif(the == math.pi/2.0):
		ey0 = (x0,y0,1.0)
	ey = (0.0,1.0,0.0)
	if(mod(ey0) != 0):
		ey = (ey0[0]/mod(ey0),ey0[1]/mod(ey0),ey0[2]/mod(ey0))
	ez = (x0,y0,z0)
	ex = cross(ey,ez)
	nearbycrds = []
	ct = 0
	for nn0 in nearbyNUM:
		nn = Grid[nn0].getXYZPos()
		nearbycrds.append((nn0,dot(nn,ey)/distancescale,-1*dot(nn,ex)/distancescale))#,-1*dot(nn,ez)/distancescale)
	return nearbycrds

def WriteLocalCRDSF(point, radius):
	nearbyNUM = NNIDW(point,radius)

	x0 = point[0]
	y0 = point[1]
	z0 = point[2]
	
	CRDS = np.load('CRDS.npy')
	X0 = CRDS[:,0]
	Y0 = CRDS[:,1]
	Z0 = CRDS[:,2]
	
	PHI = np.arcsin(Z0[:])
	THE = PHI[:]-math.pi/2.0
	North = np.array((0.0,0.0,1.0))

	if(z0 == 0):
		z0 = 0.0001
	North = (0.0,0.0,1000)
	phi = math.acos(z0)
	lat = phi*180/math.pi
	the = phi-math.pi/2.0
	if(math.cos(the) == 0):
		North = (0.0,1.0,0.0)
	else:
		North = (0.0,0.0,1.0/math.cos(the))
	ey0 = (0.0,0.0,1.0)
	if(math.sin(the) != 0 and math.cos(the) != 0):
		ey0 = cross(North,point)
	elif(the == 0):
		ey0 = (0.0,1.0,0.0)
	elif(the == math.pi/2.0):
		ey0 = (x0,y0,1.0)
	ey = (0.0,1.0,0.0)
	if(mod(ey0) != 0):
		ey = (ey0[0]/mod(ey0),ey0[1]/mod(ey0),ey0[2]/mod(ey0))
	ez = (x0,y0,z0)
	ex = cross(ey,ez)
	nearbycrds = np.zeros(shape=(len(nearbyNUM),4),dtype=np.float16)
	nearbycrds[:,0] = nearbyNUM[:]
	nearbycrds[:,1] = dot_arr1(CRDS[nearbyNUM[:]],ey)
	nearbycrds[:,2] = dot_arr1(-1*CRDS[nearbyNUM[:]],ex)
	nearbycrds[:,3] = dot_arr1(-1*CRDS[nearbyNUM[:]],ez)
#	np.save('Planes/NPY/Plane'+str(int(lat*10000))
#	nearbycrds = []
#	ct = 0
#	for nn0 in nearbyNUM:
#		nn = Grid[nn0].getXYZPos()
#		nearbycrds.append((nn0,dot(nn,ey)/distancescale,-1*dot(nn,ex)/distancescale))#,-1*dot(nn,ez)/distancescale)
	return nearbycrds

def WriteLocalCRDS2(ptID, radius):
	point = Grid[ptID].getXYZPos()
	TPSet = set([])#NNID(point,radius)
	TPSet.add(ptID)
	for i in range(int(radius)):
		for j in list(TPSet):
			for k in NNList[j].split(' '):
				if(k != '\n'):
					TPSet.add(int(k))
#	tmppts = WriteLocalCRDS(Grid[ptID].getXYZPos(),radius)
	nearbyNUM = list(TPSet)
	x0 = point[0]
	y0 = point[1]
	z0 = point[2]
	if(z0 == 0):
		z0 = 0.0001
	North = (0.0,0.0,1000)
	phi = math.acos(z0)
	the = phi-math.pi/2.0
	if(math.cos(the) == 0):
		North = (0.0,1.0,0.0)
	else:
		North = (0.0,0.0,1.0/math.cos(the))
	ey0 = (0.0,0.0,1.0)
	if(math.sin(the) != 0 and math.cos(the) != 0):
		ey0 = cross(North,point)
	elif(the == 0):
		ey0 = (0.0,1.0,0.0)
	elif(the == math.pi/2.0):
		ey0 = (x0,y0,1.0)
	ey = (0.0,1.0,0.0)
	if(mod(ey0) != 0):
		ey = (ey0[0]/mod(ey0),ey0[1]/mod(ey0),ey0[2]/mod(ey0))
	ez = (x0,y0,z0)
	ex = cross(ey,ez)
	nearbycrds = []
	ct = 0
	for nn0 in nearbyNUM:
		nn = Grid[nn0].getXYZPos()
		nearbycrds.append((nn0,dot(nn,ey)/distancescale,-1*dot(nn,ex)/distancescale))
	return nearbycrds

def WriteLocalCRDS3(point, radius):
	x0 = point[0]
	y0 = point[1]
	z0 = point[2]
	nearbyNUM = NNIDW( (-1*x0,-1*y0,-1*z0) ,radius)
	if(z0 == 0):
		z0 = 0.0001
	North = (0.0,0.0,1000)
	phi = math.acos(z0)
	the = phi-math.pi/2.0
	if(math.cos(the) == 0):
		North = (0.0,1.0,0.0)
	else:
		North = (0.0,0.0,1.0/math.cos(the))
	ey0 = (0.0,0.0,1.0)
	if(math.sin(the) != 0 and math.cos(the) != 0):
		ey0 = cross(North,point)
	elif(the == 0):
		ey0 = (0.0,1.0,0.0)
	elif(the == math.pi/2.0):
		ey0 = (x0,y0,1.0)
	ey = (0.0,1.0,0.0)
	if(mod(ey0) != 0):
		ey = (ey0[0]/mod(ey0),ey0[1]/mod(ey0),ey0[2]/mod(ey0))
	ez = (x0,y0,z0)
	ex = cross(ey,ez)
	nearbycrds = []
	ct = 0
	for nn0 in nearbyNUM:
		nn = Grid[nn0].getXYZPos()		
		nearbycrds.append((nn0,dot(nn,ey)/distancescale,-1*dot(nn,ex)/distancescale))
#		nearbycrds.append( ( nn0 , math.atan( -1*dot(nn,ey)/dot(nn,ez) )/distancescale , math.atan( dot(nn,ex)/dot(nn,ez) )/distancescale ) )
	return nearbycrds

def WriteLocalCRDSPrint():
	Ex = np.zeros(shape=(1474562,3),dtype=np.float16)
	Ey = np.zeros(shape=(1474562,3),dtype=np.float16)
	Ez = np.zeros(shape=(1474562,3),dtype=np.float16)

	for ptID in range(1474562):
		point = Grid[ptID].getXYZPos()
		x0 = point[0]
		y0 = point[1]
		z0 = point[2]
		if(z0 == 0):
			z0 = 0.0001
		North = (0.0,0.0,1000)
		phi = math.acos(z0)
		the = phi-math.pi/2.0
		if(math.cos(the) == 0):
			North = (0.0,1.0,0.0)
		else:
			North = (0.0,0.0,1.0/math.cos(the))
		ey0 = (0.0,0.0,1.0)
		if(math.sin(the) != 0 and math.cos(the) != 0):
			ey0 = cross(North,point)
		elif(the == 0):
			ey0 = (0.0,1.0,0.0)
		elif(the == math.pi/2.0):
			ey0 = (x0,y0,1.0)
		ey = (0.0,1.0,0.0)
		if(mod(ey0) != 0):
			ey = (ey0[0]/mod(ey0),ey0[1]/mod(ey0),ey0[2]/mod(ey0))
		ez = (x0,y0,z0)
		ex = cross(ey,ez)
		Ex[ptID,:] = np.array((ey[0],ey[1],ey[2]))
		Ey[ptID,:] = np.array((-1*ex[0],-1*ex[1],-1*ex[2]))
		Ez[ptID,:] = np.array((ez[0],ez[1],ez[2]))
	np.save('Ex8.npy',arr=Ex)
	np.save('Ey8.npy',arr=Ey)
	np.save('Ez8.npy',arr=Ez)

	

def LoadTexturesF(res):
	img0 = plt.imread('./imgcirc/Textures256.png')
	img = img0						# Need to resize image to res
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	img = plt.imread('./imgcirc/Textures'+str(imgsize)+'.png')
	imgarr = np.zeros(shape=(5,4,19,imgsize,imgsize,4),dtype='float32')
	imgskel = img[imgsize*0:imgsize*0+imgsize,19*imgsize:19*imgsize+imgsize]
	plt.imsave(fname='testskl.png',arr=imgskel)
	imgring = img[imgsize*1:imgsize*1+imgsize,19*imgsize:19*imgsize+imgsize]
	ct = 0
	for i in range(5):
		for j in range(3):
			for k in range(18):
				imgarr[i,j,k,:,:] = img[imgsize*ct:imgsize*ct+imgsize,k*imgsize:k*imgsize+imgsize]
			imgarr[i,j,18,:,:] = img[imgsize*14:imgsize*14+imgsize,1*imgsize:1*imgsize+imgsize]
			ct += 1
		for k in range(18):
			imgarr[i,3,k,:,:] = img[imgsize*ct:imgsize*ct+imgsize,18*imgsize:18*imgsize+imgsize]
			imgarr[i,3,18,:,:] = img[imgsize*14:imgsize*14+imgsize,1*imgsize:1*imgsize+imgsize]
	return imgskel,imgring,imgarr
	plt.imsave(fname='test.png',arr=imgskel)

def LoadTextures(res):
	img0 = plt.imread('./imgcirc/Textures256.png')
	img = img0						# Need to resize image to res
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	img = plt.imread('./imgcirc/Textures'+str(imgsize)+'.png')
	imgarr = np.zeros(shape=(5,4,19,imgsize,imgsize,4),dtype='float32')
	imgskel = img[imgsize*0:imgsize*0+imgsize,19*imgsize:19*imgsize+imgsize]
	plt.imsave(fname='testskl.png',arr=imgskel)
	imgring = img[imgsize*1:imgsize*1+imgsize,19*imgsize:19*imgsize+imgsize]
	ct = 0
	for i in range(5):
		for j in range(3):
			for k in range(18):
				imgarr[i,j,k,:,:] = img[imgsize*ct:imgsize*ct+imgsize,k*imgsize:k*imgsize+imgsize]
			imgarr[i,j,18,:,:] = img[imgsize*14:imgsize*14+imgsize,1*imgsize:1*imgsize+imgsize]
			ct += 1
		for k in range(18):
			imgarr[i,3,k,:,:] = img[imgsize*ct:imgsize*ct+imgsize,18*imgsize:18*imgsize+imgsize]
			imgarr[i,3,18,:,:] = img[imgsize*14:imgsize*14+imgsize,1*imgsize:1*imgsize+imgsize]
	return imgskel,imgring,imgarr
	plt.imsave(fname='test.png',arr=imgskel)

def LoadTextures2(res):
	ClimateCols = list(open('ClimateCols.txt'))
	IGBPCols = range(0,19,1)#list(open('IGBPCols.txt'))
	img0 = plt.imread('./imgcirc/Textures256.png')
	img = img0						# Need to resize image to res
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	img = plt.imread('./imgcirc/Textures'+str(imgsize)+'.png')
	imgarr = np.zeros(shape=(len(ClimateCols),4,len(IGBPCols),imgsize,imgsize,4),dtype='float32')
	imgskel = np.zeros(shape=(2**res,2**res,4),dtype=np.float16)
	plt.imsave(fname='testskl.png',arr=imgskel)
	imgring = np.zeros(shape=(2**res,2**res,4),dtype=np.float16)
	ct = 0
	for i in range(len(ClimateCols)):
		a = plt.imread('imgcirc/BM0/Zoom'+str(res)+'/Clim'+str(i)+'.png')
		for j in range(3):
			for k in range(len(IGBPCols)):
				for I in range(imgsize):
					for J in range(imgsize):
						imgarr[i,j,k,I,J,:] = a[imgsize*j+I,imgsize*k+J]
#						if(k == 0):
#							imgarr[i,j,k,I,J,:] = np.array((10/255.0,10/255.0,51/255.0,1.0))
#						elif(k == 1):
#							imgarr[i,j,k,I,J,:] = np.array((40/255.0,80/255.0,125/255.0,1.0))
#						else:
#							imgarr[i,j,k,I,J,:] = a[imgsize*j+I,imgsize*k+J]
	return imgskel,imgring,imgarr
	plt.imsave(fname='test.png',arr=imgskel)

def LoadTextures3(res):
	ClimateCols = list(open('ClimateCols.txt'))
	IGBPCols = range(0,19,1)#list(open('IGBPCols.txt'))
	RGNCols = list(open('RGNZones.txt'))
#	img0 = plt.imread('./imgcirc/Textures256.png')
#	img = img0						# Need to resize image to res
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
#	img = plt.imread('./imgcirc/Textures'+str(imgsize)+'.png')
	imgarr = np.ones(shape=(len(IGBPCols),len(RGNCols),len(ClimateCols),imgsize,imgsize,4),dtype='float32')
	imgskel = np.zeros(shape=(2**res,2**res,4),dtype=np.float16)
	plt.imsave(fname='testskl.png',arr=imgskel)
	imgring = np.zeros(shape=(2**res,2**res,4),dtype=np.float16)
	ct = 0
	a = plt.imread('imgcirc/NewTextures'+str(2**res)+'.png')

#	for i in range(len(IGBPTypes)):
#		for j in range(len(RGNTypes)):
#			for k in range(len(CLIMTypes)):
#				for I in range(2**res):
#					for J in range(2**res):
#						imgarr[i,j,k,I,J,:] = a[i*2**res+I,((len(RGNTypes))*k+j)*2**res+J]
#						if(k == 0):
#							imgarr[i,j,k,I,J,:] = np.array((10/255.0,10/255.0,51/255.0,1.0))
#						elif(k == 1):
#							imgarr[i,j,k,I,J,:] = np.array((10/255.0,10/255.0,51/255.0,1.0))
#						else:
#							imgarr[i,j,k,I,J,:] = a[i*2**res+I,((len(RGNTypes))*k+j)*2**res+J]

	for i in range(len(ClimateCols)):
		for j in range(len(RGNCols)):
			for k in range(len(IGBPCols)):
				for I in range(imgsize):
					for J in range(imgsize):
						imgarr[k,j,i,I,J,:] = a[imgsize*k+I,imgsize*(i*len(RGNCols)+j)+J,:]
						if(k == 0):
							imgarr[k,j,i,I,J,:] = np.array((10/255.0,10/255.0,51/255.0,1.0))
						elif(k == 1):
							imgarr[k,j,i,I,J,:] = np.array((10/255.0,10/255.0,51/255.0,1.0))
						else:
							imgarr[k,j,i,I,J,:] = a[imgsize*k+I,imgsize*(i*len(RGNCols)+j)+J,:]
	return imgskel,imgring,imgarr
	plt.imsave(fname='test.png',arr=imgskel)

#LoadTextures(8)
def AssignTile(pid,x1,y1,img0,Textures,skel,msk,res):
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	for i in range(imgsize):
		for j in range(imgsize):
			if((res == 0 or res == 1 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
				img0[y0+j+1,x0+i+1,:] = Textures[IGBP,RGN,CLIM,j,i,:]
	return img0

def AssignTile2(pid,x1,y1,img0,Textures,skel,msk,res):
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	strdat = []
	Order = (0,2,1,3)
	if(res >  0):
		strdat = Textures.split(' ')
		if(res == 1):
			Order = (0,2,1,3)
		if(res == 2):
			Order = (0,6,15,2,7,8,5,14,9,4,13,11,1,10,12,3)
			Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)		# Standard Asset Ordering
			Order = (0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15)		# Sandard  Star  Ordering
		if(res == 3):
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63)
		if(res == 4):
			Order = list(range(256))
	ct = 0
	for i in range(imgsize):
		for j in range(imgsize):
			if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
				if(res == 0):
					img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
				if(res == 1):
					img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
				if(res == 2):
					img0[y0+j+1,x0+i+1,:] = INT2Color(min(int(strdat[Order[ct]])+1,C2I( (1.0,1.0,1.0) ) ) )
#					img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]] ) )
					ct += 1
				if(res == 3):
					img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
				if(res == 4):
					img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
						
	return img0

def AssignTile3(pid,x1,y1,img0,Textures,skel0,msk,res):
	skel = 1.0*skel0
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	if(res == 2):
		skel[0,1,0:3] = 1-skel0[0,1,0:3]
		skel[0,2,0:3] = 1-skel0[0,2,0:3]
		skel[1,1,0:3] = 1-skel0[1,1,0:3]
		skel[1,2,0:3] = 1-skel0[1,2,0:3]
		skel[2,1,0:3] = 1-skel0[2,1,0:3]
		skel[2,2,0:3] = 1-skel0[2,2,0:3]
		skel[3,1,0:3] = 1-skel0[3,1,0:3]
		skel[3,2,0:3] = 1-skel0[3,2,0:3]
		skel[1,0,0:3] = 1-skel0[1,0,0:3]
		skel[2,0,0:3] = 1-skel0[2,0,0:3]
		skel[1,3,0:3] = 1-skel0[1,3,0:3]
		skel[2,3,0:3] = 1-skel0[2,3,0:3]

	for i in range(imgsize):
		for j in range(imgsize):
			if((res == 0 or res == 1 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
				img0[y0+j+1,x0+i+1,:] = Textures[IGBP,RGN,CLIM,j,i,:]
	return img0

def AssignTileH(pid,x1,y1,img0,Textures,skel,msk,lat,lon,res):
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	strdat = []
	Order = (0,2,1,3)
	if(res >  0):
		strdat = Textures.split(' ')
		if(res == 1):
			Order = (0,2,1,3)
		if(res == 2):
			Order = (0,6,15,2,7,8,5,14,9,4,13,11,1,10,12,3)
			Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)		# Standard Asset Ordering
			Order = (0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15)		# Sandard  Star  Ordering
		if(res == 3):
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63)
		if(res == 4):
			Order = list(range(256))
		if(res == 5):
			Order = list(range(1024))
	ct = 0
	lonat = Grid[pid].getLonDeg()
	if(lat >= -50 and lat <= 50):
		for i in range(imgsize):
			for j in range(imgsize):
				if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
					if(res == 0):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
					if(res == 1):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 2):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 3):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 4):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 5):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))

				ct += 1
	elif(lat < -50):
		if( (lonat-lon)%360 <= 45 or (lonat-lon)%360 > 315 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
		elif( (lonat-lon)%360 <= 135 and (lonat-lon)%360 > 45 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1				
	
		elif( (lonat-lon)%360 <= 225 and (lonat-lon)%360 > 135 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

	elif(lat > 50):
		if( (lonat-lon)%360 <= 45 or (lonat-lon)%360 > 315 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
		elif( (lonat-lon)%360 <= 135 and (lonat-lon)%360 > 45 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
	
		elif( (lonat-lon)%360 <= 225 and (lonat-lon)%360 > 135 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1		
		
	
						
	return img0
#	for i in range(imgsize):
#		for j in range(imgsize):
#			if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
#				if(res == 0):
#					img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
#				if(res == 1):
#					img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
#				if(res == 2):
#					img0[y0+j+1,x0+i+1,:] = INT2Color(min(int(strdat[Order[ct]])+1,C2I( (1.0,1.0,1.0) ) ) )
#				if(res == 3):
#					img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
#				if(res == 4):
#					img0[y0+j+1,x0+i+1,:] = ColorMapH(int(strdat[Order[ct]]))
#			ct += 1
					
						
#	return img0


def AssignTileH_3(pid,x1,y1,img0,Textures,skel,msk,lat,lon,res):
#	print(pid)
#	print(x1)
#	print(y1)
	OR225 = list(open('imgcirc/Turn16/O'+str(res)+'R225_2.txt'))
	OR45 = list(open('imgcirc/O'+str(res)+'R45.txt'))
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	strdat = []
	Order = (0,2,1,3)
	if(res >  0):
		strdat = Textures.split(' ')
		if(res == 1):
			Order = (0,2,1,3)
		if(res == 2):
			Order = (0,6,15,2,7,8,5,14,9,4,13,11,1,10,12,3)
			Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)		# Standard Asset Ordering
			Order = (0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15)		# Sandard  Star  Ordering
		if(res == 3):
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63)
		if(res == 4):
			Order = list(range(256))
		if(res == 5):
			Order = list(range(1024))
	ct = 0
	lonat = Grid[pid].getLonDeg()
	if(lat >= -10 and lat <= 10):
		for i in range(imgsize):
			for j in range(imgsize):
				if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
					if(res == 0):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
					if(res == 1):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 2):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 3):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 4):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 5):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))

				ct += 1
	elif(lat < -10):

#		NORTH
		if( (lonat-lon)%360 <= 11.25 or (lonat-lon)%360 > 348.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
#		NORTHNORTHEAST
		elif( (lonat-lon)%360 <= 33.75 and (lonat-lon)%360 > 11.25 ):

			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR225[Order[imgsize*i+j]])/imgsize)
					J = int(OR225[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5)  and (skel[I,J,0] != skel[0,0,0]) and (x0+I+1) > 0 and (x0+I+1) < XSize[res] and (y0+J+1) > 0 and (y0+J+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1				




#		NORTHEAST
		elif( (lonat-lon)%360 <= 56.25 and (lonat-lon)%360 > 33.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5)  and (skel[I,J,0] != skel[0,0,0]) and (x0+I+1) > 0 and (x0+I+1) < XSize[res] and (y0+J+1) > 0 and (y0+J+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+J+1,min(x0+I+1,1919),:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1				

#		EASTNORTHEAST	
		elif( (lonat-lon)%360 <= 78.75 and (lonat-lon)%360 > 56.25 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[int(OR225[Order[ct]])]))
					ct += 1

#		EAST	
		elif( (lonat-lon)%360 <= 101.25 and (lonat-lon)%360 > 78.75):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
#		EASTSOUTHEAST	
		elif( (lonat-lon)%360 <= 123.75 and (lonat-lon)%360 > 101.25):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR225[Order[imgsize*i+j]])/imgsize)
					J = int(OR225[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+imgsize-J) >= 0 and (x0+imgsize-J) < XSize[res] and (y0+I+1) > 0 and (y0+I+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1				

#		SOUTHEAST	
		elif( (lonat-lon)%360 <= 146.25 and (lonat-lon)%360 > 123.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+imgsize-J) >= 0 and (x0+imgsize-J) < XSize[res] and (y0+I+1) > 0 and (y0+I+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1				
#		SOUTHSOUTHEAST
		elif( (lonat-lon)%360 <= 168.75 and (lonat-lon)%360 > 146.25 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[int(OR225[Order[ct]])]))
					ct += 1
#		SOUTH	
		elif( (lonat-lon)%360 <= 191.25 and (lonat-lon)%360 > 168.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

#		SOUTHSOUTHWEST
		elif( (lonat-lon)%360 <= 213.75 and (lonat-lon)%360 > 191.25 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR225[Order[imgsize*i+j]])/imgsize)
					J = int(OR225[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+imgsize-I) > 0 and (x0+imgsize-I) < XSize[res] and (y0+imgsize-J) > 0 and (y0+imgsize-J) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1
#		SOUTHWEST
		elif( (lonat-lon)%360 <= 236.25 and (lonat-lon)%360 > 213.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+imgsize-I) > 0 and (x0+imgsize-I) < XSize[res] and (y0+imgsize-J) > 0 and (y0+imgsize-J) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1
#		WESTSOUTHWEST
		elif( (lonat-lon)%360 <= 258.75 and (lonat-lon)%360 > 236.25 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[int(OR225[Order[ct]])]))
					ct += 1
#		WEST
		elif( (lonat-lon)%360 <= 281.25 and (lonat-lon)%360 > 258.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

#		WESTNORTHWEST
		elif( (lonat-lon)%360 <= 303.75 and (lonat-lon)%360 > 281.25 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR225[Order[imgsize*i+j]])/imgsize)
					J = int(OR225[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+J+1) > 0 and (x0+J+1) < XSize[res] and (y0+imgsize-I) > 0 and (y0+imgsize-I) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1

#		NORTHWEST
		elif( (lonat-lon)%360 <= 326.25 and (lonat-lon)%360 > 303.75 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+J+1) > 0 and (x0+J+1) < XSize[res] and (y0+imgsize-I) > 0 and (y0+imgsize-I) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-I,x0+J+1,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1

#		NORTHNORTHWEST
		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[int(OR225[Order[ct]])]))

					ct += 1
	elif(lat > 10):
#		NORTH
		if( (lonat-lon)%360 <= 22.5 or (lonat-lon)%360 > 337.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
#		NORTHEAST
		elif( (lonat-lon)%360 <= 67.5 or (lonat-lon)%360 > 22.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+I+1) > 0 and (x0+I+1) < XSize[res] and (y0+J+1) > 0 and (y0+J+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+J+1,x0+I+1,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
	
					ct += 1
#		EAST
		elif( (lonat-lon)%360 <= 112.5 and (lonat-lon)%360 > 67.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

#		SOUTHEAST
		elif( (lonat-lon)%360 <= 157.5 and (lonat-lon)%360 > 112.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+J+1) > 0 and (x0+J+1) < XSize[res] and (y0+imgsize-I) > 0 and (y0+imgsize-I) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-I,x0+K+1,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1
#		SOUTH	
		elif( (lonat-lon)%360 <= 202.5 and (lonat-lon)%360 > 157.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
#		SOUTHWEST	
		elif( (lonat-lon)%360 <= 247.5 and (lonat-lon)%360 > 202.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,K,0] != skel[0,0,0]) and (x0+imgsize-I) > 0 and (x0+imgsize-I) < XSize[res] and (y0+imgsize-K) > 0 and (y0+imgsize-J) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+imgsize-J,x0+imgsize-I,:] = INT2Color(int(strdat[Order[ct]]))#INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1
#		WEST
		elif( (lonat-lon)%360 <= 292.5 and (lonat-lon)%360 > 247.5 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1	

#		NORTHWEST
		else:
			for i in range(imgsize):
				for j in range(imgsize):
					I = int(float(OR45[Order[imgsize*i+j]])/imgsize)
					J = int(OR45[Order[imgsize*i+j]])%imgsize
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or res == 5) and (skel[I,J,0] != skel[0,0,0]) and (x0+imgsize-J) >= 0 and (x0+imgsize-J) < XSize[res] and (y0+I+1) > 0 and (y0+I+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 5):
							img0[y0+I+1,x0+imgsize-J,:] = INT2Color(int(strdat[Order[ct]]))*INT2Color(int(strdat[int(OR45[Order[ct]])]))
					ct += 1		
		
	
						
	return img0

Zcols = plt.imread('ZCols.png')
def ColorMapH(hval):
	outarr = np.array((0.0,0.0,0.0,1.0))
	outarr[0:3] = Zcols[int(math.sqrt(hval)),0,0:3]
	return outarr
	if(hval == 0):
		return np.array((10/255.0,10/255.0,52/255.0,1.0))
	elif(hval < 1000):
		return np.array( (max(0.5,math.sqrt(hval/9000)) , 1.0-hval/2000 , 0.125+hval/8000 ,1.0) )
	else:
		return np.array( (max(0.5,math.sqrt(hval/9000)) , 0.5+(hval-1000)/16000 , 0.25+(hval-1000)/11000 ,1.0) )


def AssignTile4(pid,x1,y1,img0,Textures,skel,msk,lat,lon,res):
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	strdat = []
	Order = (0,2,1,3)
	if(res >  0):
		strdat = Textures.split(' ')
		if(res == 1):
			Order = (0,2,1,3)
		if(res == 2):
			Order = (0,6,15,2,7,8,5,14,9,4,13,11,1,10,12,3)
			Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
		if(res == 3):
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63)
		if(res == 4):
			Order = list(range(256))
	ct = 0
	lonat = Grid[pid].getLonDeg()
	if(lat >= -50 and lat <= 50):
		for i in range(imgsize):
			for j in range(imgsize):
				if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
					if(res == 0):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
					if(res == 1):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 2):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 3):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 4):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
				ct += 1
	elif(lat < -50):
		if( (lonat-lon)%360 <= 45 or (lonat-lon)%360 > 315 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
		elif( (lonat-lon)%360 <= 135 and (lonat-lon)%360 > 45 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1				
	
		elif( (lonat-lon)%360 <= 225 and (lonat-lon)%360 > 135 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

	elif(lat > 50):
		if( (lonat-lon)%360 <= 45 or (lonat-lon)%360 > 315 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
		elif( (lonat-lon)%360 <= 135 and (lonat-lon)%360 > 45 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
	
		elif( (lonat-lon)%360 <= 225 and (lonat-lon)%360 > 135 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1		
		
	
						
	return img0


RESF = [21,10.5,5.25,2.625,1.3125,0.65625]
def AssignINT2Color(pid,x1,y1,z1,img0,Textures,skel,msk,lat,lon,res):
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	z0 = int((z1)*(imgsize+1))#-int(imgsize/2.0)
	IGBP = int(Grid[pid].getIGBP())
	CLIM = int(Grid[pid].getCLIM())
	RGN = int(Grid[pid].getRGN())
	strdat = []
	Order = (0,2,1,3)
	if(res >  0):
		strdat = Textures.split(' ')
		if(res == 1):
			Order = (0,2,1,3)
		if(res == 2):
			Order = (0,6,15,2,7,8,5,14,9,4,13,11,1,10,12,3)
			Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
		if(res == 3):
			Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63)
		if(res == 4):
			Order = list(range(256))
	ct = 0
	lonat = Grid[pid].getLonDeg()
	if(lat >= -50 and lat <= 50):
		for i in range(imgsize):
			for j in range(imgsize):
				if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
					if(res == 0):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
					if(res == 1):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 2):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 3):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					if(res == 4):
						img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
				ct += 1
	elif(lat < -50):
		if( (lonat-lon)%360 <= 45 or (lonat-lon)%360 > 315 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
		elif( (lonat-lon)%360 <= 135 and (lonat-lon)%360 > 45 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1				
	
		elif( (lonat-lon)%360 <= 225 and (lonat-lon)%360 > 135 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

	elif(lat > 50):
		if( (lonat-lon)%360 <= 45 or (lonat-lon)%360 > 315 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
		elif( (lonat-lon)%360 <= 135 and (lonat-lon)%360 > 45 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+j+1) > 0 and (x0+j+1) < XSize[res] and (y0+imgsize-i) > 0 and (y0+imgsize-i) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-i,x0+j+1,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1
	
		elif( (lonat-lon)%360 <= 225 and (lonat-lon)%360 > 135 ):
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-i) > 0 and (x0+imgsize-i) < XSize[res] and (y0+imgsize-j) > 0 and (y0+imgsize-j) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 3):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 4):
							img0[y0+imgsize-j,x0+imgsize-i,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1

		else:
			for i in range(imgsize):
				for j in range(imgsize):
					if((res == 0 or res == 1 or res == 2 or res == 3 or res == 4 or skel[i,j,0] != skel[0,0,0]) and (x0+imgsize-j) >= 0 and (x0+imgsize-j) < XSize[res] and (y0+i+1) > 0 and (y0+i+1) < YSize[res]):
						if(res == 0):
							img0[y0+j+1,x0+i+1,:] = INT2Color(int(Textures))
						if(res == 1):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
						if(res == 2):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))			
						if(res == 3):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))	
						if(res == 4):
							img0[y0+i+1,x0+imgsize-j,:] = INT2Color(int(strdat[Order[ct]]))
					ct += 1		
		
	
						
	return img0


def AssignID(pid,x1,y1,img0,skel,res):
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)

	for i in range(imgsize):
		for j in range(imgsize):
			if((res == 0 or res == 1 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
				img0[y0+j+1,x0+i+1,:] = INT2Color(pid)
	return img0

def AssignID2(pid,x1,y1,img0,skel,res):
	if(res == 6):
		skel = plt.imread('imgcirc/skel52.png')
	imgsize = 360
	if(res >= 0):
		imgsize = 2**res
		imgsize = int(imgsize*13/16)
	x0 = int((x1)*(imgsize*16/13+1))+int(XSize[res]/2)#-int(imgsize/2.0)
	y0 = int((y1)*(imgsize*16/13+1))+int(YSize[res]/2)#-int(imgsize/2.0)

	for i in range(imgsize):
		for j in range(imgsize):
			if((res == 0 or res == 1 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
				img0[y0+j+1,x0+i+1,:] = INT2Color(pid)
	return img0

def AssignID3(pid,x1,y1,img0,skel0,res):
	imgsize = 360
	skel = 1.0*skel0
	if(res >= 0):
		imgsize = 2**res
	if(res == 2):
		skel[0,1,0:3] = 1-skel0[0,1,0:3]
		skel[0,2,0:3] = 1-skel0[0,2,0:3]
		skel[1,1,0:3] = 1-skel0[1,1,0:3]
		skel[1,2,0:3] = 1-skel0[1,2,0:3]
		skel[2,1,0:3] = 1-skel0[2,1,0:3]
		skel[2,2,0:3] = 1-skel0[2,2,0:3]
		skel[3,1,0:3] = 1-skel0[3,1,0:3]
		skel[3,2,0:3] = 1-skel0[3,2,0:3]
		skel[1,0,0:3] = 1-skel0[1,0,0:3]
		skel[2,0,0:3] = 1-skel0[2,0,0:3]
		skel[1,3,0:3] = 1-skel0[1,3,0:3]
		skel[2,3,0:3] = 1-skel0[2,3,0:3]
	x0 = int((x1)*(imgsize+1))+int(XSize[res]/2)#int(imgsize/2.0)
	y0 = int((y1)*(imgsize+1))+int(YSize[res]/2)#-int(imgsize/2.0)
	for i in range(imgsize):
		for j in range(imgsize):
			if((res == 0 or res == 1 or skel[i,j,0] != skel[0,0,0]) and (x0+i+1) > 0 and (x0+i+1) < XSize[res] and (y0+j+1) > 0 and (y0+j+1) < YSize[res]):
				img0[y0+j+1,x0+i+1,:] = INT2Color(pid)
	return img0
		
#	For Compressing Model Data into unsigned low-bit integers

def D8toCLD(uint8):
	cloudval = int(uint8/25)
	hour = uint8%25
	return cloudval, hour

def CLDtoD8(cloudval, hour):
	return np.uint8(cloudval*25+hour)

def S16toCIRV(uint16):
	climval = int(uint16/(20*16*7))
	rem16 = uint16%(20*16*7)
	igbpval = int(rem16/(16*7))
	rem16 = rem16%(16*7)
	rgnval = int(rem16/7)
	vegval = rem16%7
	return climval, igbpval, rgnval, vegval

def CIRVtoS16(climval, igbpval, rgnval, vegval):
	return np.uint16(climval*20*16*7+igbpval*16*7+rgnval*7+vegval)

#	Turns Integer ID into 4-D 24-bit color tuple
def INT2Color(pid):
	return ((pid-pid%65536)/65536/255.0 , ((pid-pid%256)/256)%256/255.0, pid%256/255.0,255/255.0)	

def INTS2Color(pid):
	outarr = np.ones(shape=(pid.shape[0],pid.shape[1],4),dtype=np.float16)
	outarr[:,:,0] = (pid-pid%65536)/65536/255.0
	outarr[:,:,1] = ((pid-pid%256)/256)%256/255.0
	outarr[:,:,2] = pid%256/255.0
	return outarr	
	
#	Inverts 24-bit color tuple or array into Integer ID
def C2I(color):
	R = int(color[0]*255.0)	
	G = int(color[1]*255.0)
	B = int(color[2]*255.0)
	return (R*256*256+G*256+B)

def CA2I(color):
	outarr = np.zeros(shape=(color.shape[0],color.shape[1]),dtype=int)
	outarr[:,:] = color[:,:,0]*255*256*256+color[:,:,1]*255*256+color[:,:,2]*255
	return outarr

def Bits2Col(bits):
	value = 0
	for i in range(len(bits)):
		value += bits[i]*2**i
	return INT2Color(value)

def Col2Bits(col):
	value = C2I(col)
	bits = []
	oldvalue = value
	for i in range(24):
		bits.append(oldvalue%2)
		oldvalue = int(oldvalue/2)
	return bits

def ColorMultiply(col1, col2):
	col3 = np.array((0.0,0.0,0.0,1.0))
	col3[0] = col1[0]*col2[0]
	col3[1] = col1[1]*col2[1]
	col3[2] = col1[2]*col2[2]
	return col3
	
#ShowNN(88,8)

#LoadTextures(0)
#*****************************************
#	For Generating Shapes and Text		
#-----------------------------------------

#	Generates text using courier 17x18 pixel characters

def GenText(textstr):
	N = len(textstr)
	A = np.ones(shape=(18,N*17,4),dtype=np.float16)
	ct = 0
	for s in textstr:
		if(s == '/'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/SLASH.png')[:,:,0:3]
		elif(s == ':'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/COLON.png')[:,:,0:3]
		elif(s == ' '):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/SPACE.png')[:,:,0:3]
		elif(s == '*'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/ASTERISK.png')[:,:,0:3]
		elif(s == '.'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/PERIOD.png')[:,:,0:3]
		elif(s == ','):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/COMMA.png')[:,:,0:3]
		elif(s == '-'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/DASH.png')[:,:,0:3]
		elif(s == '_'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/UNDERSCORE.png')[:,:,0:3]
		else:
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/'+s+'.png')[:,:,0:3]
		ct += 1
	return A

#	Col1 = Color of Text
#	Col2 = Color of Background
def GenTextCOL(textstr,col1,col2):
	N = len(textstr)
	A = np.ones(shape=(18,N*17,4),dtype=np.float16)
	ct = 0
	for s in textstr:
		if(s == '/'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/SLASH.png')[:,:,0:3]
		elif(s == ':'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/COLON.png')[:,:,0:3]
		elif(s == ' '):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/SPACE.png')[:,:,0:3]
		elif(s == '*'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/ASTERISK.png')[:,:,0:3]
		elif(s == '.'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/PERIOD.png')[:,:,0:3]
		elif(s == ','):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/COMMA.png')[:,:,0:3]
		elif(s == '-'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/DASH.png')[:,:,0:3]
		elif(s == '_'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/UNDERSCORE.png')[:,:,0:3]
		elif(s == '%'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/PCT.png')[:,:,0:3]
		elif(s == '('):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/OPENPARENTHESES.png')[:,:,0:3]
		elif(s == ')'):
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/CLOSEDPARENTHESES.png')[:,:,0:3]
		else:
			A[:,ct*17:ct*17+17,0:3] = plt.imread('imgcirc/Fonts/'+s+'.png')[:,:,0:3]
		ct += 1
	ZeroText = np.ones(shape=(18,N*17,4),dtype=np.float16)
	ZeroText[:,:,0:3] = 0.0
	B = SelectValue(A,ZeroText)
	C = np.ones(shape=(18,N*17,4),dtype=np.float16)
	C[:,:,0] = (B[:,:,0])*col2[0]
	C[:,:,1] = (B[:,:,1])*col2[1]
	C[:,:,2] = (B[:,:,2])*col2[2]
	D = np.zeros(shape=(18,N*17),dtype=int)
	D[:,:] = C2I(col1)
	
	A = INTS2Color(D)*(1-B)+C
	return A
PTFT2 = [0,0,0,0,0,0,0,0,0,0,0,13,13,0,0,0,0,0,17]

def GenTextCOL_size(textstr,col1,col2,ptft):
	ptft2 = PTFT2[ptft]
	N = len(textstr)
	A = np.ones(shape=(int(ptft),N*int(ptft2),4),dtype=np.float16)
	ct = 0
	for s in textstr:
		if(s == '/'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/SLASH.png')
		elif(s == ':'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/COLON.png')
		elif(s == ' '):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/SPACE.png')[:,:,0:3]
		elif(s == '*'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/ASTERISK.png')[:,:,0:3]
		elif(s == '.'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/PERIOD.png')[:,:,0:3]
		elif(s == ','):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/COMMA.png')[:,:,0:3]
		elif(s == '-'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/DASH.png')[:,:,0:3]
		elif(s == '_'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/UNDERSCORE.png')[:,:,0:3]
		elif(s == '%'):
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/PCT.png')[:,:,0:3]
		else:
#			print(ct)
#			print(ptft)
#			print(ptft2)
#			print(plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/'+s+'.png').shape)
#			print(A[:,ct*int(ptft):ct*int(ptft2)+int(ptft2),0:3].shape)
#			print(str(s))
			A[:,ct*int(ptft2):ct*int(ptft2)+int(ptft2),0:3] = plt.imread('imgcirc/Fonts/'+str(ptft)+'pt/'+s+'.png')[:,:,0:3]
		ct += 1
	ZeroText = np.ones(shape=(int(ptft),N*int(ptft2),4),dtype=np.float16)
	ZeroText[:,:,0:3] = 0.0
	B = SelectValue(A,ZeroText)
	C = np.ones(shape=(int(ptft),N*int(ptft2),4),dtype=np.float16)
	C[:,:,0] = (B[:,:,0])*col2[0]
	C[:,:,1] = (B[:,:,1])*col2[1]
	C[:,:,2] = (B[:,:,2])*col2[2]
	D = np.zeros(shape=(int(ptft),N*int(ptft2)),dtype=int)
	D[:,:] = C2I(col1)
	
	A = INTS2Color(D)*(1-B)+C
	return A



def GenCircle(res):
	a = np.ones(shape=(res,res,4),dtype='float32')
	for i in range(res):
		for j in range(res):
			if((i-res/2+0.5)*(i-res/2+0.5)+(j-res/2+0.5)*(j-res/2+0.5) < (res/2)*(res/2)):
				a[i,j,0:3] = np.array((0,0,0))
			else:
				a[i,j,0:3] = np.array((1,1,1))
	plt.imsave(fname='SKEL.png',arr=a)

#GenCircle(360)

#*******************************
#	Full-Image Resources
#-------------------------------

#	Select / Kronecker Delta
#	F(x,y) = 1 when x == y, F(x,y) = 0 when x != y
def SelectValue(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	return outf

def SelectValueC(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	outc = np.ones(shape=(field.shape[0],field.shape[1],4),dtype=np.float16)
	outc[:,:,0] = outf
	outc[:,:,1] = outf
	outc[:,:,2] = outf
	return outc

def SelectValue2(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	return outf

def SelectValueI(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	return outf

def SelectValue0(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = 1-np.exp(-1*np.abs(np.log( (value+1)/(field+1) ) ) )
	return outf

def SelectValuesLT(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = ((value-field)+np.abs(value-field))/(value-field)
	return outf

def SelectValuesGT(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = ((field-value)+np.abs(field-value))/(field-value)
	return outf



#	Highlight Value
#	F(x,y) = x when y == 0, F(x,y) = (color+x)/2  when y == 1
def HighlightValue(field,value,color):
	outf = np.ones(shape=field.shape,dtype=np.float16)
	selected = np.zeros(shape=field.shape,dtype=int)
	selected = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	for i in range(3):
		outf[:,:,i] = selected[:,:,i]*(field[:,:,i]+color[i])/2.0+(1-selected[:,:,i])*field[:,:,i]
	return outf

#	Highlight Selected
def HighlightSelected(field,selected,color):
	outf = np.ones(shape=field.shape,dtype=np.float16)
	for i in range(3):
		outf[:,:,i] = selected[:,:,i]*(field[:,:,i]+color[i])/2.0+(1-selected[:,:,i])*field[:,:,i]
	return outf

#	Highlight Selected (Green)
def HighlightSelected0(field,selected):
#	outf = np.ones(shape=field.shape,dtype=np.float16)
	outf = field[:,:,:]
	outf[:,:,1] = selected[:,:,1]*(field[:,:,1]+1.0)/2.0+(1-selected[:,:,1])*field[:,:,1]
#	outf[:,:,0] = selected[:,:,0]*field[:,:,0]+(1-selected[:,:,0])*field[:,:,0]
#
#	for i in range(3):
#		outf[:,:,i] = selected[:,:,i]*(field[:,:,i]+color[i])/2.0+(1-selected[:,:,i])*field[:,:,i]
	return outf

def HighlightSelected1(field,selected):
	color = np.array((0.0,1.0,0.0,1.0))
	outf = np.ones(shape=field.shape,dtype=np.float16)
	for i in range(3):
		outf[:,:,i] = selected[:,:]*(field[:,:,i]+color[i])/2.0+(1-selected[:,:])*field[:,:,i]
	return outf


def HighlightSelected2(field,selected):
	outf = field[:,:,:]
	outf[:,:,1] = selected[:,:]*(field[:,:,1]+1.0)/2.0+(1-selected[:,:])*field[:,:,1]
	return outf




#	Cloud Filter
#	F(x,y) = x when y == 0, F(x,y) = 1 when y == 1
def CloudFilter(field,value):
	return field*(1-value)+value


#	Color Filter			\\ Generalized Cloud Filter
#	F(x,y,color) = x when y == 0, F(x,y,color) = color when y == 1
def ColorFilter(field,value,color):
	outf = np.ones(shape=field.shape,dtype=np.float16)
	for i in range(3):
		outf[:,:,i] = field[:,:,i]*(1-value[:,:,i])+value[:,:,i]*color[i]
	return outf


#	Shade Filter
#	F(x,y,alpha) = alpha * x when y == 0, F(x,y) = 1 when y == 1
def ShadeFilter(field,value,alpha):
	return alpha*field*(1-value)+value


#	Concatenate Selection Arrays
#	Returns 1 if any array in arrays = 1, 0 otherwise
def CatSelection(arrays):
	outf = np.ones(shape=arrays[0].shape,dtype=np.float16)
	for a in arrays:
		outf = outf*(1-a)
	return (1-outf)

#	Select Values
def SelectValues(field,values):
	selection = []
	for value in values:
		outf = np.zeros(shape=field.shape,dtype=int)
		outf = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
		selection.append(outf)
	return CatSelection(selection)

#	Sharpen



#PIECES = ['Utility','Army','Navy','Airforce','','','','','','','','ArmyGeneral','NavyGeneral','AirforceGeneral']
#def PlacePiece(img,ptstr,imgid,pieces,res):
#	pieces = list(open('IMG/Zoom'+res+'_anc/Plane'+str(lat)+' '+str(lon)+'.txt'))
#	for p in ptstr:
#		P = p.split(' ')
#		x0 = int((P[1])*(imgsize+1))+int(XSize[res]/2)
#		y0 = int((P[1])*(imgsize+1))+int(YSize[res]/2)
#		piece = Character.getPiece(
#		
#		 
		
	