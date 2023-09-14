import math
import numpy as np
import Province
IGBPTypes = ('Ocean','Coast','EvergreenNeedleleaf','EvergreenBroadleaf','DeciduousNeedleleaf','DeciduousBroadleaf','MixedForest','ClosedShrubland','OpenShrubland','WoodedSavanna','Savanna','Grasslands','PermanentWetlands','Croplands','Urban','SparseCrops','SnowIce','Desert','Peak')
RGNTypes = list(open('RGNZones.txt'))#('Flat','Hill','Mountain','Peak')
CLIMTypes = list(open('ClimateZones.txt'))#('Grasslands','Plains','Desert','Tundra','Snow')
Provinces = list(open("Region0"))
LENS = [0,0,0,0,0,0,92162,368642,1474562,5898242]

LEVEL = 8


distancescale = math.sqrt(16.0/LENS[LEVEL])
re = 6371000.0
class Point9:
	PointList = []
	PointCount = 0
#	def __init__(self, IDnum, affiliation, xpos, ypos, zpos, gridface, isborder, terrain, nearestPlane):#xypos, projxypos, NNlist ):
#	def __init__(self, IDnum, affiliation, terrain, hgt, isborder, xpos,ypos,zpos):
	def __init__(self,IDnum, affiliation, IGBP, RGN, isborder, xpos,ypos,zpos, CLIM,POPL):
		self.ID = int(IDnum)					#Integer
		self.Province = int(affiliation)			#String
		self.XYZPos = (float(xpos),float(ypos),float(zpos))	#Vector or list
		self.XYZStr = xpos+' '+ypos+' '+zpos			#String
#		self.GridFace = int(gridface)				#Integer
		self.isBorder = int(isborder)				#Integer
#		self.Terrain = terrain 					#String
#		self.Terrain2 = hgt
		self.Population = 100
		self.Growth = 0.01
		self.Money = 1000
		self.Income = 1.0
		self.NN = []
		self.NumNN = 0
#		self.nearestPlane = nearestPlane

		self.IGBP = int(IGBP)
		self.RGN = int(RGN)
		self.HGT = 0
		self.POPL = int(POPL)
		self.CLIM = int(CLIM)
		xpos = float(xpos)
		ypos = float(ypos)
		zpos = float(zpos)
		phi = -1*math.asin(zpos/math.sqrt(xpos*xpos+ypos*ypos+zpos*zpos))
		lam = math.pi/2.0
		the = math.acos(-1*zpos)
		if(xpos	< 0 and ypos < 0):
			lam = math.atan(ypos/xpos)
		if(xpos == 0 and ypos < 0):
			lam = 1.0*math.pi/2.0
		if(xpos > 0 and ypos < 0):
			lam = math.atan(ypos/xpos)+math.pi
		if(xpos < 0 and ypos >= 0):
			lam = math.atan(ypos/xpos)
		if(xpos == 0 and ypos >= 0):
			lam = 3.0*math.pi/2.0
		if(xpos > 0 and ypos >= 0):
			lam = math.atan(ypos/xpos)+math.pi
		lam = lam-math.pi
		if(lam < 0):
			lam = lam+2*math.pi
		if(lam >= 2*math.pi):
			lam = lam-2*math.pi
		self.lat = phi
		self.lon = lam
		self.the = the
		self.North = (0.0,0.0,1000)
		if(math.cos(the) == 0):
			self.North = (0.0,1.0,0.0)
		else:
			self.North = (0.0,0.0,1.0/math.cos(the))

		ey0 = (0.0,0.0,1.0)
		if(math.sin(the) != 0 and math.cos(the) != 0):
			ey0 = (self.North[1]*zpos-self.North[2]*ypos,self.North[2]*xpos-self.North[0]*zpos,self.North[0]*ypos-self.North[1]*xpos)
		elif(the == 0):
			ey0 = (0.0,1.0,0.0)
		elif(the == math.pi/2.0):
			ey0 = (xpos,ypos,1.0)
		ey = (0.0,1.0,0.0)
		if( (ey0[0]*ey0[0]+ey0[1]*ey0[1]+ey0[2]*ey0[2]) != 0):
			ey = (ey0[0]/(ey0[0]*ey0[0]+ey0[1]*ey0[1]+ey0[2]*ey0[2]),ey0[1]/(ey0[0]*ey0[0]+ey0[1]*ey0[1]+ey0[2]*ey0[2]),ey0[2]/(ey0[0]*ey0[0]+ey0[1]*ey0[1]+ey0[2]*ey0[2]))
		ez = (xpos,ypos,zpos)
		ex = (-1*ey[1]*ez[2]+ey[2]*ez[1],-1*ey[2]*ez[0]+ey[0]*ez[2],-1*ey[0]*ez[1]+ey[1]*ez[0])

		self.ex = ey
		self.ey = ex
		self.ez = ez		

		Point9.PointList.append(self)				#Add new point to list of points
		Point9.PointCount += 1					#Append the size of point list

#		print('Adding '+str(affiliation))
		Province.addPointPRV(self,int(affiliation))

	def PrintStats(self):
		print("Point ID: "+str(self.ID))
		nnstr = ''
		for i in range(self.NumNN):
			nnstr += str(self.NN[i])+' '
		print("NN List:  "+nnstr)
		print("Province: "+str(Provinces[int(self.Province)].split(' ')[0]))
		print("XYZ Pos:  "+str(self.XYZStr))
		print("Lat:      "+str(180.0*self.lat/math.pi))
		print("Lon:      "+str(180.0*self.lon/math.pi))

		if(self.isBorder == 1):
			print("Border:   True")
		else:
			print("Border:   False")
		print("IGBP:     "+IGBPTypes[self.IGBP])
		print("Terrain:  "+RGNTypes[self.RGN][0:-1])
		print("Climate:  "+CLIMTypes[self.CLIM][0:-1])
		print(' ')
	def getID(self):
		return self.ID
	def getProvince(self):
		return self.Province
	def setProvince(self, prov):
		self.Province = prov


	def getXYZPos(self):
		return self.XYZPos
	def getXYZPosR(self):
		x8 = math.cos(self.lat)*math.sin(self.lon)
		y8 = math.cos(self.lat)*math.cos(self.lon)
		z8 = math.sin(self.lat)
		return (x8,y8,z8)

	def getXYZPosSTR(self):
		x8 = math.cos(self.lat)*math.sin(self.lon)
		y8 = math.cos(self.lat)*math.cos(self.lon)
		z8 = math.sin(self.lat)
		return str(round(x8,4))+' '+str(round(y8,4))+' '+str(round(z8,4))

	def getXYZStr(self):
		return self.XYZStr
#	def getGridFace(self):
#		return self.GridFace

	def getBorder(self):
		return self.isBorder
	def setBorder(self, border):
		self.isBorder = border
	

	def isBorder(self):
		if(self.isBorder == 1):
			return True
		else:
			return False
	def getIGBP(self):
		return self.IGBP
	def setIGBP(self, igbp):
		self.IGBP = igbp

	
	def getCLIM(self):
		return self.CLIM
	def setCLIM(self,clim):
		self.CLIM = clim

	def getRGN(self):
		return self.RGN
	def setRGN(self, rgn):
		self.RGN = rgn	

	def getHGT(self):
		return self.HGT
	def setHGT(self, hgt):
		self.HGT = hgt	

	def getPOPL(self):
		return self.POPL
	def setPOPL(self, popl):
		self.POPL = popl	

	def getTerrain(self):
		return self.Terrain
	def setTerrain(self, terrain):
		self.Terrain = terrain	

	def getTerrain2(self):
		return self.Terrain2
	def setTerrain2(self, terrain2):
		self.Terrain2 = terrain2


#	def getNearestPlane(self):
#		return self.nearestPlane

	def getLat(self):
		return self.lat
	def getLatDeg(self):
		return 180.0*self.lat/math.pi
	def getLon(self):
		return self.lon
	def getLonDeg(self):
		return 180.0*self.lon/math.pi
	def getXYPosEQR(self,yframe):
		Y = int(self.the*2160/math.pi-0.5)
		X = int(self.lon*2160/math.pi-0.5)
		if(X < 0):
			X = 0
		if(X >= 2*yframe):
			X = 2*yframe-1
		if(Y < 0):
			Y = 0
		if(Y >= yframe):
			Y = yframe-1
		return (X,Y)
	def getXYStrEQR(self,yframe):
		Y = int(self.the*2160/math.pi-0.5)
		X = int(self.lon*2160/math.pi-0.5)
		if(X < 0):
			X = 0
		if(X >= 2*yframe):
			X = 2*yframe-1
		if(Y < 0):
			Y = 0
		if(Y >= yframe):
			Y = yframe-1
		return str(X)+' '+str(Y)

#		self.Population = 100
#		self.Growth = 0.01
#		self.Money = 1000
#		self.Income = 1
#	

	def addPopulation(self, amount):
		self.Population += amount
	def getPopulation(self):
		return self.Population
	def printPopulation(self):
		print("Population: "+str(self.Populaition))
	
	def addGrowth(self,amount):
		self.Growth += amount
	def getGrowth(self):
		return self.Growth
	def printGrowth(self):
		print("Growth Rate: "+str(self.Growth)+"%")

	def addMoney(self, amount):
		self.Money += amount
	def getMoney(self, amount):
		return self.Money
	def printMoney(self):
		print("Money: "+str(self.Money))

	def addIncome(self, amount):
		self.Income += amount
	def getIncome(self):
		return self.Income
	def printIncome(self):
		print("Income: "+str(self.Income))

	def addNN(self,IDNum):
		self.NN.append(IDNum)
		self.NumNN = len(self.NN)


	def getDDX(self, field):
		vdir = self.ex
		Xval = []
		Fval = []
		Xval.append(re*(self.getXYZPos[0]*vdir[0]+self.getXYZPos[1]*vdir[1]+self.getXYZPos[2]*vdir[2]))
		Fval.append(field[self.getID()])
		for pt in self.NN:
			Xval.append(re*(Pointlist[pt].getXYZPos()[0]*vdir[0]+Pointlist[pt].getXYZPos()[1]*vdir[1]+Pointlist[pt].getXYZPos()[2]*vdir[2]))
			Fval.append(field[pt])
		XTot = 0.0
		Ftot = 0.0
		for i in range(len(Xval)):
			Xtot += Xval[i]
			Ftot += Ftot[i]
		Xavg = Xtot/len(Xval)
		Favg = Ftot/len(Fval)
		Xtot = 0.0
		Ftot = 0.0
		for i in range(len(Xval)):
			Xtot += (Xval[i]-Xavg)^2
			Ftot += (Fval[i]-Favg)*(Xval[i]-Xavg)
		return Ftot/Xtot


	def getDDY(self, field):
		vdir = self.ey
		Xval = []
		Fval = []
		Xval.append(re*(self.getXYZPos[0]*vdir[0]+self.getXYZPos[1]*vdir[1]+self.getXYZPos[2]*vdir[2]))
		Fval.append(field[self.getID()])
		for pt in self.NN:
			Xval.append(re*(Pointlist[pt].getXYZPos()[0]*vdir[0]+Pointlist[pt].getXYZPos()[1]*vdir[1]+Pointlist[pt].getXYZPos()[2]*vdir[2]))
			Fval.append(field[pt])
		XTot = 0.0
		Ftot = 0.0
		for i in range(len(Xval)):
			Xtot += Xval[i]
			Ftot += Ftot[i]
		Xavg = Xtot/len(Xval)
		Favg = Ftot/len(Fval)
		Xtot = 0.0
		Ftot = 0.0
		for i in range(len(Xval)):
			Xtot += (Xval[i]-Xavg)^2
			Ftot += (Fval[i]-Favg)*(Xval[i]-Xavg)
		return Ftot/Xtot

def getPoint(idnum):
	return Point9.PointList[idnum]

def getPointList():
	return Point9.PointList

#	Finds & returns the center of XYZ coordinates from Points list
def getCentralPoint(Points):
	xtot = 0.0
	ytot = 0.0
	ztot = 0.0
	for pt in Points:
		XYZ = pt.getXYZPos()
		xtot += XYZ[0]
		ytot += XYZ[1]
		ztot += XYZ[2]
	xNew = xtot/len(Points)
	yNew = ytot/len(Points)
	zNew = ztot/len(Points)
	size = math.sqrt(xNew*xNew+yNew*yNew+zNew*zNew)
	xNew = xNew/size
	yNew = yNew/size
	zNew = zNew/size
	return (xNew,yNew,zNew)

def getDDX2(Fval, r1):
	vdir = r1[0].ex
	Rval = []
	for pt in	 r1:
		Rval.append(re*(pt.getXYZPos[0]*vdir[0]+pt.getXYZPos[1]*vdir[1]+pt.getXYZPos[2]*vdir[2]))
	Rtot = 0.0
	Ftot = 0.0
	for i in range(len(Rval)):
		Rtot += Rval[i]
		Ftot += Ftot[i]
	Ravg = Rtot/len(Rval)
	Favg = Ftot/len(Fval)
	Rtot = 0.0
	Ftot = 0.0
	for i in range(len(Xval)):
		Rtot += (Rval[i]-Ravg)^2
		Ftot += (Fval[i]-Favg)*(Rval[i]-Ravg)
	return Ftot/Rtot

def getDDY2(Fval, r1):
	vdir = r1[1].ey
	Rval = []
	for pt in	 r1:
		Rval.append(re*(pt.getXYZPos[0]*vdir[0]+pt.getXYZPos[1]*vdir[1]+pt.getXYZPos[2]*vdir[2]))
	Rtot = 0.0
	Ftot = 0.0
	for i in range(len(Rval)):
		Rtot += Rval[i]
		Ftot += Ftot[i]
	Ravg = Rtot/len(Rval)
	Favg = Ftot/len(Fval)
	Rtot = 0.0
	Ftot = 0.0
	for i in range(len(Xval)):
		Rtot += (Rval[i]-Ravg)^2
		Ftot += (Fval[i]-Favg)*(Rval[i]-Ravg)
	return Ftot/Rtot

#	Compress/Decompress CIRV data

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

#	Loads list of points	
def LoadPoints(fname):
	fname = 'CRDS'+str(LEVEL)+'.pts'
	Points = list(open(fname))
	NNs = list(open('NN'+str(LEVEL)+'.txt'))
	ct = 0
	for pt in Points:
		tmp = pt.split(' ')
		point = Point9(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[9])
		NNtmp = NNs[ct].split(' ')
		for i in range(len(NNtmp)-2):
			point.addNN(NNtmp[i+1])
		ct += 1

def SavePoints():
	fname = 'CRDS'+str(LEVEL)+'.pts'
	Points = open(fname,'w')
	for p in Point9.PointList:
		str0 = str(p.getID())
		str1 = str(p.getProvince())
		str2 = str(p.getIGBP())
		str3 = str(p.getRGN())
		str4 = str(p.getBorder())
		str5 = str(p.getXYZStr())
		str6 = str(p.getCLIM())
		str7 = str(p.getPOPL())
		Points.write(str0+' '+str1+' '+str2+' '+str3+' '+str4+' '+str5+' '+str6+' '+str7+'\n')
	Points.close()

LoadPoints('F9.txt')
