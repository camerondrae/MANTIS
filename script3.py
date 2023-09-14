import numpy as np
import matplotlib.pyplot as plt
import math
import Builder
import os
import random
from netCDF4 import Dataset

YF = 1080
XF = 1920
Yf = 900
Xf = 1600

YSize = [YF,YF,3000,4500,6500,10500,6400,6400]
XSize = [XF,XF,3000,4500,6500,10500,6400,6400]
			
					
				

def colormultiply(col1,col2):
	col3 = np.array((0.0,0.0,0.0,1.0))
	col3[0] = col1[0]*col2[0]
	col3[1] = col1[1]*col2[1]
	col3[2] = col1[2]*col2[2]
	return col3
			

def GenCentralPtMap():
	g = open('XY2PT.txt','w')
	for i in range(-88,96,8):
		for j in range(0,360,8):
			img = plt.imread('IMG/Zoom4/Plane'+str(i)+'x'+str(j)+'.png')
			YF2 = int(Yf/2)
			XF2 = int(Xf/2)
			idnum = Builder.C2I(img[YF2,XF2])
			if(idnum == 16777215):
				nearbycols = []
				if(Builder.C2I(img[YF2-2,XF2-2]) != 16777215):
					idnum = Builder.C2I(img[YF2-2,XF2-2])
				elif(Builder.C2I(img[YF2-2,XF2-1]) != 16777215):
					idnum = Builder.C2I(img[YF2-2,XF2-1])
				elif(Builder.C2I(img[YF2-2,XF2]) != 16777215):
					idnum = Builder.C2I(img[YF2-2,XF2])
				elif(Builder.C2I(img[YF2-2,XF2+1]) != 16777215):
					idnum = Builder.C2I(img[YF2-2,XF2+1])
				elif(Builder.C2I(img[YF2-1,XF2-2]) != 16777215):
					idnum = Builder.C2I(img[YF2-1,XF2-2])
				elif(Builder.C2I(img[YF2-1,XF2-1]) != 16777215):
					idnum = Builder.C2I(img[YF2-1,XF2-1])
				elif(Builder.C2I(img[YF2-1,XF2]) != 16777215):
					idnum = Builder.C2I(img[YF2-1,XF2])
				elif(Builder.C2I(img[YF2-1,XF2+1]) != 16777215):
					idnum = Builder.C2I(img[YF2-1,XF2+1])
				elif(Builder.C2I(img[YF2,XF2-2]) != 16777215):
					idnum = Builder.C2I(img[YF2,XF2-2])
				elif(Builder.C2I(img[YF2,XF2-1]) != 16777215):
					idnum = Builder.C2I(img[YF2,XF2-1])
				elif(Builder.C2I(img[YF2,XF2]) != 16777215):
					idnum = Builder.C2I(img[YF2,XF2])
				elif(Builder.C2I(img[YF2,XF2+1]) != 16777215):
					idnum = Builder.C2I(img[YF2,XF2+1])
				elif(Builder.C2I(img[YF2+1,XF2-2]) != 16777215):
					idnum = Builder.C2I(img[YF2+1,XF2-2])
				elif(Builder.C2I(img[YF2+1,XF2-1]) != 16777215):
					idnum = Builder.C2I(img[YF2+1,XF2-1])
				elif(Builder.C2I(img[YF2+1,XF2]) != 16777215):
					idnum = Builder.C2I(img[YF2+1,XF2])
				elif(Builder.C2I(img[YF2+1,XF2+1]) != 16777215):
					idnum = Builder.C2I(img[YF2+1,XF2+1])
			g.write(str(i)+' '+str(j)+' '+str(idnum)+'\n')
			if(idnum == 16777215):
				print('Help at '+str(i)+' '+str(j))
	g.close()

#GenCentralPtMap()

def GenOnEachPlane(zoom):
	YF = YSize[zoom]
	XF = XSize[zoom]
	G = list(open('XY2PT.txt'))
	ct = 0 
	for data in G:
		print(str(ct))
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
#		print(str(y0+1-int(Yf/2))+' to '+str(y0+int(Yf/2)))
#		print(str(x0+1-int(Xf/2))+' to '+str(x0+int(Xf/2)))
		imgnew = imgpix[y0+1-int(Yf/2):y0+int(Yf/2),x0+1-int(Xf/2):x0+int(XF/2)]
		plt.imsave(fname='AssetsF/Zoom'+str(zoom)+'/Plane'+tmp[0]+'x'+tmp[1]+'.png',arr=imgnew)
		ct += 1

#GenOnEachPlane(2)

def Reorder():
	for PLANE in range(83,84,1):
		print(str(PLANE))
		if(PLANE != 80 and PLANE != 83):
			PDAT = list(open('Planes/Plane'+str(PLANE)+'.txt'))
			NDAT = np.ones(shape=(1474562,2),dtype=np.float32)*1000.0
			for i in PDAT:
				tmp = i.split(' ')
				NDAT[int(tmp[0]),0] = float(tmp[1])
				NDAT[int(tmp[0]),1] = float(tmp[2])
			g = open('Planes/P2_'+str(PLANE)+'.txt','w')
			for i in range(1474562):
				g.write(str(np.round(NDAT[i,0],4))+' '+str(np.round(NDAT[i,1],4))+'\n')
			g.close()
		elif(PLANE == 80):
			LABS = ['60','61','64','62','63']
			for LAB in LABS:
				PDAT = list(open('Planes/Plane'+str(PLANE)+LAB+'.txt'))
				NDAT = np.ones(shape=(1474562,2),dtype=np.float32)*1000.0
				for i in PDAT:
					tmp = i.split(' ')
					NDAT[int(tmp[0]),0] = float(tmp[1])
					NDAT[int(tmp[0]),1] = float(tmp[2])
				g = open('Planes/P2_'+str(PLANE)+LAB+'.txt','w')
				for i in range(1474562):
					g.write(str(np.round(NDAT[i,0],4))+' '+str(np.round(NDAT[i,1],4))+'\n')
				g.close()
		elif(PLANE == 83):
			LABS = ['74','73','72','70','71']
			for LAB in LABS:
				PDAT = list(open('Planes/Plane'+str(PLANE)+LAB+'.txt'))
				NDAT = np.ones(shape=(1474562,2),dtype=np.float32)*1000.0
				for i in PDAT:
					tmp = i.split(' ')
					NDAT[int(tmp[0]),0] = float(tmp[1])
					NDAT[int(tmp[0]),1] = float(tmp[2])
				g = open('Planes/P2_'+str(PLANE)+LAB+'.txt','w')
				for i in range(1474562):
					g.write(str(np.round(NDAT[i,0],4))+' '+str(np.round(NDAT[i,1],4))+'\n')
				g.close()
		
#Reorder()

def ExtCond():
	for PLANE in range(80,92,1):
		if(PLANE != 80 and PLANE != 83):
			PDAT = list(open('Planes/Plane'+str(PLANE)+'.txt'))
			g = open('Planes/P'+str(PLANE)+'.txt','w')
			ct = 0
			for i in range(1474562):
				tmp = '0 0 0'
				if(ct < len(PDAT)): 
					tmp = PDAT[ct].split(' ')
				if(int(tmp[0]) == i):
					g.write(str(np.round(float(tmp[1]),4))+' '+str(np.round(float(tmp[2]),4))+'\n')
					ct += 1
				else:
					g.write('1000 1000'+'\n')
			g.close()
		elif(PLANE == 80):
			LABS = ['60','61','64','62','63']
			for LAB in LABS:
				PDAT = list(open('Planes/Plane'+str(PLANE)+LAB+'.txt'))
				g = open('Planes/P'+str(PLANE)+LAB+'.txt','w')
				ct = 0
				for i in range(1474562):
					tmp = '0 0 0'
					if(ct < len(PDAT)): 
						tmp = PDAT[ct].split(' ')
					if(int(tmp[0]) == i):
						g.write(str(np.round(float(tmp[1]),4))+' '+str(np.round(float(tmp[2]),4))+'\n')
						ct += 1
					else:
						g.write('1000 1000'+'\n')
				g.close()
		elif(PLANE == 83):
			LABS = ['74','73','72','70','71']
			for LAB in LABS:
				PDAT = list(open('Planes/Plane'+str(PLANE)+LAB+'.txt'))
				g = open('Planes/P'+str(PLANE)+LAB+'.txt','w')
				ct = 0
				for i in range(1474562):
					tmp = '0 0 0'
					if(ct < len(PDAT)): 
						tmp = PDAT[ct].split(' ')
					if(int(tmp[0]) == i):
						g.write(str(np.round(float(tmp[1]),4))+' '+str(np.round(float(tmp[2]),4))+'\n')
						ct += 1
					else:
						g.write('1000 1000'+'\n')
				g.close()

def GetClims():
	A = plt.imread('KoppenClimateData.png')
	B = []
	for i in range(2160):
		for j in rangE(4320):
			B.append(Builder.C2I(A[i,j]))
	return list(set(B))

#Clims = ['Oceanic','PolarTundra','IceCap']
#ClimID = []
#ClimID.append(Builder.C2I(np.array(255,255,255)/255.0))
#ClimID.append(Builder.C2I(np.array(178,178,178)/255.0))
#ClimID.append(Builder.C2I(np.array(102,102,102)/255.0))




def ReassignHGT():
	P = plt.imread('ToXY7.png')
	Z = plt.imread('VIMG.png')
	WGT = plt.imread('WGT4320.png')
	N0 = len(Builder.Grid)
	Values = np.zeros(shape=(N0))
	N = np.zeros(shape=(N0))
	IMGDAT = np.zeros(shape=(2160,4320,4),dtype=np.float16)
	for i in range(2160):
		print(str(i))
		for j in range(4320):
			Values[Builder.C2I(P[i,j])] += Builder.C2I(Z[i,j])*Builder.C2I(WGT[i,j])
			N[Builder.C2I(P[i,j])] += Builder.C2I(WGT[i,j])
	DAT = Values/N
	G = open('V8.txt','w')
	for i in range(N0):
		G.write(str(int(DAT[i]))+'\n')
	G.close()
		 
def AvgColors(colors):
	dt = np.array((0.0,0.0,0.0,1.0))
	ct = 0
	for c in colors:
		dt += c
		ct += 1
	if(ct == 0):
		return dt
	else:
		return dt/ct

def commoncolors(colores):
	colors = []
	for c in colores:
		colors.append(Builder.C2I(c))
	
	dt = np.array((0.0,0.0,0.0,1.0))
	vals = list(set(colors))
	nums = np.zeros(shape=(len(vals)),dtype=np.float16)
	valmax = Builder.C2I(np.array((1.0,1.0,1.0,1.0)))
	nummax = 0
	ct = 0
	for v in vals:
		if(v == Builder.C2I(dt)):
			ct += 1
		else:
			c0 = 0
			for c in colors:
				if(c == v):
					c0 += 1
			if(c0 > nummax):
				c0 = nummax
				valmax = v
	return Builder.INT2Color(valmax)
		

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
		if(v == Builder.C2I(dt)):
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

def commonvalues0(colores,numans):
	colors = []
	for c in colores:
		colors.append(c)

	dt = np.array((0.0,0.0,0.0,1.0))
	vals = list(set(colors))
#	try:
#		vals.remove(Builder.C2I((1.0,1.0,1.0,1.0)))
#	except:
#		b = 1
#	try:
#		vals.remove(Builder.C2I((0.0,0.0,0.0,1.0)))
#	except:
#		b = 1
	ans = []
	anssize = 0
	for k in range(numans):
		nums = np.zeros(shape=(len(vals)),dtype=np.float16)
		valmax = 1000000
		nummax = 0
		ct = 0

#		if(len(vals) == 0):
#			ans = []
#			for k in range(numans):
#				ans.append(Builder.C2I((1.0,1.0,1.0,1.0)))
#			return ans

		for v in vals:
			if(v == Builder.C2I(dt)):
				ct += 1
			else:
				c0 = 0
				for c in colors:
					if(c == v):
						c0 += 1
				if(c0 > nummax):
					c0 = nummax
					valmax = v
		ans.append(int(valmax))
		anssize += 1
		try:
			vals.remove(valmax)
		except:
			if(len(ans) < numans):
				for k in range(len(ans), numans,1):
					ans.append(Builder.C2I((1.0,1.0,1.0,1.0)))
				return ans
			else:
				return ans
			
	return ans
		

def AssignCol():
	P = plt.imread('ToXY5400.png')
	Z = plt.imread('BigData/Topo5400.png')
#	WGT = plt.imread('WGT4320.png')
	N0 = len(Builder.Grid)
	Values = []
	for i in range(N0):
		Values.append([])
	for i in range(2700):
		print(str(i))
		for j in range(5400):
			Values[Builder.C2I(P[i,j])].append(Z[i,j])
	DAT = []
	for i in range(N0):
		DAT.append(Builder.C2I(AvgColors(Values[i])))
	G = open('T0_8.txt','w')
	for i in range(N0):
		G.write(str(int(DAT[i]))+'\n')
	G.close()

#AssignCol()


def AssignCol2():
	print('Initializing...')
	PTL = list(open('PTList/GrandXYMap8_5400.txt'))
#	P = plt.imread('ToXY7.png')
	Z = plt.imread('BigData/Topo5400.png')
#	print(Z[0,0])
#	WGT = plt.imread('WGT4320.png')
	N0 = len(Builder.Grid)
#	CornerPoints = ((1,0),(2159,1),(0,4318),(2158,4319))

	DAT = np.zeros(shape=(N0,4),dtype=int)
	print('Collecting Data...')
	for pt in range(N0):
		Ytot = 0
		Xtot = 0
		Ntot = 0
		for l in PTL[pt].split(':')[0:-1]:
			xys = l.split(' ')
			Ytot += int(xys[1])
			Xtot += int(xys[0])
			Ntot += 1
		Yavg = float(Ytot/Ntot)
		Xavg = float(Xtot/Ntot)
		NWcols = []
		SWcols = []
		NEcols = []
		SEcols = []
		for l in PTL[pt].split(':')[0:-1]:
#			print(l)
			xys = l.split(' ')
#			print(xys)
#			print(Z[int(xys[1]),int(xys[0])])
			if( int(xys[1]) < Yavg and int(xys[0]) < Xavg):
				NWcols.append(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) >= Yavg and int(xys[0]) < Xavg):
				SWcols.append(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) < Yavg and int(xys[0]) >= Xavg):
				NEcols.append(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) >= Yavg and int(xys[0]) >= Xavg):
				SEcols.append(Z[int(xys[0]),int(xys[1])])
		DAT[pt,0] = Builder.C2I(AvgColors(NWcols))
		DAT[pt,1] = Builder.C2I(AvgColors(SWcols))
		DAT[pt,2] = Builder.C2I(AvgColors(NEcols))
		DAT[pt,3] = Builder.C2I(AvgColors(SEcols))
		if(any(DAT[pt,:] == 0)):
			cols = []
			for i in range(4):
				if(DAT[pt,i] != 0):
					cols.append(Builder.INT2Color(DAT[pt,i]))
			avgcol = AvgColors(cols)
			for i in range(4):
				if(DAT[pt,i] == 0):
					DAT[pt,i] = Builder.C2I(avgcol)
	print('Writing Data...')
	G = open('T1_8.txt','w')
	for i in range(N0):
#		print(str(i)+': '+str(DAT[i,0])+' '+str(DAT[i,1])+' '+str(DAT[i,2])+' '+str(DAT[i,3]))
		G.write(str(int(DAT[i,0]))+' '+str(int(DAT[i,1]))+' '+str(int(DAT[i,2]))+' '+str(int(DAT[i,3]))+'\n')
	G.close()

def AssignCol3():
	print('Initializing...')
	PTL = list(open('PTList/GrandXYMap8_5400.txt'))
#	P = plt.imread('ToXY7.png')
	Z = plt.imread('BigData/Topo5400.png')
#	print(Z[0,0])
#	WGT = plt.imread('WGT4320.png')
	N0 = len(Builder.Grid)
#	CornerPoints = ((1,0),(2159,1),(0,4318),(2158,4319))

	DAT = np.zeros(shape=(N0,16),dtype=int)
	print('Collecting Data...')
	for pt in range(N0):
		Ytot = 0
		Xtot = 0
		Ntot = 0
		for l in PTL[pt].split(':')[0:-1]:
			xys = l.split(' ')
			Ytot += int(xys[1])
			Xtot += int(xys[0])
			Ntot += 1
		Yavg = float(Ytot/Ntot)
		Xavg = float(Xtot/Ntot)
		Yavg = Builder.toXY(Builder.Grid[pt].getXYZPos(),5400,2700)[0]
		Xavg = Builder.toXY(Builder.Grid[pt].getXYZPos(),5400,2700)[1]
		
		NWcols = []
		SWcols = []
		NEcols = []
		SEcols = []
		for l in PTL[pt].split(':')[0:-1]:
#			print(l)
			xys = l.split(' ')
#			print(xys)
#			print(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) < Yavg and int(xys[0]) < Xavg):
				NWcols.append(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) >= Yavg and int(xys[0]) < Xavg):
				SWcols.append(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) < Yavg and int(xys[0]) >= Xavg):
				NEcols.append(Z[int(xys[0]),int(xys[1])])
			if( int(xys[1]) >= Yavg and int(xys[0]) >= Xavg):
				SEcols.append(Z[int(xys[0]),int(xys[1])])	
		DAT[pt,0] = Builder.C2I(AvgColors(NWcols))
		DAT[pt,1] = Builder.C2I(AvgColors(SWcols))
		DAT[pt,2] = Builder.C2I(AvgColors(NEcols))
		DAT[pt,3] = Builder.C2I(AvgColors(SEcols))
	
#		DAT[pt,0] = Builder.C2I(Z[max(int(Xavg)-1,0),max(int(Yavg)-1,0)])
#		DAT[pt,1] = Builder.C2I(Z[max(int(Xavg)-1,0),min(int(Yavg)+2,4319)])
#		DAT[pt,2] = Builder.C2I(Z[min(int(Xavg)+2,2159),max(int(Yavg)-1,0)])
#		DAT[pt,3] = Builder.C2I(Z[min(int(Xavg)+2,2159),min(int(Yavg)+2,4319)])

		DAT[pt,4] = Builder.C2I(Z[min(int(Xavg)+1,2159),int(Yavg)])
		DAT[pt,5] = Builder.C2I(Z[int(Xavg),min(int(Yavg)+1,4319)])

		DAT[pt,5] = Builder.C2I(Z[int(Xavg),min(int(Yavg)+1,4319)])

		DAT[pt,6] = Builder.C2I(Z[max(int(Xavg)-1,0),int(Yavg)])
		DAT[pt,7] = Builder.C2I(Z[int(Xavg),max(int(Yavg)-1,0)])

		DAT[pt,8] = Builder.C2I(Z[int(Xavg),int(Yavg)])
		DAT[pt,13] = Builder.C2I(Z[min(int(Xavg)+1,2159),min(int(Yavg)+1,4319)])
		
		DAT[pt,9] = Builder.C2I(Z[min(int(Xavg)+1,2159),max(int(Yavg)-1,0)])

		DAT[pt,10] = Builder.C2I(Z[min(int(Xavg)+2,2159),int(Yavg)])
		DAT[pt,11] = Builder.C2I(Z[min(int(Xavg)+1,2159),min(int(Yavg)+1,4319)])
		DAT[pt,12] = Builder.C2I(Z[min(int(Xavg)+2,2159),min(int(Yavg)+1,4319)])

		DAT[pt,14] = Builder.C2I(Z[int(Xavg),min(int(Yavg)+2,4319)])

		DAT[pt,15] = Builder.C2I(Z[max(int(Xavg)-1,0),min(int(Yavg)+1,4319)])

		if(any(DAT[pt,:] == 0)):
			cols = []
			for i in range(9):
				if(DAT[pt,i] != 0):
					cols.append(Builder.INT2Color(DAT[pt,i]))
			avgcol = AvgColors(cols)
			for i in range(9):
				if(DAT[pt,i] == 0):
					DAT[pt,i] = Builder.C2I(avgcol)
	print('Writing Data...')
	G = open('T2_8.txt','w')
	for i in range(N0):
#		G.write(str(int(DAT[i,0]))+' '+str(int(DAT[i,6]))+' '+str(int(DAT[i,15]))+' '+str(int(DAT[i,2]))+' '+str(int(DAT[i,7]))+' '+str(int(DAT[i,8]))+' '+str(int(DAT[i,5]))+' '+str(int(DAT[i,14]))+' '+str(int(DAT[i,9]))+' '+str(int(DAT[i,4]))+' '+str(int(DAT[i,13]))+' '+str(int(DAT[i,11]))+' '+str(int(DAT[i,1]))+' '+str(int(DAT[i,10]))+' '+str(int(DAT[i,12]))+' '+str(int(DAT[i,3]))+'\n')
		G.write(str(int(DAT[i,0]))+' '+str(int(DAT[i,1]))+' '+str(int(DAT[i,2]))+' '+str(int(DAT[i,3]))+' '+str(int(DAT[i,4]))+' '+str(int(DAT[i,5]))+' '+str(int(DAT[i,6]))+' '+str(int(DAT[i,7]))+' '+str(int(DAT[i,8]))+' '+str(int(DAT[i,9]))+' '+str(int(DAT[i,10]))+' '+str(int(DAT[i,11]))+' '+str(int(DAT[i,12]))+' '+str(int(DAT[i,13]))+' '+str(int(DAT[i,14]))+' '+str(int(DAT[i,15]))+'\n')
	G.close()

#AssignCol2()
#AssignCol3()

def resizeimg(img2):
	XF = 4320
	YF = 2160
	xf = img2.shape[1]
	yf = img2.shape[0]
	col1 = img2[10,10]
	col2 = img2[7950,10]
	lasty = 0
	for i in range(7986):
		if(Builder.C2I(img2[i,10]) == Builder.C2I(col2)):
			lasty = i
	img1 = np.ones(shape=(8046,16092,4),dtype=np.float16)
	ct = 0
	for i in range(8046):
		for j in range(16092):
			if(i < 30):
				img1[i,j,0:3] = col1
			elif(i < (lasty)):
				img1[i,j,0:3] = img2[ct,j,:]
			else:
				img1[i,j,0:3] = col2
		if(i >= 30):
			ct += 1

	plt.imsave(fname='KGCCM_ER.png',arr=img1)


#resizeimg(plt.imread('BigData/KGCCM.png'))
def EnhanceImage(img):
	YF = img.shape[0]
	XF = img.shape[1]
	newimg = np.ones(shape=(YF,XF,4),dtype=np.float16)
	for i in range(YF):
		for j in range(XF):
			if(Builder.C2I(img[i,j]) == Builder.C2I((0.0,0.0,0.0,1.0)) ):
				col1 = img[max(i-1,0),max(j-1,0)]
				col2 = img[max(i-1,0),j]
				col3 = img[max(i-1,0),min(j+1,XF-1)]
				col4 = img[i,max(j-1,0)]
				col5 = img[i,j]
				col6 = img[i,min(j+1,XF-1)]
				col7 = img[min(i+1,YF-1),max(j-1,0)]
				col8 = img[min(i+1,YF-1),j]
				col9 = img[min(i+1,YF-1),min(j+1,XF-1)]
				cols = [col1,col2,col3,col4,col5,col6,col7,col8,col9]
				newcols = []
				for i in range(9):
					if(  Builder.C2I(cols[i]) != Builder.C2I((0.0,0.0,0.0,1.0)) and Builder.C2I(cols[i]) != Builder.C2I((1.0,1.0,1.0,1.0)) ):
						newcols.append(cols[i])
				if(len(newcols) != 0):
					newimg[i,j,0:3] = newcols[0]
				else:
					newimg[i,j] = np.array((0.0,0.0,0.0,1.0))
			else:
				newimg[i,j,0:3] = img[i,j]
	return newimg

#IMG = plt.imread('KGCCM4320.png')
#IMG2 = np.ones(shape=(2160,4320,4),dtype=np.float16)
#coltypes = list(open('colors.txt'))
#for i in range(2160):
#	print(str(i))
#	for j in range(4320):
#		fix = True
#		for c in coltypes:
#			strdat = c.split(' ')
#			if(Builder.C2I(IMG[i,j]) == Builder.C2I((int(strdat[0]),int(strdat[1]),int(strdat[2]))) ):
#				fix = False
#		if(fix):
#			colmix = []
#			colmix.append(IMG[(i-1)%2160,(j-1)%4320])
#			colmix.append(IMG[(i)%2160,(j-1)%4320])
#			colmix.append(IMG[(i+1)%2160,(j-1)%4320])
#			colmix.append(IMG[(i-1)%2160,(j)%4320])
#			colmix.append(IMG[(i)%2160,(j)%4320])
#			colmix.append(IMG[(i+1)%2160,(j)%4320])
#			colmix.append(IMG[(i-1)%2160,(j+1)%4320])
#			colmix.append(IMG[(i)%2160,(j+1)%4320])
#			colmix.append(IMG[(i+1)%2160,(j+1)%4320])
#			thecolor = np.array((0.0,0.0,0.0,1.0))
#			for c0 in colmix:
#				for c in coltypes:
#					strdat = c.split(' ')
#					if(Builder.C2I(c0) == Builder.C2I((int(strdat[0]),int(strdat[1]),int(strdat[2]))) ):
#						thecolor = c
#			IMG2[i,j,0:3] = thecolor[0:3]
#		else:
#			IMG2[i,j,0:3] = IMG[i,j,0:3]
#plt.imsave(fname='kcg.png',arr=IMG)



#for i in range(1):
#	IMG2 = EnhanceImage(IMG)
#	plt.imsave('KGCCM.png',arr=IMG2)

			
IGBPTypes = ('Ocean','Coast','EvergreenNeedleleaf','EvergreenBroadleaf','DeciduousNeedleleaf','DeciduousBroadleaf','MixedForest','ClosedShrubland','OpenShrubland','WoodedSavanna','Savanna','Grasslands','PermanentWetlands','Croplands','Urban','SparseCrops','SnowIce','Desert','Peak')
RGN0Types = ('Flat','Hill','Mountain','Peak')
RGNTypes = list(open('RGNZones.txt'))#('Ocean','Flat','LowHills','MedFlat','MedHills','HighFlat','HighHill','HighMtn','PeakFlat','PeakMtn')
CLIM0Types = ('Grasslands','Plains','Desert','Tundra','Snow')
CLIMTypes = range(29)

Z0 = list(open('Z0_8.txt'))
Z1 = list(open('Z1_8.txt'))
Z2 = list(open('Z2_8.txt'))

ZZ = (Z0,Z1,Z2)


def Learn(res):
	Colors = np.zeros(shape=(len(IGBPTypes),len(RGNTypes),len(CLIMTypes),4**res,4),dtype=np.float16)
	Numbrs = np.zeros(shape=(len(IGBPTypes),len(RGNTypes),len(CLIMTypes)),dtype=np.float16)
	EmptyMind = []
	EmptyThoughts = np.zeros(shape=(4),dtype=np.float16)
	IGBPS = []
	Order = (0,1)
	if(res == 1):
		Order = (0,2,1,3)
	elif(res == 2):
		Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
	ct = 0
	for p in ZZ[res]:
		print(str(ct))
		PT = Builder.Grid[ct]
		if(res == 0):
			Colors[PT.getIGBP(),PT.getRGN(),PT.getCLIM(),0,:] += Builder.INT2Color(int(p))
			Numbrs[PT.getIGBP(),PT.getRGN(),PT.getCLIM()] += 1
		else:
			tmpstr = p.split(' ')
			for k in range(4**res):
				Colors[PT.getIGBP(),PT.getRGN(),PT.getCLIM(),k,:] += Builder.INT2Color(int(tmpstr[Order[k]]))
			Numbrs[PT.getIGBP(),PT.getRGN(),PT.getCLIM()] += 1
		ct += 1
	EmptyMind = np.zeros(shape=(2**res*(len(IGBPTypes)+1),2**res*((len(RGNTypes))*len(CLIMTypes)),4),dtype=np.float16)
	if(res == 0):
		for i in range(len(IGBPTypes)):
			for j in range(len(RGNTypes)):
				for k in range(len(CLIMTypes)):
					if(Numbrs[i,j,k] != 0):
						EmptyMind[i,(len(RGNTypes))*k+j] = Colors[i,j,k,0,:]/Numbrs[i,j,k]
	else:
		for i in range(len(IGBPTypes)):
			for j in range(len(RGNTypes)):
				for k in range(len(CLIMTypes)):
					if(Numbrs[i,j,k] != 0):
						tmpct = 0
						a = np.ones(shape=(2**res,2**res,4),dtype=np.float16)
						for I in range(2**res):
							for J in range(2**res):
								EmptyMind[2**res*i+I,2**res*((len(RGNTypes))*k+j)+J] += Colors[i,j,k,tmpct,:]/Numbrs[i,j,k]
								tmpct += 1
								
				
	plt.imsave(fname='imgcirc/NewTextures'+str(2**res)+'.png',arr=EmptyMind)

#Learn(0)

def LearnB(res):
	Colors = np.zeros(shape=(len(IGBPTypes),len(RGNTypes),len(CLIMTypes),4**res,4),dtype=np.float16)
	Numbrs = np.zeros(shape=(len(IGBPTypes),len(RGNTypes),len(CLIMTypes)),dtype=np.float16)
	EmptyMind = []
	EmptyThoughts = np.zeros(shape=(4),dtype=np.float16)
	IGBPS = []
	Order = (0,1)
	if(res == 1):
		Order = (0,2,1,3)
	elif(res == 2):
		Order = (0,7,9,2,6,8,4,10,15,5,13,12,1,14,11,3)
	ct = 0
	Colors = []
	for i in range(len(IGBPTypes)):
		A = []
		for j in range(len(RGNTypes)):
			B = []
			for k in range(len(CLIMTypes)):
				B.append([])
			A.append(B)
		Colors.append(A)
	for p in ZZ[res]:
		print(str(ct))
		PT = Builder.Grid[ct]
		if(res == 0):
			Colors[PT.getIGBP()][PT.getRGN()][PT.getCLIM()].append( int(p) )
		else:
			tmpstr = p.split(' ')
			for k in range(4**res):
				Colors[PT.getIGBP()][PT.getRGN()][PT.getCLIM()].append( int(tmpstr[Order[k]]) )
		ct += 1
	EmptyMind = np.zeros(shape=(2**res*(len(IGBPTypes)+1),2**res*((len(RGNTypes))*len(CLIMTypes)),4),dtype=np.float16)
	if(res == 0):
		for i in range(len(IGBPTypes)):
			for j in range(len(RGNTypes)):
				for k in range(len(CLIMTypes)):
					if(len(Colors[i][j][k]) != 0):
						EmptyMind[i,(len(RGNTypes))*k+j] = Builder.INT2Color(commonvalues(Colors[i][j][k]))
	else:
		for i in range(len(IGBPTypes)):
			for j in range(len(RGNTypes)):
				for k in range(len(CLIMTypes)):
					tmpct = 0
					a = np.ones(shape=(2**res,2**res,4),dtype=np.float16)
					colIDs = commonvalues0(Colors[i][j][k],4**res)
					for I in range(2**res):
						for J in range(2**res):
							EmptyMind[2**res*i+I,2**res*((len(RGNTypes))*k+j)+J] = Builder.INT2Color(colIDs[tmpct])
							tmpct += 1
								
				
	plt.imsave(fname='imgcirc/NewTextures'+str(2**res)+'.png',arr=EmptyMind)

#LearnB(1)

def BuildPtList():
	print('Initializing...')
	PT = plt.imread('ToXY5400.png')
	SDAT = []
	NDAT = []
	for i in range(len(Builder.Grid)):
		Nstr = str(i)
		if(i < 1000000):
			Nstr = '0'+str(i)
			if(i < 100000):
				Nstr = '00'+str(i)
				if(i < 10000):
					Nstr = '000'+str(i)
					if(i < 1000):
						Nstr = '0000'+str(i)
						if( i < 100):
							Nstr = '00000'+str(i)
							if(i < 10):
								Nstr = '000000'+str(i)
		SDAT.append('')
		NDAT.append(Nstr)
	print('Collecting Data...')
	for i in range(2700):
		for j in range(5400):
			SDAT[Builder.C2I(PT[i,j])] += str(i)+' '+str(j)+':'
	print('Writing Data...')
	G = open('PTList/GrandXYMap8_5400.txt','w')
	for pt in range(len(Builder.Grid)):
#		G = open('PTList/'+NDAT[pt][0]+'/'+NDAT[pt][1]+'/'+NDAT[pt][2]+'/PT_'+NDAT[pt]+'.txt','w')
#		sdat = SDAT[pt].split(':')
#		for s in sdat:
#			G.write(s+'\n')
#		G.close()

		print(str(pt))
		G.write(SDAT[pt]+'\n')
	G.close()			

#BuildPtList()	

def BuildPtListBig():
	print('Initializing...')
	dirct = 0
	ct = 0
	CornerID = ('NE','NW','SE','SW')
	BX = (0,21600,0,21600)
	BY = (0,0,10800,10800)
	S0DAT = []
	for i in range(len(Builder.Grid)):
		S0DAT.append('')
	for k in range(3,4,1):
		print('Collecting Data...    '+str(k))
		SDAT = S0DAT
		ptimg = plt.imread('ToXY_'+CornerID[k]+'.png')
		LastID = 10000000
		for i in range(10800):
			print(str(i))
			for j in range(21600):
				SDAT[Builder.C2I(ptimg[i,j])] += (str(k)+' '+str(i+BY[k])+' '+str(j+BX[k])+':')
		print('Writing Data...    '+str(k))
		for i in range(len(SDAT)):
			print(str(i))
			G = open('PTList/'+str(int(i/1000))+'/'+str(i)+'.txt','a')
			
			if(len(SDAT[i]) == 0):
				G.write('')
			else:
				sdat = SDAT[i].split(':')
				for j in range(len(sdat)-1):
					G.write(sdat[j]+'\n')
			G.close()


#	PT = plt.imread('ToXY_NE.png')
#	for i in range(len(Builder.Grid)):
#		SDAT.append('')
#	print('Collecting Data... (NE)')
#	for i in range(10800):
#		for j in range(21600):
#			SDAT[Builder.C2I(PT[i,j])] += '0 '+str(i)+' '+str(j)+':'
#	print('Collecting Data... (NW)')
#	PT = plt.imread('ToXY_NW.png')
#	for i in range(10800):
#		for j in range(21600):
#			SDAT[Builder.C2I(PT[i,j])] += '1 '+str(i)+' '+str(j)+':'
#	print('Collecting Data... (SE)')
#	PT = plt.imread('ToXY_SE.png')
#	for i in range(10800):
#		for j in range(21600):
#			SDAT[Builder.C2I(PT[i,j])] += '2 '+str(i)+' '+str(j)+':'
#	print('Collecting Data... (SW)')
#	PT = plt.imread('ToXY_SW.png')
#	for i in range(10800):
#		for j in range(21600):
#			SDAT[Builder.C2I(PT[i,j])] += '3 '+str(i)+' '+str(j)+':'

#	print('Writing Data...')
#	G = open('PTList/HugeXYMap8.txt','w')
#	for pt in range(len(Builder.Grid)):


#		G = open('PTList/'+NDAT[pt][0]+'/'+NDAT[pt][1]+'/'+NDAT[pt][2]+'/PT_'+NDAT[pt]+'.txt','w')
#		sdat = SDAT[pt].split(':')
#		for s in sdat:
#			G.write(s+'\n')
#		G.close()

#		print(str(pt))
#		G.write(SDAT[pt]+'\n')
#	G.close()		

		
	
		
	
	

#ReassignHGT()
#print(str(math.cos(((0+0.5)*180/float(YF)-90)*math.pi/180.0)))

def BuildWGT(XF,YF):
	IMGDAT = np.zeros(shape=(YF,XF,4),dtype=np.float16)
	for i in range(int(YF/2)):
		print(str(i))
		for j in range(XF):
			lat = ((i+0.5)*180/float(YF)-90)*math.pi/180.0
			IMGDAT[i,j] = Builder.INT2Color(int(math.cos(lat)/math.cos(((0+0.5)*180/float(YF)-90)*math.pi/180.0)  ))
			IMGDAT[YF-1-i,j] = Builder.INT2Color(int(math.cos(lat)/math.cos(((0+0.5)*180/float(YF)-90)*math.pi/180.0)  ))
	plt.imsave(fname='WGT'+str(XF)+'.png',arr=IMGDAT)


def PointReader(PTID):
	A = list(open('PTList/'+str(int(PTID/1000))+'/'+str(PTID)+'.txt'))
	MAXY = 0
	MAXX = 0
	MINX = 100000
	MINY = 100000
	plane = 0
	for ptdat in A:
		strdat = ptdat.split(' ')
		plane = int(strdat[0])
		X = int(strdat[2])
		Y = int(strdat[1])
		if(X < MINX):
			MINX = X
		if(Y < MINY):
			MINY = Y
		if(X > MAXX):
			MAXX = X
		if(Y > MAXY):
			MAXY = Y
	DY = MAXY - MINY
	DX = MAXX - MINX
	CornerID = ('NE','NW','SE','SW')
	BX = (0,21600,0,21600)
	BY = (0,0,10800,10800)
	B = np.ones(shape=(DY+1,DX+1,4),dtype=np.float16)
#	print(str(DY)+' x '+str(DX))
#	print('MINY: '+str(MINY)+'    MINX: '+str(MINX))
	C = plt.imread('BM_'+CornerID[plane]+'.png')
	for ptdat in A:
		strdat = ptdat.split(' ')
		if(int(strdat[0]) == plane):
#			print(str(int(strdat[1])))
#			print(str(int(strdat[2])))
			B[int(strdat[1])-MINY,int(strdat[2])-MINX,:] = C[int(strdat[1])-BY[plane],int(strdat[2])-BX[plane],:]
	plt.imsave(fname='tmp.png',arr=B)
	


#BuildWGT(4320,2160)

def MLClass(res):
	IMG0 = plt.imread('BigData/BM4320.png')
	IMG1 = plt.imread('BigData/KoppenClimateData.png')
	IMG2 = plt.imread('BigData/NewIGBP.png')
	IMG3 = plt.imread('ToXY7.png')
#	img1set = []
#	img2set = []
#	print('Initializing...')
#	for i in range(2160):
#		for j in range(4320):
#			img1set.append(Builder.C2I(IMG1[i,j]))
#			img1set = list(set(img1set))
#			img2set.append(Builder.C2I(IMG2[i,j]))
#			img2set = list(set(img2set))
#	len1 = len(img1set)
#	len2 = len(img2set)
#	col1 = []
#	for i in range(len1):
#			col1.append(Builder.INT2Color(img1set[i]))
#	col2 = []
#	for i in range(len2):
#			col2.append(Builder.INT2Color(img2set[i]))
#	dict1 = {}
#	dict2 = {}
#	col1 = np.zeros(shape=(256*256*256),dtype=np.float32)
#	col2 = np.zeros(shape=(256*256*256),dtype=np.float32)
	CC = list(open('ClimateCols.txt'))
	len1 = len(CC)
	IC = list(open('IGBPCols.txt'))
	len2 = len(IC)
#	for i in range(len1):
#		dict1[img1set[i]] = i
#		dict1[i] = img1set[i]
#		col1[img1set[i]] = i
#		g.write(str(img1set[i])+'\n')
#	g.close()
#	g = open('IGBPCols.txt','w')
#	for i in range(len2):
#		dict2[img2set[i]] = i
#		dict2[i] = img2set[i]
#		col2[img2set[i]] = i
#		g.write(str(img2set[i])+'\n')
#	g.close()
#	dict2 = {}
	Tiles = np.zeros(shape=(len1,len2,4),dtype=np.float32)
	Numbs = np.zeros(shape=(len1,len2),dtype=np.float32)
	Order = (0,)
	if(res == 1):
		Order = (0,1,2,3)
	if(res == 2):
		Order = (0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
	print('Collecting Data...')
	for i in range(2160):
		print(str(i))
		for j in range(4320):
			Cct = 0
			for k in range(len(CC)):
				if(int(CC[k]) == Builder.C2I(IMG1[i,j])):
					Cct = k
			Ict = 0
			for k in range(len(IC)):
				if(int(IC[k]) == Builder.C2I(IMG2[i,j])):
					Ict = k
			
			Tiles[Cct,Ict,0:3] += IMG0[i,j,0:3]
			Numbs[Cct,Ict] += 1
	ActualTiles = np.ones(shape=(len1,len2,4),dtype=np.float32)
	print('Writing Data...')
	for i in range(len1):
		for j in range(len2):
			if(Numbs[i,j] != 0):
				ActualTiles[i,j,0:3] = Tiles[i,j,0:3]/Numbs[i,j]
	plt.imsave(fname='imgcirc/Textures'+str(res)+'.png',arr=ActualTiles)
#	N0 = len(Builder.Grid)
#	gc = open('CZ8.txt')
#	gi = open('IZ8.txt')	
#	for i in range(N0):
#		gc.write(' ')
#	gc.close()
#	gi.close()

def UpdateClimIGBP():
	IMG1 = plt.imread('BigData/KoppenClimateData.png')
	IMG2 = plt.imread('BigData/NewIGBP.png')
#	IMG3 = plt.imread('ToXY7.png')
	PTL = list(open('PTList/GrandXYMap8.txt'))
	N0 = len(Builder.Grid)
	ClimateCols = list(open('ClimateCols.txt'))
#	IGBPCols = list(open('IGBPCols.txt'))
	for pt in range(N0):
		print(str(pt))
		a = Builder.Grid[pt]
		Ccols = []
#		Icols = []
		for l in PTL[pt].split(':')[0:-1]:
#			print(str(l))
			strdat = l.split(' ')
			Ccols.append(IMG1[int(strdat[0]),int(strdat[1])])
#			Icols.append(IMG2[int(strdat[0]),int(strdat[1])])
		ccol = Builder.C2I(commoncolors(Ccols))
#		icol = Builder.C2I(commoncolors(Icols))
		Cct = 0
#		Ict = 0
		for i in range(len(ClimateCols)):
			if(int(ClimateCols[i]) == ccol):
				Cct = i
#		for i in range(len(IGBPCols)):
#			if(int(IGBPCols[i]) == icol):
#				Ict = i
		Builder.SetCLIM(pt,Cct)
#		Builder.SetIGBP(pt,Ict)
	Builder.SaveGrid()

#UpdateClimIGBP()

def MLClass2(res):
	IMG0 = plt.imread('BigData/BM4320.png')
	PTL = list(open('PTList/GrandXYMap8.txt'))
	N0 = len(Builder.Grid)
	ClimateCols = list(open('ClimateCols.txt'))
	IGBPCols = range(0,17,1)#list(open('IGBPCols.txt'))
	Cols = np.zeros(shape=(len(ClimateCols),len(IGBPCols),4,4),dtype=np.float32)
	Nums = np.zeros(shape=(len(ClimateCols),len(IGBPCols),4),dtype=np.float32)
	for pt in range(N0):
		print(str(pt))
		a = Builder.Grid[pt]
		Ccols = []
		Icols = []
		a = Builder.Grid[pt]
		for l in PTL[pt].split(':')[0:-1]:
#			print(str(l))
			strdat = l.split(' ')
			Cols[a.getCLIM(),a.getIGBP(),a.getRGN(),0:3] += IMG0[int(strdat[0]),int(strdat[1])][0:3]
			Nums[a.getCLIM(),a.getIGBP(),a.getRGN()] += 1
	for i in range(len(ClimateCols)):
		outimg = np.ones(shape=(4,len(IGBPCols),4),dtype=np.float16)
		for j in range(len(IGBPCols)):
			for k in range(4):
				if(Nums[i,j,k] != 0):
					outimg[k,j,0:3] = Cols[i,j,k,0:3]/Nums[i,j,k]
		plt.imsave('imgcirc/BM0/Zoom'+str(res)+'/Clim'+str(i)+'.png',arr=outimg)


def MLClass3(res):
	PTCOLS = list(open('Z'+str(res)+'_8.txt'))
	IMG0 = plt.imread('BigData/BM4320.png')
	PTL = list(open('PTList/GrandXYMap8.txt'))
	N0 = len(Builder.Grid)
	ClimateCols = list(open('ClimateCols.txt'))
	IGBPCols = range(0,19,1)#list(open('IGBPCols.txt'))
	Cols = np.zeros(shape=(len(ClimateCols),len(IGBPCols),4,2**res,2**res,4),dtype=np.float32)
	Nums = np.zeros(shape=(len(ClimateCols),len(IGBPCols),4),dtype=np.float32)
	for pt in range(N0):
		print(str(pt))
		ptcols = PTCOLS[pt].split(' ')
		ct = 0
		a = Builder.Grid[pt]
		for i in range(2**res):
			for j in range(2**res):
				Cols[a.getCLIM(),a.getIGBP(),a.getRGN(),i,j,:] += Builder.INT2Color(int(ptcols[ct]))
				ct += 1
		Nums[a.getCLIM(),a.getIGBP(),a.getRGN()] += 1
	for i in range(len(ClimateCols)):
		outimg = np.ones(shape=(4*2**res,len(IGBPCols)*2**res,4),dtype=np.float16)
		for j in range(len(IGBPCols)):
			for k in range(4):
				if(Nums[i,j,k] != 0):
					for I in range(2**res):
						for J in range(2**res):
							outimg[2**res*k+I,2**res*j+J,0:3] = Cols[i,j,k,I,J,0:3]/Nums[i,j,k]
		plt.imsave('imgcirc/BM0/Zoom'+str(res)+'/Clim'+str(i)+'.png',arr=outimg)

#UpdateClimIGBP()

def FillImage(plane,res):
	print('Filling Plane '+str(plane)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	OLD = plt.imread('Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float32)
	for i in range(YF):
		for j in range(XF):
			if(Builder.C2I(OLD[i,j]) == 0):
				cols = []
				if(Builder.C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
					cols.append(OLD[max(i-1,0),(j-1)%XF,:])
				if(Builder.C2I(OLD[max(i-1,0),j,:]) != 0):
					cols.append(OLD[max(i-1,0),j,:])
				if(Builder.C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
					cols.append(OLD[max(i-1,0),(j+1)%XF,:])

				if(Builder.C2I(OLD[i,(j-1)%XF,:]) != 0):
					cols.append(OLD[i,(j-1)%XF,:])
				if(Builder.C2I(OLD[i,(j+1)%XF,:]) != 0):
					cols.append(OLD[i,(j+1)%XF,:])

				if(Builder.C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
				if(Builder.C2I(OLD[min(i+1,YF-1),j,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),j,:])
				if(Builder.C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
					cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
				NEW[i,j,0:3] = AvgColors(cols)[0:3]
			else:
				NEW[i,j,0:3] = OLD[i,j,0:3]
	plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png',arr=NEW)

def EnhanceImage(plane,res):
	print('Enhancing Plane '+str(plane)+' at resolution '+str(res))
	XF = XSize[res]
	YF = YSize[res]
	OLD = plt.imread('Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png')
	NEW = np.ones(shape=(YF,XF,4),dtype=np.float32)
	for i in range(YF):
		print(str(i))
		for j in range(XF):
			cols = []
			if(Builder.C2I(OLD[max(i-1,0),(j-1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j-1)%XF,:])
			if(Builder.C2I(OLD[max(i-1,0),j,:]) != 0):
				cols.append(OLD[max(i-1,0),j,:])
			if(Builder.C2I(OLD[max(i-1,0),(j+1)%XF,:]) != 0):
				cols.append(OLD[max(i-1,0),(j+1)%XF,:])
			if(Builder.C2I(OLD[i,(j-1)%XF,:]) != 0):
				cols.append(OLD[i,(j-1)%XF,:])
			if(Builder.C2I(OLD[i,(j+1)%XF,:]) != 0):
				cols.append(OLD[i,(j+1)%XF,:])
			if(Builder.C2I(OLD[min(i+1,YF-1),(j-1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j-1)%XF,:])
			if(Builder.C2I(OLD[min(i+1,YF-1),j,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),j,:])
			if(Builder.C2I(OLD[min(i+1,YF-1),(j+1)%XF,:]) != 0):
				cols.append(OLD[min(i+1,YF-1),(j+1)%XF,:])
			
			NEW[i,j,0:3] = (len(cols)+8)*OLD[i,j,0:3]
			for k in range(len(cols)):
				NEW[i,j,0:3] = NEW[i,j,0:3]-cols[k][0:3]
			NEW[i,j,0:3] = NEW[i,j,0:3]/8.0
			for k in range(3):
				NEW[i,j,k] = max(min(NEW[i,j,k],1.0),0)	
	plt.imsave(fname='Assets32/Zoom'+str(res)+'/Plane'+str(plane)+'.png',arr=NEW)

#MLClass3(2)

def BuildHeatMap():
	print('Initializing...')
	YF = 25
	XF = 17
	DATZ = plt.imread('ZIMG.png')
	DATV = plt.imread('VIMG.png')
	N0 = len(Builder.Grid)
	HEATMAP = np.zeros(shape=(YF,XF),dtype=np.float32)
	HEATMAPV = np.ones(shape=(2160,4320,4),dtype=np.float16)
	print('Collecting Data...')
	mcol = 0
	for i in range(2160):
		for j in range(4320):
			HEATMAP[int(2*math.log(1+Builder.C2I(DATZ[i,j])/2.0)),int(math.log(1+Builder.C2I(DATV[i,j])/2.0))] += 1.0
			HEATMAPV[i,j,:] = np.array((1.0,1.0-int(Builder.C2I(DATZ[i,j])/512)/255.0,1.0-int(Builder.C2I(DATZ[i,j])/512)/255.0,1.0))
	plt.imsave(fname='ZMap2.png',arr=HEATMAPV)
	HEATIMG = np.ones(shape=(YF,XF,4),dtype=np.float16)
	print('Writing Data...')
	mcol = 100000
	for i in range(YF):
		for j in range(XF):
			print(str(HEATMAP[i,j]/mcol))
			HEATIMG[i,j,:] = np.array((1.0,max(0,1.0-(HEATMAP[i,j]/mcol)),max(0,1.0-(HEATMAP[i,j]/mcol)),1.0))
	plt.imsave(fname='HEATMAP.png',arr=HEATIMG)

#BuildHeatMap()
def UpdateHGTRGN():
	IMG1 = plt.imread('ZIMG.png')
	IMG3 = plt.imread('CRD.png')
#	IMG2 = plt.imread('PIMG.png')
	WGT = plt.imread('WGT4320.png')
	PTL = list(open('PTList/GrandXYMap8.txt'))
	N0 = len(Builder.Grid)
	ClimateCols = list(open('ClimateCols.txt'))
#	IGBPCols = list(open('IGBPCols.txt'))
	for pt in range(N0):
		print(str(pt))
		a = Builder.Grid[pt]
		Zcols = 0.0
		Vcols = 0.0
		nums = 0.0		
#		pops = 0.0
		for l in PTL[pt].split(':')[0:-1]:
#			print(str(l))
			strdat = l.split(' ')
			Zcols += 1.0*Builder.C2I(IMG1[int(strdat[0]),int(strdat[1])])*Builder.C2I(WGT[int(strdat[0]),int(strdat[1])])
			Vcols += 1.0*Builder.C2I(IMG3[int(strdat[0]),int(strdat[1])])*Builder.C2I(WGT[int(strdat[0]),int(strdat[1])])
			nums += Builder.C2I(WGT[int(strdat[0]),int(strdat[1])])
#			pops += Builder.C2I(IMG2[int(strdat[0]),int(strdat[1])])
		Zavg = 0
		Vavg = 0
#		print('    '+str(max(int(math.log(pops+1)-4),0)))
		if(nums != 0):
			Zavg = Zcols/nums
			Vavg = Vcols/nums
		Zbins = (-100000,-50,0,20,100,500,3000,8000,15000,20000,200000)
		Vbins = (-1,1,3,8,20,55,400,3000,22000,162000,100000000)
		Zints = 0
		Vints = 0
		for i in range(11):
			if(Zavg >= Zbins[i] and Zavg < Zbins[i+1]):
				Zints = i
			if(Vavg >= Vbins[i] and Vavg < Vbins[i+1]):
				Vints = i
		Act = 0
		if(Zavg < 100):
			Act = 2
		elif(Zavg < 1000):
			if(Vavg < 20):
				Act = 2
			elif(Vavg < 10000):
				Act = 3
			else:
				Act = 2
		elif(Zavg < 7000):
			if(Vavg < 20):
				Act = 2
			elif(Vavg < 100):
				Act = 3
			elif(Vavg < 10000):
				Act = 4
			else:
				Act = 2
		elif(Zavg < 35000):
			if(Vavg < 20):
				Act = 5
			elif(Vavg < 100):
				Act = 6
			elif(Vavg < 200):
				Act = 7
			elif(Vavg < 10000):
				Act = 8
			else:
				Act = 5
		elif(Zavg < 55000):
			if(Vavg < 20):
				Act = 9
			elif(Vavg < 100):
				Act = 10
			elif(Vavg < 200):
				Act = 11
			elif(Vavg < 10000):
				Act = 12
			else:
				Act = 9
		else:
			if(Vavg < 20):
				Act = 13
			elif(Vavg < 100):
				Act = 14
			elif(Vavg < 10000):
				Act = 15
			else:
				Act = 13
#		popl = max(int(math.log(pops+1)),0)		
		Builder.SetRGN(pt,Act)
#		Builder.SetPOPL(pt,popl)
	Builder.SaveGrid()

#UpdateHGTRGN()

def ExtractPopDat():
	PIMG = np.ones(shape=(4320,8640,4),dtype=np.float16)
	IMG0 = plt.imread('BigData/POPUNALTERED.png')
	POPBIN0 = (0,0.5,1,3,5,10,25,50,100,150,200,300,400,500,600,800,1000,1250,1500,1750,2000,2500,3000,4000,6000,8000,10000,25000,50000,100000)
	POPBIN = (0,0.25,0.75,2,4,7.5,17.5,37.5,75,125,175,250,350,450,550,700,900,1125,1375,1625,1875,2250,2750,3500,5000,7000,9000,17500,37500,60000)
	print(len(POPBIN))
	popcols = []
	g = open('PopCols.txt','w')
	for i in range(1245,2885,55):
		popcols.append(Builder.C2I(IMG0[i,300]))
		g.write(str(Builder.C2I(IMG0[i,300]))+'\n')
	g.close()
#	IMG0 = plt.imread('BigData/POPRAW.png')
#	ct = 0
#	I = 0
#	J = 0
#	for i in range(4320):
#		print(str(i))
#		I = int(i/2.0)
#		for j in range(8640):
#			J = int(j/2.0)
#			if(i < 119):
#				PIMG[i,j] = np.array((1.0,1.0,1.0,1.0))
#			elif(ct < 3431):
#				col = Builder.C2I(IMG0[ct,j])
#				colnum = 0
#				for k in range(len(popcols)):
#					if(col == popcols[k]):
#						colnum = k
#				PIMG[i,j] = IMG0[ct,j]#Builder.INT2Color(int(POPBIN[k]))
#			else:
#				PIMG[i,j] = np.array((1.0,1.0,1.0,1.0))
#		if(i >= 119):
#			ct += 1
#	plt.imsave(fname='PIMG.png',arr=PIMG)

#ExtractPopDat()

def EnhancePIMG():
	print('Initializing...')
	PIMG = plt.imread('PIMG0.png')
	popcols = list(open('PopCols.txt'))
	IMG0 = plt.imread('PIMG0.png')
	IMG1 = np.zeros(shape=(2160,4320),dtype=np.float32)
	POPBIN = (0,0.25,0.75,2,4,7.5,17.5,37.5,75,125,175,250,350,450,550,700,900,1125,1375,1625,1875,2250,2750,3500,5000,7000,9000,17500,37500,60000)
	WGT = plt.imread('WGT4320.png')
	maxwgt = 0
	print('Collecting Data...')
	for i in range(4320):
		print(str(i/43.2)+'%')
		I = int(i/2.0)
		for j in range(8640):
			J = int(j/2.0)	
			popcol = 0
			for k in range(len(popcols)):
				if(int(popcols[k]) == Builder.C2I(PIMG[i,j]) ):
					popcol = k
			IMG1[I,J] += POPBIN[int(popcol)]*Builder.C2I(WGT[I,J])
			if(maxwgt < Builder.C2I(WGT[I,J])):
				maxwgt = Builder.C2I(WGT[I,J])
	if(maxwgt == 0):
		print('help')
	else:
		IMG1 = IMG1/float(maxwgt)
	print(str(sum(IMG1) ) )
	IMG0 = np.ones(shape=(2160,4320,4),dtype=np.float16)
	print('Writing Data...')
	for i in range(2160):
		print(str(i/21.6)+'%')
		for j in range(4320):
			IMG0[i,j] = Builder.INT2Color(int(IMG1[i,j]))
	plt.imsave(fname='PIMG.png',arr=IMG0)
		


def NC2TIS_CLD(y0, m0, d0):
	CLDS0 = Dataset('../../../Data/MERRA/2D/CLD/MERRA2_100.tavg1_2d_rad_Nx.'+str(y0)+str(m0)+str(d0)+'.SUB.nc')
	CLDS = CLDS0.variables['CLDTOT'][:]
	latnum = 361
	lonnum = 576
	CLDOUT = np.zeros(shape=(len(Builder.Grid),24),dtype=np.uint8)
	outF = 'TISDat/CLD/CLDS'+y0+m0+d0+'.npy'
	ct = 0
#	Month = (0,31,28,31,30,31,30,31,31,30,31,30,31)
	CumMonth = np.array((0,31,59,90,120,151,181,212,243,273,304,334,365))
	if(int(y0)%4 == 0):
		CumMonth[2:] += 1
	toy = np.sum(CumMonth[0:(int(m0)-1)])+int(d0)

	SolarNoons = []
	for k in range(24):
		SolarNoons.append(Builder.FindSolarNoon(toy,k))

	for pt in Builder.Grid:
		print(str(ct))
		lat = pt.getLatDeg()
		lon = pt.getLonDeg()
		xyz = pt.getXYZPos()
		hrs = range(0,16392,683)
		for k in range(24):
			h0 = 25
			d0 = Builder.dSurf(xyz,SolarNoons[k].getXYZPos())
			for i in range(len(hrs)-1):
				if(d0*10000.0 >= hrs[i] and d0*10000.0 < hrs[i+1]):
					h0 = i
			angdat = h0
			clddat = (round(9*CLDS[k,int(360*(lat+90)/180),int(575*((lon+180)%360)/360)]))
			CLDOUT[ct,k] = (clddat*25+angdat)
		ct += 1
	np.save(outF,arr=CLDOUT)

def D8toCLD(uint8):
	cloudval = int(uint8/25)
	angle = uint8%25
	return cloudval, angle

def CLDtoD8(cloudval, angle):
	return np.uint8(cloudval*25+angle)

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
	
	

def VegIDX():
	print('hello')

def GenSmallCIRV():
	output = np.zeros(shape=(len(Builder.Grid)),dtype=np.uint16)
	ct = 0
	for pt in Builder.Grid:
		print(str(ct))
		clim = pt.getCLIM()
		igbp = pt.getIGBP()
		rgn = pt.getRGN()
		veg = max(0,pt.getPOPL()-6)
		output[ct] = CIRVtoS16(clim,igbp,rgn,veg)
		ct += 1
	np.save(file='CIRV.npy',arr=output)

def Prep4HugeT():
	for pt in range(0,len(Builder.Grid),1):
		print(str(pt))
		FILE = str(int(pt/1000))
		NAME = str(pt)
		try:
			os.system('mkdir '+FILE)
			PTS = list(open('PTList/'+FILE+'/'+NAME+'.txt'))
			PTS0 = np.zeros(shape=(len(PTS),2),dtype=np.uint16)
			ct = 0
			for p in PTS:
				sdat = p.split(' ')
				PTS0[ct,0] = int(sdat[2])
				PTS0[ct,1] = int(sdat[1])
				ct += 1
			np.save('PTList/'+FILE+'/'+NAME+'.npy',arr=PTS0)
		except:
			PTS = list(open('PTList/'+FILE+'/'+NAME+'.txt'))
			PTS0 = np.zeros(shape=(len(PTS),2),dtype=np.uint16)
			ct = 0
			for p in PTS:
				sdat = p.split(' ')
				PTS0[ct,0] = int(sdat[2])
				PTS0[ct,1] = int(sdat[1])
				ct += 1
			np.save('PTList/'+FILE+'/'+NAME+'.npy',arr=PTS0)

def Prep4Huge():
	for pt in range(1431000,len(Builder.Grid),1):
		print(str(pt))
		FILE = str(int(pt/1000))
		NAME = str(pt)
		PTS = list(open('PTList/'+FILE+'/'+NAME+'.txt'))
		PTS0 = np.zeros(shape=(len(PTS),2),dtype=np.uint16)
		ct = 0
		for p in PTS:
			sdat = p.split(' ')
			PTS0[ct,0] = int(sdat[2])
			PTS0[ct,1] = int(sdat[1])
			ct += 1
		np.save('PTList/'+FILE+'/'+NAME+'.npy',arr=PTS0)

def HugeBuild(res):
	COLS0 = []
	for i in range(16):
		col0 = []
		for j in range(2**res):
			col0.append([])
		COLS0.append(col0)
	if(res == 2):
		dat = np.zeros(shape=(len(Builder.Grid),4**res),dtype=int)
		D = np.load('BM.npy')
#		G = open('PTList/TrueTex'+str(res)+'.txt')
		for ct in range(len(Builder.Grid)):
			print(str(ct))
			if(ct == 80 or ct == 83):
				dat[ct] = Builder.C2I((1.0,0.0,0.0,1.0))
			else:
				FILE = str(int(ct/1000))
				NAME = str(ct)
				PTS = np.load('PTList/'+FILE+'/'+NAME+'.npy')
				x2 = int(np.average(PTS[:,0]))
				x1 = int(0.5*(min(PTS[:,0])+x2))
				x3 = int(0.5*(max(PTS[:,0])+x2))
				y2 = int(np.average(PTS[:,1]))
				y1 = int(0.5*(min(PTS[:,1])+y2))
				y3 = int(0.5*(max(PTS[:,1])+y2))
				A1 = []
				B1 = []
				C1 = []
				D1 = []
				A2 = []
				B2 = []
				C2 = []
				D2 = []
				A3 = []
				B3 = []
				C3 = []
				D3 = []
				A4 = []
				B4 = []
				C4 = []
				D4 = []					
				for pt in range(len(PTS[:,0])):
					if(PTS[pt,0] <= x1):
						if(PTS[pt,1] <= y1):
							A1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C1.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							D1.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x2):
						if(PTS[pt,1] <= y1):
							A2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C2.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							D2.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x3):
						if(PTS[pt,1] <= y1):
							A3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C3.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							D3.append(D[PTS[pt,1],PTS[pt,0]])
					else:
						if(PTS[pt,1] <= y1):
							A4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C4.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							D4.append(D[PTS[pt,1],PTS[pt,0]])

				dat[ct,0] = Builder.C2I(Builder.AvgColors(A1)/255.0)
				dat[ct,1] = Builder.C2I(Builder.AvgColors(B1)/255.0)
				dat[ct,2] = Builder.C2I(Builder.AvgColors(C1)/255.0)
				dat[ct,3] = Builder.C2I(Builder.AvgColors(D1)/255.0)
				dat[ct,4] = Builder.C2I(Builder.AvgColors(A2)/255.0)
				dat[ct,5] = Builder.C2I(Builder.AvgColors(B2)/255.0)
				dat[ct,6] = Builder.C2I(Builder.AvgColors(C2)/255.0)
				dat[ct,7] = Builder.C2I(Builder.AvgColors(D2)/255.0)
				dat[ct,8] = Builder.C2I(Builder.AvgColors(A3)/255.0)
				dat[ct,9] = Builder.C2I(Builder.AvgColors(B3)/255.0)
				dat[ct,10] = Builder.C2I(Builder.AvgColors(C3)/255.0)
				dat[ct,11] = Builder.C2I(Builder.AvgColors(D3)/255.0)
				dat[ct,12] = Builder.C2I(Builder.AvgColors(A4)/255.0)
				dat[ct,13] = Builder.C2I(Builder.AvgColors(B4)/255.0)
				dat[ct,14] = Builder.C2I(Builder.AvgColors(C4)/255.0)
				dat[ct,15] = Builder.C2I(Builder.AvgColors(D4)/255.0)
		np.save('imgcirc/Textures'+str(res)+'.npy',arr=dat)
	if(res == 3):
		dat = np.zeros(shape=(len(Builder.Grid),4**res),dtype=int)
		D = np.load('BM.npy')
		for ct in range(len(Builder.Grid)):
			print(str(ct))
			if(ct == 80 or ct == 83):
				dat[ct] = Builder.C2I((1.0,0.0,0.0,1.0))
			else:
				FILE = str(int(ct/1000))
				NAME = str(ct)
				PTS = np.load('PTList/'+FILE+'/'+NAME+'.npy')

				x4 = int(np.average(PTS[:,0]))
				x2 = int(0.5*(min(PTS[:,0])+x4))
				x1 = int(0.5*(min(PTS[:,0])+x2))
				x3 = int(0.5*(x2+x4))
				x6 = int(0.5*(max(PTS[:,0])+x4))
				x8 = int(0.5*(max(PTS[:,0])+x6))
				x5 = int(0.5*(x4+x6))
				x7 = int(0.5*(x6+x8))

				y4 = int(np.average(PTS[:,1]))
				y2 = int(0.5*(min(PTS[:,1])+y4))
				y1 = int(0.5*(min(PTS[:,1])+y2))
				y3 = int(0.5*(y2+y4))
				y6 = int(0.5*(max(PTS[:,1])+y4))
				y8 = int(0.5*(max(PTS[:,1])+y6))
				y5 = int(0.5*(y4+y6))
				y7 = int(0.5*(y6+y8))


				A1 = []
				B1 = []
				C1 = []
				D1 = []
				E1 = []
				F1 = []
				G1 = []
				H1 = []

				A2 = []
				B2 = []
				C2 = []
				D2 = []
				E2 = []
				F2 = []
				G2 = []
				H2 = []

				A3 = []
				B3 = []
				C3 = []
				D3 = []
				E3 = []
				F3 = []
				G3 = []
				H3 = []

				A4 = []
				B4 = []
				C4 = []
				D4 = []
				E4 = []
				F4 = []
				G4 = []
				H4 = []

				A5 = []
				B5 = []
				C5 = []
				D5 = []
				E5 = []
				F5 = []
				G5 = []
				H5 = []

				A6 = []
				B6 = []
				C6 = []
				D6 = []
				E6 = []
				F6 = []
				G6 = []
				H6 = []

				A7 = []
				B7 = []
				C7 = []
				D7 = []
				E7 = []
				F7 = []
				G7 = []
				H7 = []

				A8 = []
				B8 = []
				C8 = []
				D8 = []
				E8 = []
				F8 = []
				G8 = []
				H8 = []

				
				for pt in range(len(PTS[:,0])):
					if(PTS[pt,0] <= x1):
						if(PTS[pt,1] <= y1):
							A1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G1.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H1.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x2):
						if(PTS[pt,1] <= y1):
							A2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G2.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H2.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x3):
						if(PTS[pt,1] <= y1):
							A3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G3.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H3.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x4):
						if(PTS[pt,1] <= y1):
							A4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G4.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H4.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x5):
						if(PTS[pt,1] <= y1):
							A5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G5.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H5.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x6):
						if(PTS[pt,1] <= y1):
							A6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G6.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H6.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x7):
						if(PTS[pt,1] <= y1):
							A7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G7.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H7.append(D[PTS[pt,1],PTS[pt,0]])
					else:
						if(PTS[pt,1] <= y1):
							A8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G8.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							H8.append(D[PTS[pt,1],PTS[pt,0]])


				dat[ct,0] = Builder.C2I(Builder.AvgColors(A1)/255.0)
				dat[ct,1] = Builder.C2I(Builder.AvgColors(B1)/255.0)
				dat[ct,2] = Builder.C2I(Builder.AvgColors(C1)/255.0)
				dat[ct,3] = Builder.C2I(Builder.AvgColors(D1)/255.0)
				dat[ct,4] = Builder.C2I(Builder.AvgColors(E1)/255.0)
				dat[ct,5] = Builder.C2I(Builder.AvgColors(F1)/255.0)
				dat[ct,6] = Builder.C2I(Builder.AvgColors(G1)/255.0)
				dat[ct,7] = Builder.C2I(Builder.AvgColors(H1)/255.0)

				dat[ct,8] = Builder.C2I(Builder.AvgColors(A2)/255.0)
				dat[ct,9] = Builder.C2I(Builder.AvgColors(B2)/255.0)
				dat[ct,10] = Builder.C2I(Builder.AvgColors(C2)/255.0)
				dat[ct,11] = Builder.C2I(Builder.AvgColors(D2)/255.0)
				dat[ct,12] = Builder.C2I(Builder.AvgColors(E2)/255.0)
				dat[ct,13] = Builder.C2I(Builder.AvgColors(F2)/255.0)
				dat[ct,14] = Builder.C2I(Builder.AvgColors(G2)/255.0)
				dat[ct,15] = Builder.C2I(Builder.AvgColors(H2)/255.0)

				dat[ct,16] = Builder.C2I(Builder.AvgColors(A3)/255.0)
				dat[ct,17] = Builder.C2I(Builder.AvgColors(B3)/255.0)
				dat[ct,18] = Builder.C2I(Builder.AvgColors(C3)/255.0)
				dat[ct,19] = Builder.C2I(Builder.AvgColors(D3)/255.0)
				dat[ct,20] = Builder.C2I(Builder.AvgColors(E3)/255.0)
				dat[ct,21] = Builder.C2I(Builder.AvgColors(F3)/255.0)
				dat[ct,22] = Builder.C2I(Builder.AvgColors(G3)/255.0)
				dat[ct,23] = Builder.C2I(Builder.AvgColors(H3)/255.0)

				dat[ct,24] = Builder.C2I(Builder.AvgColors(A4)/255.0)
				dat[ct,25] = Builder.C2I(Builder.AvgColors(B4)/255.0)
				dat[ct,26] = Builder.C2I(Builder.AvgColors(C4)/255.0)
				dat[ct,27] = Builder.C2I(Builder.AvgColors(D4)/255.0)
				dat[ct,28] = Builder.C2I(Builder.AvgColors(E4)/255.0)
				dat[ct,29] = Builder.C2I(Builder.AvgColors(F4)/255.0)
				dat[ct,30] = Builder.C2I(Builder.AvgColors(G4)/255.0)
				dat[ct,31] = Builder.C2I(Builder.AvgColors(H4)/255.0)	

				dat[ct,32] = Builder.C2I(Builder.AvgColors(A5)/255.0)
				dat[ct,33] = Builder.C2I(Builder.AvgColors(B5)/255.0)
				dat[ct,34] = Builder.C2I(Builder.AvgColors(C5)/255.0)
				dat[ct,35] = Builder.C2I(Builder.AvgColors(D5)/255.0)
				dat[ct,36] = Builder.C2I(Builder.AvgColors(E5)/255.0)
				dat[ct,37] = Builder.C2I(Builder.AvgColors(F5)/255.0)
				dat[ct,38] = Builder.C2I(Builder.AvgColors(G5)/255.0)
				dat[ct,39] = Builder.C2I(Builder.AvgColors(H5)/255.0)	

				dat[ct,40] = Builder.C2I(Builder.AvgColors(A6)/255.0)
				dat[ct,41] = Builder.C2I(Builder.AvgColors(B6)/255.0)
				dat[ct,42] = Builder.C2I(Builder.AvgColors(C6)/255.0)
				dat[ct,43] = Builder.C2I(Builder.AvgColors(D6)/255.0)
				dat[ct,44] = Builder.C2I(Builder.AvgColors(E6)/255.0)
				dat[ct,45] = Builder.C2I(Builder.AvgColors(F6)/255.0)
				dat[ct,46] = Builder.C2I(Builder.AvgColors(G6)/255.0)
				dat[ct,47] = Builder.C2I(Builder.AvgColors(H6)/255.0)	

				dat[ct,48] = Builder.C2I(Builder.AvgColors(A7)/255.0)
				dat[ct,49] = Builder.C2I(Builder.AvgColors(B7)/255.0)
				dat[ct,50] = Builder.C2I(Builder.AvgColors(C7)/255.0)
				dat[ct,51] = Builder.C2I(Builder.AvgColors(D7)/255.0)
				dat[ct,52] = Builder.C2I(Builder.AvgColors(E7)/255.0)
				dat[ct,53] = Builder.C2I(Builder.AvgColors(F7)/255.0)
				dat[ct,54] = Builder.C2I(Builder.AvgColors(G7)/255.0)
				dat[ct,55] = Builder.C2I(Builder.AvgColors(H7)/255.0)	

				dat[ct,56] = Builder.C2I(Builder.AvgColors(A8)/255.0)
				dat[ct,57] = Builder.C2I(Builder.AvgColors(B8)/255.0)
				dat[ct,58] = Builder.C2I(Builder.AvgColors(C8)/255.0)
				dat[ct,59] = Builder.C2I(Builder.AvgColors(D8)/255.0)
				dat[ct,60] = Builder.C2I(Builder.AvgColors(E8)/255.0)
				dat[ct,61] = Builder.C2I(Builder.AvgColors(F8)/255.0)
				dat[ct,62] = Builder.C2I(Builder.AvgColors(G8)/255.0)
				dat[ct,63] = Builder.C2I(Builder.AvgColors(H8)/255.0)

		np.save('imgcirc/Textures'+str(res)+'.npy',arr=dat)

	if(res == 4):

		dat = np.zeros(shape=(len(Builder.Grid),4**res),dtype=int)
#		dat = np.zeros(shape=(100,4**res),dtype=int)
		D = np.load('BM.npy')
		for ct in range(len(Builder.Grid)):
			print(str(ct))
			if(ct == 80 or ct == 83):
				dat[ct] = Builder.C2I((1.0,0.0,0.0,1.0))
			else:
				FILE = str(int(ct/1000))
				NAME = str(ct)
				PTS = np.load('PTList/'+FILE+'/'+NAME+'.npy')

				xs = np.linspace(min(PTS[:,0]),max(PTS[:,0]),num=17)
				xs[0] = xs[0] - 0.1
				ys = np.linspace(min(PTS[:,1]),max(PTS[:,1]),num=17)
				ys[0] = ys[0] - 0.1
				x1 = xs[1]
				x2 = xs[2]
				x3 = xs[3]
				x4 = xs[4]
				x5 = xs[5]
				x6 = xs[6]
				x7 = xs[7]
				x8 = xs[8]
				x9 = xs[9]
				x10 = xs[10]
				x11 = xs[11]
				x12 = xs[12]
				x13 = xs[13]
				x14 = xs[14]
				x15 = xs[15]
				x16 = xs[16]
				y1 = ys[1]
				y2 = ys[2]
				y3 = ys[3]
				y4 = ys[4]
				y5 = ys[5]
				y6 = ys[6]
				y7 = ys[7]
				y8 = ys[8]
				y9 = ys[9]
				y10 = ys[10]
				y11 = ys[11]
				y12 = ys[12]
				y13 = ys[13]
				y14 = ys[14]
				y15 = ys[15]
				y16 = ys[16]
			


				A1 = []
				B1 = []
				C1 = []
				D1 = []
				E1 = []
				F1 = []
				G1 = []
				H1 = []
				I1 = []
				J1 = []
				K1 = []
				L1 = []
				M1 = []
				N1 = []
				O1 = []
				P1 = []

				A2 = []
				B2 = []
				C2 = []
				D2 = []
				E2 = []
				F2 = []
				G2 = []
				H2 = []
				I2 = []
				J2 = []
				K2 = []
				L2 = []
				M2 = []
				N2 = []
				O2 = []
				P2 = []

				A3 = []
				B3 = []
				C3 = []
				D3 = []
				E3 = []
				F3 = []
				G3 = []
				H3 = []
				I3 = []
				J3 = []
				K3 = []
				L3 = []
				M3 = []
				N3 = []
				O3 = []
				P3 = []

				A4 = []
				B4 = []
				C4 = []
				D4 = []
				E4 = []
				F4 = []
				G4 = []
				H4 = []
				I4 = []
				J4 = []
				K4 = []
				L4 = []
				M4 = []
				N4 = []
				O4 = []
				P4 = []

				A5 = []
				B5 = []
				C5 = []
				D5 = []
				E5 = []
				F5 = []
				G5 = []
				H5 = []
				I5 = []
				J5 = []
				K5 = []
				L5 = []
				M5 = []
				N5 = []
				O5 = []
				P5 = []

				A6 = []
				B6 = []
				C6 = []
				D6 = []
				E6 = []
				F6 = []
				G6 = []
				H6 = []
				I6 = []
				J6 = []
				K6 = []
				L6 = []
				M6 = []
				N6 = []
				O6 = []
				P6 = []
				
				A7 = []
				B7 = []
				C7 = []
				D7 = []
				E7 = []
				F7 = []
				G7 = []
				H7 = []
				I7 = []
				J7 = []
				K7 = []
				L7 = []
				M7 = []
				N7 = []
				O7 = []
				P7 = []

				A8 = []
				B8 = []
				C8 = []
				D8 = []
				E8 = []
				F8 = []
				G8 = []
				H8 = []
				I8 = []
				J8 = []
				K8 = []
				L8 = []
				M8 = []
				N8 = []
				O8 = []
				P8 = []

				A9 = []
				B9 = []
				C9 = []
				D9 = []
				E9 = []
				F9 = []
				G9 = []
				H9 = []
				I9 = []
				J9 = []
				K9 = []
				L9 = []
				M9 = []
				N9 = []
				O9 = []
				P9 = []

				A10 = []
				B10 = []
				C10 = []
				D10 = []
				E10 = []
				F10 = []
				G10 = []
				H10 = []
				I10 = []
				J10 = []
				K10 = []
				L10 = []
				M10 = []
				N10 = []
				O10 = []
				P10 = []

				A11 = []
				B11 = []
				C11 = []
				D11 = []
				E11 = []
				F11 = []
				G11 = []
				H11 = []
				I11 = []
				J11 = []
				K11 = []
				L11 = []
				M11 = []
				N11 = []
				O11 = []
				P11 = []

				A12 = []
				B12 = []
				C12 = []
				D12 = []
				E12 = []
				F12 = []
				G12 = []
				H12 = []
				I12 = []
				J12 = []
				K12 = []
				L12 = []
				M12 = []
				N12 = []
				O12 = []
				P12 = []

				A13 = []
				B13 = []
				C13 = []
				D13 = []
				E13 = []
				F13 = []
				G13 = []
				H13 = []
				I13 = []
				J13 = []
				K13 = []
				L13 = []
				M13 = []
				N13 = []
				O13 = []
				P13 = []

				A14 = []
				B14 = []
				C14 = []
				D14 = []
				E14 = []
				F14 = []
				G14 = []
				H14 = []
				I14 = []
				J14 = []
				K14 = []
				L14 = []
				M14 = []
				N14 = []
				O14 = []
				P14 = []

				A15 = []
				B15 = []
				C15 = []
				D15 = []
				E15 = []
				F15 = []
				G15 = []
				H15 = []
				I15 = []
				J15 = []
				K15 = []
				L15 = []
				M15 = []
				N15 = []
				O15 = []
				P15 = []

				A16 = []
				B16 = []
				C16 = []
				D16 = []
				E16 = []
				F16 = []
				G16 = []
				H16 = []
				I16 = []
				J16 = []
				K16 = []
				L16 = []
				M16 = []
				N16 = []
				O16 = []
				P16 = []

				for pt in range(len(PTS[:,0])):
					if(PTS[pt,0] <= x1):
						if(PTS[pt,1] <= y1):
							A1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N1.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O1.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P1.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x2):
						if(PTS[pt,1] <= y1):
							A2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N2.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O2.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P2.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x3):
						if(PTS[pt,1] <= y1):
							A3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N3.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O3.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P3.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x4):
						if(PTS[pt,1] <= y1):
							A4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N4.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O4.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P4.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x5):
						if(PTS[pt,1] <= y1):
							A5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N5.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O5.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P5.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x6):
						if(PTS[pt,1] <= y1):
							A6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N6.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O6.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P6.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x7):
						if(PTS[pt,1] <= y1):
							A7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N7.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O7.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P7.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x8):
						if(PTS[pt,1] <= y1):
							A8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N8.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O8.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P8.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x9):
						if(PTS[pt,1] <= y1):
							A9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N9.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O9.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P9.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x10):
						if(PTS[pt,1] <= y1):
							A10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N10.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O10.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P10.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x11):
						if(PTS[pt,1] <= y1):
							A11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N11.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O11.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P11.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x12):
						if(PTS[pt,1] <= y1):
							A12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N12.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O12.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P12.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x13):
						if(PTS[pt,1] <= y1):
							A13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N13.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O13.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P13.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x14):
						if(PTS[pt,1] <= y1):
							A14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N14.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O14.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P14.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x15):
						if(PTS[pt,1] <= y1):
							A15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N15.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O15.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P15.append(D[PTS[pt,1],PTS[pt,0]])
					elif(PTS[pt,0] <= x16):
						if(PTS[pt,1] <= y1):
							A16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y2):
							B16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y3):
							C16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y4):
							D16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y5):
							E16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y6):
							F16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y7):
							G16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y8):
							H16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y9):
							I16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y10):
							J16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y11):
							K16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y12):
							L16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y13):
							M16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y14):
							N16.append(D[PTS[pt,1],PTS[pt,0]])
						elif(PTS[pt,1] <= y15):
							O16.append(D[PTS[pt,1],PTS[pt,0]])
						else:
							P16.append(D[PTS[pt,1],PTS[pt,0]])
				dat[ct,0] = Builder.C2I(Builder.AvgColors(A1)/255.0)
				dat[ct,1] = Builder.C2I(Builder.AvgColors(B1)/255.0)
				dat[ct,2] = Builder.C2I(Builder.AvgColors(C1)/255.0)
				dat[ct,3] = Builder.C2I(Builder.AvgColors(D1)/255.0)
				dat[ct,4] = Builder.C2I(Builder.AvgColors(E1)/255.0)
				dat[ct,5] = Builder.C2I(Builder.AvgColors(F1)/255.0)
				dat[ct,6] = Builder.C2I(Builder.AvgColors(G1)/255.0)
				dat[ct,7] = Builder.C2I(Builder.AvgColors(H1)/255.0)
				dat[ct,8] = Builder.C2I(Builder.AvgColors(I1)/255.0)
				dat[ct,9] = Builder.C2I(Builder.AvgColors(J1)/255.0)
				dat[ct,10] = Builder.C2I(Builder.AvgColors(K1)/255.0)
				dat[ct,11] = Builder.C2I(Builder.AvgColors(L1)/255.0)
				dat[ct,12] = Builder.C2I(Builder.AvgColors(M1)/255.0)
				dat[ct,13] = Builder.C2I(Builder.AvgColors(N1)/255.0)
				dat[ct,14] = Builder.C2I(Builder.AvgColors(O1)/255.0)
				dat[ct,15] = Builder.C2I(Builder.AvgColors(P1)/255.0)


				dat[ct,16] = Builder.C2I(Builder.AvgColors(A2)/255.0)
				dat[ct,17] = Builder.C2I(Builder.AvgColors(B2)/255.0)
				dat[ct,18] = Builder.C2I(Builder.AvgColors(C2)/255.0)
				dat[ct,19] = Builder.C2I(Builder.AvgColors(D2)/255.0)
				dat[ct,20] = Builder.C2I(Builder.AvgColors(E2)/255.0)
				dat[ct,21] = Builder.C2I(Builder.AvgColors(F2)/255.0)
				dat[ct,22] = Builder.C2I(Builder.AvgColors(G2)/255.0)
				dat[ct,23] = Builder.C2I(Builder.AvgColors(H2)/255.0)
				dat[ct,24] = Builder.C2I(Builder.AvgColors(I2)/255.0)
				dat[ct,25] = Builder.C2I(Builder.AvgColors(J2)/255.0)
				dat[ct,26] = Builder.C2I(Builder.AvgColors(K2)/255.0)
				dat[ct,27] = Builder.C2I(Builder.AvgColors(L2)/255.0)
				dat[ct,28] = Builder.C2I(Builder.AvgColors(M2)/255.0)
				dat[ct,29] = Builder.C2I(Builder.AvgColors(N2)/255.0)
				dat[ct,30] = Builder.C2I(Builder.AvgColors(O2)/255.0)
				dat[ct,31] = Builder.C2I(Builder.AvgColors(P2)/255.0)	

				dat[ct,32] = Builder.C2I(Builder.AvgColors(A3)/255.0)
				dat[ct,33] = Builder.C2I(Builder.AvgColors(B3)/255.0)
				dat[ct,34] = Builder.C2I(Builder.AvgColors(C3)/255.0)
				dat[ct,35] = Builder.C2I(Builder.AvgColors(D3)/255.0)
				dat[ct,36] = Builder.C2I(Builder.AvgColors(E3)/255.0)
				dat[ct,37] = Builder.C2I(Builder.AvgColors(F3)/255.0)
				dat[ct,38] = Builder.C2I(Builder.AvgColors(G3)/255.0)
				dat[ct,39] = Builder.C2I(Builder.AvgColors(H3)/255.0)	
				dat[ct,40] = Builder.C2I(Builder.AvgColors(I3)/255.0)
				dat[ct,41] = Builder.C2I(Builder.AvgColors(J3)/255.0)
				dat[ct,42] = Builder.C2I(Builder.AvgColors(K3)/255.0)
				dat[ct,43] = Builder.C2I(Builder.AvgColors(L3)/255.0)
				dat[ct,44] = Builder.C2I(Builder.AvgColors(M3)/255.0)
				dat[ct,45] = Builder.C2I(Builder.AvgColors(N3)/255.0)
				dat[ct,46] = Builder.C2I(Builder.AvgColors(O3)/255.0)
				dat[ct,47] = Builder.C2I(Builder.AvgColors(P3)/255.0)	

				dat[ct,48] = Builder.C2I(Builder.AvgColors(A4)/255.0)
				dat[ct,49] = Builder.C2I(Builder.AvgColors(B4)/255.0)
				dat[ct,50] = Builder.C2I(Builder.AvgColors(C4)/255.0)
				dat[ct,51] = Builder.C2I(Builder.AvgColors(D4)/255.0)
				dat[ct,52] = Builder.C2I(Builder.AvgColors(E4)/255.0)
				dat[ct,53] = Builder.C2I(Builder.AvgColors(F4)/255.0)
				dat[ct,54] = Builder.C2I(Builder.AvgColors(G4)/255.0)
				dat[ct,55] = Builder.C2I(Builder.AvgColors(H4)/255.0)
				dat[ct,56] = Builder.C2I(Builder.AvgColors(I4)/255.0)
				dat[ct,57] = Builder.C2I(Builder.AvgColors(J4)/255.0)
				dat[ct,58] = Builder.C2I(Builder.AvgColors(K4)/255.0)
				dat[ct,59] = Builder.C2I(Builder.AvgColors(L4)/255.0)
				dat[ct,60] = Builder.C2I(Builder.AvgColors(M4)/255.0)
				dat[ct,61] = Builder.C2I(Builder.AvgColors(N4)/255.0)
				dat[ct,62] = Builder.C2I(Builder.AvgColors(O4)/255.0)
				dat[ct,63] = Builder.C2I(Builder.AvgColors(P4)/255.0)

				dat[ct,64] = Builder.C2I(Builder.AvgColors(A5)/255.0)
				dat[ct,65] = Builder.C2I(Builder.AvgColors(B5)/255.0)
				dat[ct,66] = Builder.C2I(Builder.AvgColors(C5)/255.0)
				dat[ct,67] = Builder.C2I(Builder.AvgColors(D5)/255.0)
				dat[ct,68] = Builder.C2I(Builder.AvgColors(E5)/255.0)
				dat[ct,69] = Builder.C2I(Builder.AvgColors(F5)/255.0)
				dat[ct,70] = Builder.C2I(Builder.AvgColors(G5)/255.0)
				dat[ct,71] = Builder.C2I(Builder.AvgColors(H5)/255.0)
				dat[ct,72] = Builder.C2I(Builder.AvgColors(I5)/255.0)
				dat[ct,73] = Builder.C2I(Builder.AvgColors(J5)/255.0)
				dat[ct,74] = Builder.C2I(Builder.AvgColors(K5)/255.0)
				dat[ct,75] = Builder.C2I(Builder.AvgColors(L5)/255.0)
				dat[ct,76] = Builder.C2I(Builder.AvgColors(M5)/255.0)
				dat[ct,77] = Builder.C2I(Builder.AvgColors(N5)/255.0)
				dat[ct,78] = Builder.C2I(Builder.AvgColors(O5)/255.0)
				dat[ct,79] = Builder.C2I(Builder.AvgColors(P5)/255.0)

				dat[ct,80] = Builder.C2I(Builder.AvgColors(A6)/255.0)
				dat[ct,81] = Builder.C2I(Builder.AvgColors(B6)/255.0)
				dat[ct,82] = Builder.C2I(Builder.AvgColors(C6)/255.0)
				dat[ct,83] = Builder.C2I(Builder.AvgColors(D6)/255.0)
				dat[ct,84] = Builder.C2I(Builder.AvgColors(E6)/255.0)
				dat[ct,85] = Builder.C2I(Builder.AvgColors(F6)/255.0)
				dat[ct,86] = Builder.C2I(Builder.AvgColors(G6)/255.0)
				dat[ct,87] = Builder.C2I(Builder.AvgColors(H6)/255.0)
				dat[ct,88] = Builder.C2I(Builder.AvgColors(I6)/255.0)
				dat[ct,89] = Builder.C2I(Builder.AvgColors(J6)/255.0)
				dat[ct,90] = Builder.C2I(Builder.AvgColors(K6)/255.0)
				dat[ct,91] = Builder.C2I(Builder.AvgColors(L6)/255.0)
				dat[ct,92] = Builder.C2I(Builder.AvgColors(M6)/255.0)
				dat[ct,93] = Builder.C2I(Builder.AvgColors(N6)/255.0)
				dat[ct,94] = Builder.C2I(Builder.AvgColors(O6)/255.0)
				dat[ct,95] = Builder.C2I(Builder.AvgColors(P6)/255.0)

				dat[ct,96] = Builder.C2I(Builder.AvgColors(A7)/255.0)
				dat[ct,97] = Builder.C2I(Builder.AvgColors(B7)/255.0)
				dat[ct,98] = Builder.C2I(Builder.AvgColors(C7)/255.0)
				dat[ct,99] = Builder.C2I(Builder.AvgColors(D7)/255.0)
				dat[ct,100] = Builder.C2I(Builder.AvgColors(E7)/255.0)
				dat[ct,101] = Builder.C2I(Builder.AvgColors(F7)/255.0)
				dat[ct,102] = Builder.C2I(Builder.AvgColors(G7)/255.0)
				dat[ct,103] = Builder.C2I(Builder.AvgColors(H7)/255.0)
				dat[ct,104] = Builder.C2I(Builder.AvgColors(I7)/255.0)
				dat[ct,105] = Builder.C2I(Builder.AvgColors(J7)/255.0)
				dat[ct,106] = Builder.C2I(Builder.AvgColors(K7)/255.0)
				dat[ct,107] = Builder.C2I(Builder.AvgColors(L7)/255.0)
				dat[ct,108] = Builder.C2I(Builder.AvgColors(M7)/255.0)
				dat[ct,109] = Builder.C2I(Builder.AvgColors(N7)/255.0)
				dat[ct,110] = Builder.C2I(Builder.AvgColors(O7)/255.0)
				dat[ct,111] = Builder.C2I(Builder.AvgColors(P7)/255.0)

				dat[ct,112] = Builder.C2I(Builder.AvgColors(A8)/255.0)
				dat[ct,113] = Builder.C2I(Builder.AvgColors(B8)/255.0)
				dat[ct,114] = Builder.C2I(Builder.AvgColors(C8)/255.0)
				dat[ct,115] = Builder.C2I(Builder.AvgColors(D8)/255.0)
				dat[ct,116] = Builder.C2I(Builder.AvgColors(E8)/255.0)
				dat[ct,117] = Builder.C2I(Builder.AvgColors(F8)/255.0)
				dat[ct,118] = Builder.C2I(Builder.AvgColors(G8)/255.0)
				dat[ct,119] = Builder.C2I(Builder.AvgColors(H8)/255.0)
				dat[ct,120] = Builder.C2I(Builder.AvgColors(I8)/255.0)
				dat[ct,121] = Builder.C2I(Builder.AvgColors(J8)/255.0)
				dat[ct,122] = Builder.C2I(Builder.AvgColors(K8)/255.0)
				dat[ct,123] = Builder.C2I(Builder.AvgColors(L8)/255.0)
				dat[ct,124] = Builder.C2I(Builder.AvgColors(M8)/255.0)
				dat[ct,125] = Builder.C2I(Builder.AvgColors(N8)/255.0)
				dat[ct,126] = Builder.C2I(Builder.AvgColors(O8)/255.0)
				dat[ct,127] = Builder.C2I(Builder.AvgColors(P8)/255.0)

				dat[ct,128] = Builder.C2I(Builder.AvgColors(A9)/255.0)
				dat[ct,129] = Builder.C2I(Builder.AvgColors(B9)/255.0)
				dat[ct,130] = Builder.C2I(Builder.AvgColors(C9)/255.0)
				dat[ct,131] = Builder.C2I(Builder.AvgColors(D9)/255.0)
				dat[ct,132] = Builder.C2I(Builder.AvgColors(E9)/255.0)
				dat[ct,133] = Builder.C2I(Builder.AvgColors(F9)/255.0)
				dat[ct,134] = Builder.C2I(Builder.AvgColors(G9)/255.0)
				dat[ct,135] = Builder.C2I(Builder.AvgColors(H9)/255.0)
				dat[ct,136] = Builder.C2I(Builder.AvgColors(I9)/255.0)
				dat[ct,137] = Builder.C2I(Builder.AvgColors(J9)/255.0)
				dat[ct,138] = Builder.C2I(Builder.AvgColors(K9)/255.0)
				dat[ct,139] = Builder.C2I(Builder.AvgColors(L9)/255.0)
				dat[ct,140] = Builder.C2I(Builder.AvgColors(M9)/255.0)
				dat[ct,141] = Builder.C2I(Builder.AvgColors(N9)/255.0)
				dat[ct,142] = Builder.C2I(Builder.AvgColors(O9)/255.0)
				dat[ct,143] = Builder.C2I(Builder.AvgColors(P9)/255.0)

				dat[ct,144] = Builder.C2I(Builder.AvgColors(A10)/255.0)
				dat[ct,145] = Builder.C2I(Builder.AvgColors(B10)/255.0)
				dat[ct,146] = Builder.C2I(Builder.AvgColors(C10)/255.0)
				dat[ct,147] = Builder.C2I(Builder.AvgColors(D10)/255.0)
				dat[ct,148] = Builder.C2I(Builder.AvgColors(E10)/255.0)
				dat[ct,149] = Builder.C2I(Builder.AvgColors(F10)/255.0)
				dat[ct,150] = Builder.C2I(Builder.AvgColors(G10)/255.0)
				dat[ct,151] = Builder.C2I(Builder.AvgColors(H10)/255.0)
				dat[ct,152] = Builder.C2I(Builder.AvgColors(I10)/255.0)
				dat[ct,153] = Builder.C2I(Builder.AvgColors(J10)/255.0)
				dat[ct,154] = Builder.C2I(Builder.AvgColors(K10)/255.0)
				dat[ct,155] = Builder.C2I(Builder.AvgColors(L10)/255.0)
				dat[ct,156] = Builder.C2I(Builder.AvgColors(M10)/255.0)
				dat[ct,157] = Builder.C2I(Builder.AvgColors(N10)/255.0)
				dat[ct,158] = Builder.C2I(Builder.AvgColors(O10)/255.0)
				dat[ct,159] = Builder.C2I(Builder.AvgColors(P10)/255.0)

				dat[ct,160] = Builder.C2I(Builder.AvgColors(A11)/255.0)
				dat[ct,161] = Builder.C2I(Builder.AvgColors(B11)/255.0)
				dat[ct,162] = Builder.C2I(Builder.AvgColors(C11)/255.0)
				dat[ct,163] = Builder.C2I(Builder.AvgColors(D11)/255.0)
				dat[ct,164] = Builder.C2I(Builder.AvgColors(E11)/255.0)
				dat[ct,165] = Builder.C2I(Builder.AvgColors(F11)/255.0)
				dat[ct,166] = Builder.C2I(Builder.AvgColors(G11)/255.0)
				dat[ct,167] = Builder.C2I(Builder.AvgColors(H11)/255.0)
				dat[ct,168] = Builder.C2I(Builder.AvgColors(I11)/255.0)
				dat[ct,169] = Builder.C2I(Builder.AvgColors(J11)/255.0)
				dat[ct,170] = Builder.C2I(Builder.AvgColors(K11)/255.0)
				dat[ct,171] = Builder.C2I(Builder.AvgColors(L11)/255.0)
				dat[ct,172] = Builder.C2I(Builder.AvgColors(M11)/255.0)
				dat[ct,173] = Builder.C2I(Builder.AvgColors(N11)/255.0)
				dat[ct,174] = Builder.C2I(Builder.AvgColors(O11)/255.0)
				dat[ct,175] = Builder.C2I(Builder.AvgColors(P11)/255.0)


				dat[ct,176] = Builder.C2I(Builder.AvgColors(A12)/255.0)
				dat[ct,177] = Builder.C2I(Builder.AvgColors(B12)/255.0)
				dat[ct,178] = Builder.C2I(Builder.AvgColors(C12)/255.0)
				dat[ct,179] = Builder.C2I(Builder.AvgColors(D12)/255.0)
				dat[ct,180] = Builder.C2I(Builder.AvgColors(E12)/255.0)
				dat[ct,181] = Builder.C2I(Builder.AvgColors(F12)/255.0)
				dat[ct,182] = Builder.C2I(Builder.AvgColors(G12)/255.0)
				dat[ct,183] = Builder.C2I(Builder.AvgColors(H12)/255.0)
				dat[ct,184] = Builder.C2I(Builder.AvgColors(I12)/255.0)
				dat[ct,185] = Builder.C2I(Builder.AvgColors(J12)/255.0)
				dat[ct,186] = Builder.C2I(Builder.AvgColors(K12)/255.0)
				dat[ct,187] = Builder.C2I(Builder.AvgColors(L12)/255.0)
				dat[ct,188] = Builder.C2I(Builder.AvgColors(M12)/255.0)
				dat[ct,189] = Builder.C2I(Builder.AvgColors(N12)/255.0)
				dat[ct,190] = Builder.C2I(Builder.AvgColors(O12)/255.0)
				dat[ct,191] = Builder.C2I(Builder.AvgColors(P12)/255.0)

				dat[ct,192] = Builder.C2I(Builder.AvgColors(A13)/255.0)
				dat[ct,193] = Builder.C2I(Builder.AvgColors(B13)/255.0)
				dat[ct,194] = Builder.C2I(Builder.AvgColors(C13)/255.0)
				dat[ct,195] = Builder.C2I(Builder.AvgColors(D13)/255.0)
				dat[ct,196] = Builder.C2I(Builder.AvgColors(E13)/255.0)
				dat[ct,197] = Builder.C2I(Builder.AvgColors(F13)/255.0)
				dat[ct,198] = Builder.C2I(Builder.AvgColors(G13)/255.0)
				dat[ct,199] = Builder.C2I(Builder.AvgColors(H13)/255.0)
				dat[ct,200] = Builder.C2I(Builder.AvgColors(I13)/255.0)
				dat[ct,201] = Builder.C2I(Builder.AvgColors(J13)/255.0)
				dat[ct,202] = Builder.C2I(Builder.AvgColors(K13)/255.0)
				dat[ct,203] = Builder.C2I(Builder.AvgColors(L13)/255.0)
				dat[ct,204] = Builder.C2I(Builder.AvgColors(M13)/255.0)
				dat[ct,205] = Builder.C2I(Builder.AvgColors(N13)/255.0)
				dat[ct,206] = Builder.C2I(Builder.AvgColors(O13)/255.0)
				dat[ct,207] = Builder.C2I(Builder.AvgColors(P13)/255.0)

				dat[ct,208] = Builder.C2I(Builder.AvgColors(A14)/255.0)
				dat[ct,209] = Builder.C2I(Builder.AvgColors(B14)/255.0)
				dat[ct,210] = Builder.C2I(Builder.AvgColors(C14)/255.0)
				dat[ct,211] = Builder.C2I(Builder.AvgColors(D14)/255.0)
				dat[ct,212] = Builder.C2I(Builder.AvgColors(E14)/255.0)
				dat[ct,213] = Builder.C2I(Builder.AvgColors(F14)/255.0)
				dat[ct,214] = Builder.C2I(Builder.AvgColors(G14)/255.0)
				dat[ct,215] = Builder.C2I(Builder.AvgColors(H14)/255.0)
				dat[ct,216] = Builder.C2I(Builder.AvgColors(I14)/255.0)
				dat[ct,217] = Builder.C2I(Builder.AvgColors(J14)/255.0)
				dat[ct,218] = Builder.C2I(Builder.AvgColors(K14)/255.0)
				dat[ct,219] = Builder.C2I(Builder.AvgColors(L14)/255.0)
				dat[ct,220] = Builder.C2I(Builder.AvgColors(M14)/255.0)
				dat[ct,221] = Builder.C2I(Builder.AvgColors(N14)/255.0)
				dat[ct,222] = Builder.C2I(Builder.AvgColors(O14)/255.0)
				dat[ct,223] = Builder.C2I(Builder.AvgColors(P14)/255.0)

				dat[ct,224] = Builder.C2I(Builder.AvgColors(A15)/255.0)
				dat[ct,225] = Builder.C2I(Builder.AvgColors(B15)/255.0)
				dat[ct,226] = Builder.C2I(Builder.AvgColors(C15)/255.0)
				dat[ct,227] = Builder.C2I(Builder.AvgColors(D15)/255.0)
				dat[ct,228] = Builder.C2I(Builder.AvgColors(E15)/255.0)
				dat[ct,229] = Builder.C2I(Builder.AvgColors(F15)/255.0)
				dat[ct,230] = Builder.C2I(Builder.AvgColors(G15)/255.0)
				dat[ct,231] = Builder.C2I(Builder.AvgColors(H15)/255.0)
				dat[ct,232] = Builder.C2I(Builder.AvgColors(I15)/255.0)
				dat[ct,233] = Builder.C2I(Builder.AvgColors(J15)/255.0)
				dat[ct,234] = Builder.C2I(Builder.AvgColors(K15)/255.0)
				dat[ct,235] = Builder.C2I(Builder.AvgColors(L15)/255.0)
				dat[ct,236] = Builder.C2I(Builder.AvgColors(M15)/255.0)
				dat[ct,237] = Builder.C2I(Builder.AvgColors(N15)/255.0)
				dat[ct,238] = Builder.C2I(Builder.AvgColors(O15)/255.0)
				dat[ct,239] = Builder.C2I(Builder.AvgColors(P15)/255.0)

				dat[ct,240] = Builder.C2I(Builder.AvgColors(A16)/255.0)
				dat[ct,241] = Builder.C2I(Builder.AvgColors(B16)/255.0)
				dat[ct,242] = Builder.C2I(Builder.AvgColors(C16)/255.0)
				dat[ct,243] = Builder.C2I(Builder.AvgColors(D16)/255.0)
				dat[ct,244] = Builder.C2I(Builder.AvgColors(E16)/255.0)
				dat[ct,245] = Builder.C2I(Builder.AvgColors(F16)/255.0)
				dat[ct,246] = Builder.C2I(Builder.AvgColors(G16)/255.0)
				dat[ct,247] = Builder.C2I(Builder.AvgColors(H16)/255.0)
				dat[ct,248] = Builder.C2I(Builder.AvgColors(I16)/255.0)
				dat[ct,249] = Builder.C2I(Builder.AvgColors(J16)/255.0)
				dat[ct,250] = Builder.C2I(Builder.AvgColors(K16)/255.0)
				dat[ct,251] = Builder.C2I(Builder.AvgColors(L16)/255.0)
				dat[ct,252] = Builder.C2I(Builder.AvgColors(M16)/255.0)
				dat[ct,253] = Builder.C2I(Builder.AvgColors(N16)/255.0)
				dat[ct,254] = Builder.C2I(Builder.AvgColors(O16)/255.0)
				dat[ct,255] = Builder.C2I(Builder.AvgColors(P16)/255.0)


		np.save('imgcirc/Textures'+str(res)+'.npy',arr=dat)
							

def NPYtoZ(res):
	NPY = np.load('imgcirc/Textures'+str(res)+'.npy')
	G = open('imgcirc/Z'+str(res)+'_8.txt','w')
	for i in range(len(NPY[:,0])):
		print(str(i))
		str0 = ''
		for j in range(len(NPY[0,:])):
			str0 += str(NPY[i,j])+' '
		G.write(str0+'\n')
	G.close()		


def Buildit(res):
	for i in range(-720,765,45):
		I = i/10
		if(i%10 == 0):
			I = int(i/10)
		else:
			I = i/10
		print(str(I))
		for j in range(0,3600,45):
			J = j/10
			if(j%10 == 0):
				J = int(j/10)
			else:
				J = j/10
			Builder.FillImageLL0(I,J,res)

def Bits2Col(bits):
	value = 0
	for i in range(len(bits)):
		value += bits[i]*2**i
	return Builder.INT2Color(value)

def Col2Bits(col):
	value = Builder.C2I(col)
	bits = []
	oldvalue = value
	for i in range(24):
		bits.append(oldvalue%2)
		oldvalue = int(oldvalue/2)
	return bits
def NC2PNG_CLD(y0, m0, d0,X,Y):
	if(X*Y < len(Builder.Grid)):
		print('Array Not Big Enough for Data!')
		return np.ones(shape=(Y,X,4),dtype=np.float16)
	CLDS0 = Dataset('../../../Data/MERRA/CLD/goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/MERRA2_400.tavg1_2d_rad_Nx.'+str(y0)+str(m0)+str(d0)+'.SUB.nc')
	CLDS = CLDS0.variables['CLDTOT'][:]
	latnum = 361
	lonnum = 576
	CLDOUT = np.ones(shape=(X,Y,4),dtype=np.float16)
	outF = 'TISDat/CLD/CLDS'+y0+m0+d0+'.png'
	ct = 0
#	Month = (0,31,28,31,30,31,30,31,31,30,31,30,31)
	CumMonth = np.array((0,31,59,90,120,151,181,212,243,273,304,334,365))
	if(int(y0)%4 == 0):
		CumMonth[2:] += 1
	toy = np.sum(CumMonth[0:(int(m0)-1)])+int(d0)

#	SolarNoons = []
#	for k in range(24):
#		SolarNoons.append(Builder.FindSolarNoon(toy,k))

	for pt in Builder.Grid:
#		print(str(ct))
		lat = pt.getLatDeg()
		lon = pt.getLonDeg()
#		xyz = pt.getXYZPos()
#		hrs = range(0,16392,683)
		bits = []
		for k in range(24):
#			h0 = 25
			value = 0
#			d0 = Builder.dSurf(xyz,SolarNoons[k].getXYZPos())
#			for i in range(len(hrs)-1):
#				if(d0*10000.0 >= hrs[i] and d0*10000.0 < hrs[i+1]):
#					h0 = i
			if(random.random() < CLDS[k,int(360*(lat+90)/180),int(575*((lon+180)%360)/360)]):
				bits.append(1)
			else:
				bits.append(0)
#		print(bits)
#		print(Bits2Col(bits))
		CLDOUT[int(ct/X),int(ct%X),:] = Bits2Col(bits)
#		print(str(int(ct/X))+' '+str(ct%X))
#
#			h0 = i
#			angdat = h0
#			clddat = (round(9*CLDS[k,int(360*(lat+90)/180),int(575*((lon+180)%360)/360)]))
#			CLDOUT[ct,k] = (clddat*25+angdat)
		ct += 1
#	np.save(outF,arr=CLDOUT)
	plt.imsave(fname='TISDat/CLD/CLDS'+y0+m0+d0+'.png',arr=CLDOUT)


def NC2PNG_ANG(m0, d0, X, Y):
	print('ANG Data for '+m0+'/'+d0)
#	CLDS0 = Dataset('../../../Data/MERRA/2D/CLD/MERRA2_100.tavg1_2d_rad_Nx.'+str(y0)+str(m0)+str(d0)+'.SUB.nc')
#	CLDS = CLDS0.variables['CLDTOT'][:]
#	latnum = 361
#	lonnum = 576
	CLDOUT = np.ones(shape=(X,Y,4),dtype=np.float16)
#	CLDOUT = np.zeros(shape=(len(Builder.Grid),24),dtype=np.uint8)
	outF = 'TISDat/ANG/ANGS'+m0+d0+'.png'
	ct = 0
#	Month = (0,31,28,31,30,31,30,31,31,30,31,30,31)
	CumMonth = np.array((0,31,59,90,120,151,181,212,243,273,304,334,365))

	toy = CumMonth[int(m0)-1]+int(d0)
#	print(str(toy))
	SolarNoons = []
	for k in range(24):
		SolarNoons.append(Builder.FindSolarNoon(toy,k))
	for pt in Builder.Grid:
#		print(str(ct))
		lat = pt.getLatDeg()
		lon = pt.getLonDeg()
		xyz = pt.getXYZPos()
		hrs = range(0,16392,683)
		bits = []
		for k in range(24):
			h0 = 1
			d0 = Builder.dSurf(xyz,SolarNoons[k].getXYZPos())
			if(d0 > math.pi/2.0):
				h0 = 0
			bits.append(h0)
		CLDOUT[int(ct/X),int(ct%X),:] = Bits2Col(bits)
		ct += 1
	plt.imsave(fname=outF,arr=CLDOUT)

def NC2PNG_ANG_2(m0, d0, h0, lat,zoom):
	A = plt.imread('IMG/Zoom'+str(zoom)+'/Plane'+str(lat*10000)+'x0.png')

	CLDOUT = np.zeros(shape=(1080,1920,4),dtype=np.float16)
	outF = 'TISDat/ANG/Zoom'+str(zoom)+'/ANGS'+m0+d0+h0+'_'+str(lat*10000)+'.png'
	ct = 0
	CumMonth = np.array((0,31,59,90,120,151,181,212,243,273,304,334,365))

	toy = CumMonth[int(m0)-1]+int(d0)
	SolarNoon = Builder.FindSolarNoon(toy,int(h0))

	for i in range(1080):
#		print(str(i))
		for j in range(1920):
			ptid = Builder.C2I(A[i,j])
			if(ptid < 1474562):
				xyz = Builder.Grid[ptid].getXYZPos()
				d0 = Builder.dSurf(xyz,SolarNoon.getXYZPos())
				if(d0 > math.pi/2.0):
					CLDOUT[i,j,:] = np.array((0.2,0.2,0.2,1.0))
				else:
					CLDOUT[i,j,:] = np.array((1.0,1.0,1.0,1.0))
			else:
				CLDOUT[i,j,:] = np.array((0.0,0.0,0.0,1.0))
	plt.imsave(fname=outF,arr=CLDOUT)

IGBPTypes = ('Ocean','Coast','EvergreenNeedleleaf','EvergreenBroadleaf','DeciduousNeedleleaf','DeciduousBroadleaf','MixedForest','ClosedShrubland','OpenShrubland','WoodedSavanna','Savanna','Grasslands','PermanentWetlands','Croplands','Urban','SparseCrops','SnowIce','Desert','Peak')
print(IGBPTypes[0])
print(IGBPTypes[1])
print(IGBPTypes[16])
def LearnHeights():
	B = plt.imread('BigData/Topo5400.png')
	Z = Builder.CA2I(plt.imread('Z5400.png'))
	PTS = Builder.CA2I(plt.imread('ToXY5400.png'))

	cols = []
	for i in range(100):
		cols.append([])
	
	Z0 = np.zeros(shape=(100,4),dtype=np.float16)
	N0 = np.zeros(shape=100,dtype=int)
	for i in range(2700):
		print(str(i))
		for j in range(5400):
			typ = Builder.Grid[PTS[i,j]].getIGBP()
			if(typ != 0 and typ != 1 and typ != 16):
				cols[1+int(math.sqrt(Z[i,j]))].append(B[i,j])
				N0[1+int(math.sqrt(Z[i,j]))] += 1

	C = np.ones(shape=(100,1,4),dtype=np.float16)
#	C[:,0,0] = (0.001+Z0[:,0])/(0.001+N0[:])
#	C[:,0,1] = (0.001+Z0[:,1])/(0.001+N0[:])
#	C[:,0,2] = (0.001+Z0[:,2])/(0.001+N0[:])
	C[0,0,0] = 10/255.0
	C[0,0,1] = 10/255.0
	C[0,0,2] = 52/255.0
	for i in range(1,100,1):
		C[i,0,:] = Builder.commoncolors(cols[i])

	plt.imsave(fname='ZCols.png',arr=C)
				
	



#Buildit(4)
#HugeBuild(4)
#NPYtoZ(4)
#GenSmallCIRV()




#Prep4Huge()
HR = ('00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23')
DN = [31,28,31,30,31,30,31,31,30,31,30,31]
DYS = ('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31')


#def GenANGS(M,D,res):
#	for l in range(-90,105,15):
#		for h in HR:
#	print(D+'/'+M+': '+str(l)+' '+h+':00')
#NC2PNG_ANG('01','01', 1215, 1215)
#			NC2PNG_ANG_2(M,D,h,l,res)	

#BuildWGT(8640,4320)
#for i in range(4,DN[0],1):
#	print(DYS[i]+'/01')
#	GenANGS('01',DYS[i],0)

#NC2TIS_CLD('1980','01','02')

#def GenPlanes24(res):
#	for i in range(-90,105,15):
#		for j in range(0,360,15):
#			Builder.GenIMG2(i,j,res)
#			Builder.GenIMG4(i,j,res)
#			Builder.FillImageLL(i*10000,j*10000,res)


#GenANGS('01','01',0)
#for y in range(2011,2020,1):
#	for m in range(0,12,1):
#		for i in range(0,DN[m],1):
#			NC2PNG_CLD(str(y),DYS[m],DYS[i],1215,1215)


#GenPlanes24(0)

#EnhancePIMG()
#BuildPtListBig()
#PointReader(73)


#BuildPtList()	


#AssignCol3()
LearnHeights()
#Learn(2)