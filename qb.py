import numpy as np
import math
import matplotlib.pyplot as plt


Yc0 = 420
Yc1 = 1500
S0 = [21,10.5,5.25,2.625,1.3175,0.60875,0.304375]
#D0 = [25484,12742,6371,3186,1593]
#D0 = [384400,192200,96100,48050,24025,12012,6006]
D0 = [40000,26667,15238,8127,4195,2131,1065,533,266,133]
XF = 1920
YF = 1920
XOUT = 8500
YOUT = 8500
#re = 6371000.0
#re = 1
#re = np.sqrt(16/1474562)
re = 0.0655822
Ex = np.load('Ex8.npy')
Ey = np.load('Ey8.npy')
Ez = np.load('Ez8.npy')
NNs = np.load('NN8.npy')
Ns = np.load('CRDS.npy')
shadr = np.zeros(shape=(1920,1920,4),dtype=np.float16)
shadr[:,:,0:3] = 0.983
shadr[:,:,3] = 1.00
Xcrds = np.zeros(shape=(1920,1920),dtype=int)
Ycrds = np.zeros(shape=(1920,1920),dtype=int)
for i in range(1920):
	Xcrds[:,i] = i
	Ycrds[i,:] = i

y11 = np.zeros(shape=(3840,3840),dtype=int)
x11 = np.zeros(shape=(3840,3840),dtype=int)

for i in range(3840):
	y11[i,:] = i
	x11[:,i] = i



#	Color Distance Metric is Euclidean distance / 1+Correlation Coefficient (R^2)

def ColD(col1,col2):
	SG1 = np.sqrt(np.average((col1[0:3]-np.average(col1[0:3]))*(col1[0:3]-np.average(col1[0:3]))))
	CV12 = np.average((col1[0:3]-np.average(col1[0:3]))*(col2[0:3]-np.average(col2[0:3])))
	SG2 = np.sqrt(np.average((col2[0:3]-np.average(col2[0:3]))*(col2[0:3]-np.average(col2[0:3]))))
	ED12 = np.average((col1[0:3]-col2[0:3])**2)
	return ED12/(1+CV12**2/(SG1*SG2)**2)

def ColDA(col1,col2):
	colD = np.ones(shape=col1.shape,dtype=np.float16)
	SG1 = np.zeros(shape=(col1.shape[0],col1.shape[1]),dtype=float)
	SG2 = np.zeros(shape=(col1.shape[0],col1.shape[1]),dtype=float)
	CV = np.zeros(shape=(col1.shape[0],col1.shape[1]),dtype=float)
	C1A = ((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3)
	SG1[:,:]  = ( ( col1[:,:,0]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) )**2 + ( col1[:,:,1]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) )**2 + ( col1[:,:,2]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) )**2 )/3
	SG2[:,:]  = ( ( col2[:,:,0]-((col2[:,:,0]+col2[:,:,1]+col2[:,:,2])/3) )**2 + ( col2[:,:,1]-((col2[:,:,0]+col2[:,:,1]+col2[:,:,2])/3) )**2 + ( col2[:,:,2]-((col2[:,:,0]+col2[:,:,1]+col2[:,:,2])/3) )**2 )/3
	CV[:,:] = ( ( ( col1[:,:,0]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) ) - ( col2[:,:,0]-((col2[:,:,0]+col2[:,:,1]+col2[:,:,2])/3) ) )**2 + ( ( col1[:,:,1]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) ) - ( col2[:,:,1]-((col2[:,:,0]+col2[:,:,1]+col2[:,:,2])/3) ) )**2 + ( ( col1[:,:,0]-((col1[:,:,2]+col1[:,:,1]+col1[:,:,2])/3) ) - ( col2[:,:,2]-((col2[:,:,0]+col2[:,:,1]+col2[:,:,2])/3) ) )**2 )/3

	ED12 = np.sqrt( (col1[:,:,0]-col2[:,:,0])**2 + (col1[:,:,1]-col2[:,:,1])**2 + (col1[:,:,2]-col2[:,:,2])**2 )
	return ED12/(1+CV**2/(SG1*SG2)**2)

def ColDA0(col1,col2):
#	colD = np.ones(shape=col1.shape,dtype=np.float16)
#	SG1 = np.zeros(shape=(col1.shape[0],col1.shape[1]),dtype=float)
#	CV = np.zeros(shape=(col1.shape[0],col1.shape[1]),dtype=float)
#	C1A = ((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3)
#	SG1[:,:]  = ( ( col1[:,:,0]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) )**2 + ( col1[:,:,1]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) )**2 + ( col1[:,:,2]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) )**2 )/3
#	SG2  = ( ( col2[0]-((col2[0]+col2[1]+col2[2])/3) )**2 + ( col2[1]-((col2[0]+col2[1]+col2[2])/3) )**2 + ( col2[2]-((col2[0]+col2[1]+col2[2])/3) )**2 )/3
#	CV[:,:] = ( ( ( col1[:,:,0]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) ) - ( col2[0]-((col2[0]+col2[1]+col2[2])/3) ) )**2 + ( ( col1[:,:,1]-((col1[:,:,0]+col1[:,:,1]+col1[:,:,2])/3) ) - ( col2[1]-((col2[0]+col2[1]+col2[2])/3) ) )**2 + ( ( col1[:,:,0]-((col1[:,:,2]+col1[:,:,1]+col1[:,:,2])/3) ) - ( col2[2]-((col2[0]+col2[1]+col2[2])/3) ) )**2 )/3

	ED12 = np.sqrt( (col1[:,:,0]-col2[0])**2 + (col1[:,:,1]-col2[1])**2 + (col1[:,:,2]-col2[2])**2 )
	return ED12#/(1+CV**2/(0.0001+SG1*SG2)**2)


def Cluster2Colors(fname):
	A = list(open(fname))
	B = np.ones(shape=(len(A)-1,4),dtype=np.float16)
	ct = 0
	for i in range(len(A)):
		if( i != 0 ):
			B[i-1,0] = float(A[i].split(',')[1])/255.0
			B[i-1,1] = float(A[i].split(',')[2])/255.0
			B[i-1,2] = float(A[i].split(',')[3])/255.0
	return B



	
def ReplaceColors0(img,colors):
	imgnew = np.ones(shape=(img.shape[0],img.shape[1],4),dtype=np.float16)
	for i in range(img.shape[0]):
		print(str(i))
		for j in range(img.shape[1]):
			dmax = 10.0
			col0 = np.ones(shape=4,dtype=np.float16)
			for k in colors:
				dcol = ColD(img[i,j],k)
				if(dcol < dmax):
					dmax = dcol
					col0 = k
			imgnew[i,j,:] = col0
	return imgnew

def ReplaceColors1(img,colors):
	imgnew = np.zeros(shape=(img.shape[0],img.shape[1]),dtype=int)
	dmax = 2550.0*np.ones(shape=(img.shape[0],img.shape[1]),dtype=float)
	dcol = np.zeros(shape=(img.shape[0],img.shape[1]),dtype=float)
	ct = 0
	for k in colors:
#		print(str(ct))
		dcol[:,:] = ColDA0(img,k)
		cmask = SelectValuesLT(dcol[:,:],dmax[:,:]-0.0000001)
		dmax[:,:] = dcol[:,:]*cmask+dmax[:,:]*(1-cmask)
#		print(dmax)
		imgnew[:,:] = C2I(k)*cmask+imgnew[:,:]*(1-cmask)
		ct += 1
	return INTS2Color(imgnew)


def DDX(INarr):
	outarr = np.zeros(shape=INarr.shape,dtype=INarr.dtype)
	Xmax = INarr.shape[1]
	outarr[:,0:Xmax-1] = (INarr[:,1:Xmax]-INarr[:,0:Xmax-1])/re
	return outarr

def DDX4(INarr):
	outarr = np.zeros(shape=(INarr.shape[0],INarr.shape[1],4),dtype=INarr.dtype)
	Xmax = INarr.shape[1]
	outarr[:,0:Xmax-1,0] = INarr[:,0:Xmax-1]-INarr[:,1:Xmax]
	outarr[:,:,1] = outarr[:,:,0]
	outarr[:,:,2] = outarr[:,:,0]
	outarr[:,:,3] = outarr[:,:,0]
	return outarr

def DDY(INarr):
	outarr = np.zeros(shape=INarr.shape,dtype=INarr.dtype)
	Ymax = INarr.shape[0]
	outarr[0:Ymax-1,:] = (INarr[1:Ymax,:]-INarr[0:Ymax-1,:])/re
	return outarr

def DDY4(INarr):
	outarr = np.zeros(shape=(INarr.shape[0],INarr.shape[1],4),dtype=INarr.dtype)
	Ymax = INarr.shape[0]
	outarr[0:Ymax-1,:,0] = INarr[1:Ymax,:]-INarr[0:Ymax-1,:]
	outarr[:,:,1] = outarr[:,:,0]
	outarr[:,:,2] = outarr[:,:,0]
	outarr[:,:,3] = outarr[:,:,0]
	return outarr

def UFF(img):
	CMASK = SelectValueC(CA2I(img),0)
	img[0:img.shape[0]-2,:] = img[0:img.shape[0]-2,:]*(1-CMASK[0:img.shape[0]-2,:])+CMASK[0:img.shape[0]-2,:]*img[1:img.shape[0]-1,:]
#	img[:,0:img.shape[1]-2] = img[:,0:img.shape[1]-2]*(1-CMASK[:,0:img.shape[1]-2])+CMASK[:,0:img.shape[1]-2]*img[:,1:img.shape[1]-1]
	return img

def PositiveOrZerosINT(X):
	outarr = np.zeros(shape=X.shape,dtype=int)
	outarr[:,:] = (X[:,:]+np.abs(X[:,:]))/2
	return outarr


def Lift(img,hgt,res,angle):

	z0 = GenSphereDepth(res)
	if(angle != 0):
		outimg = np.zeros(shape=(img.shape[0],img.shape[1],4),dtype=np.float16)
		outimg[:,:,3] = 1.0
		cmask = SelectValue(img[:,:],0)
		H = np.zeros(shape=(img.shape[0],img.shape[1]),dtype=int)
		if(res >= 2):
			H[:,:] = 2**(res)*(hgt[:,:]/1000.0+z0[:,:])*(math.sin(math.pi*angle/180.0))/14
		else:
			H[:,:] = 2**(res)*(z0[:,:])*(math.sin(math.pi*angle/180.0))/14
		curvmask = cmask[1919-PositiveOrZerosINT(1919-(Ycrds[:,:]+H[:,:])) , Xcrds[:,:] ]
		outimg[:,:] = img[ 1919-PositiveOrZerosINT(1919-(Ycrds[:,:]+H[:,:])) , Xcrds[:,:] ]#*(1-cmask)
		return outimg
	else:	
		return img



def PositiveOrZeros(X):
	return (X+np.abs(X))/(2)

#	Stores arrays for Rotate

yn = np.zeros(shape=(YF,XF),dtype=int)
xn = np.zeros(shape=(YF,XF),dtype=int)
for i in range(XF):
	if(i < YF):
		yn[i,:] = i
	xn[:,i] = i

y0 = np.zeros(shape=(Yc1,XF),dtype=float)
x0 = np.zeros(shape=(Yc1,XF),dtype=float)
yn = np.zeros(shape=(Yc1,XF),dtype=float)
xn = np.zeros(shape=(Yc1,XF),dtype=float)
for i in range(XF):
	if(i < Yc1):
		yn[i,:] = Yc1-i
	xn[:,i] = i-int(XF/2)


y1 = np.zeros(shape=(Yc1,XF),dtype=int)
x1 = np.zeros(shape=(Yc1,XF),dtype=int)



#	Returns Brightness (Intensity) of Lambertion Reflectance off MANTIS Grid (NL == grid number where sun is overhead)

def LambertReflect(NL):
	ds = 22445.0		# Distance to the Sun
	B = np.zeros(shape=(1474562),dtype=np.float16)
	L = np.zeros(shape=(1474562,3),dtype=np.float16)
	L[:,0] = Ns[NL,0]-Ns[:,0]/ds
	L[:,1] = Ns[NL,1]-Ns[:,1]/ds
	L[:,2] = Ns[NL,2]-Ns[:,2]/ds
#	dl = np.sqrt(L[:,0]**2+L[:,1]**2+L[:,2]**2)
#	L[:,0] = L[:,0]/dl					# Normalze
#	L[:,1] = L[:,1]/dl
#	L[:,2] = L[:,2]/dl
	return Ns[:,0]*L[:,0]+Ns[:,1]*L[:,1]+Ns[:,2]*L[:,2]		# N dot L = B

def SpecularReflect(NL,m):
	cosi = LambertReflect(NL)
	Rs = 0.5*( ((cosi-np.sqrt(m**2+cosi**2-1))/(cosi+np.sqrt(m**2+cosi**2-1)))**2+((m**2*cosi-np.sqrt(m**2+cosi**2-1))/(m**2*cosi+np.sqrt(m**2+cosi**2-1)))**2)
	return Rs

def SpecularReflectPhong(NL,m):

	cosi = LambertReflect(NL)
	Rs = 0.5*( ((cosi-np.sqrt(m**2+cosi**2-1))/(cosi+np.sqrt(m**2+cosi**2-1)))**2+((m**2*cosi-np.sqrt(m**2+cosi**2-1))/(m**2*cosi+np.sqrt(m**2+cosi**2-1)))**2)
	return Rs

def SpecularReflectBP(NL,lat,lon,n):
	ds = 22445.0		# Distance to the Sun
	if(lat > 90 or lat < 90):
		if(lon > 360):
			lat = lat/10000.0
			lon = lon/10000.0
	lat = math.pi*lat/180.0
	lon = math.pi*lon/180.0
	xv = math.cos(lat)*math.sin(lon)
	yv = math.cos(lat)*math.cos(lon)
	zv = math.sin(lat)
	L = np.zeros(shape=(1474562,3),dtype=np.float16)
	L[:,0] = Ns[NL,0]-Ns[:,0]/ds
	L[:,1] = Ns[NL,1]-Ns[:,1]/ds
	L[:,2] = Ns[NL,2]-Ns[:,2]/ds
	H = np.zeros(shape=(1474562,3),dtype=np.float16)
	H[:,0] = 0.5*(L[:,0]+xv)
	H[:,1] = 0.5*(L[:,1]+yv)
	H[:,2] = 0.5*(L[:,2]+zv)
	return (Ns[:,0]*H[:,0]+Ns[:,1]*H[:,1]+Ns[:,2]*H[:,2])**n		# N dot H = B

def SelectValuesLT(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = ((value-field)+np.abs(value-field))/(value-field)/2
	return outf

def SelectValuesLT3(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = ((value-field)+np.abs(value-field))/(value-field)/2
	return outf

def SelectValuesGT(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = ((field-value)+np.abs(field-value))/(field-value)/2
	return outf

#	Turns Integer ID into 4-D 24-bit color tuple
def INT2Color(pid):
	return ((pid-pid%65536)/65536/255.0 , ((pid-pid%256)/256)%256/255.0, pid%256/255.0,255/255.0)	

def INTS2Color(pid):
	outarr = np.ones(shape=(pid.shape[0],pid.shape[1],4),dtype=np.float16)
	outarr[:,:,0] = (pid-pid%65536)/65536/255.0
	outarr[:,:,1] = ((pid-pid%256)/256)%256/255.0
	outarr[:,:,2] = pid%256/255.0
	return outarr	

def INTS2Color1(pid):
	outarr = np.ones(shape=(pid.shape[0],4),dtype=np.float16)
	outarr[:,0] = (pid-pid%65536)/65536/255.0
	outarr[:,1] = ((pid-pid%256)/256)%256/255.0
	outarr[:,2] = pid%256/255.0
	return outarr	
	
#	Inverts 24-bit color tuple or array into Integer ID
def C2I(color):
	R = int(color[0]*255.0+0.5)	
	G = int(color[1]*255.0+0.5)
	B = int(color[2]*255.0+0.5)
	return (R*256*256+G*256+B)

def CA2I(color):
	outarr = np.zeros(shape=(color.shape[0],color.shape[1]),dtype=int)
	outarr[:,:] = color[:,:,0]*255*256*256+color[:,:,1]*255*256+color[:,:,2]*255
	return outarr

def SelectValuesList(field,values):
	outf = np.zeros(shape=(shape.field[0],shape.field[1]),dtype=float)

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
	return outc*1.0

def SelectValueH3(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:,:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	return outf

def SelectValue1(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	return outf

def SelectValue2(field,value):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = np.exp(-1*np.abs(np.log((value+1)/(field+1)) ))
	return outf

def SelectValues(field,values):
	ourf = np.zeros(shape=(len(values),field.shape),dtype=int)



def SelectValue3(field,value1,value2):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = 1-(1-np.exp(-1*np.abs(np.log((value1+1)/(field+1)) )))*(1-np.exp(-1*np.abs(np.log((value2+1)/(field+1)) )))
	return outf

def DescColors252(color):
	R = color[0]
	G = color[1]
	B = color[2]
	R0 = 0.99
	G0 = 0.99
	B0 = 0.99
	#	Red Colors
	if(R < 0.1):
		R0 = 0.04167
	elif(R < 0.2):
		R0 = 0.15
	elif(R < 0.3):
		R0 = 0.25
	elif(R <  0.55):
		R0 = 0.35
	elif(R < 0.85):
		R0 = 0.65
	else:
		R0 = 0.99

	#	Green Colors
	if(G < 0.1):
		G0 = 0.04167
	elif(G < 0.3):
		G0 = 0.2
	elif(G <  0.375):
		G0 = 0.35
	elif(G < 0.5):
		G0 = 0.4
	elif(G < 0.7):
		G0 = 0.6
	elif(G < 0.85):
		G0 = 0.75
	else:
		G0 = 0.99

	if(B < 0.075):
		B0 = 0.04167
	elif(B < 0.15):
		B0 = 0.12
	elif(B < 0.25):
		B0 = 0.2
	elif(B <  0.45):
		B0 = 0.35
	elif(B < 0.85):
		B0 = 0.55
	else:
		B0 = 0.99

	return np.array((R0,G0,B0,1.0))


#*********************
#	Misc Tools	
#---------------------

def InjectiveResize(imgarr,y1,x1):
	outarr = np.zeros(shape=(y1,x1,4),dtype=np.float16)
	outarr[:,:,3] = 1.0
	for i in range(imgarr.shape[0]):
		for j in range(imgarr.shape[1]):
			outarr[int(i*(y1-1)/(imgarr.shape[0]-1)),int(j*(x1-1)/(imgarr.shape[1]-1)),0:3] = imgarr[i,j,0:3]
	return outarr

def AbsHumidity(tempF,relhumid):
	tempI = int(tempF)
	tempD = tempF-tempI
	Vals = list(open('../../../SaturatedAir32Fto80F.txt'))[0].split(',')
	if(tempF < 32):
		print('Temperature is Below 32')
	elif(tempF > 80):
		print('Temperature is Above 80')
	else:
		return (float(Vals[tempI-32])*relhumid/100)**(1-tempD)*(float(Vals[tempI-31])*relhumid/100)**tempD


def Hours2Dry(initialwater,T1,H1,T2,H2):
	inval = AbsHumidity(T1,H1)
	outval = AbsHumidity(T2,H2)
	return initialwater/((outval-intval)*25)

#**********************
#	Array Tools
#----------------------


def SelectBorders(field,rad):
	out1 = np.zeros(shape=(field.shape[0]+rad,field.shape[1]+rad),dtype=int)
	out2 = np.zeros(shape=(field.shape[0]+rad,field.shape[1]+rad),dtype=int)
	out3 = np.zeros(shape=(field.shape[0]+rad,field.shape[1]+rad),dtype=int)
	out4 = np.zeros(shape=(field.shape[0]+rad,field.shape[1]+rad),dtype=int)
	out1[rad:field.shape[0]+rad,0:field.shape[1]] = field
	out1[0:field.shape[0]+0,0:field.shape[1]] = field
	out2[0:field.shape[0]+0,0:field.shape[1]] = field
	out2[rad:field.shape[0]+rad,0:field.shape[1]] = field
	out3[0:field.shape[0]+0,rad:field.shape[1]+rad] = field
	out3[0:field.shape[0]+0,0:field.shape[1]+0] = field
	out4[0:field.shape[0]+0,0:field.shape[1]+0] = field
	out4[0:field.shape[0]+0,rad:field.shape[1]+rad] = field
	map2 = SelectValue2(out1[0:field.shape[0],0:field.shape[1]],out2[0:field.shape[0],0:field.shape[1]])
	map3 = SelectValue2(out3[0:field.shape[0],0:field.shape[1]],out4[0:field.shape[0],0:field.shape[1]])
	return (1-map2*map3)
	
	
	
	

def FastFillAVG(OLD):
	YF0 = OLD.shape[0]
	XF0 = OLD.shape[1]
	OLD = ChangeColorI(OLD,16777215,16777214)
	val = 0
	CMASK = SelectValue2( CA2I(OLD) , val)
	MASK = SelectValueC( CA2I(OLD) , val)
	NEW = np.zeros(shape=(YF0+2,XF0+2,4),dtype=np.float32)
	TOT = np.zeros(shape=(YF0+2,XF0+2),dtype=int)	
	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(YF0,XF0,4),dtype=np.float16)
	for i in range(3):
		for j in range(3):
			NEW[i:YF0+i,j:XF0+j,0:3] += OLD[:,:,0:3]
			TOT[i:YF0+i,j:XF0+j] += (1-CMASK[:,:])	
	NEW2[:,:,0] = (NEW[1:YF0+1,1:XF0+1,0]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
	NEW2[:,:,1] = (NEW[1:YF0+1,1:XF0+1,1]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
	NEW2[:,:,2] = (NEW[1:YF0+1,1:XF0+1,2]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
	OLD[:,:,0:3] += NEW2[:,:,0:3]*MASK[:,:,0:3]
	MASK = SelectValueC(CA2I(OLD),16777215)
	OLD[:,:,0:3] = OLD[:,:,0:3]*(1-MASK[:,:,0:3])
	return OLD

def FastFillHGT(OLD):
	YF0 = OLD.shape[0]
	XF0 = OLD.shape[1]
	val = 0
	CMASK = SelectValue2( OLD , val)
#	MASK = SelectValueC( OLD , val)
	NEW = np.zeros(shape=(YF0+2,XF0+2),dtype=float)
	TOT = np.zeros(shape=(YF0+2,XF0+2),dtype=int)	
#	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(YF0,XF0),dtype=int)
	for i in range(3):
		for j in range(3):
			NEW[i:YF0+i,j:XF0+j] += OLD[:,:]
			TOT[i:YF0+i,j:XF0+j] += (1-CMASK[:,:])	
	NEW2[:,:] = (NEW[1:YF0+1,1:XF0+1]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
#	NEW2[:,:,1] = (NEW[1:YF0+1,1:XF0+1,1]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
#	NEW2[:,:,2] = (NEW[1:YF0+1,1:XF0+1,2]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
	OLD[:,:] += NEW2[:,:]*CMASK[:,:]
	MASK = SelectValue2(OLD,1)
	OLD[:,:] = OLD[:,:]*(1-MASK[:,:])
	return OLD


def FastFillHGT3(OLD):
	YF0 = OLD.shape[1]
	XF0 = OLD.shape[2]
	val = 0
	CMASK = SelectValueH3( OLD , val)
#	MASK = SelectValueC( OLD , val)
	NEW = np.zeros(shape=(1474562,YF0+2,XF0+2),dtype=float)
	TOT = np.zeros(shape=(1474562,YF0+2,XF0+2),dtype=int)	
#	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(1474562,YF0,XF0),dtype=int)
	for i in range(3):
		for j in range(3):
			NEW[:,i:YF0+i,j:XF0+j] += OLD[:,:,:]
			TOT[:,i:YF0+i,j:XF0+j] += (1-CMASK[:,:,:])	
	NEW2[:,:,:] = (NEW[1:YF0+1,1:XF0+1]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
#	NEW2[:,:,1] = (NEW[1:YF0+1,1:XF0+1,1]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
#	NEW2[:,:,2] = (NEW[1:YF0+1,1:XF0+1,2]+0.001)/(TOT[1:YF0+1,1:XF0+1]+0.001)
	OLD[:,:,:] += NEW2[:,:,:]*CMASK[:,:,:]
	MASK = SelectValueH3(OLD,1)
	OLD[:,:,:] = OLD[:,:,:]*(1-MASK[:,:,:])
	return OLD

def FastFillBAY(OLD,i,j):
	YF0 = OLD.shape[0]
	XF0 = OLD.shape[1]
	val = 0
	CMASK = SelectValue2( CA2I(OLD) , val)
	NEW = np.zeros(shape=(YF0+2,XF0+2,4),dtype=np.float16)
	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(YF0,XF0,4),dtype=np.float16)
	NEW[i:YF0+i,j:XF0+j,0:3] = OLD[:,:,0:3]
	NEW2[:,:,0] = OLD[:,:,0]*(1-CMASK[:,:])+NEW[1:YF0+1,1:XF0+1,0]*CMASK[:,:]
	NEW2[:,:,1] = OLD[:,:,1]*(1-CMASK[:,:])+NEW[1:YF0+1,1:XF0+1,1]*CMASK[:,:]
	NEW2[:,:,2] = OLD[:,:,2]*(1-CMASK[:,:])+NEW[1:YF0+1,1:XF0+1,2]*CMASK[:,:]	
	return NEW2

def FFB(OLD):
	NEW0 = FastFillBAY(OLD,0,1)
	NEW0 = FastFillBAY(NEW0,2,1)
	NEW0 = FastFillBAY(NEW0,1,0)
	NEW0 = FastFillBAY(NEW0,1,2)
	return NEW0

def ChangeColor(OLD,colold,colnew):
	MASK = SelectValue2( CA2I(OLD) , C2I(colold))
	COL = CA2I(OLD)*(1-MASK)+C2I(colnew)*MASK
	return INTS2Color(COL)


def ChangeColorI(OLD,colold,colnew):
	MASK = SelectValue2( CA2I(OLD) , colold)
	COL = CA2I(OLD)*(1-MASK)+colnew*MASK
	return INTS2Color(COL)



def DDYM(Field):
#	RNs = np.zeros(shape=(1474562,7,3),dtype=Field.dtype)
	conserve = np.sum(Field)
	FNs = Field[NNs[:,:]]
	RNs = Ns[NNs[:,:],:]
	RLocs = np.zeros(shape=(1474562,7),dtype=float)
	Ravg = np.zeros(shape=(1474562),dtype=float)
	Favg = np.zeros(shape=(1474562),dtype=float)
	for i in range(7):
		RLocs[:,i] = re*(RNs[:,i,0]*Ey[:,0]+RNs[:,i,1]*Ey[:,1]+RNs[:,i,2]*Ey[:,2])
		Ravg[:] += RLocs[:,i]
		Favg[:] += FNs[:,i]
	Ravg = Ravg/7
	Favg = Favg/7
	Fsqr = np.zeros(shape=(1474562),dtype=float)
	Nsqr = np.zeros(shape=(1474562),dtype=float)
	for i in range(7):
		Fsqr[:] += (FNs[:,i]-Favg[:])*(RLocs[:,i]-Ravg[:])
		Nsqr[:] += (RLocs[:,i]-Ravg[:]+0.000001)**2
#	return np.sum((Field[:]-np.average(Field))*(Rloc[:]-np.average(Rloc)))/npsum((Rloc[:]-np.average(Rloc))**2)
	return Fsqr/Nsqr#*conserve/np.sum(Fsqr/Nsqr)	
		

def DDXM(Field):
#	RNs = np.zeros(shape=(1474562,7,3),dtype=Field.dtype)
	conserve = np.sum(Field)
	FNs = Field[NNs[:,:]]
	RNs = Ns[NNs[:,:],:]
	RLocs = np.zeros(shape=(1474562,7),dtype=float)
	Ravg = np.zeros(shape=(1474562),dtype=float)
	Favg = np.zeros(shape=(1474562),dtype=float)
	for i in range(7):
		RLocs[:,i] = re*(RNs[:,i,0]*Ex[:,0]+RNs[:,i,1]*Ex[:,1]+RNs[:,i,2]*Ex[:,2])
		Ravg[:] += RLocs[:,i]
		Favg[:] += FNs[:,i]
	Ravg = Ravg/7
	Favg = Favg/7
	Fsqr = np.zeros(shape=(1474562),dtype=float)
	Nsqr = np.zeros(shape=(1474562),dtype=float)
	for i in range(7):
		Fsqr[:] += (FNs[:,i]-Favg[:])*(RLocs[:,i]-Ravg[:])
		Nsqr[:] += (RLocs[:,i]-Ravg[:]+0.000001)**2
#	return np.sum((Field[:]-np.average(Field))*(Rloc[:]-np.average(Rloc)))/npsum((Rloc[:]-np.average(Rloc))**2)
	return (Fsqr/Nsqr)#*conserve/np.sum(Fsqr/Nsqr)
	


def WaveEvolve(Field,FieldT,c,dt):
	FieldT2 = c**2*(DDYM(DDYM(Field[0:1474562]))+DDXM(DDXM(Field[0:1474562])))
	Field = Field[0:1474562]+FieldT[0:1474562]*dt+0.5*FieldT2[0:1474562]
	FieldT20 = c**2*(DDYM(DDYM(Field[0:1474562]))+DDXM(DDXM(Field[0:1474562])))
	FieldT = FieldT[0:1474562]+0.5*(FieldT2[0:1474562]+FieldT20[0:1474562])*dt
	return Field, FieldT

def GenSphereDepth(zoom):
	HGTS = list(open('H_Obs.txt'))
	hgt0 = float(HGTS[zoom])
#	phi = lat*180.0/math.pi
#	lam = lat*180.0/math.pi
#	Zdepth = np.ones(shape=(1920,1920),dtype=float)
	re = 6371.0
	Xi = np.zeros(shape=(1920,1920),dtype=float)
	Yi = np.zeros(shape=(1920,1920),dtype=float)
	RL0 = 3830/hgt0
	for i in range(1920):
		Xi[:,i] = i-959.5
		Yi[i,:] = i-959.5
	if(zoom >= 2):
		return np.sqrt((re)**2-(Xi/RL0)**2-(Yi/RL0)**2)-np.max(np.sqrt((re)**2-(Xi/RL0)**2-(Yi/RL0)**2))
	else:
		ans0 = (re)**2-(Xi/RL0)**2-(Yi/RL0)**2
		return np.sqrt((ans0+np.abs(ans0))/2)-np.max(np.sqrt((ans0+np.abs(ans0))/2))
#		return np.sqrt((ans0+np.abs(ans0))/2)

#		return (re/RL0)-1*np.sqrt(SelectValuesLT((re/L0)**2-(Xi/RL0)**2-(Yi/RL0)**2,0)*(re/(RL0)**2-(Xi/RL0)**2-(Yi/RL0)**2,0)+np.zeros(shape=(1920,1920),dtype=float)*(1-qb.SelectValuesLT(re/(RL0)**2-(Xi/RL0)**2-(Yi/RL0)**2,0)))


	

def Rotate(img,angle,res):
	if(angle == 0):
		return img[Yc0:Yc1,:]
	if(res < 2):
		outimg = np.zeros(shape=(1920,1920,4),dtype=np.float16)
		outimg[:,:,3] = 1.0
		Ssize = [38,153,335]
		dpp = 197/7/1920.0
		if(res == 0 or res == 1):
			if(angle == 15):
				outimg[int(15/dpp):1920,:] = img[0:1920-int(15/dpp),:]
				return outimg
			else:
				return outimg
			

	direction = 'N'
	yf = img.shape[0]
	xf = img.shape[1]
	outimg = np.zeros(shape=(xf,yf,4),dtype=np.float16)
	outimg[:,:,3] = 1.0
	if(xf == yf):
	
		imgold = img[0:Yc1,:]
		img = imgold[:,:]		

#	YF = img.shape[0]
#	XF = img.shape[1]	


#	YOUT = int(1.2*YF)
#	XOUT = int(1.2*XF)
	
	newimg = np.zeros(shape=(YOUT,XOUT,4),dtype=np.float16)
	newimg[:,:,3] = 1.0

	if(direction == 'N'):
		y0[:,:] = yn[:,:]
		x0[:,:] = xn[:,:]
	if(direction == 'E'):
		y0[:,:] = xn[:,:]
		x0[:,:] = -1*yn[:,:]
	if(direction == 'S'):
		y0[:,:] = -1*yn[:,:]
		x0[:,:] = -1*xn[:,:]
	if(direction == 'W'):
		y0[:,:] = -1*xn[:,:]
		x0[:,:] = yn[:,:]
#	z0 = GenSphereDepth(res)
	y1[:,:] = y0[:,:]*math.cos(math.pi*angle/180.0)#+z0[0:Yc1,:]*math.sin(math.pi*angle/180.0)
	x1[:,:] = x0[:,:]*((y0[:,:]*math.sin(math.pi*angle/180.0)+D0[res])/D0[res])



	newimg[(YOUT-1-y1[:,:]),int(XOUT/2.0)+x1[:,:],0:3] = img[:,:,0:3]
	if(angle >= 60):
		outimg[Yc0-300:Yc1-300,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0-300:Yc1-300,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg
	else:
		outimg[Yc0:Yc1,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg

#	return newimg[Yc0:Yc1,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),:]


def RotateLarge(img,angle,res):
	Yc0 = int(img.shape[0]/2)-540
	Yc1 = int(img.shape[0]/2)+540
	Xc0 = int(img.shape[1]/2)-960
	Xc1 = int(img.shape[1]/2)+960
	if(angle == 0):
		return img[Yc0:Yc1,:]
	if(res < 2):
		outimg = np.zeros(shape=(3840,3840,4),dtype=np.float16)
		outimg[:,:,3] = 1.0
		Ssize = [38,153,335]
		dpp = 197/7/1920.0
		if(res == 0 or res == 1):
			if(angle == 15):
				outimg[int(15/dpp):3840,:] = img[0:3840-int(15/dpp),:]
				return outimg
			else:
				return outimg
			

	direction = 'N'
	yf = img.shape[0]
	xf = img.shape[1]
	outimg = np.zeros(shape=(xf,yf,4),dtype=np.float16)
	outimg[:,:,3] = 1.0
	if(xf == yf):
	
		imgold = img[0:Yc1,:]
		img = imgold[:,:]		


	
	newimg = np.zeros(shape=(YOUT,XOUT,4),dtype=np.float16)
	newimg[:,:,3] = 1.0


	y11 = np.zeros(shape=(Yc1,3840),dtype=int)
	x11 = np.zeros(shape=(Yc1,3840),dtype=int)
	for i in range(3840):
		if(i < Yc1):
			y11[i,:] = Yc1-i
		x11[:,i] = i-int(3840/2)






#	newimg[(YOUT-1-y1[:,:]),int(XOUT/2.0)+x1[:,:],0:3] = img[:,:,0:3]

	y1 = np.zeros(shape=(2460,3840),dtype=int)
	x1 = np.zeros(shape=(2460,3840),dtype=int)
#	z0 = GenSphereDepth(res)
	y1[:,:] = y11[:,:]*math.cos(math.pi*angle/180.0)#+z0[0:Yc1,:]*math.sin(math.pi*angle/180.0)
	x1[:,:] = x11[:,:]*((y11[:,:]*math.sin(math.pi*angle/180.0)+D0[res])/D0[res])


#	newimg[YOUT-1-y1[:,:],x1[:,:],0:3] = img[:,:,0:3]
	newimg[(YOUT-1-y1[:,:]),int(XOUT/2.0)+x1[:,:],0:3] = img[:,:,0:3]
	if(angle >= 90):
		outimg[Yc0-300:Yc1-300,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
		return outimg
	else:
#		outimg[Yc0:Yc1,Xc0:Xc1,0:3] = newimg[YOUT-2460:YOUT,int((XOUT-3840)/2.0):XOUT-int((XOUT-3840)/2.0),0:3]
#		outimg[:,:,0:3] = newimg[YOUT-2460:YOUT,0:3840,0:3]
#		return newimg[YOUT-2460:YOUT,0:3840,0:3]
		return newimg[YOUT-2460:YOUT,int((XOUT-3840)/2.0):XOUT-int((XOUT-3840)/2.0),0:3]




def RotateAny(img,angle,res):
	if(angle == 0):
		return img
	if(res < 2):
		outimg = np.zeros(shape=(img.shape[0],img.shape[1],4),dtype=np.float16)
		outimg[:,:,3] = 1.0
		Ssize = [38,153,335]
		dpp = 197/7/1920.0
		if(res == 0):
			if(angle == 15):
				outimg[int(15/dpp):img.shape[0],:] = img[0:img.shape[0]-int(15/dpp),:]
				return outimg
			else:
				return outimg

	direction = 'N'
	yf = img.shape[0]
	xf = img.shape[1]
	Yc0 = int(img.shape[0]/2)-540
	Yc1 = int(img.shape[0]/2)+540
	Xc0 = int(img.shape[1]/2)-960 
	Xc1 = int(img.shape[1]/2)+960 
	outimg = np.zeros(shape=(yf,xf,4),dtype=np.float16)
	outimg[:,:,3] = 1.0
#	if(xf == yf):
#	
#		imgold = img[:,:]
#		img = imgold[:,:]		
	
	newimg = np.zeros(shape=(YOUT,XOUT,4),dtype=np.float16)
	newimg[:,:,3] = 1.0

	if(direction == 'N'):
		y0[:,:] = yn[:,:]
		x0[:,:] = xn[:,:]
	if(direction == 'E'):
		y0[:,:] = xn[:,:]
		x0[:,:] = -1*yn[:,:]
	if(direction == 'S'):
		y0[:,:] = -1*yn[:,:]
		x0[:,:] = -1*xn[:,:]
	if(direction == 'W'):
		y0[:,:] = -1*xn[:,:]
		x0[:,:] = yn[:,:]

	y12 = np.zeros(shape=(img.shape[0],img.shape[1]),dtype=int)
	x12 = np.zeros(shape=(img.shape[0],img.shape[1]),dtype=int)

	y12[:,:] = y11[:,:]*math.cos(math.pi*angle/180.0)#+z0[0:Yc1,:]*math.sin(math.pi*angle/180.0)
	x12[:,:] = x11[:,:]*((y11[:,:]*math.sin(math.pi*angle/180.0)+D0[res])/D0[res])


	newimg[(YOUT-1-y12[:,:]),int(XOUT/4.0)+x12[:,:],0:3] = img[:,:,0:3]
	if(angle >= 60):
		outimg[Yc0-300:Yc1-300,:,0:3] = newimg[YOUT-img.shape[0]:YOUT,int((XOUT-img.shape[1])/2.0):XOUT-int((XOUT-img.shape[1])/2.0),0:3]
		return outimg
	else:
		outimg[(yf-1-y11[:,:]),x11[:,:],0:3] = newimg[YOUT-img.shape[0]:YOUT,int((XOUT-img.shape[1])/2.0):XOUT-int((XOUT-img.shape[1])/2.0),0:3]
		return outimg



def RotateLift(img,hgt,angle,res):
	if(angle == 0):
		return img[Yc0:Yc1,:]
	direction = 'N'
	yf = img.shape[0]
	xf = img.shape[1]
	outimg = np.zeros(shape=(xf,yf,4),dtype=np.float16)
	outimg[:,:,3] = 1.0
	if(xf == yf):
	
		imgold = img[0:Yc1,:]
		img = imgold[:,:]		

#	YF = img.shape[0]
#	XF = img.shape[1]	


#	YOUT = int(1.2*YF)
#	XOUT = int(1.2*XF)
	
	newimg = np.zeros(shape=(YOUT,XOUT,4),dtype=np.float16)
	newimg[:,:,3] = 1.0

	if(direction == 'N'):
		y0[:,:] = yn[:,:]
		x0[:,:] = xn[:,:]
	if(direction == 'E'):
		y0[:,:] = xn[:,:]
		x0[:,:] = -1*yn[:,:]
	if(direction == 'S'):
		y0[:,:] = -1*yn[:,:]
		x0[:,:] = -1*xn[:,:]
	if(direction == 'W'):
		y0[:,:] = -1*xn[:,:]
		x0[:,:] = yn[:,:]
#	z0 = hgt+GenSphereDepth(res)
	y1[:,:] = y0[:,:]*math.cos(math.pi*angle/180.0)#+z0[:,:]*math.sin(math.pi*angle/180.0)
	x1[:,:] = x0[:,:]*((y0[:,:]*math.sin(math.pi*angle/180.0)+D0[res])/D0[res])

	

	newimg[(YOUT-1-y1[:,:]),int(XOUT/2.0)+x1[:,:],0:3] = img[:,:,0:3]
	if(angle == 60):
		outimg[Yc0-300:Yc1-300,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0-300:Yc1-300,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg
	else:
		outimg[Yc0:Yc1,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg

#	return newimg[Yc0:Yc1,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),:]


def Rotate0(img,angle,res,E2):

	if(angle == 0):
		return img[Yc0:Yc1,:]
	direction = 'N'
	yf = img.shape[0]
	xf = img.shape[1]
	outimg = np.zeros(shape=(xf,yf,4),dtype=np.float16)
	outimg[:,:,3] = 1.0
	if(xf == yf):
	
		imgold = img[0:Yc1,:]
		img = imgold[:,:]		

#	YF = img.shape[0]
#	XF = img.shape[1]	


#	YOUT = int(1.2*YF)
#	XOUT = int(1.2*XF)
	
	newimg = np.zeros(shape=(YOUT,XOUT,4),dtype=np.float16)
	newimg[:,:,3] = 1.0

	if(direction == 'N'):
		y0[:,:] = yn[:,:]
		x0[:,:] = xn[:,:]
	if(direction == 'E'):
		y0[:,:] = xn[:,:]
		x0[:,:] = -1*yn[:,:]
	if(direction == 'S'):
		y0[:,:] = -1*yn[:,:]
		x0[:,:] = -1*xn[:,:]
	if(direction == 'W'):
		y0[:,:] = -1*xn[:,:]
		x0[:,:] = yn[:,:]

	y1[:,:] = y0[:,:]*math.cos(math.pi*angle/180.0)#+z0[:,:]*math.sin(math.pi*angle/180.0)
	x1[:,:] = x0[:,:]*((y0[:,:]*math.sin(math.pi*angle/180.0)+D0[res])/D0[res])
	newimg[(YOUT-1-y1[:,:]),int(XOUT/2.0)+x1[:,:],0:3] = img[:,:,0:3]
	if(angle >= 60):
		outimg[Yc0-300:Yc1-300,:,0:3] = FFB(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		
		return outimg
	else:
		outimg[Yc0:Yc1,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],1,0)[:,:,0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],1,2)[:,:,0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],0,1)[:,:,0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],2,1)[:,:,0:3]
		return outimg

#	return newimg[Yc0:Yc1,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),:], E2


def RotateAnticlockwise90(img):
	YF = img.shape[0]
	XF = img.shape[1]
	
	if(YF != XF):
		print('Base image is not a square!')
		XF = 1/0
		YF = 0/0
	yn = np.zeros(shape=(YF,XF),dtype=int)
	xn = np.zeros(shape=(YF,XF),dtype=int)
	for i in range(XF):
		if(i < YF):
			yn[i,:] = YF-i
		xn[:,i] = i#-int(XF/2)
	return img[YF-1-xn[:,:],-1*yn[:,:]]

def RotateClockwise90(img):
	YF = img.shape[0]
	XF = img.shape[1]
	
	if(YF != XF):
		print('Base image is not a square!')
		XF = 1/0
		YF = 0/0
	yn = np.zeros(shape=(YF,XF),dtype=int)
	xn = np.zeros(shape=(YF,XF),dtype=int)
	for i in range(XF):
		if(i < YF):
			yn[i,:] = YF-i-1
		xn[:,i] = i#-int(XF/2)
	return img[-xn[:,:],yn[:,:]]
	

#Y3 = np.zeros(shape=(3600,7200),dtype=int)
#X3 = np.zeros(shape=(3600,7200),dtype=int)

#for i in range(7200):
#	X3[:,i] = 3*i
#	if(i < 3600):
#		Y3[i,:] = 3*i



def AssetQuan(res):
	C64 = Cluster2Colors('BMCols_4320.csv')
	Dres = [15000,15000,15000,7500,3750,1875]
	for i in range(-90000,91875,Dres[res]):
		for j in range(0,360000,Dres[res]):
			print(str(i)+' '+str(j))
			A = plt.imread('Assets/Zoom'+str(res)+'/Plane'+str(i*10)+'x'+str(j*10)+'.png')
			B = ReplaceColors1(A,C64)
			plt.imsave(fname='Assets64/Zoom'+str(res)+'/Plane'+str(i*10)+'x'+str(j*10)+'.png',arr=B)


#A = plt.imread('BM_SW.png')

#B = plt.imread('BigData/BM7200_6.png')
#B0 = np.zeros(shape=(10800,21600,4),dtype=np.float16)
#for i in range(3):
#	for j in range(3):
#		B0[Y3[:,:]+i,X3[:,:]+j] = B[:,:]
#B = np.ones(shape=(10800,21600),dtype=int)

#C64 = Cluster2Colors('BMCols_4320.csv')
#for i in range(A.shape[0]):
#	A[i,:] = ReplaceColors1(A[i,:],C64)

#C = np.ones(shape=A.shape,dtype=A.dtype)\
#xjump = int(A.shape[1]/10)
#yjump = int(A.shape[0]/10)
#for i in range(10):
#	for j in range(10):
#		print(str((10*i+j))+'%')

#		cmask51 = SelectValue2(CA2I(B0[i*1080:i*1080+1080,j*2160:j*2160+2160]),C2I( (10/255.0,10/255.0,51/255.0)) )
#		cmask52 = SelectValue2(CA2I(B0[i*1080:i*1080+1080,j*2160:j*2160+2160]),C2I( (10/255.0,10/255.0,52/255.0)) )
#		B[i*1080:i*1080+1080,j*2160:j*2160+2160] = C2I((10/255.0,10/255.0,51/255.0))*cmask51+C2I((10/255.0,10/255.0,51/255.0))*cmask52+CA2I(A[i*yjump:i*yjump+yjump,j*2160:j*2160+2160])*(1-cmask51-cmask52)

#		A[i*yjump:i*yjump+yjump,j*xjump:j*xjump+xjump,0:3] = ReplaceColors1(A[i*yjump:i*yjump+yjump,j*xjump:j*xjump+xjump],C64)[:,:,0:3]


#plt.imsave(fname='BM_SW_SML.png',arr=A)


	