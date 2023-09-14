import numpy as np
import math


Yc0 = 420
Yc1 = 1500
S0 = [21,10.5,5.25,2.625,1.3175,0.60875]
#D0 = [25484,12742,6371,3186,1593]
D0 = [384400,192200,96100,48050,24025,12012]
XF = 1920
YF = 1920
XOUT = 7500
YOUT = 7500
Ns = np.load('CRDS.npy')


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

def SelectValue3(field,value1,value2):
	outf = np.zeros(shape=field.shape,dtype=int)
	outf[:,:] = 1-(1-np.exp(-1*np.abs(np.log((value1+1)/(field+1)) )))*(1-np.exp(-1*np.abs(np.log((value2+1)/(field+1)) )))
	return outf

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
	NEW = np.zeros(shape=(YF0+2,XF0+2,4),dtype=np.float16)
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

def FastFillBAY(OLD,i,j):
	YF0 = OLD.shape[0]
	XF0 = OLD.shape[1]
	val = 0
	CMASK = SelectValue2( CA2I(OLD) , val)
	NEW = np.zeros(shape=(YF0+2,XF0+2,4),dtype=np.float16)
	NEW[:,:,3] = 1.0
	NEW2 = np.ones(shape=(YF0,XF0,4),dtype=np.float16)
	NEW[i:YF0+i,j:XF0+j,0:3] += OLD[:,:,0:3]
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


def Rotate(img,angle,res):
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
	if(angle == 60):
		outimg[Yc0-300:Yc1-300,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0-300:Yc1-300,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg
	else:
		outimg[Yc0:Yc1,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg

#	return newimg[Yc0:Yc1,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),:]

def Rotate0(img,angle,res):

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
	if(angle == 60):
		outimg[Yc0-300:Yc1-300,:,0:3] = FastFillAVG(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3])[:,:,0:3]
		return outimg
	else:
		outimg[Yc0:Yc1,:,0:3] = newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],1,0)[:,:,0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],1,2)[:,:,0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],0,1)[:,:,0:3]
#		outimg[Yc0:Yc1,:,0:3] = FastFillBAY(newimg[YOUT-1080:YOUT,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),0:3],2,1)[:,:,0:3]
		return outimg

#	return newimg[Yc0:Yc1,int((XOUT-1920)/2.0):XOUT-int((XOUT-1920)/2.0),:]


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
	


	