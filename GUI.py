import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QWidget, QShortcut, QMainWindow, QApplication, QLabel, QFileDialog, QAction
from PyQt5.QtGui import QPixmap, QKeySequence, QImage, QPainter
from PyQt5.QtCore import Qt
from netCDF4 import Dataset
import Point9
import Builder
import Province
import math
import qb
import random
import Faction

PList = Province.getProvinceList()

YF = 1920
XF = 1920

Yc0 = 420#330
Yc1 = 1500#+330

IGBPTypes = ('Ocean','Coast','EvergreenNeedleleaf','EvergreenBroadleaf','DeciduousNeedleleaf','DeciduousBroadleaf','MixedForest','ClosedShrubland','OpenShrubland','WoodedSavanna','Savanna','Grasslands','PermanentWetlands','Croplands','Urban','SparseCrops','SnowIce','Desert','Peak')
RGNTypes = list(open('RGNZones.txt'))#  ('Flat','Hill','Mountain','Peak')
CLIMTypes = list(open('ClimateZones.txt'))# ('Grasslands','Plains','Desert','Tundra','Snow')
Provinces = list(open("Region0"))
PRVL = np.load('PRV.npy')
MonthStarts = np.array((0,31,59,90,120,151,181,212,243,273,304,334,365))
CLDSI = np.zeros(shape=(1215,1215),dtype=int)
CLDSD = np.zeros(shape=(24,1215,1215),dtype=np.float16)
ANGSI = np.zeros(shape=(1215,1215),dtype=int)
ANGSD = np.zeros(shape=(24,1215,1215),dtype=np.float16)
cr1 = np.zeros(shape=(YF,XF),dtype=int)
cr2 = np.zeros(shape=(YF,XF),dtype=int)
COL4 = np.zeros(shape=(YF,XF),dtype=np.float16)
COL5 = np.ones(shape=(YF,XF,4),dtype=np.float16)
HGTS = plt.imread('TISDat/Z_HGT.png')
Ys = np.zeros(shape=(YF,XF),dtype=int)
Xs = np.zeros(shape=(YF,XF),dtype=int)
Y2 = np.zeros(shape=(int(YF/2),int(XF/2)),dtype=int)
X2 = np.zeros(shape=(int(YF/2),int(XF/2)),dtype=int)
for i in range(YF):
	Ys[i,:] = i
	Xs[:,i] = i
	if(i < int(YF/2)):
		Y2[i,:] = 2*i
		X2[:,i] = 2*i

CHARS = []
SETTS = []
for F in Faction.getFactionList():
	for c in F.GetCharacters():
		CHARS.append(c)
	for s in F.GetSettlements():
		SETTS.append(s)

ALPHABET = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

Hr = ('00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23')
Mn = ('01','02','03','04','05','06','07','08','09','10','11','12')
Dn = ('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31')
Ml = (31,28,31,30,31,30,31,31,30,31,30,31)
Ly  = (1,0,0,0)

MONTH = np.zeros(shape=(365),dtype=np.uint8)
ct = 0
mct = 0
for i in range(365):
	MONTH[i] = mct
	ct += 1
	if(ct == Ml[mct]):
		ct = 0
		mct += 1

Grid = Point9.getPointList()
#g = open('Field.txt','w')
#for i in Grid:
#	z = i.getXYZPos()[2]
#	g.write(str(z)+'\n')

A = plt.imread('Assets/Zoom0/Plane0x0.png')
plt.imsave(fname='tmp.png',arr=A)
C = np.ones(shape=A.shape,dtype=np.float16)
C[:,:,0:3] = 0.0
D0 = qb.SelectValueC(qb.CA2I(A),0)
#E0 = np.zeros(shape=(1920,1920),dtype=np.uint8)
#E0[:,:] = (1-qb.SelectValue2(qb.CA2I(plt.imread('Assets/Zoom0/Plane0x0.png')),0))
#print(np.sum(E0))
A = plt.imread('Assets/Zoom1/Plane0x0.png')
D1 = qb.SelectValueC(qb.CA2I(A),0)
D = np.ones(shape=A.shape,dtype=np.uint8)
D[:,:,:] = D0[:,:,:]
E0 = (1-D0)
E0[:,:,3] = 1.0

E1 = (1-D1)
E1[:,:,3] = 1.0
E2 = np.zeros(shape=(1920,1920),dtype=np.uint8)


class MainWindow(QMainWindow):

	def __init__(self, parent = None):
		self.title = 'MANTIS'
		super(MainWindow, self).__init__(parent)
	
		self.setMouseTracking(True)
		self.VDisp = True
		self.FDisp = False
		self.RDisp = False	
		self.ODisp = False
		self.Place = 'Assets'
		self.GeodesicSetReady = False

		self.fldmap = np.loadtxt('Field.txt',dtype=np.float32)
		self.maxABS = max(abs(self.fldmap))

		self.lat = 0
		self.lon = 0
		self.latdeg = 0
		self.londeg = 0
		self.ycursorpos = 0
		self.lastypos = 0
		self.xcursorpos = 0
		self.lastxpos = 0
		self.zoom = 0
		self.zoom0 = 0
		self.scale = 10000
#		self.plane = list(open('Planes/TanXY'+str(self.lat)+'x'+str(self.lon)+'.txt'))
		self.zmlvl = [0,1,2,4,8,16,32,64,128,256]
#		self.jplvl = [15,15,15,7.5,3.75,1.875,0.9375]*self.scale
		self.jplvl = [150000,150000,150000,75000,37500,18750,18750,9375]
		self.selected = 88
		self.selectedProvince = 0
		self.lastselected = 88
		self.tobeedited = []
		self.isPole = 74
		self.lastPole = 74
		self.nearbyplane = 88


#		Time Information

		self.sol = False
		self.TODD = False
		self.CLDD = False
		self.GADD = False
		self.PRSD = False
		self.PTSD = False
		self.HUDD = False
		self.BLDD = False
		self.GRID = False
		self.TOPD = False

		self.year = 1980
		self.month = 0
		self.day = 0
		self.hour = 18
		self.img = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.hgt = plt.imread('HGT/Zoom2/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.img00 = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.angle = 0
		self.orientation = 'N'

		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))
#		self.cld = np.load('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.npy')

		self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
		CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255
		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')
#		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
		ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255
		for i in range(24):							#	Updates global CLDSD info which
			CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
			ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	


		
#		self.textures = plt.imread('imgcirc/New2/Textures'+str(self.zmlvl[self.zoom])+'.png')
		menubar = self.menuBar()
		fileMenu = menubar.addMenu('File')
		editMenu = menubar.addMenu('Edit')
		self.resize(500, 500)

		
#		editIGBP = menubar.addMenu(

		openAction = QAction('Open Image', self)  
		openAction.triggered.connect(self.openImage) 
		fileMenu.addAction(openAction)

		closeAction = QAction('Exit', self)  
		closeAction.triggered.connect(self.close) 
		fileMenu.addAction(closeAction)

#		Edit Scroll

		editCell = editMenu.addMenu('Edit Cell')
		editPRVM = editMenu.addMenu('Edit Province')
		editFLDV = editMenu.addMenu('Edit Field View')
		editDISP = editMenu.addMenu('Edit Display Field')

#		Edit Cell

		editPRV = editCell.addMenu('Edit Province')
		editIGBP = editCell.addMenu('Edit IGBP')
		editRGN = editCell.addMenu('Edit Terrain')
		editBRD = editCell.addMenu('Edit BorderInfo')
		editCLIM = editCell.addMenu('Edit Climate Zone')

		addPRV = editPRV.addMenu('Add Province')
#		addAfrica = addPRV.addMenu('Africa')
#		addAntarctica = addPRV.addMenu('Antarctica')
#		addAsia = addPRV.addMenu('Asia')
#		addAustralia = addPRV.addMenu('Australia')
#		addEurope = addPRV.addMenu('Europe')
		remPRV = editPRV.addMenu('Remove Province')
		
		self.field = np.load('PRV.npy')
		self.fieldt = np.zeros(shape=(1474562),dtype=self.field.dtype)

		Letter1 = []
		Letter2 = []
		for i in range(26):
			Letter1.append(addPRV.addMenu(ALPHABET[i]))
			Letter2.append([])
			for j in range(26):
				Letter2[i].append(Letter1[i].addMenu(ALPHABET[j]))



#		for i in range(26):
#			addPRV.addMenu(ALPHABET[i])
				

		Action1 = []
		Setter1 = []
		for i in range(len(Provinces)):
			Name = Provinces[i].split(' ')
#			PList[i].RefreshPRV()
			Action1.append(QAction('Find '+str(Name[0]),self))
			Setter1.append(0)
			for j in range(26):
				if(Name[0][0] == ALPHABET[j] or Name[0][0] == alphabet[j]):
					for k in range(26):
						if(Name[0][1] == ALPHABET[k] or Name[0][1] == alphabet[k]):
							Setter1[i] = Letter2[j][k].addMenu(str(Name[0]))
							
#							setPRV = QAction('Set '+str(Name[0]),self)
#							findPRV = QAction('Find '+str(Name[0]),self)
#							idxPRV = QAction(str(Name[0])+': '+str(i),self)
#							setPRV.triggered.connect(lambda: self.FindPRV(i,str(Name[0])))
							Action1[i].triggered.connect(lambda: self.FindPRV(i,str(Name[0])))
#							idxPRV.triggered.connect(lambda: self.FindPRV(i,str(Name[0])))
#							prvstat.addAction(setPRV)
							Setter1[i].addAction(Action1[i])
#							prvstat.addAction(idxPRV)

#		For Setting IGBP

		setIGBPOcean = QAction('Set to Ocean', self)
		setIGBPOcean.triggered.connect(lambda: self.SetIGBP(0))
		editIGBP.addAction(setIGBPOcean)
		setIGBPCoast = QAction('Set to Coast', self)
		setIGBPCoast.triggered.connect(lambda: self.SetIGBP(1))
		editIGBP.addAction(setIGBPCoast)
		setIGBPEvergreenNeedleleaf = QAction('Set to EvergreenNeedleleaf', self)
		setIGBPEvergreenNeedleleaf.triggered.connect(lambda: self.SetIGBP(2))
		editIGBP.addAction(setIGBPEvergreenNeedleleaf)
		setIGBPEvergreenBroadleaf = QAction('Set to EvergreenBroadleaf', self)
		setIGBPEvergreenBroadleaf.triggered.connect(lambda: self.SetIGBP(3))
		editIGBP.addAction(setIGBPEvergreenBroadleaf)
		setIGBPDeciduousNeedleleaf = QAction('Set to DeciduousNeedleleaf', self)
		setIGBPDeciduousNeedleleaf.triggered.connect(lambda: self.SetIGBP(4))
		editIGBP.addAction(setIGBPDeciduousNeedleleaf)
		setIGBPDeciduousBroadleaf = QAction('Set to DeciduousBroadleaf', self)
		setIGBPDeciduousBroadleaf.triggered.connect(lambda: self.SetIGBP(5))
		editIGBP.addAction(setIGBPDeciduousBroadleaf)
		setIGBPMixedForest = QAction('Set to MixedForest', self)
		setIGBPMixedForest.triggered.connect(lambda: self.SetIGBP(6))
		editIGBP.addAction(setIGBPMixedForest)
		setIGBPClosedShrubland = QAction('Set to ClosedShrubland', self)
		setIGBPClosedShrubland.triggered.connect(lambda: self.SetIGBP(7))
		editIGBP.addAction(setIGBPClosedShrubland)
		setIGBPOpenShrubland = QAction('Set to OpenShrubland', self)
		setIGBPOpenShrubland.triggered.connect(lambda: self.SetIGBP(8))
		editIGBP.addAction(setIGBPOpenShrubland)
		setIGBPWoodedSavanna = QAction('Set to WoodedSavanna', self)
		setIGBPWoodedSavanna.triggered.connect(lambda: self.SetIGBP(9))
		editIGBP.addAction(setIGBPWoodedSavanna)
		setIGBPSavanna = QAction('Set to Savanna', self)
		setIGBPSavanna.triggered.connect(lambda: self.SetIGBP(10))
		editIGBP.addAction(setIGBPSavanna)
		setIGBPGrasslands = QAction('Set to Grasslands', self)
		setIGBPGrasslands.triggered.connect(lambda: self.SetIGBP(11))
		editIGBP.addAction(setIGBPGrasslands)
		setIGBPPermanentWetlands = QAction('Set to PermanentWetlands', self)
		setIGBPPermanentWetlands.triggered.connect(lambda: self.SetIGBP(12))
		editIGBP.addAction(setIGBPPermanentWetlands)
		setIGBPCroplands = QAction('Set to Croplands', self)
		setIGBPCroplands.triggered.connect(lambda: self.SetIGBP(13))
		editIGBP.addAction(setIGBPCroplands)
		setIGBPUrban = QAction('Set to Urban', self)
		setIGBPUrban.triggered.connect(lambda: self.SetIGBP(14))
		editIGBP.addAction(setIGBPUrban)
		setIGBPSparseCrops = QAction('Set to SparseCrops', self)
		setIGBPSparseCrops.triggered.connect(lambda: self.SetIGBP(15))
		editIGBP.addAction(setIGBPSparseCrops)
		setIGBPSnowIce = QAction('Set to SnowIce', self)
		setIGBPSnowIce.triggered.connect(lambda: self.SetIGBP(16))
		editIGBP.addAction(setIGBPSnowIce)
		setIGBPDesert = QAction('Set to Desert', self)
		setIGBPDesert.triggered.connect(lambda: self.SetIGBP(17))
		editIGBP.addAction(setIGBPDesert)
		setIGBPPeak = QAction('Set to Peak', self)
		setIGBPPeak.triggered.connect(lambda: self.SetIGBP(18))
		editIGBP.addAction(setIGBPPeak)

#		For Setting Terrain

		setRGNFlat = QAction('Set to Flat', self)
		setRGNFlat.triggered.connect(lambda: self.SetRGN(0))
		editRGN.addAction(setRGNFlat)
		setRGNHill = QAction('Set to Hill', self)
		setRGNHill.triggered.connect(lambda: self.SetRGN(1))
		editRGN.addAction(setRGNHill)
		setRGNMountain = QAction('Set to Mountain', self)
		setRGNMountain.triggered.connect(lambda: self.SetRGN(2))
		editRGN.addAction(setRGNMountain)
		setRGNPeak = QAction('Set to Peak', self)
		setRGNPeak.triggered.connect(lambda: self.SetRGN(3))
		editRGN.addAction(setRGNPeak)

#		For Setting Climate

		setCLIMGrasslands = QAction('Set to Grasslands', self)
		setCLIMGrasslands.triggered.connect(lambda: self.SetCLIM(0))
		editCLIM.addAction(setCLIMGrasslands)
		setCLIMPlains = QAction('Set to Plains', self)
		setCLIMPlains.triggered.connect(lambda: self.SetCLIM(1))
		editCLIM.addAction(setCLIMPlains)
		setCLIMDesert = QAction('Set to Desert', self)
		setCLIMDesert.triggered.connect(lambda: self.SetCLIM(2))
		editCLIM.addAction(setCLIMDesert)
		setCLIMTundra = QAction('Set to Tundra', self)
		setCLIMTundra.triggered.connect(lambda: self.SetCLIM(3))
		editCLIM.addAction(setCLIMTundra)
		setCLIMSnow = QAction('Set to Snow', self)
		setCLIMSnow.triggered.connect(lambda: self.SetCLIM(4))
		editCLIM.addAction(setCLIMSnow)
		

		setBRD = QAction('Set to No Boarder',self)
		setBRD.triggered.connect(lambda: self.SetBRD(0))
		editBRD.addAction(setBRD)

		setBRD = QAction('Set to Boarder',self)
		setBRD.triggered.connect(lambda: self.SetBRD(1))
		editBRD.addAction(setBRD)


#		Edit Province

		editOwner = editPRVM.addMenu('Change Owner')
		editContinent = editPRVM.addMenu('Change Continent')
		editWatershed = editPRVM.addMenu('Edit Watershed')
		editAddSettlement = editPRVM.addMenu('Add Settlement')
		editRemSettlement = editPRVM.addMenu('Remove Settlement')
		editTile = editPRVM.addMenu('Edit Tiles')

#		Edit Field

		editField = editFLDV.addMenu('Edit Field')

		setFieldPRV = QAction('Provinces',self)
		setFieldPRV.triggered.connect(lambda: self.ChangeField(np.load('PRV.npy')))
		editField.addAction(setFieldPRV)
		setFieldHGT = QAction('Heights Z',self)
		setFieldHGT.triggered.connect(lambda: self.ChangeField(np.load('H0.npy')[:,0] ))
		editField.addAction(setFieldHGT)
		setFieldRGN = QAction('Roughness Z_sigma',self)
		setFieldRGN.triggered.connect(lambda: self.ChangeField(np.load('V0.npy')[:,0] ))
		editField.addAction(setFieldRGN)
		setFieldVIS = QAction('Blue Marble Visuals',self)
		setFieldVIS.triggered.connect(lambda: self.ChangeField(np.load('Z0.npy')))
		editField.addAction(setFieldVIS)
		setFieldRAW = QAction('Raw Visuals',self)
		setFieldRAW.triggered.connect(lambda: self.ChangeField(np.array(list(range(1474562)))))
		editField.addAction(setFieldRAW)

		setFieldZEROS = QAction('Zeros',self)
		setFieldZEROS.triggered.connect(lambda: self.ChangeField(np.zeros(shape=(1474562),dtype=int)))
		editField.addAction(setFieldZEROS)

#		Edit Display

		setDispEarthBM = QAction('Earth (Blue Marble)',self)
		setDispEarthBM.triggered.connect(lambda: self.OtherDisplay('Assets'))
		editDISP.addAction(setDispEarthBM)

		setDispMoon = QAction('Moon (Lunar)',self)
		setDispMoon.triggered.connect(lambda: self.OtherDisplay('Lunar'))
		editDISP.addAction(setDispMoon)

		setDispMars = QAction('Mars',self)
		setDispMars.triggered.connect(lambda: self.OtherDisplay('Mars'))
		editDISP.addAction(setDispMars)

		
		

		editColorBar = editFLDV.addMenu('Edit ColorBar')
		editRange = editFLDV.addMenu('Edit Range')
			
		
		
	
		
		




		self.shortcut = QShortcut(QKeySequence("W"), self)
		self.shortcut.activated.connect(self.moveNorth)

		self.shortcut = QShortcut(QKeySequence("A"), self)
		self.shortcut.activated.connect(self.moveWest)

		self.shortcut = QShortcut(QKeySequence("S"), self)
		self.shortcut.activated.connect(self.moveSouth)

		self.shortcut = QShortcut(QKeySequence("D"), self)
		self.shortcut.activated.connect(self.moveEast)

		self.shortcut = QShortcut(QKeySequence("C"),self)
		self.shortcut.activated.connect(self.CldDisplay)

		self.shortcut = QShortcut(QKeySequence("Z"),self)
		self.shortcut.activated.connect(self.RotateAnticlockwise)

		self.shortcut = QShortcut(QKeySequence("X"),self)
		self.shortcut.activated.connect(self.RotateClockwise)


		self.shortcut = QShortcut(QKeySequence("Q"), self)
		self.shortcut.activated.connect(self.zoomOut)

		self.shortcut = QShortcut(QKeySequence("E"), self)
		self.shortcut.activated.connect(self.zoomIn)

		self.shortcut = QShortcut(QKeySequence("F"), self)
		self.shortcut.activated.connect(self.FieldDisplay)

		self.shortcut = QShortcut(QKeySequence("K"), self)
		self.shortcut.activated.connect(self.DecreaseAngle)

		self.shortcut = QShortcut(QKeySequence("V"), self)
		self.shortcut.activated.connect(self.VisualDisplay)

		self.shortcut = QShortcut(QKeySequence("G"), self)
		self.shortcut.activated.connect(self.HighlightGrid)

		self.shortcut = QShortcut(QKeySequence("R"), self)
		self.shortcut.activated.connect(self.RawDisplay)

		self.shortcut = QShortcut(QKeySequence("I"), self)
		self.shortcut.activated.connect(self.IncreaseAngle)

		self.shortcut = QShortcut(QKeySequence("T"), self)
		self.shortcut.activated.connect(self.TimeOfDayDisplay)

		self.shortcut = QShortcut(QKeySequence("O"), self)
		self.shortcut.activated.connect(self.IncreaseTime)

		self.shortcut = QShortcut(QKeySequence("Shift+O"), self)
		self.shortcut.activated.connect(self.IncreaseTimeDay)

		self.shortcut = QShortcut(QKeySequence("Ctrl+O"), self)
		self.shortcut.activated.connect(self.IncreaseTimeMonth)

		self.shortcut = QShortcut(QKeySequence("Alt+O"), self)
		self.shortcut.activated.connect(self.IncreaseTimeYear)

		self.shortcut = QShortcut(QKeySequence("L"), self)
		self.shortcut.activated.connect(self.DecreaseTime)

		self.shortcut = QShortcut(QKeySequence("Shift+L"), self)
		self.shortcut.activated.connect(self.DecreaseTimeDay)

		self.shortcut = QShortcut(QKeySequence("Ctrl+L"), self)
		self.shortcut.activated.connect(self.DecreaseTimeMonth)

		self.shortcut = QShortcut(QKeySequence("Alt+L"), self)
		self.shortcut.activated.connect(self.DecreaseTimeYear)


		self.shortcut = QShortcut(QKeySequence("H"), self)
		self.shortcut.activated.connect(self.LowerHUDDisplay)

		self.shortcut = QShortcut(QKeySequence("U"), self)
		self.shortcut.activated.connect(self.Update)

		self.shortcut = QShortcut(QKeySequence("Ctrl+U"), self)
		self.shortcut.activated.connect(self.Update2)


		self.shortcut = QShortcut(QKeySequence("B"), self)
		self.shortcut.activated.connect(self.BuildingDisplay)

		self.shortcut = QShortcut(QKeySequence("P"), self)
		self.shortcut.activated.connect(self.GalacticDisplay)

		self.shortcut = QShortcut(QKeySequence("Y"), self)
		self.shortcut.activated.connect(self.TopoDisplay)

		self.shortcut = QShortcut(QKeySequence("F5"), self)
		self.shortcut.activated.connect(self.Refresh)

		self.shortcut = QShortcut(QKeySequence("Ctrl+F5"), self)
		self.shortcut.activated.connect(self.TotalRefresh)

		self.shortcut = QShortcut(QKeySequence("M"), self)
		self.shortcut.activated.connect(self.AddMiniMap)

		self.shortcut = QShortcut(QKeySequence("Ctrl+P"), self)
		self.shortcut.activated.connect(self.take_screenshot)

		self.shortcut = QShortcut(QKeySequence("Ctrl+S"), self)
		self.shortcut.activated.connect(self.SaveGrid)
		self.shortcut = QShortcut(QKeySequence("Ctrl+E"), self)
		self.shortcut.activated.connect(self.EngageEditor)

		self.shortcut = QShortcut(QKeySequence("Ctrl+X"), self)
		self.shortcut.activated.connect(self.EngageEditor)

		self.label = QLabel()
		self.setCentralWidget(self.label)
		
#		img = QImage(1920,1080,QImage.Format_RGB888)
#		pixmap = QPixmap(img)
		pixmap = QPixmap('BigData/Cover.png')
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()

	def SetPRV(self,prv):
#		print('Relocating ...')
		Builder.SetPRV(self.selected,prv)
		self.tobeedited.append(self.selectedProvince)

	def FindPRV(self,prv,name):
		self.lat = int(Builder.Grid[PList[prv].GetCapital()].getLatDeg()*10000/self.jplvl[self.zoom])*self.jplvl[self.zoom]
		self.lon = int(Builder.Grid[PList[prv].GetCapital()].getLonDeg()*10000/self.jplvl[self.zoom])*self.jplvl[self.zoom]
		print('For '+str(name))
		print('Nearest Plane to Capital '+str(Provinces[prv])+' : '+str(self.lat)+' x '+str(self.lon))
		self.StandardProcedure()
#		Builder.SetPRV(self.selected,prv)
#		self.tobeedited.append(self.selectedProvince)
		
		

	def SetIGBP(self,igbp):
		Builder.SetIGBP(self.selected,igbp)
		print('Setting IGBP of point '+str(self.selected)+' to type '+IGBPTypes[igbp])
		if(igbp == 17):
			Builder.SetCLIM(self.selected,2)
			print('Setting climate of point '+str(self.selected)+' to type Desert')
		self.tobeedited.append(self.selected)

	def SetRGN(self,rgn):
		Builder.SetRGN(self.selected,rgn)
		print('Setting terrain of point '+str(self.selected)+' to type '+RGNTypes[rgn])
		self.tobeedited.append(self.selected)

	def SetBoarder(self,brd):
		Builder.SetBoarder(self.selected,brd)
		self.tobeedited.append(self.selected)
	
	def SetCLIM(self,clim):
		Builder.SetCLIM(self.selected,clim)
		print('Setting climate of point '+str(self.selected)+' to type '+CLIMTypes[clim])
		self.tobeedited.append(self.selected)

#*************************************
#	Adds Captions to Image
#-------------------------------------

	def AddDTCaption(self):
		self.img[Yc0:Yc0+30,:,:] = Builder.INTS2Color(np.zeros(shape=(30,1920),dtype=int))
		self.img[Yc0:Yc0+30,0:30,0:3] = plt.imread('imgmisc/Compass/'+self.orientation+'.png')[:,:,0:3]
		self.img[Yc0+28:Yc0+30,:,:] = Builder.INTS2Color(np.ones(shape=(2,1920),dtype=int)*16777215)
		self.img[Yc0+7:Yc0+25,45:45+17*10,:] = Builder.GenText(Dn[self.day]+'/'+Mn[self.month]+'/'+str(self.year))	
		self.img[Yc0+7:Yc0+25,230:230+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')
		self.AddLowerHUD()			
#		if(self.HUDD == True):
#			STR1 = "Point ID:    "+str(self.selected)
#			STR2 = "Local Clock: "+Hr[int(self.hour+Builder.Grid[self.selected].getLonDeg()/15)%24]+":00"
#			STR3 = "IGBP:        "+IGBPTypes[Builder.Grid[self.selected].getIGBP()]
#			STR4 = "Climate:     "+CLIMTypes[Builder.Grid[self.selected].getCLIM()][0:-1]
#			STR5 = "Terrain:     "+RGNTypes[Builder.Grid[self.selected].getRGN()][0:-1]
#			STR6 = "Height:      "+str(Builder.C2I(HGTS[int(self.selected/1215),self.selected%1215]))+"m"


#			self.img[867:885,20:20+17*len(STR1),:] = Builder.GenText(STR1)	
#			self.img[892:910,20:20+17*len(STR2),:] = Builder.GenText(STR2)	
#			self.img[917:935,20:20+17*len(STR3),:] = Builder.GenText(STR3)	
#			self.img[942:960,20:20+17*len(STR4),:] = Builder.GenText(STR4)	
#			self.img[967:985,20:20+17*len(STR5),:] = Builder.GenText(STR5)	
#			self.img[992:1010,20:20+17*len(STR6),:] = Builder.GenText(STR6)	

		self.AddSettlementHUD()



	def AddDTCaption0(self,A):

		A[7:25,20:20+17*10,:] = Builder.GenText(Dn[self.day]+'/'+Mn[self.month]+'/'+str(self.year))	
		A[7:25,200:200+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')
		plt.imsave(fname='tmp.png',arr=A)


	def LowerHUDDisplay(self):
		self.img = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		if(self.HUDD == True):
			self.HUDD = False
			print('HUD Deactivated...')
			self.HighlightProvince()
			if(self.CLDD == True):
				self.AddCloudView()
			if(self.TODD == True):
				self.AddDayView()
			if(self.GADD == True):
				self.AddGalaxyView()
			self.AddDTCaption()
			plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
			pixmap = QPixmap('tmp.png')	
			self.label.setPixmap(pixmap)
			self.resize(pixmap.size())
			self.adjustSize()
		else:
			self.HUDD = True
			print('HUD Activated...')
			self.HighlightProvince()
			if(self.CLDD == True):
				self.AddCloudView()
			if(self.TODD == True):
				self.AddDayView()
			if(self.GADD == True):
				self.AddGalaxyView()
			self.AddDTCaption()
			plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
			pixmap = QPixmap('tmp.png')	
			self.label.setPixmap(pixmap)
			self.resize(pixmap.size())
			self.adjustSize()

	def BuildingDisplay(self):
		self.img = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		if(self.BLDD == True):
			self.BLDD = False
			print('Building Display Deactivated...')
			self.HighlightProvince()
			if(self.CLDD == True):
				self.AddCloudView()
			if(self.TODD == True):
				self.AddDayView()
			if(self.GADD == True):
				self.AddGalaxyView()
			self.AddDTCaption()
			plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
			pixmap = QPixmap('tmp.png')	
			self.label.setPixmap(pixmap)
			self.resize(pixmap.size())
			self.adjustSize()
		else:
			self.BLDD = True
			print('Building Display Activated...')
			self.HighlightProvince()
			if(self.CLDD == True):
				self.AddCloudView()
			if(self.TODD == True):
				self.AddDayView()
			if(self.GADD == True):
				self.AddGalaxyView()
			self.AddDTCaption()
			plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
			pixmap = QPixmap('tmp.png')	
			self.label.setPixmap(pixmap)
			self.resize(pixmap.size())
			self.adjustSize()

	def AddLowerHUD(self):
		col1 = np.array((0.5,0.25,0.25,1.0))
#		col2 = np.array((1.0,0.5,0.5,1.0))
		col2 = np.array((0.0,0.0,0.0,1.0))
		col3 = np.array((0.0,0.0,0.0,1.0))
		if(self.angle == 60):
			self.HUDD = True
		if(self.HUDD == True):
#			self.img[Yc0+780:Yc0+810,:,0] = (self.img[Yc0+780:Yc0+810,:,0]+col2[0])/2.0
#			self.img[Yc0+780:Yc0+810,:,1] = (self.img[Yc0+780:Yc0+810,:,1]+col2[1])/2.0
#			self.img[Yc0+780:Yc0+810,:,2] = (self.img[Yc0+780:Yc0+810,:,2]+col2[2])/2.0
			self.img[Yc0+780:Yc0+810,:,:] = Builder.INTS2Color(np.ones(shape=(30,1920),dtype=int)*Builder.C2I(col2))
			self.img[Yc0+780:Yc0+782,:,:] = Builder.INTS2Color(np.ones(shape=(2,1920),dtype=int)*16777215)
			self.img[Yc0+808:Yc0+810,:,:] = Builder.INTS2Color(np.ones(shape=(2,1920),dtype=int)*16777215)
			self.img[Yc0+810:Yc0+1080,:,:] = Builder.INTS2Color(np.ones(shape=(270,1920),dtype=int)*Builder.C2I(col2))
#			self.img[Yc0+810:Yc0+1080,:,0] = (self.img[Yc0+810:Yc0+1080,:,0]+col3[0])/2.0
#			self.img[Yc0+810:Yc0+1080,:,1] = (self.img[Yc0+810:Yc0+1080,:,1]+col3[1])/2.0
#			self.img[Yc0+810:Yc0+1080,:,2] = (self.img[Yc0+810:Yc0+1080,:,2]+col3[2])/2.0
#			NL1 = len(Provinces[self.selectedProvince].split(' ')[0])
#			self.img[Yc0+787:Yc0+805,20:20+17*NL1,:] = Builder.GenTextCOL(Provinces[self.selectedProvince].split(' ')[0],np.array((1.0,1.0,1.0,1.0)),col2)
			self.AddMiniMap()
			self.AddGuideView()
			self.img[Yc0+810:Yc0+1080,480:481,:] = Builder.INTS2Color(np.ones(shape=(270,1),dtype=int)*16777215)
			self.img[Yc0+810:Yc0+1080,509:510,:] = Builder.INTS2Color(np.ones(shape=(270,1),dtype=int)*16777215)


#	Guide View
	def AddGuideView(self):
		L1 = 'MODEL ATLAS (MANTIS)'		# Title String
		NL1 = len(L1)
		col2 = np.array((0.0,0.0,0.0,1.0))
		self.img[Yc0+787:Yc0+805,20:20+17*NL1,:] = Builder.GenTextCOL(L1,np.array((1.0,1.0,1.0,1.0)),col2)
		tmprv = Provinces[self.selectedProvince].split(' ')[0]
		if(tmprv == 'Ocean' and self.zoom < 2):
			tmprv = 'Earth'
		if(self.PRSD == True):
			TXT = np.array((1.0,1.0,1.0,1.0))
			BKG = np.array((0.0,0.0,0.0,1.0))
			tmprv = 'Region: '+tmprv
			LN = 'Income:     '+str(int(PList[self.selectedProvince].GetIncome()))
			self.img[Yc0+817:Yc0+835,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Order:      '+str(PList[self.selectedProvince].GetPublicOrder())+'%'
			self.img[Yc0+847:Yc0+865,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Population: '+str(PList[self.selectedProvince].GetPopulation())
			self.img[Yc0+877:Yc0+895,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Growth:     '+str(PList[self.selectedProvince].GetPopGrowth())+'%'
			self.img[Yc0+907:Yc0+925,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)


		if(self.PTSD == True):
			tmprv = 'Point '+str(self.selected)+' in region '+tmprv
			TXT = np.array((1.0,1.0,1.0,1.0))
			BKG = np.array((0.0,0.0,0.0,1.0))


#				print('Color:    NullSpace')
#				print('Point ID: NullSpace')
#				print('NN List:  NullSpace')
#				print('Province: NullSpace')
#				print('XYZ Pos:  NullSpace')
#				print('Lat:      NullSpace')
#				print('Lon:      NullSpace')
#				print('Border:   NullSpace')
#				print('IGBP:     NullSpace')
#				print('Terrain:  NullSpace')
#				print('Climate:  NullSpace')

			LN = 'NNList:      '+list(open('NN8.txt'))[self.selected][0:-1]
			self.img[Yc0+817:Yc0+835,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'XYZPos:      '+str(Builder.Grid[self.selected].getXYZPosSTR())
			self.img[Yc0+847:Yc0+865,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Latitude:    '+str(round(Builder.Grid[self.selected].getLatDeg(),6))
			self.img[Yc0+877:Yc0+895,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Longitude:   '+str(round(Builder.Grid[self.selected].getLonDeg(),6))
			self.img[Yc0+907:Yc0+925,530:530+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)


			LN = 'IGBP:    '+IGBPTypes[Builder.Grid[self.selected].getIGBP()]
			self.img[Yc0+847:Yc0+865,1120:1120+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Terrain: '+RGNTypes[Builder.Grid[self.selected].getRGN()][0:-1]
			self.img[Yc0+877:Yc0+895,1120:1120+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)

			LN = 'Climate: '+CLIMTypes[Builder.Grid[self.selected].getCLIM()][0:-1]
			if(len(LN) < 47):
				self.img[Yc0+907:Yc0+925,1120:1120+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)
			else:
				self.img[Yc0+907:Yc0+925,1120:1120+17*47,:] = Builder.GenTextCOL(LN[0:47],TXT,BKG)

			
		NL1 = len(tmprv)
		self.img[Yc0+787:Yc0+805,510+20:510+20+17*NL1,:] = Builder.GenTextCOL(tmprv,np.array((1.0,1.0,1.0,1.0)),col2)
		
				

#	Add Minimap

	def AddMiniMap(self):

		lathere = int(self.lat/self.jplvl[0])*self.jplvl[0]
		lonhere = int(self.lon/self.jplvl[0])*self.jplvl[0]
		colM = Builder.C2I(np.array((22/255.0,134/255.0,128/255.0,1.0)))

		rev1920 = list(range(1919,-1,-1))
		MASK0 = Builder.SelectValueC(Builder.CA2I(plt.imread('imgmisc/Shell02.png')),colM)
		if(self.GADD == True):
			lat0 = int(self.lat/self.jplvl[2])*self.jplvl[2]
			lon0 = int(self.lon/self.jplvl[2])*self.jplvl[2]
	
			A = plt.imread('Assets/STAR0/Plane'+str(lat0)+'x'+str((-1*lon0-10000*int(self.SolarOffset)-150000*(self.hour-12) )%3600000 )+'.png')[:,rev1920[:],:]
			self.img[Yc0+810:Yc0+1080,0:480,0:3] = A[960-135:960+135,960-215:960+265,0:3]

		colM = Builder.C2I(np.array((1.0,127/255.0,39/255.0,1.0)))
#		HUD0 = plt.imread('imgmisc/MTW2/Screenshots/MiniMap2.png')
#		MASK0 = Builder.SelectValueC(Builder.CA2I(HUD0),colM)

#		HUD0 = plt.imread('imgmisc/MTW2/Screenshots/MiniMap2.png')
#		MASK0 = Builder.SelectValueC(Builder.CA2I(HUD0),colM)

#		self.img[Yc0:Yc1,:,0:3] = self.img[Yc0:Yc1,:,0:3]*(MASK0[:,:,0:3])
		HUD1 = plt.imread('../CRDS6/Assets/Zoom0/Plane'+str(lathere)+'x'+str(lonhere)+'.png')
		self.img[Yc0+810:Yc0+1080,215-175:215+175,0:3] = self.img[Yc0+810:Yc0+1080,215-175:215+175,0:3]*MASK0[540-135:540+135,960-175:960+175,0:3]+HUD1[540-135:540+135,960-175:960+175,0:3]*(1-MASK0[540-135:540+135,960-175:960+175,0:3])

		X0 = np.zeros(shape=(6,4),dtype=int)
		Y0 = np.zeros(shape=(6,4),dtype=int)
		DX = np.zeros(shape=(6,4),dtype=int)
		DY = np.zeros(shape=(6,4),dtype=int)

		n0 = 0

#		if(self.orientation == 'E'):
#			n0 = 1
#		elif(self.orientation == 'S'):
#			n0 = 2
#		elif(self.orientation == 'W'):
#			n0 = 3
		

		X0[:,0] = [0,56,120,165,188,200]
		Y0[:,0] = [820,860,893,917,928,937]

		DX[:,0] = [480,328,190,108,58,32]
		DY[:,0] = [188,126,76,44,22,14]
	
		DX[:,1] = DY[:,0]
		DY[:,1] = DX[:,0]
	
		DX[:,2] = DX[:,0]
		DY[:,2] = DY[:,0]
	
		DX[:,3] = DY[:,0]
		DY[:,3] = DX[:,0]

		for i in range(len(X0)):
	
			X0[i] = X0[i]+2*int((self.lon-lonhere)/10000.0)
			Y0[i] = Y0[i]-2*int((self.lat-lathere)/10000.0)

		xc = np.zeros(shape=6,dtype=int)
		yc = np.zeros(shape=6,dtype=int)
		xc[:] = X0[:,0]+DX[:,0]/2
		yc[:] = Y0[:,0]+DY[:,0]/2
		for i in range(1,4,1):
			for j in range(0,6,1):
				X0[j,i] = max(xc[j]-DX[j,i]/2,0)
				Y0[j,i] = max(yc[j]-DY[j,i]/2,0)
			
 
		reds = np.zeros(shape=(1080,1920),dtype=int)
		reds[:,:] = 256*256*256-256*256
		REDS = Builder.INTS2Color(reds)
		self.img[Yc0+Y0[self.zoom,n0]:Yc0+Y0[self.zoom,n0]+DY[self.zoom,n0],X0[self.zoom,n0],0:3] = REDS[Y0[self.zoom,n0]:Y0[self.zoom,n0]+DY[self.zoom,n0],X0[self.zoom,n0],0:3]
		self.img[Yc0+Y0[self.zoom,n0]:Yc0+Y0[self.zoom,n0]+DY[self.zoom,n0],X0[self.zoom,n0]+DX[self.zoom,n0],0:3] = REDS[Y0[self.zoom,n0]:Y0[self.zoom,n0]+DY[self.zoom,n0],X0[self.zoom,n0]+DX[self.zoom,n0],0:3]
		self.img[Yc0+Y0[self.zoom,n0],X0[self.zoom,n0]:X0[self.zoom,n0]+DX[self.zoom,n0],0:3] = REDS[Y0[self.zoom,n0],X0[self.zoom,n0]:X0[self.zoom,n0]+DX[self.zoom,n0],0:3]
		self.img[Yc0+Y0[self.zoom,n0]+DY[self.zoom,n0],X0[self.zoom,n0]:X0[self.zoom,n0]+DX[self.zoom,n0],0:3] = REDS[Y0[self.zoom,n0]+DY[self.zoom,n0],X0[self.zoom,n0]:X0[self.zoom,n0]+DX[self.zoom,n0],0:3]



#***********************
#	MTW2 HUD
#-----------------------

	def AddLowerHUDD(self):

		col1 = np.array((0.5,0.25,0.25,1.0))
		col2 = np.array((1.0,0.5,0.5,1.0))
		col3 = np.array((0.0,0.0,0.0,1.0))
		
		if(self.HUDD == True):

			colM = Builder.C2I(np.array((1.0,127/255.0,39/255.0,1.0)))
			HUD0 = plt.imread('imgmisc/MTW2/Screenshots/HUD2.png')
			MASK = Builder.SelectValueC(Builder.CA2I(HUD0),colM)
			self.img[Yc0:Yc1,:,0:3] = self.img[Yc0:Yc1,:,0:3]*(MASK[:,:,0:3])+HUD0[:,:,0:3]*(1-MASK[:,:,0:3])

			lathere = int(self.lat/self.jplvl[0])*self.jplvl[0]
			lonhere = int(self.lon/self.jplvl[0])*self.jplvl[0]
			colM = Builder.C2I(np.array((1.0,127/255.0,39/255.0,1.0)))
			HUD0 = plt.imread('imgmisc/MTW2/Screenshots/MiniMap2.png')
			MASK0 = Builder.SelectValueC(Builder.CA2I(HUD0),colM)
			self.img[Yc0:Yc1,:,0:3] = self.img[Yc0:Yc1,:,0:3]*(MASK0[:,:,0:3])
			HUD1 = plt.imread('../CRDS6/Assets/Zoom0/Plane'+str(lathere)+'x'+str(lonhere)+'.png')
			self.img[Yc0+896:Yc0+1073,117:117+312,0:3] = HUD1[540-88:540+89,960-156:960+156,0:3]


			X0 = [63,110,175,222,244,258]
			Y0 = [896,896,930,952,968,975]
			for i in range(len(X0)):
				X0[i] = X0[i]+2*int((self.lon-lonhere)/10000.0)
				Y0[i] = Y0[i]-2*int((self.lat-lathere)/10000.0)
	
			DX = [408,328,190,107,61,32]
			DY = [176,176,88,50,27,15]		

			reds = np.zeros(shape=(1080,1920),dtype=int)
			reds[:,:] = 256*256*256-256*256
			REDS = Builder.INTS2Color(reds)
			self.img[Yc0+Y0[self.zoom]:Yc0+Y0[self.zoom]+DY[self.zoom],X0[self.zoom],0:3] = REDS[Y0[self.zoom]:Y0[self.zoom]+DY[self.zoom],X0[self.zoom],0:3]
			self.img[Yc0+Y0[self.zoom]:Yc0+Y0[self.zoom]+DY[self.zoom],X0[self.zoom]+DX[self.zoom],0:3] = REDS[Y0[self.zoom]:Y0[self.zoom]+DY[self.zoom],X0[self.zoom]+DX[self.zoom],0:3]
			self.img[Yc0+Y0[self.zoom],X0[self.zoom]:X0[self.zoom]+DX[self.zoom],0:3] = REDS[Y0[self.zoom],X0[self.zoom]:X0[self.zoom]+DX[self.zoom],0:3]
			self.img[Yc0+Y0[self.zoom]+DY[self.zoom],X0[self.zoom]:X0[self.zoom]+DX[self.zoom],0:3] = REDS[Y0[self.zoom]+DY[self.zoom],X0[self.zoom]:X0[self.zoom]+DX[self.zoom],0:3]



#			self.img[195:195+156,908:908+156,0:3] = HUD1[883:883+156,464:464+156,0:3]

#			if(self.zoom0 == 0):
#				HUD1 = plt.imread('../CRDS5/Assets/Zoom0/Plane'+str(lathere)+'x'+str(lonhere)+'.png')
#				self.img[908:908+156,195:195+156,0:3] = HUD1[464:464+156,883:883+156,0:3]
#				self.img[896:1073,117:117+312,0:3] = HUD1[540-88:540+89,960-156:960+156,0:3]
#			if(self.zoom0 == 1):
#				HUD1 = plt.imread('../CRDS6/Assets/Zoom0/Plane'+str(lathere)+'x'+str(lonhere)+'.png')
#				self.img[896:1073,117:117+312,0:3] = HUD1[540-88:540+89,960-156:960+156,0:3]

#			self.img[830:860,:,0] = (self.img[830:860,:,0]+col2[0])/2.0
#			self.img[830:860,:,1] = (self.img[830:860,:,1]+col2[1])/2.0
#			self.img[830:860,:,2] = (self.img[830:860,:,2]+col2[2])/2.0
#			self.img[860:1080,:,0] = (self.img[860:1080,:,0]+col3[0])/2.0
#			self.img[860:1080,:,1] = (self.img[860:1080,:,1]+col3[1])/2.0
#			self.img[860:1080,:,2] = (self.img[860:1080,:,2]+col3[2])/2.0
#			NL1 = len(Provinces[self.selectedProvince].split(' ')[0])
#			self.img[837:855,20:20+17*NL1,:] = Builder.GenTextCOL(Provinces[self.selectedProvince].split(' ')[0],np.array((1.0,1.0,1.0,1.0)),col1)
#



#			self.img[902:820,20:20+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')

	def AddSettlementHUDD(self):
		col1 = np.array((0.5,0.25,0.25,1.0))
		col2 = np.array((1.0,0.5,0.5,1.0))
		col3 = np.array((1.0,1.0,1.0,1.0))
		if(self.BLDD == True):
			self.img[Yc0+15:Yc0+45,1120:,0] = col1[0]#(self.img[15:45,1120:,0]+col2[0])/2.0
			self.img[Yc0+15:Yc0+45,1120:,1] = col1[1]#(self.img[15:45,1120:,1]+col2[1])/2.0
			self.img[Yc0+15:Yc0+45,1120:,2] = col1[2]#(self.img[15:45,1120:,2]+col2[2])/2.0
			self.img[Yc0+45:Yc0+815,1120:,0] = (self.img[Yc0+45:Yc0+815,1120:,0]+2*col3[0])/3.0
			self.img[Yc0+45:Yc0+815,1120:,1] = (self.img[Yc0+45:Yc0+815,1120:,1]+2*col3[1])/3.0
			self.img[Yc0+45:Yc0+815,1120:,2] = (self.img[Yc0+45:Yc0+815,1120:,2]+2*col3[2])/3.0
			outstr = Provinces[self.selectedProvince].split(' ')[0]+' - Settlement Buildings'
			NL1 = len(outstr)
			self.img[Yc0+22:Yc0+40,1140:1140+17*NL1,:] = Builder.GenTextCOL(outstr,np.array((1.0,1.0,1.0,1.0)),col1)
			ct = 0
			for i in list(open('UI/bestfiles.txt')):
				if(ct < 120):
					print(i[0:-1])
					self.img[Yc0+50+64*int(ct/10):Yc0+112+64*int(ct/10),1121+79*(ct%10):1199+79*(ct%10),:] = plt.imread('UI/png/'+i[0:-1]+'.png')
				ct += 1


#			self.img[Yc0+902:Yc0+820,20:20+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')
			


	def AddSettlementHUD(self):

		if(self.BLDD == True):

			colM = Builder.C2I(np.array((1.0,127/255.0,39/255.0,1.0)))
			HUD0 = plt.imread('imgmisc/MTW2/Screenshots/SettlementHUD2.png')
			MASK = Builder.SelectValueC(Builder.CA2I(HUD0),colM)
			self.img[Yc0:Yc1,:,0:3] = self.img[Yc0:Yc1,:,0:3]*(MASK[:,:,0:3])+HUD0[:,:,0:3]*(1-MASK[:,:,0:3])
			TXT = np.array((153/255.0,141/255.0,118/255.0,1.0))
			BKG = np.array((228/255.0,209/255.0,185/255.0,1.0))

#			if(self.selectedProvince != 0):
#				LN = 'Income           '
#				self.img[290:308,1025:1025+17*len(LN),:] = Builder.GenTextCOL(LN,TXT,BKG)
			LN = str(int(PList[self.selectedProvince].GetIncome()))+'    '
			self.img[Yc0+291:Yc0+303,1295:1295+13*len(LN),:] = Builder.GenTextCOL_size(LN,TXT,BKG,12)

			LN = str(PList[self.selectedProvince].GetPublicOrder())+'%    '
			self.img[Yc0+314:Yc0+326,1295:1295+13*len(LN),:] = Builder.GenTextCOL_size(LN,TXT,BKG,12)

			LN = str(PList[self.selectedProvince].GetPopulation())+'  '
			self.img[Yc0+336:Yc0+348,1295:1295+13*len(LN),:] = Builder.GenTextCOL_size(LN,TXT,BKG,12)

			LN = str(PList[self.selectedProvince].GetPopGrowth())+'%    '
			self.img[Yc0+359:Yc0+371,1295:1295+13*len(LN),:] = Builder.GenTextCOL_size(LN,TXT,BKG,12)


#			self.img[Yc0+290:Yc0+308,1295:1295+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')


#			ct = 0
#			for i in list(open('UI/bestfiles.txt')):
#				if(ct < 120):
#					print(i[0:-1])
#					self.img[Yc0+50+64*int(ct/10):Yc0+112+64*int(ct/10),1121+79*(ct%10):1199+79*(ct%10),:] = plt.imread('UI/png/'+i[0:-1]+'.png')
#				ct += 1


#			self.img[902:820,20:20+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')
		
	


	def SaveGrid(self):
		Builder.SaveGrid()

	def Refresh(self):
		pixmap = QPixmap('MANTIS/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		if(len(self.tobeedited) != 0):
			pixmap = QPixmap('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
			if(self.FDisp == True):
				plane = str(self.nearbyplane)+self.isPole
#				Builder.GenIMG7(list(set(self.tobeedited)),[plane],self.zoom)
				pixmap = QPixmap('MANTIS/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')

			else:
				Builder.GenIMG2(self.lat,self.lon,self.zoom)
				pixmap = QPixmap('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()

	def TotalRefresh(self):
		pixmap = QPixmap('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		if(self.FDisp == True):
			Builder.GenIMG5(self.nearbyplane,self.zoom)
			pixmap = QPixmap('MANTIS/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		else:
			Builder.GenIMG2(self.lat,self.lon,self.zoom)
			pixmap = QPixmap('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()

	def ChangeField(self,newfield):
		if(np.average(newfield) == 0 and np.var(newfield) == 0):
#			newfield[83] = 100
#			outarr = np.zeros(shape=(1476225,3),dtype=float)
#			numarr = np.zeros(shape=(1476225),dtype=int)
			newarr = plt.imread('BigData/KoppenClimateData.png')
			krnl = qb.CA2I(plt.imread('ToXY7.png'))
#			outarr[krnl[:,:],0] += newarr[:,:,0]
#			outarr[krnl[:,:],1] += newarr[:,:,1]
#			outarr[krnl[:,:],2] += newarr[:,:,2]
#			numarr[krnl[:,:]] = np.ones(shape=(2700,5400),dtype=int)
			self.field[krnl[:,:]] = qb.CA2I(newarr[:,:])
#			self.field[:] = outarr[:,0]/(0.01+numarr[:])*255*256*256+outarr[:,1]/(0.01+numarr[:])*255*256+outarr[:,2]/(0.01+numarr[:])*255
		else:
			self.field[0:1474562] = newfield[0:1474562]
		self.StandardProcedure()

	def SettlementView(self):
		for c in SETTS:
			cmask = qb.SelectValue2(qb.CA2I(self.img0),c.getPos())
			colmask = qb.SelectValueC(qb.CA2I(self.img0),c.getPos())
			denom = np.sum(cmask[:,:])
			x0 = 0
			y0 = 0
			if(denom != 0):
				x0 = int(np.sum(cmask[:,:]*Xs[:,:])/denom)
				y0 = int(np.sum(cmask[:,:]*Ys[:,:])/denom)
			cimg = c.getIMG(self.zoom)
			self.img[y0-2**max(self.zoom-1,0):y0+2**max(self.zoom-1,0),x0-2**max(self.zoom-1,0):x0+2**max(self.zoom-1,0),:] = self.img[y0-2**max(self.zoom-1,0):y0+2**max(self.zoom-1,0),x0-2**max(self.zoom-1,0):x0+2**max(self.zoom-1,0),:]*qb.SelectValueC(qb.CA2I(cimg[0:2**max(self.zoom,1),0:2**max(self.zoom,1),:]),16777215)+cimg[0:2**max(self.zoom,1),0:2**max(self.zoom,1),:]*(1-qb.SelectValueC(qb.CA2I(cimg[0:2**max(self.zoom,1),0:2**max(self.zoom,1),:]),16777215))


	def CharacterView(self):
		for c in CHARS:
			cmask = qb.SelectValue2(qb.CA2I(self.img0),c.getPos())
			colmask = qb.SelectValueC(qb.CA2I(self.img0),c.getPos())
			denom = np.sum(cmask[:,:])
			x0 = 0
			y0 = 0
			if(denom != 0):
				x0 = int(np.sum(cmask[:,:]*Xs[:,:])/denom)
				y0 = int(np.sum(cmask[:,:]*Ys[:,:])/denom)
			cimg = c.getIMG(self.zoom)
			self.img[y0-2**max(self.zoom-1,0):y0+2**max(self.zoom-1,0),x0-2**max(self.zoom-1,0):x0+2**max(self.zoom-1,0),:] = self.img[y0-2**max(self.zoom-1,0):y0+2**max(self.zoom-1,0),x0-2**max(self.zoom-1,0):x0+2**max(self.zoom-1,0),:]*qb.SelectValueC(qb.CA2I(cimg[0:2**max(self.zoom,1),0:2**max(self.zoom,1),:]),16777215)+cimg[0:2**max(self.zoom,1),0:2**max(self.zoom,1),:]*(1-qb.SelectValueC(qb.CA2I(cimg[0:2**max(self.zoom,1),0:2**max(self.zoom,1),:]),16777215))
			


	def Colorbar(self,field):
		outarr = np.zeros(shape=field.shape,dtype=int)
		outarr[:,:] = field
#		outarr[:,:] = outarr[:,:]*qb.SelectValuesLT(outarr[:,:],16777215)
#		outarr[:,:] = outarr[:,:]*qb.SelectValuesGT(outarr[:,:],0)
		return qb.INTS2Color(outarr)	


	def Colorbar0(self,field):
		outarr = np.zeros(shape=field.shape,dtype=float)
		
#		outarr[:,:] = 16*(field-np.average(field))/np.sqrt(1+np.var(field))
		outarr[:,:] = (field-np.min(field))/(np.max(field)-np.min(field))
		realoutarr = np.ones(shape=(field.shape[0],field.shape[1],4),dtype=np.float16)
		realoutarr[:,:,0] = 1-outarr[:,:]/2
		realoutarr[:,:,1] = 1-outarr[:,:]
		realoutarr[:,:,2] = 1-outarr[:,:]
#		return qb.INTS2Color(2*outarr)	
		return realoutarr
	def StandardProcedure(self):
		Place = ''

		if(self.ODisp == True):
			Place = self.Place
			if(self.zoom < 6):
				self.img = plt.imread(Place+'/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
				self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')

		if(self.VDisp == True):
			Place = 'Assets'
			if(self.zoom < 6):
				self.img = plt.imread(Place+'/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
				self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')

		if(self.FDisp == True):
			Place = 'HGT'
			self.orientation = 'N'
			if(self.zoom >= 0):
				self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
				self.img = self.Colorbar(self.field[qb.CA2I(self.img0[:,:])])
				plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
				pixmap = QPixmap('tmp.png')	
				self.label.setPixmap(pixmap)
				self.resize(pixmap.size())
				self.adjustSize()	
				return
			else:
				self.img = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
				self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
#				plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
#				pixmap = QPixmap('tmp.png')	
#				self.label.setPixmap(pixmap)
#				self.resize(pixmap.size())
#				self.adjustSize()	
#				return
		if(self.RDisp == True):
			Place = 'IMG'
			self.img = plt.imread(Place+'/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
			self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')

#			self.orientation = 'N'
#			if(self.zoom >= 0):
#				self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
#				self.img = self.Colorbar(self.field[qb.CA2I(self.img0[:,:])])
#				plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
#				pixmap = QPixmap('tmp.png')	
#				self.label.setPixmap(pixmap)
#				self.resize(pixmap.size())
#				self.adjustSize()	
#				return

		if(self.zoom == 6):
			lat5 = int(self.lat/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
			lon5 = int(self.lon/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
			if(True):
				if(True):
					centerimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(YF/4):int(3*YF/4),int(XF/4):int(3*XF/4)]
					centerimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(YF/4):int(3*YF/4),int(XF/4):int(3*XF/4)]
					self.img[Y2[:,:],X2[:,:]] = centerimg[:,:]
					self.img0[Y2[:,:],X2[:,:]] = centerimg0[:,:]
					self.img[Y2[:,:]+1,X2[:,:]] = centerimg[:,:]
					self.img0[Y2[:,:]+1,X2[:,:]] = centerimg0[:,:]
					self.img[Y2[:,:],X2[:,:]+1] = centerimg[:,:]
					self.img0[Y2[:,:],X2[:,:]+1] = centerimg0[:,:]
					self.img[Y2[:,:]+1,X2[:,:]+1] = centerimg[:,:]
					self.img0[Y2[:,:]+1,X2[:,:]+1] = centerimg0[:,:]
#				else:
#					leftimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(YF/4):int(3*YF/4),int(3*XF/4):XF]
#					leftimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(YF/4):int(3*YF/4),int(3*XF/4):XF]
#					rightimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[int(YF/4):int(3*YF/4),0:int(1*XF/4)]
#					rightimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[int(YF/4):int(3*YF/4),0:int(1*XF/4)]
#
#
#					self.img[Y2[:,0:480],X2[:,0:480]] = leftimg[:,:]
#					self.img0[Y2[:,0:480],X2[:,0:480]] = leftimg0[:,:]	
#					self.img[Y2[:,0:480],X2[:,0:480]+1] = leftimg[:,:]
#					self.img0[Y2[:,0:480],X2[:,0:480]+1] = leftimg0[:,:]	
#					self.img[Y2[:,0:480]+1,X2[:,0:480]] = leftimg[:,:]
#					self.img0[Y2[:,0:480]+1,X2[:,0:480]] = leftimg0[:,:]	
#					self.img[Y2[:,0:480]+1,X2[:,0:480]+1] = leftimg[:,:]
#					self.img0[Y2[:,0:480]+1,X2[:,0:480]+1] = leftimg0[:,:]	
#
#					self.img[Y2[:,480:960],X2[:,480:960]] = rightimg[:,:]
#					self.img0[Y2[:,480:960],X2[:,480:960]] = rightimg0[:,:]
#					self.img[Y2[:,480:960],X2[:,480:960]+1] = rightimg[:,:]
#					self.img0[Y2[:,480:960],X2[:,480:960]+1] = rightimg0[:,:]
#					self.img[Y2[:,480:960]+1,X2[:,480:960]] = rightimg[:,:]
#					self.img0[Y2[:,480:960]+1,X2[:,480:960]] = rightimg0[:,:]
#					self.img[Y2[:,480:960]+1,X2[:,480:960]+1] = rightimg[:,:]
#					self.img0[Y2[:,480:960]+1,X2[:,480:960]+1] = rightimg0[:,:]
#			elif(self.latdeg != 90):
#				if(self.lon == lon5):
#					bottomimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[0:int(1*YF/4),int(1*XF/4):int(3*XF/4)]
#					bottomimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[0:int(1*YF/4),int(1*XF/4):int(3*XF/4)]
#					topimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5+self.jplvl[5])+'x'+str(lon5)+'.png')[int(3*YF/4):YF,int(1*XF/4):int(3*XF/4)]
#					topimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5+self.jplvl[5])+'x'+str(lon5)+'.png')[int(3*YF/4):YF,int(1*XF/4):int(3*XF/4)]
#
#
#					self.img[Y2[0:480,:],X2[0:480,:]] = topimg[:,:]
#					self.img0[Y2[0:480,:],X2[0:480,:]] = topimg0[:,:]
#					self.img[Y2[0:480,:],X2[0:480,:]+1] = topimg[:,:]
#					self.img0[Y2[0:480,:],X2[0:480,:]+1] = topimg0[:,:]
#					self.img[Y2[0:480,:]+1,X2[0:480,:]] = topimg[:,:]
#					self.img0[Y2[0:480,:]+1,X2[0:480,:]] = topimg0[:,:]
#					self.img[Y2[0:480,:]+1,X2[0:480,:]+1] = topimg[:,:]
#					self.img0[Y2[0:480,:]+1,X2[0:480,:]+1] = topimg0[:,:]
#
#					self.img[Y2[480:960,:],X2[480:960,:]] = bottomimg[:,:]
#					self.img0[Y2[480:960,:],X2[480:960,:]] = bottomimg0[:,:]
#					self.img[Y2[480:960,:],X2[480:960,:]+1] = bottomimg[:,:]
#					self.img0[Y2[480:960,:],X2[480:960,:]+1] = bottomimg0[:,:]
#					self.img[Y2[480:960,:]+1,X2[480:960,:]] = bottomimg[:,:]
#					self.img0[Y2[480:960,:]+1,X2[480:960,:]] = bottomimg0[:,:]
#					self.img[Y2[480:960,:]+1,X2[480:960,:]+1] = bottomimg[:,:]
#					self.img0[Y2[480:960,:]+1,X2[480:960,:]+1] = bottomimg0[:,:]
#				else:
#					bottomleftimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[0:int(1*YF/4),int(3*XF/4):XF]
#					bottomleftimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[0:int(1*YF/4),int(3*XF/4):XF]
#					bottomrightimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[0:int(1*YF/4),0:int(1*XF/4)]
#					bottomrightimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[0:int(1*YF/4),0:int(1*XF/4)]
#					topleftimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(3*YF/4):YF,int(3*XF/4):XF]
#					topleftimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(3*YF/4):YF,int(3*XF/4):XF]
#					toprightimg = plt.imread(Place+'/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[int(3*YF/4):YF,0:int(1*XF/4)]
#					toprightimg0 = plt.imread('IMG/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[int(3*YF/4):YF,0:int(1*XF/4)]
#
#
#					self.img[Y2[0:480,0:480],X2[0:480,0:480]] = topleftimg[:,:]
#					self.img0[Y2[0:480,0:480],X2[0:480,0:480]] = topleftimg0[:,:]
#					self.img[Y2[0:480,0:480],X2[0:480,0:480]+1] = topleftimg[:,:]
#					self.img0[Y2[0:480,0:480],X2[0:480,0:480]+1] = topleftimg0[:,:]
#					self.img[Y2[0:480,0:480]+1,X2[0:480,0:480]] = topleftimg[:,:]
#					self.img0[Y2[0:480,0:480]+1,X2[0:480,0:480]] = topleftimg0[:,:]
#					self.img[Y2[0:480,0:480]+1,X2[0:480,0:480]+1] = topleftimg[:,:]
#					self.img0[Y2[0:480,0:480]+1,X2[0:480,0:480]+1] = topleftimg0[:,:]
#
#					self.img[Y2[0:480,480:960],X2[0:480,480:960]] = toprightimg[:,:]
#					self.img0[Y2[0:480,480:960],X2[0:480,480:960]] = toprightimg0[:,:]
#					self.img[Y2[0:480,480:960],X2[0:480,480:960]+1] = toprightimg[:,:]
#					self.img0[Y2[0:480,480:960],X2[0:480,480:960]+1] = toprightimg0[:,:]
#					self.img[Y2[0:480,480:960]+1,X2[0:480,480:960]] = toprightimg[:,:]
#					self.img0[Y2[0:480,480:960]+1,X2[0:480,480:960]] = toprightimg0[:,:]
#					self.img[Y2[0:480,480:960]+1,X2[0:480,480:960]+1] = toprightimg[:,:]
#					self.img0[Y2[0:480,480:960]+1,X2[0:480,480:960]+1] = toprightimg0[:,:]
#
#					self.img[Y2[480:960,0:480],X2[480:960,0:480]] = bottomleftimg[:,:]
#					self.img0[Y2[480:960,0:480],X2[480:960,0:480]] = bottomleftimg0[:,:]
#					self.img[Y2[480:960,0:480],X2[480:960,0:480]+1] = bottomleftimg[:,:]
#					self.img0[Y2[480:960,0:480],X2[480:960,0:480]+1] = bottomleftimg0[:,:]
#					self.img[Y2[480:960,0:480]+1,X2[480:960,0:480]] = bottomleftimg[:,:]
#					self.img0[Y2[480:960,0:480]+1,X2[480:960,0:480]] = bottomleftimg0[:,:]
#					self.img[Y2[480:960,0:480]+1,X2[480:960,0:480]+1] = bottomleftimg[:,:]
#					self.img0[Y2[480:960,0:480]+1,X2[480:960,0:480]+1] = bottomleftimg0[:,:]
#
#					self.img[Y2[480:960,480:960],X2[480:960,480:960]] = bottomrightimg[:,:]
#					self.img0[Y2[480:960,480:960],X2[480:960,480:960]] = bottomrightimg0[:,:]
#					self.img[Y2[480:960,480:960],X2[480:960,480:960]+1] = bottomrightimg[:,:]
#					self.img0[Y2[480:960,480:960],X2[480:960,480:960]+1] = bottomrightimg0[:,:]
#					self.img[Y2[480:960,480:960]+1,X2[480:960,480:960]] = bottomrightimg[:,:]
#					self.img0[Y2[480:960,480:960]+1,X2[480:960,480:960]] = bottomrightimg0[:,:]
#					self.img[Y2[480:960,480:960]+1,X2[480:960,480:960]+1] = bottomrightimg[:,:]
#					self.img0[Y2[480:960,480:960]+1,X2[480:960,480:960]+1] = bottomrightimg0[:,:]
#		
#		
#			
		if(self.zoom >= 2 and self.zoom <= 5):
			self.hgt = plt.imread('HGT/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		if(self.zoom == 6):
			lat5 = int(self.lat/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
			lon5 = int(self.lon/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
			if(True):
				if(True):
					centerimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(YF/4):int(3*YF/4),int(XF/4):int(3*XF/4)]
					self.hgt[Y2[:,:],X2[:,:]] = centerimg0[:,:]
					self.hgt[Y2[:,:]+1,X2[:,:]] = centerimg0[:,:]
					self.hgt[Y2[:,:],X2[:,:]+1] = centerimg0[:,:]
					self.hgt[Y2[:,:]+1,X2[:,:]+1] = centerimg0[:,:]
#				else:
#					leftimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(YF/4):int(3*YF/4),int(3*XF/4):XF]
#					rightimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[int(YF/4):int(3*YF/4),0:int(1*XF/4)]
#
#					self.hgt[Y2[:,0:480],X2[:,0:480]] = leftimg0[:,:]	
#					self.hgt[Y2[:,0:480],X2[:,0:480]+1] = leftimg0[:,:]	
#					self.hgt[Y2[:,0:480]+1,X2[:,0:480]] = leftimg0[:,:]	
#					self.hgt[Y2[:,0:480]+1,X2[:,0:480]+1] = leftimg0[:,:]	
#
#					self.hgt[Y2[:,480:960],X2[:,480:960]] = rightimg0[:,:]
#					self.hgt[Y2[:,480:960],X2[:,480:960]+1] = rightimg0[:,:]
#					self.hgt[Y2[:,480:960]+1,X2[:,480:960]] = rightimg0[:,:]
#					self.hgt[Y2[:,480:960]+1,X2[:,480:960]+1] = rightimg0[:,:]
#			elif(self.latdeg != 90):
#				if(self.lon == lon5):
#					bottomimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[0:int(1*YF/4),int(1*XF/4):int(3*XF/4)]
#					topimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5+self.jplvl[5])+'x'+str(lon5)+'.png')[int(3*YF/4):YF,int(1*XF/4):int(3*XF/4)]
#
#
#					self.hgt[Y2[0:480,:],X2[0:480,:]] = topimg0[:,:]
#					self.hgt[Y2[0:480,:],X2[0:480,:]+1] = topimg0[:,:]
#					self.hgt[Y2[0:480,:]+1,X2[0:480,:]] = topimg0[:,:]
#					self.hgt[Y2[0:480,:]+1,X2[0:480,:]+1] = topimg0[:,:]
#
#					self.hgt[Y2[480:960,:],X2[480:960,:]] = bottomimg0[:,:]
#					self.hgt[Y2[480:960,:],X2[480:960,:]+1] = bottomimg0[:,:]
#					self.hgt[Y2[480:960,:]+1,X2[480:960,:]] = bottomimg0[:,:]
#					self.hgt[Y2[480:960,:]+1,X2[480:960,:]+1] = bottomimg0[:,:]
#				else:
#					bottomleftimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[0:int(1*YF/4),int(3*XF/4):XF]
#					bottomrightimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[0:int(1*YF/4),0:int(1*XF/4)]
#					topleftimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5)+'.png')[int(3*YF/4):YF,int(3*XF/4):XF]
#					toprightimg0 = plt.imread('HGT/Zoom5/Plane'+str(lat5)+'x'+str(lon5+self.jplvl[5])+'.png')[int(3*YF/4):YF,0:int(1*XF/4)]
#
#
#					self.hgt[Y2[0:480,0:480],X2[0:480,0:480]] = topleftimg0[:,:]
#					self.hgt[Y2[0:480,0:480],X2[0:480,0:480]+1] = topleftimg0[:,:]
#					self.hgt[Y2[0:480,0:480]+1,X2[0:480,0:480]] = topleftimg0[:,:]
#					self.hgt[Y2[0:480,0:480]+1,X2[0:480,0:480]+1] = topleftimg0[:,:]
#
#					self.hgt[Y2[0:480,480:960],X2[0:480,480:960]] = toprightimg0[:,:]
#					self.hgt[Y2[0:480,480:960],X2[0:480,480:960]+1] = toprightimg0[:,:]
#					self.hgt[Y2[0:480,480:960]+1,X2[0:480,480:960]] = toprightimg0[:,:]
#					self.hgt[Y2[0:480,480:960]+1,X2[0:480,480:960]+1] = toprightimg0[:,:]
#
#					self.img0[Y2[480:960,0:480],X2[480:960,0:480]] = bottomleftimg0[:,:]
#					self.img0[Y2[480:960,0:480],X2[480:960,0:480]+1] = bottomleftimg0[:,:]
#					self.img0[Y2[480:960,0:480]+1,X2[480:960,0:480]] = bottomleftimg0[:,:]
#					self.img0[Y2[480:960,0:480]+1,X2[480:960,0:480]+1] = bottomleftimg0[:,:]
#
#					self.hgt[Y2[480:960,480:960],X2[480:960,480:960]] = bottomrightimg0[:,:]
#					self.hgt[Y2[480:960,480:960],X2[480:960,480:960]+1] = bottomrightimg0[:,:]
#					self.hgt[Y2[480:960,480:960]+1,X2[480:960,480:960]] = bottomrightimg0[:,:]
#					self.hgt[Y2[480:960,480:960]+1,X2[480:960,480:960]+1] = bottomrightimg0[:,:]
			

		if(self.orientation == 'E'):
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.hgt = qb.RotateAnticlockwise90(self.hgt)
			self.hgt = qb.RotateAnticlockwise90(self.hgt)
			self.hgt = qb.RotateAnticlockwise90(self.hgt)
		elif(self.orientation == 'S'):
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.hgt = qb.RotateAnticlockwise90(self.hgt)
			self.hgt = qb.RotateAnticlockwise90(self.hgt)
		elif(self.orientation == 'W'):
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.hgt = qb.RotateAnticlockwise90(self.hgt)

			self.img[Yc0:Yc1,:] = qb.Rotate(self.img,self.angle,self.zoom)
			selfimg0 = Builder.ChangeColorI(self.img0,0,1474562)
			self.img0[:,:] = qb.Rotate(self.img0,self.angle,self.zoom)
			self.img0 = qb.FastFillBay(self.img0,1,0)
			self.img0 = qb.FastFillBay(self.img0,1,2)
			self.img0 = qb.FastFillBay(self.img0,0,1)
			self.img0 = qb.FastFillBay(self.img0,2,1)
			self.img0 = Builder.ChangeColorI(self.img0,0,1476224)
			self.img0 = Builder.ChangeColorI(self.img0,1474562,0)

#		self.img0 = Builder.ChangeColorI(self.img0,1476225,1476224)

		if(self.PTSD == True):		
			self.HighlightPoint()
		if(self.PRSD == True):
			self.HighlightProvince()
		if(self.CLDD == True):
			self.AddCloudView()
		if(self.TODD == True):
			self.AddDayView()			
		if(self.zoom >= 0 and self.TOPD == True):
			self.img[:,:] = self.img[:,:]*qb.shadr[:,:]**(qb.PositiveOrZeros(qb.DDX4(Builder.CA2I(self.hgt[:,:])))**2 + qb.PositiveOrZeros(qb.DDY4(Builder.CA2I(self.hgt[:,:])))**2)**(0.3)
			self.img = qb.Lift(self.img,qb.CA2I(self.hgt),self.zoom,self.angle)
		if(self.angle != 0 and self.zoom >= 2):
#			fullimg = np.zeros(shape=(3840,3840,4),dtype=np.float16)
#			tmpimg = plt.imread('Assets/Zoom'+str(self.zoom-1)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
#			fullimg[Ys[:,:]*2,Xs[:,:]] = tmpimg[:,:]
#			fullimg[Ys[:,:]*2+1,Xs[:,:]] = tmpimg[:,:]
#			fullimg[Ys[:,:]*2,Xs[:,:]+1] = tmpimg[:,:]
#			fullimg[Ys[:,:]*2+1,Xs[:,:]+1] = tmpimg[:,:]

			self.img[:,:] = qb.Rotate(self.img,self.angle,self.zoom)
#			self.img[:,:,0:3] = qb.RotateLarge(fullimg,self.angle,self.zoom)[2460-1920:2460,0:1920,0:3]
			if(self.TOPD == True):
				self.img[:,:,0] = self.img[:,:,0]*E2[:,:]
				self.img[:,:,1] = self.img[:,:,1]*E2[:,:]
				self.img[:,:,2] = self.img[:,:,2]*E2[:,:]
			self.img0[:,:] = qb.Rotate0(self.img0,self.angle,self.zoom,E2)
		if(self.GADD == True):
			self.AddGalaxyView()
#		self.img = qb.UFF(self.img)
#		self.img = qb.UFF(self.img)
#		self.SettlementView()
#		self.CharacterView()

		self.AddDTCaption()
#		plt.imsave(fname='tmp.png',arr=self.img)
		plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])

#		pixmap = QPixmap(Place+'/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		pixmap = QPixmap('tmp.png')	
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()

	def WindowedProcedure(self):
		Place = ''
		if(self.VDisp == True):
			Place = 'Assets'
		if(self.FDisp == True):
			Place = 'HGT'
		if(self.RDisp == True):
			Place = 'IMG'
		self.img = plt.imread(Place+'/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')

		if(self.orientation == 'E'):
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
		elif(self.orientation == 'S'):
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)
		elif(self.orientation == 'W'):
			self.img = qb.RotateAnticlockwise90(self.img)
			self.img0 = qb.RotateAnticlockwise90(self.img0)

		self.img0 = Builder.ChangeColorI(self.img0,1476225,1476224)

		if(self.PTSD == True):		
			self.HighlightPoint()
		if(self.PRSD == True):
			self.HighlightProvince()
		if(self.CLDD == True):
			self.AddCloudView()
		if(self.TODD == True):
			self.AddDayView()			
		if(self.GADD == True):
			self.AddGalaxyView()
		if(self.angle != 0 and self.zoom >= 0):
			self.img[:,:] = qb.Rotate(self.img,self.angle,self.zoom)

		self.AddDTCaption()

		plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
		pixmap = QPixmap('tmp.png')	
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()

	def moveNorth(self):		

		if(self.orientation == 'N'):
			self.lat = min(self.lat+self.jplvl[self.zoom],900000)
			self.latdeg = self.lat/self.scale
		elif(self.orientation == 'E'):
			self.lon = (self.lon+self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		elif(self.orientation == 'S'):
			self.lat = max(self.lat-self.jplvl[self.zoom],-900000)
			self.latdeg = self.lat/self.scale
		elif(self.orientation == 'W'):
			self.lon = (self.lon-self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		if(self.zoom >= 5):
			self.StandardProcedure()
		else:
			self.StandardProcedure()

	def moveWest(self):
		if(self.orientation == 'N'):
			self.lon = (self.lon-self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		elif(self.orientation == 'E'):
			self.lat = min(self.lat+self.jplvl[self.zoom],900000)
			self.latdeg = self.lat/self.scale
		elif(self.orientation == 'S'):
			self.lon = (self.lon+self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		elif(self.orientation == 'W'):
			self.lat = max(self.lat-self.jplvl[self.zoom],-900000)
			self.latdeg = self.lat/self.scale

		if(self.zoom >= 5):
			self.StandardProcedure()
		else:
			self.StandardProcedure()
	def moveSouth(self):

		if(self.orientation == 'N'):
			self.lat = max(self.lat-self.jplvl[self.zoom],-900000)
			self.latdeg = self.lat/self.scale
		elif(self.orientation == 'E'):
			self.lon = (self.lon-self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		elif(self.orientation == 'S'):
			self.lat = min(self.lat+self.jplvl[self.zoom],900000)
			self.latdeg = self.lat/self.scale
		elif(self.orientation == 'W'):
			self.lon = (self.lon+self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24

		if(self.zoom >= 5):
			self.StandardProcedure()
		else:
			self.StandardProcedure()
	def moveEast(self):
		if(self.orientation == 'N'):
			self.lon = (self.lon+self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		elif(self.orientation == 'E'):
			self.lat = max(self.lat-self.jplvl[self.zoom],-900000)
			self.latdeg = self.lat/self.scale
		elif(self.orientation == 'S'):
			self.lon = (self.lon-self.jplvl[self.zoom])%3600000
			self.londeg = self.lon/self.scale
			self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		elif(self.orientation == 'W'):
			self.lat = min(self.lat+self.jplvl[self.zoom],900000)
			self.latdeg = self.lat/self.scale

		if(self.zoom >= 5):
			self.StandardProcedure()
		else:
			self.StandardProcedure()

	def IncreaseTime(self):
		if(self.FDisp == True):
#			for i in range(10):
#				print(str(i))
			print('Incrementing one timestep')
			self.field, self.fieldt = qb.WaveEvolve(self.field,self.fieldt,1,0.01)
			self.StandardProcedure()
			return
		self.hour += 1
		oldday = self.day
		if(self.hour == 24):
			self.hour = 0
			self.day += 1
			if(self.day >= Ml[self.month]):
				if(self.month == 1 and self.year%4 == 0):
					self.day += 1
				else:
					self.day = 0
					self.month += 1
					if(self.month == 12):
						self.year += 1
			self.dayofyear = MonthStarts[self.month]+self.day

		if(self.day != oldday):
			self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')

#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
			print('Effective Time of Year: '+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ])
			CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
			ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
			for i in range(24):							#	Updates global CLDSD info which
				CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
				ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	


		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))

		self.StandardProcedure()

	def IncreaseTimeDay(self):
		oldday = self.day
		if(self.day == 27 and self.month == 1 and self.year%4 == 0):
			self.day += 1
		elif(self.day >= ( Ml[self.month]-1) ):
			self.day = 0
			self.month += 1
			if(self.month == 12):
				self.year += 1
				self.month = 0 
		else:
			self.day += 1

		self.dayofyear = MonthStarts[self.month]+self.day

		if(self.day != oldday):
			self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')

#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
			print('Effective Time of Year: '+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ])
			CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
			ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
			for i in range(24):							#	Updates global CLDSD info which
				CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
				ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	

		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))
			
		self.StandardProcedure()

	def IncreaseTimeMonth(self):
		self.month += 1
		self.day = min(Ml[self.month]-1,self.day)
		if(self.month == 12):
			self.year += 1
			self.month = 0

		self.dayofyear = MonthStarts[self.month]+self.day
		self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')

#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
		print('Effective Time of Year: '+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ])
		CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
		ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
		for i in range(24):							#	Updates global CLDSD info which
			CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
			ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	

		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))
			
		self.StandardProcedure()

	def IncreaseTimeYear(self):
		if(self.month == 1 and self.year%4 == 0 and self.day == 28 and self.year != 2020):
			self.day = self.day-1
		if(self.year < 2020):
			self.year += 1
		else:
			print('Cannot increase year any further!')
			
		self.dayofyear = MonthStarts[self.month]+self.day
		self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')

#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
		print('Effective Time of Year: '+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ])
		CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
		ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
		for i in range(24):							#	Updates global CLDSD info which
			CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
			ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	

		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))

		self.StandardProcedure()

		
	def DecreaseTime(self):
		oldday = self.day
		if(self.hour == 0 and self.day == 0 and self.month == 0 and self.year == 1980):
			print('Cannot go back in time any further...')
		else:
			self.hour = self.hour -1
			if(self.hour < 0):
				self.hour = 23
				self.day = self.day-1
				if(self.day < 0):
					if(self.month == 2 and self.year%4 == 0):
						self.day = 29
						self.month = 1
					else:
						self.day = Ml[(self.month-1)%12]
						self.month = self.month-1
						if(self.month < 0):
							self.year = self.year -1
							self.month = 11
							self.day = 30
		if(self.day != oldday):
			self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')
#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
			CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
			ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
			for i in range(24):							#	Updates global CLDSD info which
				CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
				ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	


		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))

		self.StandardProcedure()

	def DecreaseTimeDay(self):
		oldday = self.day
		if(self.day == 0 and self.month == 0 and self.year == 1980):
			print('Cannot go back in time any further...')
		else:
			self.day = self.day -1
			if(self.day < 0):
				if(self.month == 2 and self.year%4 == 0):
					self.day = 28
					self.month = 1
				else:
					self.day = Ml[(self.month-1)%12]
					self.month = self.month-1
					if(self.month < 0):
						self.year = self.year-1
						self.month = 11
						self.day = 30
		self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')
#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
		CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
		ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
		for i in range(24):							#	Updates global CLDSD info which
			CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
			ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	
		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))

		self.StandardProcedure()

	def DecreaseTimeMonth(self):
		oldday = self.day
		if(self.month == 0 and self.year == 1980):
			print('Cannot go back in time any further...')
		else:
			self.month = self.month -1
			if(self.month < 0):
				self.year = self.year -1
				self.month = 11
			self.day = min(Ml[self.month]-1,self.day)
		self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')
#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
		CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
		ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
		for i in range(24):							#	Updates global CLDSD info which
			CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
			ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	
		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))

		self.StandardProcedure()

	def DecreaseTimeYear(self):
		oldday = self.day
		if(self.month == 1 and self.year%4 == 0 and self.day == 28 and self.year != 1980):
			self.day = self.day-1

		if(self.year == 1980):
			print('Cannot go back in time any further...')
		else:
			self.year = self.year -1
	
		self.cld = plt.imread('TISDat/CLD/CLDS'+str(self.year)+Mn[self.month]+Dn[self.day]+'.png')
		self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[self.month]+Dn[self.day]+'.png')
#			self.ang = plt.imread('TISDat/ANG/ANGS'+Mn[MONTH[(self.dayofyear + 182)%365]]+Dn[(self.dayofyear + 182)%365-MonthStarts[MONTH[(self.dayofyear + 182)%365]] ]+'.png')
		CLDSI[:,:] = self.cld[:,:,0]*255*256*256+self.cld[:,:,1]*255*256+self.cld[:,:,2]*255	
		ANGSI[:,:] = self.ang[:,:,0]*255*256*256+self.ang[:,:,1]*255*256+self.ang[:,:,2]*255	
		for i in range(24):							#	Updates global CLDSD info which
			CLDSD[i,:,:] = (CLDSI%(2**(i+1))-CLDSI%(2**i))/(2**i)		#	should be done daily
			ANGSD[i,:,:] = (ANGSI%(2**(i+1))-ANGSI%(2**i))/(2**i)	
		self.localtime = (self.hour+int(self.lon/self.jplvl[0]))%24
		self.dayofyear = MonthStarts[self.month]+self.day
		self.SolarOffset = (self.dayofyear - (183.0625-7))

		self.StandardProcedure()


	def FieldDisplay(self):
		self.FDisp = True
		self.VDisp = False
		self.RDisp = False

		self.StandardProcedure()

	def VisualDisplay(self):
		self.FDisp = False
		self.VDisp = True
		self.RDisp = False

		self.StandardProcedure()

	def RawDisplay(self):
		self.FDisp = False
		self.VDisp = False
		self.RDisp = True
		self.StandardProcedure()

	def OtherDisplay(self,Name):
		self.Place = Name
		self.ODisp = True
		self.RDisp = False
		self.VDisp = False
		self.FDisp = False
		self.StandardProcedure()


	def CldDisplay(self):
		if(self.CLDD == False):
			self.CLDD = True
			print('Cloud Display Activated...')
		else:
			self.CLDD = False
			print('Cloud Display Deactivated...')
		self.StandardProcedure()			


	def TopoDisplay(self):
		if(self.TOPD == False):
			self.TOPD = True	
			print('Topographic Display Activated...')			
		else:
			print('Topographic Display Deactivated...')
			self.TOPD = False
		self.StandardProcedure()
	

	def GalacticDisplay(self):
		if(self.GADD == False):
			self.GADD = True	
			print('Galactic Display Activated...')			
		else:
			print('Galactic Display Deactivated...')
			self.GADD = False
		self.StandardProcedure()

#	Run After DayView	
	def AddGalaxyView(self):
#		if(self.zoom <= 1):
			A = plt.imread('Assets/STAR2/Plane'+str(self.lat)+'x'+str((-1*self.lon-10000*int(self.SolarOffset)-150000*(self.hour-12) )%3600000 )+'.png')
#			rev1920 = list(range(1919,-1,-1))

			if(self.orientation == 'E'):
				A = qb.RotateAnticlockwise90(A)
				A = qb.RotateAnticlockwise90(A)
				A = qb.RotateAnticlockwise90(A)
			elif(self.orientation == 'S'):
				A = qb.RotateAnticlockwise90(A)
				A = qb.RotateAnticlockwise90(A)
			elif(self.orientation == 'W'):
				A = qb.RotateAnticlockwise90(A)

			self.img[:,:,0] = A[:,:,0]*(1-E2[:,:])+self.img[:,:,0]*E2[:,:]
			self.img[:,:,1] = A[:,:,1]*(1-E2[:,:])+self.img[:,:,1]*E2[:,:]
			self.img[:,:,2] = A[:,:,2]*(1-E2[:,:])+self.img[:,:,2]*E2[:,:]

#			if(self.zoom == 0):
#				self.img[:,:,:] = A[:,rev1920[:],:]*D0[:,:,:]+self.img[:,:,:]*(1-D0[:,:,:])
#				self.img[:,:,:] = A[:,:,:]*D0[:,:,:]+self.img[:,:,:]*(1-D0[:,:,:])

#				self.img[:,:,0] = A[:,:,0]*(1-E0[:,:])+self.img[:,:,0]*E0[:,:]
#				self.img[:,:,1] = A[:,:,1]*(1-E0[:,:])+self.img[:,:,1]*E0[:,:]
#				self.img[:,:,2] = A[:,:,2]*(1-E0[:,:])+self.img[:,:,2]*E0[:,:]

#			elif(self.zoom == 1):
#				self.img[:,:,:] = A[:,rev1920[:],:]*D1[:,:,:]+self.img[:,:,:]*(1-D1[:,:,:])
#				self.img[:,:,:] = A[:,:,:]*D1[:,:,:]+self.img[:,:,:]*(1-D1[:,:,:])
#			else:
#				self.img[:,:,0] = A[:,:,0]*(1-E2[:,:])+self.img[:,:,0]*E2[:,:]
#				self.img[:,:,1] = A[:,:,1]*(1-E2[:,:])+self.img[:,:,1]*E2[:,:]
#				self.img[:,:,2] = A[:,:,2]*(1-E2[:,:])+self.img[:,:,2]*E2[:,:]
#		plt.imsave(fname='tmp.png',arr=C)

	def Update(self):
		self.StandardProcedure()
		print('Updated GUI Display')		

	def Update2(self):
		if(self.zoom >= 2):
#			Builder.GenIMG2(self.lat/10000.0,self.lon/10000,self.zoom)
			Builder.AssetFill(self.lat,self.lon,self.zoom)
			Builder.IMGFill(self.lat,self.lon,self.zoom)
		self.StandardProcedure()
		print('Updated GUI Display')		
		

#	Run 1st
	def AddCloudView(self):
#		print('Adding Cloud View')
		cldcol = np.array((0.9,0.9,0.9,1.0))
		cldcol = np.array((1.0,1.0,1.0,1.0))

		cr1[:,:] = CA2I(self.img0[:,:])/1215
		cr2[:,:] = CA2I(self.img0[:,:])%1215
		COL4[:,:] = CLDSD[self.hour,cr1[:,:],cr2[:,:]]
		COL5[:,:,0] = COL4[:,:]
		COL5[:,:,1] = COL4[:,:]
		COL5[:,:,2] = COL4[:,:]						
#		B = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.img[:,:,0] = (4*(self.img[:,:,0]*(1-COL5[:,:,0])+COL5[:,:,0]*cldcol[0])+self.img[:,:,0])/5.0
		self.img[:,:,1] = (4*(self.img[:,:,1]*(1-COL5[:,:,1])+COL5[:,:,1]*cldcol[1])+self.img[:,:,1])/5.0
		self.img[:,:,2] = (4*(self.img[:,:,2]*(1-COL5[:,:,2])+COL5[:,:,2]*cldcol[2])+self.img[:,:,2])/5.0


#	Run After CloudView
	def AddDayView(self):
		cr1[:,:] = CA2I(self.img0)/1215
		cr2[:,:] = CA2I(self.img0)%1215
		COL4[:,:] = ANGSD[self.hour,cr1[:,:],cr2[:,:]]
		COL5[:,:,0] = COL4[:,:]
		COL5[:,:,1] = COL4[:,:]
		COL5[:,:,2] = COL4[:,:]
		if(self.zoom == 0):
			self.img[:,:,:] = self.img[:,:,:]*(0.25+(COL5[:,:,:]))*0.8*(E0[:,:,:])
		elif(self.zoom == 1):
			self.img[:,:,:] = self.img[:,:,:]*(0.25+(COL5[:,:,:]))*0.8*(E1[:,:,:])
		else:
			self.img[:,:,:] = self.img[:,:,:]*(0.25+(COL5[:,:,:]))*0.8
			

	def HighlightGrid(self):
		if(self.GRID == False):
			self.GRID = True
			print('Grid Display Activated...')
		else:
			self.GRID = False
			print('Grid Display Deactivated...')
		self.StandardProcedure()
			
	def HighlightProvince(self):
		if(self.GRID == True):
#			cr1[:,:] = CA2I(self.img0)/1215
#			cr2[:,:] = CA2I(self.img0)%1215
#			brd = plt.imread('Mappings/BRD.png')[cr1[:,:],cr2[:,:]]
			self.img[:,:,0] = self.img[:,:,0]*(1-qb.SelectBorders(PRVL[Builder.CA2I(self.img0)],2))
			self.img[:,:,1] = self.img[:,:,1]*(1-qb.SelectBorders(PRVL[Builder.CA2I(self.img0)],2))
			self.img[:,:,2] = self.img[:,:,2]*(1-qb.SelectBorders(PRVL[Builder.CA2I(self.img0)],2))
		if(self.PRSD == True and self.selectedProvince != 0 ):
			PT = int(Grid[self.selected].getProvince())
			print('Highlighting Province '+Provinces[PT].split(' ')[0])
#		PRVID = np.zeros(shape=(1080,1920),dtype=int)
			PRVID = PRVL[CA2I(self.img0)]
			self.img = Builder.HighlightSelected1( self.img , Builder.SelectValue2(PRVID,PT) )

	def HighlightPoint(self):
		if(self.GRID == True and self.zoom >= 2):
			self.img[:,:,0] = self.img[:,:,0]*(1-qb.SelectBorders(Builder.CA2I(self.img0),1))
			self.img[:,:,1] = self.img[:,:,1]*(1-qb.SelectBorders(Builder.CA2I(self.img0),1))
			self.img[:,:,2] = self.img[:,:,2]*(1-qb.SelectBorders(Builder.CA2I(self.img0),1))
#			if(self.PRSD == True):
			self.img[:,:,0] = self.img[:,:,0]*(1-qb.SelectBorders(PRVL[Builder.CA2I(self.img0)],2))
			self.img[:,:,1] = self.img[:,:,1]*(1-qb.SelectBorders(PRVL[Builder.CA2I(self.img0)],2))
			self.img[:,:,2] = self.img[:,:,2]*(1-qb.SelectBorders(PRVL[Builder.CA2I(self.img0)],2))
		if(self.PTSD == True):	
			if(self.selected < 1474562):#len(Grid)):
				print('Highlighting Point '+str(self.selected))
				self.selectedProvince = int(Grid[self.selected].getProvince())
				print('NearbyPl: '+str(Builder.NearestPlane(self.lat,self.lon)))
				tmpcol = Builder.INT2Color(self.selected)
				print('Color:    '+str(int(tmpcol[0]*255.0))+' '+str(int(tmpcol[1]*255.0))+' '+str(int(tmpcol[2]*255.0)))
				print('DateTime: '+Mn[self.month]+'/'+Dn[self.day]+'/'+str(self.year)+' at '+Hr[self.hour]+':00 GMT')
				print('Clock:    '+Hr[self.localtime]+':00')
				Grid[self.selected].PrintStats()
				PTSID = CA2I(self.img0)
				self.img = Builder.HighlightSelected1( self.img , Builder.SelectValue2(PTSID , self.selected ) )
			else:
				print('NearbyPl: '+str(Builder.NearestPlane(self.lat,self.lon)))
				print('Color:    NullSpace')
				print('Point ID: NullSpace')
				print('NN List:  NullSpace')
				print('Province: NullSpace')
				print('XYZ Pos:  NullSpace')
				print('Lat:      NullSpace')
				print('Lon:      NullSpace')
				print('Border:   NullSpace')
				print('IGBP:     NullSpace')
				print('Terrain:  NullSpace')
				print('Climate:  NullSpace')
				print(' ')
			

		

	def RemoveGalaxyView(self):
		C = np.zeros(shape=A.shape,dtype=np.float16)
		B = plt.imread('tmp.png')
		  
		C[:,:,:] = B*(1-D)
		plt.imsave(fname='tmp.png',arr=C)	
	

	def TimeOfDayDisplay(self):

		if(self.TODD == False):
			self.TODD = True
			print('Time of Day Display Activated...')
		else:
			self.TODD = False
			print('Time of Day Display Deactivated...')
		self.StandardProcedure()

	def IncreaseAngle(self):
		self.angle = min(self.angle+15,90)
		E2[:,:] = plt.imread('Assets/EShade/Zoom'+str(self.zoom)+'/A'+str(self.angle)+'.png')[:,:,0]
		if(self.angle == 90):
			self.HUDD = True
		self.StandardProcedure()

	def DecreaseAngle(self):
		self.angle = max(self.angle-15,0)
		if(self.angle == 0):
			E2[:,:] = np.ones(shape=(1920,1920),dtype=int)
		else:
			E2[:,:] = plt.imread('Assets/EShade/Zoom'+str(self.zoom)+'/A'+str(self.angle)+'.png')[:,:,0]
		self.StandardProcedure()


	def zoomIn(self):
		if(self.zoom >= 6):
			print("Can't zoom any further!")
			if(self.angle == 75):
				print("Can't zoom any further!")
			self.angle = min(self.angle+15,75)	
			self.StandardProcedure()
		else:
			self.zoom = self.zoom+1
			if(self.zoom == 1):
				D[:,:,:] = D1[:,:,:]
				print('Changing star mask to 1')
			if(self.zoom >= 2):
				E2[:,:] = plt.imread('Assets/EShade/Zoom'+str(self.zoom)+'/A'+str(self.angle)+'.png')[:,:,0]		
			self.StandardProcedure()

	def zoomOut(self):
		if(self.zmlvl[self.zoom] == 0):
			print("Can't zoom out any further!")
#		elif(self.zoom >= 6):
#			if(self.angle != 0):
#				self.angle = self.angle - 15
#			else:
#				self.lat = int(self.lat/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
#				self.lon = int(self.lon/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
#				self.zoom = self.zoom-1
		else:
			self.lat = int(self.lat/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
			self.lon = int(self.lon/self.jplvl[self.zoom-1])*self.jplvl[self.zoom-1]
			self.zoom = self.zoom-1
			if(self.zoom == 1):
				D[:,:,:] = D1[:,:,:]
				print('Changing StarMask to 1')
			if(self.zoom == 0):
				D[:,:,:] = D0[:,:,:]
				print('Changing StarMask to 0')
			if(self.zoom >= 2):
#				D[:,:,:] = D0[:,:,:]
				E2[:,:] = plt.imread('Assets/EShade/Zoom'+str(self.zoom)+'/A'+str(self.angle)+'.png')[:,:,0]		
				print('Changing StarMask to '+str(self.zoom))
		self.StandardProcedure()

	def RotateClockwise(self):
		if(self.orientation == 'N'):
			self.orientation = 'E'
		elif(self.orientation == 'E'):
			self.orientation = 'S'
		elif(self.orientation == 'S'):
			self.orientation = 'W'
		elif(self.orientation == 'W'):
			self.orientation = 'N'
		self.StandardProcedure()

	def RotateAnticlockwise(self):
		if(self.orientation == 'N'):
			self.orientation = 'W'
		elif(self.orientation == 'E'):
			self.orientation = 'N'
		elif(self.orientation == 'S'):
			self.orientation = 'E'
		elif(self.orientation == 'W'):
			self.orientation = 'S'
		self.StandardProcedure()

	def take_screenshot(self):
		widget = QWidget()
		p = QPixmap.grabWindow(widget.winID())
		p.save('Screenshot.jpg','jpg')

	def take_screenshot(self):
		widget = QWidget()
		p = QPixmap.grabWindow(widget.winID())
		p.save('Screenshot.jpg','jpg')
		

	def openImage(self):
		imagePath, _ = QFileDialog.getOpenFileName()
		pixmap = QPixmap(imagePath)
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()




	def EngageEditor(self):
		print('Engaging Editor...')
		ptID = int(Builder.NearestPlane(self.lat,self.lon))
		point = Point9.getPoint(ptID)
		self.lat = point.getLatDeg()
		self.lon = point.getLonDeg()
		pixmap = QPixmap('Assets8/Zoom5/Plane'+str(ptID)+'.png')
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()
#		self.lat = 

	def EngageEditorX(self):
		self.selected
		print('Engaging Editor...')
		ptID = int(Builder.NearestPlane(self.lat,self.lon))
		point = Point9.getPoint(ptID)
		self.lat = point.getLatDeg()
		self.lon = point.getLonDeg()
		pixmap = QPixmap('Assets8/Zoom5/Plane'+str(ptID)+'.png')
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()
#		self.lat = 

	def HighlightNN(self):
		pt = self.selected
		NN = Builder.FindNN(self.selected)
		
		mask1 = np.ones(shape=(900,1600),dtype='float16')
		mask2 = np.ones(shape=(900,1600),dtype='float16')
		mask3 = np.ones(shape=(900,1600),dtype='float16')
		mask4 = np.ones(shape=(900,1600),dtype='float16')
		NNIMG = np.zeros(shape=(len(NN),3),dtype='float16')
		ct = 0
		mask0 = np.ones(shape=(900,1600),dtype='float16')
		imgtmp = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		for i in NN:
			colsel = Builder.INT2Color(i)
			mask1 = np.isin(imgtmp[:,:,0],colsel[0])
			mask2 = np.isin(imgtmp[:,:,1],colsel[1])
			mask3 = np.isin(imgtmp[:,:,2],colsel[2])
		
	
#			mask1 = 1-(1-mask1)*(1-np.floor(1-(imgtmp[:,:,0]-colsel[0])*(imgtmp[:,:,0]-colsel[0])))
#			mask2 = 1-(1-mask2)*(1-np.floor(1-(imgtmp[:,:,0]-colsel[0])*(imgtmp[:,:,0]-colsel[0])))
#			mask3 = 1-(1-mask3)*(1-np.floor(1-(imgtmp[:,:,0]-colsel[0])*(imgtmp[:,:,0]-colsel[0])))
		masktot = (mask1 and mask2) and mask3
		imgnew = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')

		imgnew0 = imgnew[~masktot[:,:],:]
		imgnew0[22:40,20:20+17*10,:] = Builder.GenText(Dn[self.day]+'/'+Mn[self.month]+'/'+str(self.year))	
		imgnew0[52:70,20:20+17*10,:] = Builder.GenText(Hr[self.hour]+':00  GMT')	
		plt.imsave(fname='tmp.png',arr=imgnew0)
		pixmap = QPixmap('tmp.png')
		self.label.setPixmap(pixmap)
		self.resize(pixmap.size())
		self.adjustSize()

		
	def mouseMoveEvent(self, event):
		qtRectangle = self.frameGeometry()
		imgtmp = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
		self.ycursorpos = event.y()
		self.xcursorpos = event.x()
#		print(str(imgtmp[event.y(),event.x()]))	
#		coldat = imgtmp[event.y(),event.x()]*255
#		if(event.button() == Qt.MidButton):
#			print(str(round(coldat[0])*256*256+round(coldat[1])*256+round(coldat[2])))		
#		print(str(round(coldat[0])*256*256+round(coldat[1])*256+round(coldat[2])))
		
	def HighlightGeodesic(self):
		offsetY = -16
		coldat = self.img0[self.ycursorpos+offsetY,self.xcursorpos]*255
		IDnum = round(coldat[0])*256*256+round(coldat[1])*256+round(coldat[2])

		if(self.GeodesicSetReady and IDnum < len(Grid)):
			GPath = Builder.GeodesicPathID(self.TMPPoint, IDnum)


			print('Highlighting '+str(len(GPath))+' Points...')
			PRVID = np.zeros(shape=(1080,1920),dtype=int)
			PTSID = CA2I(self.img0)
			selecteds = []
			for IDnum0 in GPath:
				print(str(IDnum0))
				selecteds.append( Builder.SelectValue2(PTSID,IDnum0) )
			StandardProcedure()
			self.img = Builder.HighlightSelected2( self.img , Builder.CatSelection(selecteds) )

			plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
			pixmap = QPixmap('tmp.png')
			self.label.setPixmap(pixmap)
			self.resize(pixmap.size())
			self.adjustSize()
			self.GeodesicSetReady = False
		else:
			if(IDnum < len(Grid)):
				self.TMPPoint = IDnum
				self.lastypos = self.ycursorpos
				self.lastxpos = self.xcursorpos
				self.GeodesicSetReady = True
				print('Geodesic Point Set')

	def mousePressEvent(self,event):
		offsetY = -16
		self.ycursorpos = event.y()
		self.xcursorpos = event.x()
		if(self.PRSD == True):
			self.PRSD = False
		self.selectedProvince = 0
#		self.img = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
#		if(self.FDisp == True):
#			self.img = plt.imread('MANTIS/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')


		if(event.button() == Qt.LeftButton and (self.zoom <= 6)):
			qtRectangle = self.frameGeometry()
			imgtmp = self.img0[Yc0:Yc1,:]#plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
			if(self.angle != 0 and self.zoom >= 2):
				imgtmp0 = self.img0[:,:,:]
				imgtmp = qb.Rotate0(imgtmp0,self.angle,self.zoom)[Yc0:Yc1,:]
			coldat = imgtmp[event.y()+offsetY,event.x()]*255
			IDnum = round(coldat[0])*256*256+round(coldat[1])*256+round(coldat[2])
			if(self.HUDD == True):
				if(self.ycursorpos >= 850+Yc0 and self.ycursorpos  < 882+Yc0 and self.xcursorpos >= 53 and self.xcursorpos < 94):
					print('Zooming in mini map...')
					self.zoom0 = min(self.zoom0+1,1)
				elif(self.ycursorpos >= 850+Yc0 and self.ycursorpos  < 882+Yc0 and self.xcursorpos >= 119 and self.xcursorpos < 160):
					print('Zooming out mini map...')
					self.zoom0 = max(self.zoom0-1,0)

			if(self.selected != IDnum):
				self.selected = IDnum
				self.PTSD = True
			else:
				self.selected = len(Builder.Grid)
				self.PTSD = False
			if(self.selected <= len(Builder.Grid) ):
				self.StandardProcedure()
			else:
				self.WindowedProcedure()


			
#				self.HighlightProvince()
#				self.img0 = plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
#				plt.imsave(fname='tmp.png',arr=self.img[int(self.img.shape[0]/2)-540:int(self.img.shape[0]/2)+540,:])
#				
#				pixmap = QPixmap('tmp.png')
#				self.label.setPixmap(pixmap)
#				self.resize(pixmap.size())
#				self.adjustSize()




		if(event.button() == Qt.RightButton and (self.zoom <= 6)):
			if(self.PTSD == True):
				self.PTSD = False
			if(self.PRSD == False):
				self.PRSD = True
			imgtmp = self.img0[Yc0:Yc1,:]#plt.imread('IMG/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
			if(self.angle != 0 and self.zoom >= 2):
				imgtmp0 = self.img0[:,:,:]
				imgtmp = qb.Rotate0(imgtmp0,self.angle,self.zoom)[Yc0:Yc1,:]
			coldat = imgtmp[event.y()+offsetY,event.x()]*255

#			coldat = self.img0[event.y()+offsetY,event.x()]*255
			IDnum = round(coldat[0])*256*256+round(coldat[1])*256+round(coldat[2])


#			self.img = plt.imread('Assets/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
#			if(self.FDisp == True):
#				self.img = plt.imread('MANTIS/Zoom'+str(self.zoom)+'/Plane'+str(self.lat)+'x'+str(self.lon)+'.png')
			sameasbefore = 1*self.selectedProvince
#			print(sameasbefore)

			if(self.selected != IDnum and IDnum < len(Grid) ):
				self.selected = IDnum
			else:
				self.PRSD = False
			

			if(IDnum < len(Grid) ):
				self.selectedProvince = int(Grid[IDnum].getProvince())	
			print(self.selectedProvince)
			if(self.PRSD == True and self.selectedProvince == sameasbefore):
				self.PRSD = False
				self.selectedProvince = 0
				print('Deactivated')

			self.StandardProcedure()




						

	def hoverMoveEvent(self,event):
		point = event.pos().toPoint()
		print(point)
		
#	def paintEvent(self,event):
#		qp = QPainter(self)
#		qp.setPen(Qt.black)
#		if(event.button() == Qt.RightButton):
#			qp.drawPoint(event.x(),event.y(y))

#def addTextures(fname):
	

def main():
	app = QApplication(sys.argv)
	win = MainWindow()
#	ex = MouseTracker()
	win.show()
	return app.exec_()

#	A filter that brightens an RGB color vector, good for visualizing the sun-facing portions of the grid
def Brighten(color):
	R = color[0]
	G = color[1]
	B = color[2]
	color1 = np.ones(shape=4,dtype='float16')
	color1[0] = R
	color1[1] = G
	color1[2] = B
	color1[0] = min((1.1*R),((1+R)/2.0))
	color1[1] = min((1.1*G),((1+G)/2.0))
	color1[2] = min((1.1*B),((1+B)/2.0))
	return color1

def QBrighten(colors):
	return (colors + np.ones(shape=colors.shape,dtype=np.float16))/2.0
def CA2I(colors):
	sizes = colors.shape
	INTS = np.zeros(shape=(sizes[0],sizes[1]),dtype=int)
	INTS[:,:] = colors[:,:,0]*255*256*256+colors[:,:,1]*255*256+colors[:,:,2]*255
#	for i in range(sizes[0]):
#		for j in range(sizes[1]):
#			if(INTS[i,j] >= 1474562):
#				print('help at '+str(i)+' '+str(j))
	return INTS

def MapCB(x,y):
	return x*(1-y)+y


def LoadTextures(n):
	base = plt.imread('imgcirc/BigBlank'+str(n)+'.png')
	outarr = np.array(shape=(20,4,5),dtype='object')
	for i in range(20):
		for j in range(4):
			for k in range(5):
				tmpimg = base[n*i:n*(i+1),n*(4*k+j):n*(4*k+j+1),:].imsave('tmptexture.png')
				ourarr[i,j,k] = QImage('tmptexture.png')			
	return outarr
if __name__ == '__main__':
	sys.exit(main()) 