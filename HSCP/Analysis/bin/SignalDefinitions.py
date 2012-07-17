class SignalSample:
  def __init__(self,name,mass,crossSection,ptCut,iasCut,massCut):
    self.name = name
    self.mass = mass
    self.crossSection = crossSection
    self.ptCut = ptCut
    self.iasCut = iasCut
    self.massCut = massCut


# name, crossSectionPb,ptCut,iasCut,massCut
# Gluinos
Gluino300 = SignalSample('Gluino300',300,'65.8',60,'0.400',180)
Gluino400 = SignalSample('Gluino400',400,'11.2',80,'0.400',230)
Gluino500 = SignalSample('Gluino500',500,'2.54',50,'0.375',290)
Gluino600 = SignalSample('Gluino600',600,'0.693',50,'0.325',360)
Gluino700 = SignalSample('Gluino700',700,'0.214',50,'0.300',410)
Gluino800 = SignalSample('Gluino800',800,'0.0725',50,'0.275',460)
Gluino900 = SignalSample('Gluino900',900,'0.0262',100,'0.250',500)
Gluino1000 = SignalSample('Gluino1000',1000,'0.00987',120,'0.225',540)
Gluino1100 = SignalSample('Gluino1100',1100,'0.00386',120,'0.225',570)
Gluino1200 = SignalSample('Gluino1200',1200,'0.00154',155,'0.200',600)
# GluinoN
GluinoN300 = SignalSample('Gluino300N',300,'65.8',60,'0.400',180)
GluinoN400 = SignalSample('Gluino400N',400,'11.2',85,'0.400',230)
GluinoN500 = SignalSample('Gluino500N',500,'2.54',50,'0.375',290)
GluinoN600 = SignalSample('Gluino600N',600,'0.693',50,'0.350',330)
GluinoN700 = SignalSample('Gluino700N',700,'0.214',50,'0.325',270)
GluinoN800 = SignalSample('Gluino800N',800,'0.0725',50,'0.300',400)
GluinoN900 = SignalSample('Gluino900N',900,'0.0262',95,'0.300',380)
GluinoN1000 = SignalSample('Gluino1000N',1000,'0.00987',50,'0.300',410)
GluinoN1100 = SignalSample('Gluino1100N',1100,'0.00386',50,'0.275',460)
GluinoN1200 = SignalSample('Gluino1200N',1200,'0.00154',50,'0.300',400)
# Stop
Stop130 = SignalSample('Stop130',130,'120.0',70,'0.400',70)
Stop200 = SignalSample('Stop200',200,'13.0',60,'0.400',130)
Stop300 = SignalSample('Stop300',300,'1.31',70,'0.400',190)
Stop400 = SignalSample('Stop400',400,'0.218',65,'0.375',250)
Stop500 = SignalSample('Stop500',500,'0.0478',50,'0.350',310)
Stop600 = SignalSample('Stop600',600,'0.0125',50,'0.325',360)
Stop700 = SignalSample('Stop700',700,'0.00356',50,'0.300',410)
Stop800 = SignalSample('Stop800',800,'0.00114',50,'0.275',450)
# StopN
StopN130 = SignalSample('Stop130N',130,'120.0',70,'0.400',70)
StopN200 = SignalSample('Stop200N',200,'13.0',70,'0.400',120)
StopN300 = SignalSample('Stop300N',300,'1.31',70,'0.400',170)
StopN400 = SignalSample('Stop400N',400,'0.218',60,'0.400',240)
StopN500 = SignalSample('Stop500N',500,'0.0478',50,'0.375',280)
StopN600 = SignalSample('Stop600N',600,'0.0125',50,'0.350',320)
StopN700 = SignalSample('Stop700N',700,'0.00356',50,'0.325',350)
StopN800 = SignalSample('Stop800N',800,'0.00114',50,'0.325',370)
# GMSB Stau
GMStau100 = SignalSample('GMStau100',100,'1.3398',65,'0.400',20)
GMStau126 = SignalSample('GMStau126',126,'0.274591',70,'0.400',50)
GMStau156 = SignalSample('GMStau156',156,'0.0645953',85,'0.400',80)
GMStau200 = SignalSample('GMStau200',200,'0.0118093',70,'0.400',120)
GMStau247 = SignalSample('GMStau247',247,'0.00342512',60,'0.400',150)
GMStau308 = SignalSample('GMStau308',308,'0.00098447',70,'0.400',190)
GMStau370 = SignalSample('GMStau370',370,'0.000353388',60,'0.400',230)
GMStau432 = SignalSample('GMStau432',432,'0.000141817',55,'0.375',260)
GMStau494 = SignalSample('GMStau494',494,'0.00006177',65,'0.350',300)
# PP Stau
PPStau100 = SignalSample('PPStau100',100,'0.0382',70,'0.400',40)
PPStau126 = SignalSample('PPStau126',126,'0.0161',80,'0.400',60)
PPStau156 = SignalSample('PPStau156',156,'0.00704',85,'0.400',80)
PPStau200 = SignalSample('PPStau200',200,'0.00247',120,'0.400',120)
PPStau247 = SignalSample('PPStau247',247,'0.00101',110,'0.400',150)
PPStau308 = SignalSample('PPStau308',308,'0.000353',70,'0.400',190)
## DiChamps - TODO


modelList = []
modelList.extend([Gluino300,Gluino400,Gluino500,Gluino600,Gluino700,Gluino800,Gluino900,Gluino1000,Gluino1100,Gluino1200])
modelList.extend([Stop200,Stop300,Stop400,Stop500,Stop600,Stop700,Stop800]) # remove Stop130
modelList.extend([GMStau200,GMStau247,GMStau308,GMStau370,GMStau432,GMStau494]) # remove GMStau100, GMStau126, GMStau156
modelList.extend([GluinoN300,GluinoN400,GluinoN500,GluinoN600,GluinoN700,GluinoN800,GluinoN900,GluinoN1000,GluinoN1100,GluinoN1200])
modelList.extend([StopN200,StopN300,StopN400,StopN500,StopN600,StopN700,StopN800]) # remove StopN130

#modelList.extend([PPStau100,PPStau126,PPStau156,PPStau200,PPStau247,PPStau308])


SignalBasePath = "    'file:/hdfs/cms/user/cooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Nov11/"

