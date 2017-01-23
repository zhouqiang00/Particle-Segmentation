#!/home/zhouqiang/softwares/EMAN2/extlib/bin/python
from EMAN2 import *
import os,sys
import random
from optparse import OptionParser
import multiprocessing as mul

usage="""prog [options]
A programm to generate simulated micrographs.
%s --model1name ref.mrc  --NrMgr 2 --NrMgr_start 1 --particle_diameter
"""%(sys.argv[0])
optParser= OptionParser(usage)
optParser.add_option("--model1name",action="store",type="str",dest="model1name",default="",help="Name of volume 1 to generate particles. [default: %default]",metavar="test1.mrc")
optParser.add_option("--model2name",action="store",type="str",dest="model2name",default="",help="Name of volume 2 to generate particles with flexible component. [default: %default]",metavar="test2.mrc")
optParser.add_option("--NrMgr",action="store",type="int",dest="NrMgr",default=1,help="Number of micrographs to be generated. [default: %default]",metavar="Number of micrographs")
optParser.add_option("--NrMgr_start",action="store",type="int",dest="NrMgr_start",default=1,help="Starting No. of simulation micrographs. [default: %default]",metavar="Start number of micrographs")
optParser.add_option("--Mgrroot",action="store",type="str",dest="Mgrroot",default="Micrographs/simmgr",help="root name of micrographs including the folder. [default: %default]",metavar="Root name of micrographs")
optParser.add_option("--mgrsize",action="store",type="int",dest="mgrsize",default=4096,help="Size of micrographs. [default: %default]",metavar="Size of micrographs")
optParser.add_option("--NrPtcl",action="store",type="int",dest="NrPtcl",default=150,help="Number of particles in each micrograph. [default: %default]",metavar="150")
optParser.add_option("--particle_diameter",action="store",type="float",dest="particle_diameter",default=160.,help="Particles diameter in pixel, also used as the minimum distance between particles. [default: %default]",metavar="particle diameter")
optParser.add_option("--TrMean",action="store",type="float",dest="TrMean",default=0.0,help="Mean of Gaussian distribution of particles translation in pixels. [default: %default]",metavar="Translation mean")
optParser.add_option("--TrSigma",action="store",type="float",dest="TrSigma",default=2.0,help="Sigma of Gaussian distribution of particles translation in pixel. [default: %default]",metavar="Translation sigma")
optParser.add_option("--RotMean",action="store",type="float",dest="RotMean",default=0.0,help="Mean of Gaussian distribution of rotation of flexible part of particles in degrees. [default: %default]",metavar="Rotation mean")
optParser.add_option("--RotSigma",action="store",type="float",dest="RotSigma",default=10.0,help="Sigma of Gaussian distribution of rotation of flexible part of particles in degrees. [default: %default]",metavar="Rotation sigma")
optParser.add_option("--apix",action="store",type="float",dest="apix",default=1.32,help="Pixel size. [default: %default]",metavar="Pixel size")
optParser.add_option("--cs",action="store",type="float",dest="cs",default=2.7,help="Coefficient of spherical in mm. [default: %default]",metavar="Cs")
optParser.add_option("--HT",action="store",type="float",dest="HT",default=300.0,help="High tension in kV. [default: %default]",metavar="HT")
optParser.add_option("--ac",action="store",type="float",dest="ac",default=10.0,help="Amplitude constant of CTF in EMAN definition. [default: %default]",metavar="amplitude constant")
optParser.add_option("--bfactorMean",action="store",type="float",dest="bfactorMean",default=50.0,help="Mean of Gaussian distribution of B factor (positive). [default: %default]",metavar="Mean B factor")
optParser.add_option("--bfactorSigma",action="store",type="float",dest="bfactorSigma",default=2.0,help="Sigma of Gaussian distribution of B factor. [default: %default]",metavar="Sigma of B factor")
optParser.add_option("--mag",action="store",type="float",dest="mag",default=37879.,help="Magnification. [default: %default]",metavar="Magnification")
optParser.add_option("--Detector_pixel_size",action="store",type="float",dest="DP",default=5.,help="Pixel size in um of detector. [default: %default]",metavar="Detector pixel size")
optParser.add_option("--ctfmerit",action="store",type="float",dest="ctfmerit",default=0.3,help="Ctf merit of Ctffind3. [default: %default]",metavar="Ctf merit")
optParser.add_option("--dfmin",action="store",type="float",dest="dfmin",default=1.0,help="Minimum defocus in um (positive). [default: %default]",metavar="Minimum defocus")
optParser.add_option("--dfmax",action="store",type="float",dest="dfmax",default=3.0,help="Maximum defocus in um (positive). [default: %default]",metavar="Maximum defocus")
optParser.add_option("--dfdiffMean",action="store",type="float",dest="dfdiffMean",default=0.1,help="Mean of Gaussian distribution of defocus astigmatism in um. [default: %default]",metavar="Mean of amplitude of astigmatism in um")
optParser.add_option("--dfdiffSigma",action="store",type="float",dest="dfdiffSigma",default=0.005,help="Sigma of Gaussian distribution of defocus astigmatism in um. [default: %default]",metavar="Sigma of Amplitude of astigmatism in um")
optParser.add_option("--dfangMean",action="store",type="float",dest="dfangMean",default=50.0,help="Mean of Gaussian distribution of defocus astigmatism angle in degrees. [default: %default]",metavar="Mean of angle of astigmatism in degrees")
optParser.add_option("--dfangSigma",action="store",type="float",dest="dfangSigma",default=2.0,help="Sigma of Gaussian distribution of defocus astigmatism angle in degrees. [default: %default]",metavar="Sigma of angle of astigmatism in degrees")
optParser.add_option("--projector",action="store",type="str",dest="projector",default="standard",help="Projector. [default: %default]",metavar="standard or fourier_gridding")
optParser.add_option("--skip_translate",action="store",type="int",dest="skip_translate",default=0,help="Skip translate . [default: %default]",metavar="0 or 1")
optParser.add_option("--skip_applyCTFandNoise",action="store",type="int",dest="skip_applyCTFandNoise",default=0,help="Skip applying CTF and CTF-dependent noise to micrographs. [default: %default]",metavar="0 or 1")
optParser.add_option("--skip_applyCTFFreeNoise",action="store",type="int",dest="skip_applyCTFFreeNoise",default=0,help="Skip applying CTF-free noise to micrographs. [default: %default]",metavar="0 or 1")
optParser.add_option("--pad",action="store",type="float",dest="pad",default=1.0,help="Padding factor during applying CTF. [default: %default]",metavar="Padding factor")
optParser.add_option("--noisegaussiansigma",action="store",type="float",dest="noisegaussiansigma",default=10.24,help="Sigma of fourier gaussian noise. Higher value means sharper lowpass filter effect.  [default: %default]",metavar="Sigma of gaussian")
optParser.add_option("--water2Damp",action="store",type="str",dest="water2Damp",default="%s/water_amplitude_rotavg.mrc"%(os.path.dirname(sys.argv[0])),help="Micrograph of 2D rotationary average of water amplitude. [default: %default]",metavar="water_amplitude_rotavg")
optParser.add_option("--icenoiselevel",action="store",type="float",dest="icenoiselevel",default=1.0,help="Ice thickness in pixel . [default: %default]",metavar="CTF dependent noise level")
optParser.add_option("--noiselevel",action="store",type="float",dest="noiselevel",default=1.0,help="Noise level of micrographs. [default: %default]",metavar="Noise level")
optParser.add_option("--noisemode",action="store",type="str",dest="noisemode",default="k2",help="Noise mode used to add to micrographs. [default: %default]",metavar="white or ccd")
optParser.add_option("--multicpu",action="store",type="int",dest="multicpu",default=4,help="Number of cpu used for parallel processing. [default: %default]",metavar="2")
optParser.add_option("--verbose",action="store",type="int",dest="verbose",default=1,help="Verbose level. [default: %default]",metavar="Verbose level")

(options,args)=optParser.parse_args()

def getMinDist2(x0,y0,coords):
	"""
	Get the minimum distance between a coordinate and the other coordinates.

	"""
	if len(coords) > 0:
		distances=[ ( x0 - x ) * (x0 - x) + ( y0 -y ) * ( y0 - y ) for (x,y) in coords]
		return min(distances)
	else:
		return -1

def putParticleIntoMgr(mgr,ptcl,emanx,emany):
	"""
	Put a particle into a micrograph at coordinate (emanx,emany) in EMAN format.
	
	Maybe similar to EMData().insert_clip
	"""
	for j in range(ptcl.get_ysize()):
		for i in range(ptcl.get_xsize()):
			orivalue=mgr.get_value_at(i+emanx,j+emany)
			mgr.set_value_at(i+emanx,j+emany,orivalue + ptcl.get_value_at(i,j))


def applyCTF2Img(image,ctfmodel,pad=1.0):
	"""
	Apply CTF model to image.
	
	"""
	boxsize=image.get_xsize()
	if abs(pad - 1.) > 0.0001:
		b=image.process("xform.scale",{"scale":1.0,"clip":int(boxsize * pad)})
	else:
		b=image
	bfft=b.do_fft()
	bfftamp=bfft.copy()
	ctfamp=ctfmodel.CtfType.CTF_AMP
	ctfmodel.compute_2d_complex(bfftamp,ctfamp)
	bfft.mult(bfftamp)
	bctf=bfft.do_ift()
	if abs(pad - 1.) > 0.0001:
		bctf.process_inplace("xform.scale",{"scale":1.0,"clip":boxsize})
	return bctf

def getRandomAlign(Mean,Sigma):
	x=random.gauss(Mean,Sigma)
	y=random.gauss(Mean,Sigma)
	z=random.gauss(Mean,Sigma)
	return (x,y,z)
def getRandomEulerAngle(sym="c1"):
	if sym == "c1" or sym == "C1":
		phirange=359.99
	else:
		phirange=359.99
	phi=random.uniform(0,phirange)
	theta=random.uniform(0,179.99)
	psi=random.uniform(0,359.99)
	return (phi,theta,psi)
def getRandomCTF(bfactorMean=200.,bfactorSigma=10.,dfmin=1.0,dfmax=3.0,dfdiffMean=0.1,dfdiffSigma=0.005,dfangMean=50,dfangSigma=2.5):
	df=random.uniform(dfmin,dfmax)
	bfactor=random.gauss(bfactorMean,bfactorSigma)
	while abs(bfactor - bfactorMean) > 3 * bfactorSigma:
		bfactor=random.gauss(bfactorMean,bfactorSigma)
	dfdiff=random.gauss(dfdiffMean,dfdiffSigma)
	while abs(dfdiff - dfdiffMean) > 3 * dfdiffSigma:
		dfdiff=random.gauss(dfdiffMean,dfdiffSigma)
	dfang=random.gauss(dfangMean,dfangSigma)
	while abs(dfang - dfangMean) > 3 * dfangSigma:
		dfang=random.gauss(dfangMean,dfangSigma)
	return bfactor,df,dfdiff,dfang
def checkApix(map,apix):
	model_apix=map.get_attr("apix_x")
	if abs(model_apix - apix) > 0.01:
		return False
	else:
		return True
def checkmag(detect,apix,mag):
	if abs(detect * 10000. / apix - mag) > 1. :
		return False
	else:
		return True
def write_star_header(simstar,startype="mgr"):
	fsimstar=open(simstar,"w")
	fsimstar.write("\ndata_\n\nloop_\n")
	if startype == "mgr":
		fsimstar.write("_rlnMicrographName #1\n" + 
			"_rlnDefocusU #2\n" + 
			"_rlnDefocusV #3\n" + 
			"_rlnDefocusAngle #4\n" + 
			"_rlnVoltage #5\n" + 
			"_rlnSphericalAberration #6\n" + 
			"_rlnAmplitudeContrast #7\n" + 
			"_rlnMagnification #8\n" + 
			"_rlnDetectorPixelSize #9\n" + 
			"_rlnCtfFigureOfMerit #10\n")
	if startype == "ptcls":
		fsimstar.write("_rlnMicrographName #1\n" + 
			"_rlnCoordinateX #2\n" + 
			"_rlnCoordinateY #3\n" + 
			"_rlnImageName #4\n" + 
			"_rlnDefocusU #5\n" + 
			"_rlnDefocusV #6\n" + 
			"_rlnDefocusAngle #7\n" + 
			"_rlnVoltage #8\n" + 
			"_rlnSphericalAberration #9\n" + 
			"_rlnAmplitudeContrast #10\n" + 
			"_rlnMagnification #11\n" + 
			"_rlnDetectorPixelSize #12\n" + 
			"_rlnCtfFigureOfMerit #13\n")
	if startype == "ptcls_align":
		fsimstar.write("_rlnMicrographName #1\n" + 
			"_rlnCoordinateX #2\n" + 
			"_rlnCoordinateY #3\n" + 
			"_rlnImageName #4\n" + 
			"_rlnOriginX #5\n" +  
			"_rlnOriginY #6\n" + 
			"_rlnAngleRot #7\n" + 
			"_rlnAngleTilt #8\n" +
			"_rlnAnglePsi #9\n" + 
			"_rlnDefocusU #10\n" + 
			"_rlnDefocusV #11\n" + 
			"_rlnDefocusAngle #12\n" + 
			"_rlnVoltage #13\n" + 
			"_rlnSphericalAberration #14\n" + 
			"_rlnAmplitudeContrast #15\n" + 
			"_rlnMagnification #16\n" + 
			"_rlnDetectorPixelSize #17\n" + 
			"_rlnCtfFigureOfMerit #18\n")
	fsimstar.close()
def append_star(simstar,lines):
	fsimstar=open(simstar,"a")
	for line in lines:
		fsimstar.write(line)
	fsimstar.close()	

def sim2DWaterLikeNoise(fn_avgamp2D,target_apix=1.32,target_bs=4096):
        ice_amp_rot_i=EMData(fn_avgamp2D)
        myboxsize=ice_amp_rot_i.get_xsize()
        myapix=ice_amp_rot_i['apix_x']
	if abs(myapix - target_apix) > 0.001 :
		ice_amp_rot_i.process_inplace("xform.scale",{"scale":( myboxsize * myapix)/(target_bs * target_apix ),"clip":target_bs})
		myboxsize = target_bs
		myapix = target_apix
        ice_amp_rot_i.clip_inplace(Region(myboxsize / 2 - 1,0,myboxsize / 2 + 1,myboxsize))
        ice_amp_rot_i = ice_amp_rot_i.real2complex()
        wl_noise = test_image(1,size=(myboxsize,myboxsize))
        wl_noise.do_fft_inplace()
        wl_noise.fft_shuffle()
        wl_noise *= ice_amp_rot_i
        wl_noise.fft_shuffle()
        wl_noise = wl_noise.do_ift()
        wl_noise.process_inplace("normalize")
        wl_noise.update()
        return wl_noise


# Read in options.model 3D images
model1=EMData()
model1.read_image(options.model1name)
if not checkApix(model1,options.apix):
	print "Different apix between model1 and options"
	sys.exit(-1)
boxsize=model1.get_xsize()
if options.model2name:
	model2=EMData()
	model2.read_image(options.model2name)
	boxsize2=model2.get_xsize()
	if boxsize2 != boxsize:
		print "Error, model 2 have different size with model 1."
		exit(-1)
	if not checkApix(model2,options.apix):
		print "Different apix between model2 and options"
		sys.exit(-1)
# Check magnification, apix and detector size
if not checkmag(options.DP,options.apix,options.mag):
	options.mag=options.DP * 10000. / options.apix
	print "Warning: different magnification. change it to %f"%(options.mag)
# Check folder for micrographs
if not os.path.exists(os.path.split(options.Mgrroot)[0]):
	print "Warning, creat a folder %s"%(os.path.split(options.Mgrroot)[0])
	os.mkdir(os.path.split(options.Mgrroot)[0])
else:
	if not os.path.isdir(os.path.split(options.Mgrroot)[0]):
		print "%s exists and is not a folder. Exit."%(os.path.split(options.Mgrroot)[0])
		exit(-1)


# Read in options.model 3D images
class simMicrograph(object):
	def __init__(self,options,Nr_start):
		super(simMicrograph,self).__init__()
		self.options = options
		self.Nr_start = Nr_start
	def initialize(self):
		self.model1=EMData()
		self.model1.read_image(self.options.model1name)
		self.boxsize=self.model1.get_xsize()
		if self.options.model2name:
			self.model2=EMData()
			self.model2.read_image(self.options.model2name)
			self.boxsize2=self.model2.get_xsize()
		self.Mgrnum=self.options.Mgrroot + "_%04d"%(self.Nr_start)
		self.myMgr=self.Mgrnum + ".mrc"
		self.myMicrograph=EMData(self.options.mgrsize,self.options.mgrsize)
		self.get_Ctf()
		self.get_BoxCoords()
		self.get_alignment_1()
		if self.options.model2name:
			self.get_alignment_2()
	def __call__(self):
		self.run()
	def run(self):
		self.initialize()
		self.insertPtcls2Mgr()
		if not self.options.skip_applyCTFandNoise:
			self.applyWLNoise2Mgr()
		self.applyCTFtomgr()
		if not self.options.skip_applyCTFFreeNoise:
			self.applyCTFFreeNoise()
		self.output()
	def output(self):
		self.myMicrograph.write_image(self.myMgr)
		self.write_Mgrstar()
		self.write_PtclsStar()
		self.write_Ptclstar_1()
		if self.options.model2name:
			self.write_Ptclstar_2()
		self.write_Ctflog()
		self.write_Boxes()
	def get_Ctf(self):
		self.ctfmgr=EMAN2Ctf()
		bfactor,df,dfdiff,dfang= getRandomCTF(self.options.bfactorMean,self.options.bfactorSigma,self.options.dfmin,self.options.dfmax,self.options.dfdiffMean,self.options.dfdiffSigma,self.options.dfangMean,self.options.dfangSigma)
		self.ctfmgr.from_dict({"ampcont":self.options.ac,"apix":self.options.apix,"bfactor":bfactor,"cs":self.options.cs,"defocus":df,"dfang":dfang,"dfdiff":dfdiff,"voltage":self.options.HT})
		self.ctf_etc_line="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%((df + dfdiff/2) * 10000.,(df - dfdiff/2) * 10000.,dfang,self.options.HT,self.options.cs,self.options.ac/100.,self.options.mag,self.options.DP,self.options.ctfmerit)
	def get_BoxCoords(self):
		self.coords=[]
		for j in range(self.options.NrPtcl):
			if self.options.verbose > 0 :
				print "Try to get random coordinate for particle %06d@%s"%(j,self.myMgr)
			max_edge=self.options.mgrsize - self.boxsize
			x=random.randint(0,max_edge)
			y=random.randint(0,max_edge)
			min_distance2 = getMinDist2(x,y,self.coords)
			while min_distance2 >= 0. and min_distance2 < self.options.particle_diameter * self.options.particle_diameter:
				if self.options.verbose > 0 :
					print "Retry to get random coordinate for particle %06d@%s"%(j,self.myMgr)
				x=random.randint(0,max_edge)
				y=random.randint(0,max_edge)
				min_distance2 = getMinDist2(x,y,self.coords)
			self.coords.append((x,y))
	def get_alignment_1(self):
		self.alignment_1=[]
		for j in range(self.options.NrPtcl):
			tx1,ty1,tz1=getRandomAlign(self.options.TrMean,self.options.TrSigma)
			if self.options.skip_translate :
				tx1, ty1 = 0., 0. 
			phi1,theta1,psi1=getRandomEulerAngle("c1")
			self.alignment_1.append((tx1,ty1,tz1,phi1,theta1,psi1))
	def get_alignment_2(self):
		self.alignment_2=[]
		for j in range(self.options.NrPtcl):
			tx2,ty2,tz2=getRandomAlign(self.options.TrMean,self.options.TrSigma)
			if self.options.skip_translate :
				tx2, ty2 = 0., 0. 
			phi1,theta1,psi1=self.alignment_1[j][3:6]
			dphi,dtheta,dpsi=getRandomAlign(self.options.RotMean,self.options.RotSigma)
			phi2,theta2,psi2 = phi1 + dphi, theta1 + dtheta, psi1 + dpsi 
			self.alignment_2.append((tx2,ty2,tz2,phi2,theta2,psi2))
	def write_Mgrstar(self):
		myMgrstar=self.Mgrnum + ".star"
		write_star_header(myMgrstar,startype="mgr")
		mgrline="%s\t%s\n"%(self.myMgr,self.ctf_etc_line)
		append_star(myMgrstar,[mgrline])
	
	def write_PtclsStar(self):
		myParticles=[]
		myParticlestar=self.Mgrnum + "_ptcls"  + ".star"
		write_star_header(myParticlestar,startype="ptcls")
		for j in range(self.options.NrPtcl):
			ptclsline  ="%s\t%s\t%s\t%06d@%s\t%s\n"%(
			self.myMgr,self.coords[j][0],self.coords[j][1],j+1, "Particles/" + self.Mgrnum + "_particles.mrcs",self.ctf_etc_line)
			myParticles.append(ptclsline)
		append_star(myParticlestar,myParticles)
	def write_Ptclstar_1(self):
		myParticlestar_1=self.Mgrnum + "_ptcls1"  + ".star"
		write_star_header(myParticlestar_1,startype="ptcls_align")
		myParticles_1=[]
		for j in range(self.options.NrPtcl):
			ptclsline_1="%s\t%s\t%s\t%06d@%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(
				self.myMgr,self.coords[j][0],self.coords[j][1],j+1, "Particles/" + 
				self.Mgrnum + "_particles.mrcs",-self.alignment_1[j][0],
				-self.alignment_1[j][1],self.alignment_1[j][3],
				self.alignment_1[j][4],self.alignment_1[j][5],self.ctf_etc_line)
			myParticles_1.append(ptclsline_1)
		append_star(myParticlestar_1,myParticles_1)
	def write_Ptclstar_2(self):
		myParticlestar_2=self.Mgrnum + "_ptcls2"  + ".star"
		write_star_header(myParticlestar_2,startype="ptcls_align")
		myParticles_2=[]
		for j in range(self.options.NrPtcl):
			ptclsline_2="%s\t%s\t%s\t%06d@%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(
				self.myMgr,self.coords[j][0],self.coords[j][1],j+1, "Particles/" + 
				self.Mgrnum + "_particles.mrcs",-self.alignment_2[j][0],
				-self.alignment_2[j][1],self.alignment_2[j][3],
				self.alignment_2[j][4],self.alignment_2[j][5],self.ctf_etc_line)
			myParticles_2.append(ptclsline_2)
		append_star(myParticlestar_2,myParticles_2)
	def write_Ctflog(self):	
		myMgrctflog=self.Mgrnum + "_ctffind3.log"
		fctflog=open(myMgrctflog,"w")
		fctflog.write("CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]\n")
		fctflog.write("%f %f %f %f %f\n"%(self.options.cs,self.options.HT,self.options.ac / 100.,self.options.mag,self.options.DP))
		fctflog.write("%f %f %f %f Final Values\n"%((self.ctfmgr.defocus + self.ctfmgr.dfdiff / 2) * 10000., (self.ctfmgr.defocus - self.ctfmgr.dfdiff / 2) * 10000., self.ctfmgr.dfang, self.options.ctfmerit))
		fctflog.close()
	def write_Boxes(self):
		# Write out box coordinates
		myMgrbox=self.Mgrnum + ".box"
		fbox=open(myMgrbox,"w")
		for x,y in self.coords:
			fbox.write("%d\t%d\t%d\t%d\n"%(x,y,self.boxsize,self.boxsize))
		fbox.close()

	def insertPtcls2Mgr(self):
		for j in range(self.options.NrPtcl):
			tx1,ty1,tz1,phi1,theta1,psi1=self.alignment_1[j]
			t1=Transform({"type":"spider","psi":psi1,"theta":theta1,"phi":phi1})
			proj1=self.model1.project(self.options.projector,t1)
			if not self.options.skip_translate:
				proj1.translate(int(tx1),int(ty1),0)
			# For model 2 if specified
			if self.options.model2name :
				tx2,ty2,tz2,phi2,theta2,psi2=self.alignment_2[j]
				t2=Transform({"type":"spider","psi":psi2,"theta":theta2,"phi":phi2})
				proj2=self.model2.project(self.options.projector,t2)
				if not self.options.skip_translate:
					proj2.translate(int(tx2),int(ty2),0)
				proj1 += proj2
			# Insert projection into micrograph
			putParticleIntoMgr(self.myMicrograph,proj1,self.coords[j][0],self.coords[j][1])

	def applyWLNoise2Mgr(self):
		icenoise = sim2DWaterLikeNoise(self.options.water2Damp,self.options.apix,self.options.mgrsize)
		# 0.024: coefficient between ice thickness in pixel and  ice sigma 
		# 0.228: coefficient between ice thickness in pixel and  ice mean
		sigma = self.options.icenoiselevel * 0.024
		mean = self.options.icenoiselevel * 0.228
		icenoise = icenoise * sigma
		icenoise = icenoise + mean
		# 4. :   coefficient between protein projection value and protein thickness in pixel
		thick_frac = self.myMicrograph * 4. / self.options.icenoiselevel
		mgr_one=EMData(self.options.mgrsize,self.options.mgrsize)
		mgr_one.to_one()
		water_diff_frac = mgr_one - thick_frac 
		#water_diff_frac.process_inplace("threshold.belowtozero",{"minval":0.0})
		water_diff = icenoise * water_diff_frac 
		self.myMicrograph += water_diff
		self.myMicrograph.process_inplace("normalize")
	def applyCTFtomgr(self):
		# Apply CTF to micrograph
		self.myMicrograph=applyCTF2Img(self.myMicrograph,self.ctfmgr,self.options.pad)
	def applyCTFFreeNoise(self):
		if self.options.noisemode == "DDD" or self.options.noisemode == "ddd" or self.options.noisemode == "Ddd" or self.options.noisemode == "K2" or self.options.noisemode == "k2":
			noise = test_image(1,size=(self.options.mgrsize,self.options.mgrsize))
			if abs(self.options.noiselevel - 1.)  > 0.001 :
				noise *= self.options.noiselevel
			self.myMicrograph += noise
		self.myMicrograph.process_inplace("normalize")
		self.myMicrograph.set_ctf(self.ctfmgr)

def runSimMgr(s):
	s.run()
myPool=mul.Pool(options.multicpu)
for i in range(options.NrMgr_start,options.NrMgr_start + options.NrMgr):
	simmgr=simMicrograph(options,i)
	myPool.apply_async(runSimMgr,(simmgr,))
myPool.close()
myPool.join()


