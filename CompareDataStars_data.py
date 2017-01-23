#!/home/zhouqiang/softwares/EMAN2/extlib/bin/python
import os,sys
from metadata import *
import random
import math
from EMAN2 import *
from optparse import OptionParser


def getEulersIndexes(labels):
	i_rot=labels["_rlnAngleRot"] - 1
	i_tilt=labels["_rlnAngleTilt"] - 1
	i_psi=labels["_rlnAnglePsi"] - 1
	return i_rot,i_tilt,i_psi
def getCoordOriginIndexes(labels):
	i_x=labels["_rlnCoordinateX"] - 1
	i_y=labels["_rlnCoordinateY"] - 1
	i_sx=labels["_rlnOriginX"] - 1
	i_sy=labels["_rlnOriginY"] - 1
	return i_x,i_y,i_sx,i_sy

def getEulersFromRelionDataRecord(record,i_rot,i_tilt,i_psi):
	return float(record[i_rot]),float(record[i_tilt]),float(record[i_psi])
def getCoordOriginFromRelionDataRecord(record,i_x,i_y,i_sx,i_sy):
	return float(record[i_x]),float(record[i_y]),float(record[i_sx]),float(record[i_sy])

def CLIP(v,max=0.999999,min=-0.999999):
	if v > max:
		v = max
	if v < min:
		v = min
	return v
def RegularAnglePhiPsi(angle):
	while True:
		if angle >= 360.:
			angle -= 360.
		if angle < 360.:
			break
	while True:
		if angle < 0.0:
			angle += 360.
		if angle >=0.0:
			break
	return angle

def calculateEulerDistance(e1,e2,e3,e4,e5,e6,type="spider",symmetry="c1"):
	if type == "spider":
		syms=Symmetries.get(symmetry)
		diff_angs_min = 3000.
		ts1=Transform({"type":"spider","phi":e1,"theta":e2,"psi":e3})
		ts2=Transform({"type":"spider","phi":e4,"theta":e5,"psi":e6})
		# The following line does not work for some particles, at least for c6 symmetry with unknown reason
		#for sym in syms.get_syms():
		matrix1=ts1.get_matrix()
		for n in range(syms.get_nsym()):
			if symmetry == "c1" or symmetry == "C1":
				matrix2=ts2.get_matrix()
			else:
				# Must be use this reduce method for c6 symmetry.
				ts2_reduced=syms.reduce(ts2,n)
				matrix2=ts2_reduced.get_matrix()
			v1_1=Vec3f(matrix1[0],matrix1[4],matrix1[8])
			v1_2=Vec3f(matrix1[1],matrix1[5],matrix1[9])
			v1_3=Vec3f(matrix1[2],matrix1[6],matrix1[10])
	
			v2_1=Vec3f(matrix2[0],matrix2[4],matrix2[8])
			v2_2=Vec3f(matrix2[1],matrix2[5],matrix2[9])
			v2_3=Vec3f(matrix2[2],matrix2[6],matrix2[10])

			diff_angs = degrees(acos(CLIP(v1_1.dot(v2_1))))
			diff_angs += degrees(acos(CLIP(v1_2.dot(v2_2))))
			diff_angs += degrees(acos(CLIP(v1_3.dot(v2_3))))

			if diff_angs < diff_angs_min:
				diff_angs_min = diff_angs
		return diff_angs_min  / 3
def calculateTranslationDistance(x1,y1,sx1,sy1,b1,x2,y2,sx2,sy2,b2):
	return	((x2 - (sx2 * b2) - (x1 - sx1 * b1) ) ** 2 + ( y2 - ( sy2 * b2) - (y1 - sy1 * b1) ) ** 2 ) ** 0.5

def calculateEulerTransformation(e1,e2,e3,e4,e5,e6,type="spider",symmetry="c1"):
	if type == "spider":
		syms=Symmetries.get(symmetry)
		ts1=Transform({"type":"spider","phi":e1,"theta":e2,"psi":e3})
		ts2=Transform({"type":"spider","phi":e4,"theta":e5,"psi":e6})
		# Find Euler transformation of volume 1 with Euler 1 relative to volume 2 with Euler 2
		# ts1: 30S, ts2: 50S. For ratched motion, 30S rotated against 50S by transform_1relativeto2.
		# And then 30S and 50S rotate together by ts2, getting the final position ts1 of 30S.
		# ts1 = ts2 * transform_1relativeto2 
		ts2.invert()
		transform_1relativeto2 = ts2 * ts1 
		Eulers = transform_1relativeto2.get_rotation("spider")
		return Eulers["phi"],Eulers["theta"],Eulers["psi"]
		
def main():
	prog_name=os.path.basename(sys.argv[0])
	usage="""
	# Compare Euler angle differences between different segmented parts of particles:
	# Corresponding segmented particles have same Image Name
	# data2.star is hierachied by Micrograph.
	{prog} < --input_data_star1 data1.star> < --input_data_star2 data2.star  > < -o output star file > < --sameImageName >

	# substitue imageName of data1.star with imageName of data2.star. Particles have same coordinates @ same micrograph:
	# data2.star is hierachied by Micrograph.
	{prog} < --input_data_star1 data1.star> < --input_data_star2 data2.star  > < -o output star file > <--changename> < --sameCoordinate>
	# Find rotating angle of volume of data1.star against to volume of data2.star:
	{prog} < --input_data_star1 data1.star> < --input_data_star2 data2.star  > < -o output star file >  --calculateEulerTransform --sameImageName
	#Transform euler angle of particles against to a rotated or flipped ref:
	{prog} < --input_data_star1 data1.star> < -o output star file >  < --TransformEuler > < [ --target_rot --target_tilt --target_psi  ] | [ --target_emanaz --target_emanalt --target_emanphi ] | [ --flipX | --flipY | --flipZ ]>
	""".format(prog=prog_name)
	optParser= OptionParser(usage)
	optParser.add_option("--input_data_star1",action="store",type="str",dest="input_data_star1",default="",help="Input relion data.star file #1. [default: %default]")
	optParser.add_option("--input_data_star2",action="store",type="str",dest="input_data_star2",default="",help="Input relion data.star file #2, usually a reference data.star with pre-determined values. [default: %default]")
	optParser.add_option("--sameImageName",action="store_true",dest="sameImageName",default=False,help="Two input data.star files have same image names for same particles. [default: %default]")
	optParser.add_option("--sameCoordinate",action="store_true",dest="sameCoordinate",default=False,help="Two input data.star files have same coordinates and same micrographs for same particles. [default: %default]")
	optParser.add_option("-o","--outputstar",action="store",type="str",dest="outputstar",default="",help="Output relion data.star file. [default: %default]")
	optParser.add_option("--calculateTransError",action="store_true",dest="calculateTransError",default=False,help="Calculate Translation differences between two input data.star. [default: %default]")
	optParser.add_option("--binningfactor1",action="store",type="float",dest="binningfactor1",default=1.0,help="Binning factor of input_data_star #1. [default: %default]")
	optParser.add_option("--binningfactor2",action="store",type="float",dest="binningfactor2",default=1.0,help="Binning factor of input_data_star #2. [default: %default]")
	optParser.add_option("--calculateEulerError",action="store_true",dest="calculateEulerError",default=False,help="Calculate Euler angles differences between two input data.star. [default: %default]")
	optParser.add_option("--flipX",action="store_true",dest="flipX",default=False,help="The references for the two input data.star files are flipped along X axis. [default: %default]")
	optParser.add_option("--flipY",action="store_true",dest="flipY",default=False,help="The references for the two input data.star files are flipped along Y axis. [default: %default]")
	optParser.add_option("--flipZ",action="store_true",dest="flipZ",default=False,help="The references for the two input data.star files are flipped along Z axis. [default: %default]")
	optParser.add_option("--noflip",action="store_true",dest="noflip",default=False,help="The references for the two input data.star files are not flipped each other. [default: %default]")
	optParser.add_option("--calculateEulerTransform",action="store_true",dest="calculateEulerTransform",default=False,help="Calculate Euler angle transformation of volume of input data.star 1 relative to volume of input data.star 2. [default: %default]")
	optParser.add_option("--TransformEuler",action="store_true",dest="TransformEuler",default=False,help="Transform Euler angle of particles against a rotated or flipped  reference. [default: %default]")
	optParser.add_option("--target_rot",action="store",type="float",dest="target_rot",default=0.0,help="Euler angle rot of the rotated reference relative to original ref. Overwrite eman angle. [default: %default]")
	optParser.add_option("--target_tilt",action="store",type="float",dest="target_tilt",default=0.0,help="Euler angle tilt of the rotated reference relative to original ref. Overwrite eman angle. [default: %default]")
	optParser.add_option("--target_psi",action="store",type="float",dest="target_psi",default=0.0,help="Euler angle psi of the rotated reference relative to original ref.. Overwrite eman angle. [default: %default]")
	optParser.add_option("--target_emanaz",action="store",type="float",dest="target_emanaz",default=0.0,help="Euler angle EMAN az of the rotated reference relative to original ref. [default: %default]")
	optParser.add_option("--target_emanalt",action="store",type="float",dest="target_emanalt",default=0.0,help="Euler angle EMAN alt of the rotated reference relative to original ref. [default: %default]")
	optParser.add_option("--target_emanphi",action="store",type="float",dest="target_emanphi",default=0.0,help="Euler angle EMAN phi of the rotated reference relative to original ref. [default: %default]")
	optParser.add_option("--TransformEulerbySym",action="store_true",dest="TransformEulerbySym",default=False,help="Transform Euler angle by symmetry. Output a series of files with suffix '_sym?'. [default: %default]")
	optParser.add_option("--RandomEulerbySym",action="store_true",dest="RandomEulerbySym",default=False,help="Random Euler angle by symmetry. (Random choose a symmetry-related Euler angle for particles)'. [default: %default]")
	optParser.add_option("--symmetry",action="store",type="str",dest="symmetry",default="C1",help="The symmetry used in the 3D reconstruction. [default: %default]")
	#optParser.add_option("--Euleroff_Phi",action="store",type="float",dest="Euleroff_Phi",default=0.0,help="PHI of Euler angle of reference of input data.star 1 relative to that of input 2. [default: %default]")
	optParser.add_option("--changename",action="store_true",dest="changename",default=False,help="Substitute ImageName of input data.star 1 with that of data.star 2 . [default: %default]")
	optParser.add_option("--filter_NumberatStack",action="store",dest="filter_NumberatStack",default="",help="File used for filter input data.star 1 with ImageName of a format of number@stack. [default: %default]")
	(options,args)=optParser.parse_args()
	
	data1=data_meta(options.input_data_star1)
	
	# find indexes
	i_MicrographName_1 = data1.data_["labels"]["_rlnMicrographName"] - 1 
	i_ImageName_1 = data1.data_["labels"]["_rlnImageName"] - 1
	try:
		i_CoordinateX_1 = data1.data_["labels"]["_rlnCoordinateX"] - 1 
		i_CoordinateY_1 = data1.data_["labels"]["_rlnCoordinateY"] - 1
	except:
		pass
	if options.input_data_star2:
		data2=data_meta(options.input_data_star2)
		# find indexes
		i_MicrographName_2 = data2.data_["labels"]["_rlnMicrographName"] - 1 
		i_ImageName_2 = data2.data_["labels"]["_rlnImageName"] - 1
		i_CoordinateX_2 = data2.data_["labels"]["_rlnCoordinateX"] - 1 
		i_CoordinateY_2 = data2.data_["labels"]["_rlnCoordinateY"] - 1
		# Hierachy second data.star by micrographs 
		allmics2={}
		for i in range(len(data2.data_["datas"])):
			mic = data2.data_["datas"][i][i_MicrographName_2]
			if allmics2.has_key(mic):
				allmics2[mic].append(data2.data_["datas"][i])
			else:
				allmics2[mic]=[]
				allmics2[mic].append(data2.data_["datas"][i])
	
	if options.changename:
		for record1 in data1.data_["datas"]:
			mic = record1[i_MicrographName_1]
			if allmics2.has_key(mic):
				for record2 in allmics2[mic]:
					if options.sameCoordinate:
						if record1[i_CoordinateX_1]  == record2[i_CoordinateX_2] and record1[i_CoordinateY_1]  == record2[i_CoordinateY_2]:
							record1[i_ImageName_1]  = record[i_ImageName_2]
							break
		if options.outputstar:
			data1.write(options.outputstar)
	

	if options. calculateTransError:
		i_x1,i_y1,i_sx1,i_sy1=getCoordOriginIndexes(data1.data_["labels"])
		i_x2,i_y2,i_sx2,i_sy2=getCoordOriginIndexes(data2.data_["labels"])
		for record1 in data1.data_["datas"]:
			mic = record1[i_MicrographName_1]
			if allmics2.has_key(mic):
				for record2 in allmics2[mic]:
					if options.sameCoordinate:
						pass
					if options.sameImageName:
						if record1[i_ImageName_1]  == record2[i_ImageName_2] :
							x1,y1,sx1,sy1=getCoordOriginFromRelionDataRecord(record1,i_x1,i_y1,i_sx1,i_sy1)
							x2,y2,sx2,sy2=getCoordOriginFromRelionDataRecord(record2,i_x2,i_y2,i_sx2,i_sy2)
							diff_translation = calculateTranslationDistance(x1,y1,sx1,sy1,options.binningfactor1,x2,y2,sx2,sy2,options.binningfactor2)
							if options.outputstar:
								record1.append(str(diff_translation))
							break
		if options.outputstar:
			data1.write(options.outputstar)

	if options. calculateEulerError:
		diff_eulers = 0.0
		count = 0
		i_rot1,i_tilt1,i_psi1=getEulersIndexes(data1.data_["labels"])
		i_rot2,i_tilt2,i_psi2=getEulersIndexes(data2.data_["labels"])
		outlines=[]
		for record1 in data1.data_["datas"]:
			mic = record1[i_MicrographName_1]
			if allmics2.has_key(mic):
				for record2 in allmics2[mic]:
					if options.sameCoordinate:
						pass
					if options.sameImageName:
						if record1[i_ImageName_1]  == record2[i_ImageName_2] :
							e1,e2,e3=getEulersFromRelionDataRecord(record1,i_rot1,i_tilt1,i_psi1)
							e4,e5,e6=getEulersFromRelionDataRecord(record2,i_rot2,i_tilt2,i_psi2)
							if options.flipX:
								e4 *= -1
								e5 = 180. - e5
							elif options.flipY:
								e4 = 180. - e4
								e5 = 180. - e5
							elif options.flipZ:
								e5 *= -1
							#e4 += options.Euleroff_Phi
							diff_euler = calculateEulerDistance(e1,e2,e3,e4,e5,e6,type="spider",symmetry=options.symmetry)
							if diff_euler <  1000. :
								count += 1
								diff_eulers += diff_euler
								if options.outputstar:
									if options.flipZ or options.flipX or options.flipY or options.noflip:
										record1.append(str(diff_euler))
									else :
										outlines.append("%s %f %f %f %f %f %f %f %f %f %f\n"%(record1[i_ImageName_1],e1,e2,e3,e4,e5,e6,diff_euler,RegularAnglePhiPsi(e4-e1),e5-e2,e6-e3))
								break
		if options.outputstar:
			if options.flipZ or options.flipX or options.flipY or options.noflip :
				data1.write(options.outputstar)
			else:
				f = open(options.outputstar,"w")
				for outline in outlines:
					f.write(outline)
				f.close()
				print "The match records between two star file is %d particles"%(count)
				print "The average Euler angle difference between two input star files is: %f degrees"%(diff_eulers / count)
				print "The fields in output file are:"
				print "ImageName,Phi1,Theta1,Psi1,Phi2,Theta2,Psi2,deltaEuler,deltaPhi,deltaTheta,deltaPsi"
	
	if options.TransformEulerbySym:
		i_rot1,i_tilt1,i_psi1=getEulersIndexes(data1.data_["labels"])
		allEulers=[]
		for record1 in data1.data_["datas"]:
			e1,e2,e3=getEulersFromRelionDataRecord(record1,i_rot1,i_tilt1,i_psi1)
			allEulers.append((e1,e2,e3))
		allTrans=[Transform({"type":"spider","phi":e[0],"theta":e[1],"psi":e[2]}) for e in allEulers]
		allEulers=None
		syms=Symmetries.get(options.symmetry)
		for index_sym in range(syms.get_nsym()):
			dirname=os.path.dirname(options.input_data_star1)
			basename=os.path.splitext(os.path.basename(options.input_data_star1))[0]
			symed_star = os.path.join(dirname,basename+"_sym%d.star"%(index_sym))
			my_current_sym = syms.get_sym(index_sym)
			n=0
			for record1 in data1.data_["datas"]:
				current_trans=allTrans[n]
				trans_by_sym = current_trans * my_current_sym
				Eulers_by_sym = trans_by_sym.get_rotation("spider")
				record1[i_rot1]=str(Eulers_by_sym["phi"])
				record1[i_tilt1]=str(Eulers_by_sym["theta"])
				record1[i_psi1]=str(Eulers_by_sym["psi"])
				n += 1
			data1.write(symed_star)
	if options.TransformEuler:
		TE=Transform()
		do_flip=False
		if (options.target_rot != 0 or options.target_tilt != 0 or options.target_psi != 0):
			TE.set_rotation({'type':spider,'phi':options.target_rot,'theta':options.target_tilt,'psi':options.target_psi })
		elif (options.target_emanaz != 0 or options.target_emanalt != 0 or options.target_emanphi != 0) :
			TE.set_rotation({'type':"eman",'az':options.target_emanaz,'alt':options.target_emanalt,'phi':options.target_emanphi })
		elif (options.flipX or options.flipY or options.flipZ):
			do_flip=True
		else:
			print "Nothing to do. Exit"
			exit()
		TE.invert()
		i_rot1,i_tilt1,i_psi1=getEulersIndexes(data1.data_["labels"])
		for record1 in data1.data_["datas"]:
			e1,e2,e3=getEulersFromRelionDataRecord(record1,i_rot1,i_tilt1,i_psi1)
			if not do_flip:
				myTrans=Transform({"type":"spider","phi":e1,"theta":e2,"psi":e3})
				trans_by_TE = myTrans * TE
				Eulers_by_TE = trans_by_TE.get_rotation("spider")
				record1[i_rot1]=str(Eulers_by_TE["phi"])
				record1[i_tilt1]=str(Eulers_by_TE["theta"])
				record1[i_psi1]=str(Eulers_by_TE["psi"])
			else:
				if options.flipX:
					e1 *= -1
					e2 = 180. - e2
				elif options.flipY:
					e1 = 180. - e1
					e2 = 180. - e2
				elif options.flipZ:
					e2 *= -1
				record1[i_rot1]=str(e1)
				record1[i_tilt1]=str(e2)
				record1[i_psi1]=str(e3)
		data1.write(options.outputstar)
	if options.RandomEulerbySym:
		from random import randrange
		i_rot1,i_tilt1,i_psi1=getEulersIndexes(data1.data_["labels"])
		syms=Symmetries.get(options.symmetry)
		nsyms =syms.get_nsym()
		for record1 in data1.data_["datas"]:
			e1,e2,e3 = getEulersFromRelionDataRecord(record1,i_rot1,i_tilt1,i_psi1)
			current_trans = Transform({"type":"spider","phi":e1,"theta":e2,"psi":e3})
			# Random choose a symmetry index
			my_chosen_sym = syms.get_sym(randrange(nsyms))
			trans_by_sym = current_trans * my_chosen_sym
			Eulers_by_sym = trans_by_sym.get_rotation("spider")
			record1[i_rot1]=str(Eulers_by_sym["phi"])
			record1[i_tilt1]=str(Eulers_by_sym["theta"])
			record1[i_psi1]=str(Eulers_by_sym["psi"])
		dirname = os.path.dirname(options.input_data_star1)
		basename = os.path.splitext(os.path.basename(options.input_data_star1))[0]
		randomsymed_star = os.path.join(dirname,basename+"_randomsym.star")
		data1.write(randomsymed_star)
	if options.calculateEulerTransform:
		i_rot1,i_tilt1,i_psi1=getEulersIndexes(data1.data_["labels"])
		i_rot2,i_tilt2,i_psi2=getEulersIndexes(data2.data_["labels"])
		if options.outputstar:
			f = open(options.outputstar,"w")
		for record1 in data1.data_["datas"]:
			mic = record1[i_MicrographName_1]
			if allmics2.has_key(mic):
				for record2 in allmics2[mic]:
					if options.sameCoordinate:
						pass
					if options.sameImageName:
						if record1[i_ImageName_1]  == record2[i_ImageName_2] :
							e1,e2,e3=getEulersFromRelionDataRecord(record1,i_rot1,i_tilt1,i_psi1)
							e4,e5,e6=getEulersFromRelionDataRecord(record2,i_rot2,i_tilt2,i_psi2)
							phi,theta,psi = calculateEulerTransformation(e1,e2,e3,e4,e5,e6,type="spider",symmetry=options.symmetry)
							f.write("%s %f %f %f %f %f %f %f\n"%(record1[i_ImageName_1],phi,theta,psi,RegularAnglePhiPsi(phi + psi),e4 - e1,e5 - e2 ,e6 - e3))
							break
		if options.outputstar:
			f.close()
		#print "The match records between two star file is %d particles"%(count)
		#print "The average Euler angle difference between two input star files is: %f degrees"%(diff_eulers / count)
		print "The fields in output file are:"
		print "ImageName, PhiET, ThetaET, PsiET, PhiET+PsiET, deltaPhi, deltaTheta, deltaPsi"
	if options.filter_NumberatStack :
		# hierachy filter_NumberatStack files by image stack name and number
		filter=open(options.filter_NumberatStack)
		filter_record={}
		for record in filter:
			splrec=record.split()
			number,img=splrec[0].split("@")
			if filter_record.has_key(img):
				filter_record[img].append(number)
			else:
				filter_record[img]=[number,]
		filter.close()
	
		dummy=[]
		for record1 in data1.data_["datas"]:
			number,img = record1[i_ImageName_1].split("@")
			if filter_record.has_key(img):
				if number in filter_record[img]:
					dummy.append(record1)

		data1.data_["datas"]=dummy
		data1.write(options.outputstar)
		
if __name__ == "__main__":
	main()
	
	
