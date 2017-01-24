#!/usr/bin/env python
#
# Author: Qiang Zhou (zhouqiang00@tsinghua.org.cn)
# School of Life Sciences, Tsinghua University, Beijing, China
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA

import re
class Metadata(object):
	def __init__(self,fn_meta):
		self.read(fn_meta)
	def __call__(self):
		self.run()
	def run(self):
		pass
	def read(self,fn_meta):
		pat_sector=re.compile(r"^(data_.*)$")
		pat_loop=re.compile(r"^loop_\s*$")
		pat_labels=re.compile(r"^(_rln\S+)\s+#*(\S+)\s*$")
		pat_datas=re.compile(r"^\s*(\S+)(\s+\S+)*\s*$")
		f=open(fn_meta,"r")
		for line in f:
			if line.startswith("#"):
				continue
			lspt=line.split()
			if len(lspt) >= 1:
				sectorm = pat_sector.match(line)
				loopm = pat_loop.match(line)
				labelsm = pat_labels.match(line)
				datasm = pat_datas.match(line)
				if sectorm:
					mysector=sectorm.groups()[0]
					setattr(self,mysector,{})
					continue
				elif loopm:
					continue
				elif labelsm:
					labelname,labelvalue=labelsm.groups()
					try :
						labelvalue = int(labelvalue)
					except:
						pass
					if not getattr(self,mysector).has_key("labels"):
						getattr(self,mysector)["labels"]={}
					getattr(self,mysector)["labels"][labelname]=labelvalue
					continue
				elif datasm:
					if not getattr(self,mysector).has_key("datas"):
						getattr(self,mysector)["datas"]=[]
					getattr(self,mysector)["datas"].append(lspt)
					continue
		f.close()
	def write(self,fn_meta):
		# Ok for data.star, but inferior for other type of metadata
		# The order of sectors can be controled for the other types
		f=open(fn_meta,"w")
		for sector,aux in self.__dict__.items():
			f.write("\n")
			f.write("%s\n"%(sector))
			f.write("\n")
			if aux.has_key("datas"):
				f.write("loop_\n")
				sorted_labels=sorted(aux["labels"].items(),key=lambda e:e[1],reverse=False)
				# Output _rlnLabel
				for label in sorted_labels :
				        f.write("%s\t#%d\n"%(label[0],label[1]))
			        # Output data records
				for record in aux["datas"]:
		        	        f.write(" ".join(record) + "\n")
			else:
				# only output _rlnLabel
				for label in aux["labels"].items():
					if isinstance(label[1],int):
					        f.write("%s\t%d\n"%(label[0],label[1]))
					else:
					        f.write("%s\t%s\n"%(label[0],label[1]))
		f.close()

					
# sampling is a subclass of Metadata
class sampling_meta(Metadata):
	def __init__(self,fn_sampling_meta):
		super(sampling_meta,self).__init__(fn_sampling_meta)
# sampling is a subclass of Metadata
class data_meta(Metadata):
	def __init__(self,fn_data_meta):
		super(data_meta,self).__init__(fn_data_meta)
# sampling is a subclass of Metadata
class optimiser_meta(Metadata):
	def __init__(self,fn_opti_meta):
		super(optimiser_meta,self).__init__(fn_opti_meta)


#my_sampling=sampling_meta("run1_it028_sampling.star")
#print my_sampling.data_sampling_general["labels"]
#print my_sampling.data_sampling_directions["labels"]
#print my_sampling.data_sampling_directions["datas"][:5]
#my_data=data_meta("run1_data.star")
#print my_data.data_images["datas"][:5]


					
		
