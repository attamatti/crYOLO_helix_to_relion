#!/usr/bin/env python

import sys
import glob
import os

def read_parts_file(partsfile,boxdir):
	n=0
	pd = open(partsfile,'r').readlines()
	labels,header,data = {},[],[]
	indata = False
	inplabs = False
	for i in pd:	
		if indata == False:
			header.append(i)
			if i.split() == ['data_particles']:
				inplabs = True
			if inplabs == True and '_rln' in i:
				labels[i.split()[0]] = int(i.split()[1].replace('#',''))-1	
			if '_rln' in i and '_rln' not in pd[n+1] and inplabs == True:
				indata = True			
		elif indata == True and len(i.split()) == len(labels) and len(i.split()) != 0:
			data.append(i.split())
		n+=1
	mics = {} 		# {micrograph:[particle,particle,paticle]}

### get all the micrographs/particles that need to be analysed	
	for i in data:
		mic = i[labels['_rlnMicrographName']]
		boxfilename ='{0}{1}'.format(boxdir,mic.split('/')[-1].replace('.mrc','.box'))		
		partid ='{0}{1}'.format(i[labels['_rlnCoordinateX']].split('.')[0],i[labels['_rlnCoordinateY']].split('.')[0])
		try:
			mics[mic].append(partid)
		except:
			mics[mic] = [partid]

### make the dictionary of filaments
	boxdic = {}		#{filename:{{micrograph:part:fil, part:fil}}
	for i in mics:
		micbox = i.split('/')[-1].replace('.mrc','.box')
		boxfile = open('{0}{1}'.format(boxdir,micbox),'r').readlines()
		boxdic[i] = {}
		fil = 0
		partcount = 0
		for line in boxfile:
			if line.split()[0] == '#helix:':
				fil+=1 
				box = float(line.split(',')[-1])
			elif '#' not in line:
				partid ='{0:0.0f}{1:0.0f}'.format(float(line.split()[0])+(0.5*box),float(line.split()[1])+(0.5*box))	
				boxdic[i][partid] = fil
				partcount+=1
		print(' '.join([i.split('/')[-1],'filament count:',str(fil),'particle count:',str(partcount)]))
	return(labels,header,data,boxdic)	
## program
errmsg = '\nUSAGE: rln3p1_crYOLO_add_filaments <particles file> <boxfiles directory>'
try:
	boxdir = sys.argv[2]
except:
	sys.exit('\nERROR no particles directoty {0}'.format(errmsg))
if os.path.isdir(boxdir) == False:
	sys.exit('\nERROR: {0} is not a vaild boxfiles directory{1}'.format(sys.argv[2],errmsg))
if os.path.isfile(sys.argv[1]) == False:
	sys.exit('\nERROR reading particles starfile\n{0} is not a valid star file{1}'.format(sys.argv[1],errmsg))
labels,header,data,boxdic = read_parts_file(sys.argv[1],sys.argv[2])
output = open('crYOLO_helix_parts.star','w')
for i in header:
	output.write(i)
output.write('_rlnHelicalTubeID #{0}'.format(int(header[-1].split('#')[-1])+1))
for i in data:
	mic = i[labels['_rlnMicrographName']]
	partid ='{0}{1}'.format(i[labels['_rlnCoordinateX']].split('.')[0],i[labels['_rlnCoordinateY']].split('.')[0])	
	output.write('\n{0}  {1}'.format('  '.join(i),boxdic[mic][partid]))
