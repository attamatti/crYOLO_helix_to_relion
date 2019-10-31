#!/usr/bin/env python

import sys
import glob
import os
import numpy as np

def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return(np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))))

 
def fit_line(points,position):
	'''verticies of right triangle = xo,yo, x0,y+1, xn,y+1 '''
	xs = [x[0] for x in points]
	ys = [x[1] for x in points]
	poly = np.polyfit(xs,ys,1)
	ox,oy = points[position][0],points[position][1]
	
	if poly[0] >=-0.01:
		trivert2 = (((oy+1)-poly[1])/poly[0],oy+1)
		trivert1 = (ox,oy+1)
	else:
		trivert2 = (((oy-1)-poly[1])/poly[0],oy-1)
		trivert1 = (ox,oy-1)
	print('point: ({0},{1})'.format(ox,oy))
	print('position: {0}'.format(position))
	print('right triangle verticies ({0},{1}){2} {3}'.format(ox,oy,trivert1,trivert2))
	print('fit = y={0}x+{1}'.format(poly[0],poly[1]))
	opplength = np.abs(trivert1[0]-trivert2[0])+np.abs(trivert1[1]-trivert2[1])
	hyplength = np.abs(ox-trivert2[0])+np.abs(oy-trivert2[1])
	sin_alpha = opplength/hyplength
	alpha = np.degrees(np.arcsin(sin_alpha))  
	print('oppo length: {0}'.format(opplength))
	print('hypot length: {0}'.format(hyplength))
	
	if poly[0] >= 0:
		alpha = (90-alpha)
	elif poly[0] < 0:
		alpha = (alpha-90)
	print('alpha: {0:.4f}'.format(alpha))
	return(alpha)

	
def read_parts_file(partsfile,boxdir):
### read the star file
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



### get all the micrographs/particles that need to be analysed	assign particles to micrographs
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
				boxdic[i][partid] = [fil]
				partcount+=1
		print(' '.join([i.split('/')[-1],'filament count:',str(fil),'particle count:',str(partcount)]))

	## assign priors
	### put the particles in groups by filament
	fils = {}		#{filament:[[dataline],[dataline]]}
	for i in data:
		partid ='{0:0.0f}{1:0.0f}'.format(float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']]))	
		try:
			fils[boxdic[i[labels['_rlnMicrographName']]][partid][0]].append(i)
		except:
			fils[boxdic[i[labels['_rlnMicrographName']]][partid][0]]= [i]
	for fil in fils:
		dat = fils[fil]
		## do 0th through -4th segments
		n=1
		for i in dat[1:-1]:

			partid ='{0:0.0f}{1:0.0f}'.format(float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']]))
			partxy = np.array([float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']])])
			npxy = np.array([float(dat[n+1][labels['_rlnCoordinateX']]),float(dat[n+1][labels['_rlnCoordinateY']])])
			ppxy = np.array([float(dat[n-1][labels['_rlnCoordinateX']]),float(dat[n-1][labels['_rlnCoordinateY']])])
			angtonext = angle_between(partxy,npxy)
			if (partxy-npxy)[1] < 0:
				segang = -1*(angle_between(partxy-npxy,np.array([1.0,0.0])))
				psegang = -1*(angle_between(ppxy-partxy,np.array([1.0,0.0])))
			else:
				segang = angle_between(partxy-npxy,np.array([1.0,0.0]))
				psegang = angle_between(ppxy-partxy,np.array([1.0,0.0]))
			bend = abs(abs(segang)-abs(psegang))	
			boxdic[mic][partid].append(90.0)
			fsegang = np.mean([segang,psegang])
			boxdic[mic][partid].append(fsegang)
			print('\nfil position: {0}'.format(n))
			print('filament {0}'.format(fil))
			print('part xy, nextpart xy',partxy,npxy)
			print('segment vector',partxy-npxy)
			print('segment angle',segang)
			print('prev seg angle',psegang)
			print('final segment angle',fsegang)
			print('bend',bend)
			n+=1







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
output.write('_rlnHelicalTubeID #{0}\n'.format(int(header[-1].split('#')[-1])+1))
output.write('_rlnAngleTiltPrior #{0}\n'.format(int(header[-1].split('#')[-1])+2))
output.write('_rlnAnglePsiPrior #{0}\n'.format(int(header[-1].split('#')[-1])+3))

for i in data:
	mic = i[labels['_rlnMicrographName']]
	partid ='{0}{1}'.format(i[labels['_rlnCoordinateX']].split('.')[0],i[labels['_rlnCoordinateY']].split('.')[0])	
	output.write('\n{0}  {1}'.format('  '.join(i),'  '.join([str(x) for x in boxdic[mic][partid]])))


