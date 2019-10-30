#!/usr/bin/env python

import sys
import glob
import os
import numpy as np

 
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

def get_group(partslist,labels):
	groupxs = [float(x[labels['_rlnCoordinateX']]) for x in partslist]
	groupys = [float(x[labels['_rlnCoordinateY']]) for x in partslist]
	coords = list(zip(groupxs,groupys))
	return(coords)

def get_alpha(data,labels,pos,n,rangelow,rangehigh):
	partid ='{0:0.0f}{1:0.0f}'.format(float(data[n][labels['_rlnCoordinateX']]),float(data[n][labels['_rlnCoordinateY']]))
	group = get_group(data[rangelow:rangehigh],labels)
	alpha = fit_line(group,pos)		
	return('{0:.4f}'.format(alpha),partid)
	
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
				boxdic[i][partid] = [fil]
				partcount+=1
		print(' '.join([i.split('/')[-1],'filament count:',str(fil),'particle count:',str(partcount)]))
	## assign priors
	### put the particles in groups by filament
	fils = {}		#{filament:[[dataline],[dataline]]}
	for i in data:
		partid ='{0:0.0f}{1:0.0f}'.format(float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']]))	
		try:
			#print(partid,boxdic[i[labels['_rlnMicrographName']]][partid][0])
			fils[boxdic[i[labels['_rlnMicrographName']]][partid][0]].append(i)
		except:
			#print(partid,boxdic[i[labels['_rlnMicrographName']]][partid][0])
			fils[boxdic[i[labels['_rlnMicrographName']]][partid][0]]= [i]
	## do each filament
	filsums = {}			##{filament:[value,value,value,value]}
	for fil in fils:
		dat = fils[fil]
	## do 0th segment:
		print('\nfil position: 0')
		print('filament {0}'.format(fil))
		alpha,partid = get_alpha(dat,labels,0,0,0,3)
		boxdic[mic][partid].append(90.0)
		filsums[fil]= [float(alpha)]

		## do 1st segment
		print('\nfil position: 1')
		print('filament {0}'.format(fil))
		alpha,partid = get_alpha(dat,labels,1,1,0,4)
		boxdic[mic][partid].append(90.0)
		filsums[fil].append(float(alpha))

		## do middle segments
		n=2
		for i in dat[2:-2]:
			print('\nfil position: {0}'.format(n))
			print('filament {0}'.format(fil))
			alpha,partid = get_alpha(dat,labels,2,n,n-2,n+3)
			boxdic[mic][partid].append(90.0)
			filsums[fil].append(float(alpha))

			n+=1
		## do -2th segment
		print('\nfil position: -2')
		print('filament {0}'.format(fil))
		alpha,partid = get_alpha(dat,labels,-2,-2,-4,-1)
		boxdic[mic][partid].append(90.0)
		filsums[fil].append(float(alpha))

		## do last segment
		print('\nfil position: -1')
		print('filament {0}'.format(fil))
		alpha,partid = get_alpha(dat,labels,-1,-1,-3,-1)
		boxdic[mic][partid].append(90.0)
		filsums[fil].append(float(alpha))

	
	return(labels,header,data,boxdic,filsums)

def chunk_list(seq,num):
    avg = len(seq)/float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return(out)


## program
errmsg = '\nUSAGE: rln3p1_crYOLO_add_filaments <particles file> <boxfiles directory> <number of segments for averaging>'
try:
	boxdir = sys.argv[2]
except:
	sys.exit('\nERROR no particles directoty {0}'.format(errmsg))
if os.path.isdir(boxdir) == False:
	sys.exit('\nERROR: {0} is not a vaild boxfiles directory{1}'.format(sys.argv[2],errmsg))
if os.path.isfile(sys.argv[1]) == False:
	sys.exit('\nERROR reading particles starfile\n{0} is not a valid star file{1}'.format(sys.argv[1],errmsg))
try:
	numsegs = int(sys.argv[3])
except:
	sys.exit('\nERROR number of segments per fibril not specified{0}'.format(errmsg))




labels,header,data,boxdic,filsums = read_parts_file(sys.argv[1],sys.argv[2])
filmeans = {}
print('\nfibril\tsubsection\tmean\t\tstd\toutliers')
for i in filsums:
	sscount = 0
	filmeans[i] = []
	subsets = chunk_list(filsums[i],int(len(filsums[i])/numsegs))
	for ss in subsets:
		calcgroup = []
		mean,std = np.mean(ss),np.std(ss)
		sscount+=1
		for val in ss:
			if np.abs(mean-val) < std:
				calcgroup.append(val)
			elif val == -0.0:
				calcgroup.append(val)				
		for val in ss:
			filmeans[i].append(np.mean(calcgroup))
		cmean,cstd = np.mean(calcgroup),np.std(calcgroup)
		print('{0:03.0f}\t{1:03.0f}\t\t{2:0.3f}\t\t{3:0.3f}\t{4}/{5}'.format(i,sscount,cmean,cstd,len(ss)-len(calcgroup),len(ss)))
		#print('\t\t\t\t\t\t      {0}\n'.format(calcgroup))
	print('\n')
output = open('crYOLO_helix_parts.star','w')
for i in header:
	output.write(i)
output.write('_rlnHelicalTubeID #{0}\n'.format(int(header[-1].split('#')[-1])+1))
output.write('_rlnAngleTiltPrior #{0}\n'.format(int(header[-1].split('#')[-1])+2))
output.write('_rlnAnglePsiPrior #{0}\n'.format(int(header[-1].split('#')[-1])+3))

filn = 0
dn = 0
for i in data:
	mic = i[labels['_rlnMicrographName']]
	partid ='{0}{1}'.format(i[labels['_rlnCoordinateX']].split('.')[0],i[labels['_rlnCoordinateY']].split('.')[0])	
	output.write('\n{0}  {1} {2}'.format('  '.join(i),'  '.join([str(x) for x in boxdic[mic][partid]]),filmeans[boxdic[mic][partid][0]][filn]))
	filn+=1
	try:
		if boxdic[mic][partid][0] != boxdic[mic]['{0}{1}'.format(data[dn+1][labels['_rlnCoordinateX']].split('.')[0],data[dn+1][labels['_rlnCoordinateY']].split('.')[0])][0]:
			filn=0
	except:
		pass
	dn+=1

