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

 
def read_parts_file(partsfile,boxdir,overlap):
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
    mics = {}       # {micrograph:[particle,particle,paticle]}



### get all the micrographs/particles that need to be analysed  assign particles to micrographs
    for i in data:
        mic = i[labels['_rlnMicrographName']]
        boxfilename ='{0}{1}'.format(boxdir,mic.split('/')[-1].replace('.mrc','.box'))      
        partid ='{0}{1}'.format(i[labels['_rlnCoordinateX']].split('.')[0],i[labels['_rlnCoordinateY']].split('.')[0])
        try:
            mics[mic].append(partid)
        except:
            mics[mic] = [partid]

### make the dictionary of filaments
    boxdic = {}     #{filename:{{micrograph:part:fil, part:fil}}    
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
                partid ='{0:0.0f}{1:0.0f}'.format(float(line.split()[0])+(0.5*boxsize),float(line.split()[1])+(0.5*boxsize))    
                boxdic[i][partid] = [fil]
                partcount+=1
        print(' '.join([i.split('/')[-1],'filament count:',str(fil),'particle count:',str(partcount)]))
    for i in boxdic:
        print (i,boxdic[i])
    ## assign priors
    ### put the particles in groups by filament
    fils = {}       #{filament:[[dataline],[dataline]]}
    for i in data:
        mic = i[labels['_rlnMicrographName']]
        partid ='{0:0.0f}{1:0.0f}'.format(float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']]))  
        filid = '{0}-{1}'.format(mic,boxdic[i[labels['_rlnMicrographName']]][partid][0])
        try:
            fils[filid].append(i)
        except:
            fils[filid]= [i]

    for fil in fils:
        plength = 0
        dat = fils[fil]
        ## do first segment
        partid ='{0:0.0f}{1:0.0f}'.format(float(dat[0][labels['_rlnCoordinateX']]),float(dat[0][labels['_rlnCoordinateY']]))
        mic = dat[0][labels['_rlnMicrographName']]
        partxy = np.array([float(dat[0][labels['_rlnCoordinateX']]),float(dat[0][labels['_rlnCoordinateY']])])
        npxy = np.array([float(dat[1][labels['_rlnCoordinateX']]),float(dat[1][labels['_rlnCoordinateY']])])
        if (partxy-npxy)[1] < 0:
            segang = -180-(-1*(angle_between(partxy-npxy,np.array([1.0,0.0]))))
        else:
            segang = -(angle_between(partxy-npxy,np.array([1.0,0.0])))
        boxdic[mic][partid].append(90.0)
        boxdic[mic][partid].append(0.500)
        boxdic[mic][partid].append(0.500)
        boxdic[mic][partid].append(plength)
        fsegang = (segang)
        boxdic[mic][partid].append(segang)
        print('\nfilament {0}'.format(fil))
        print('fil position: 0 {0}'.format(mic))
        print('part xy, nextpart xy',partxy,npxy)
        print('segment vector',partxy-npxy)
        print('segment angle',segang)
        print('final segment angle',fsegang)

            
        ## do 2nd through -1th segments
        n=1
        for i in dat[1:-1]:
            plength += overlap
            mic = i[labels['_rlnMicrographName']]
            partid ='{0:0.0f}{1:0.0f}'.format(float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']]))
            partxy = np.array([float(i[labels['_rlnCoordinateX']]),float(i[labels['_rlnCoordinateY']])])
            npxy = np.array([float(dat[n+1][labels['_rlnCoordinateX']]),float(dat[n+1][labels['_rlnCoordinateY']])])
            ppxy = np.array([float(dat[n-1][labels['_rlnCoordinateX']]),float(dat[n-1][labels['_rlnCoordinateY']])])
            angtonext = angle_between(partxy,npxy)
            if (partxy-npxy)[1] < 0:
                segang = -180-(-1*(angle_between(partxy-npxy,np.array([1.0,0.0]))))
                psegang = -180-(-1*(angle_between(ppxy-partxy,np.array([1.0,0.0]))))
            else:
                segang = -180-(angle_between(partxy-npxy,np.array([1.0,0.0])))
                psegang = -180-(angle_between(ppxy-partxy,np.array([1.0,0.0])))
            bend = abs(abs(segang)-abs(psegang))    
            boxdic[mic][partid].append(90.0)
            boxdic[mic][partid].append(0.500)
            boxdic[mic][partid].append(0.500)
            boxdic[mic][partid].append(plength)
            fsegang = np.mean([segang,psegang])
            boxdic[mic][partid].append(fsegang)
            print('\nfilament {0}'.format(fil))
            print('fil position: {0} {1}'.format(n,mic))
            print('part xy, nextpart xy',partxy,npxy)
            print('segment vector',partxy-npxy)
            print('segment angle',segang)
            print('prev seg angle',psegang)
            print('final segment angle',fsegang)
            print('bend',bend)
            n+=1

        ## do last segment
        plength += overlap
        mic = dat[-1][labels['_rlnMicrographName']]
        partid ='{0:0.0f}{1:0.0f}'.format(float(dat[-1][labels['_rlnCoordinateX']]),float(dat[-1][labels['_rlnCoordinateY']]))
        partxy = np.array([float(dat[-1][labels['_rlnCoordinateX']]),float(dat[-1][labels['_rlnCoordinateY']])])
        npxy = np.array([float(dat[-2][labels['_rlnCoordinateX']]),float(dat[-2][labels['_rlnCoordinateY']])])
        if (partxy-npxy)[1] < 0:
            segang = -(-1*angle_between(partxy-npxy,np.array([1.0,0.0])))
        else:
            segang = -(angle_between(partxy-npxy,np.array([1.0,0.0])))
        boxdic[mic][partid].append(90.0)
        boxdic[mic][partid].append(0.500)
        boxdic[mic][partid].append(0.500)
        boxdic[mic][partid].append(plength)
        fsegang = (segang)
        boxdic[mic][partid].append(segang)
        print('\nfilament {0}'.format(fil))
        print('fil position: {0} {1}'.format(n,mic))
        print('part xy, prevpart xy',partxy,npxy)
        print('segment vector',partxy-npxy)
        print('segment angle',-1*segang,'** last segment reversed')
        print('final segment angle',fsegang)
    return(labels,header,data,boxdic)




## program
errmsg = '\nUSAGE: rln3p1_crYOLO_add_filaments <particles file> <boxfiles directory> <overlap in pixels> <A/pix>'
try:
    boxdir = sys.argv[2]
except:
    sys.exit('\nERROR no particles directoty {0}'.format(errmsg))
if os.path.isdir(boxdir) == False:
    sys.exit('\nERROR: {0} is not a vaild boxfiles directory{1}'.format(sys.argv[2],errmsg))
if os.path.isfile(sys.argv[1]) == False:
    sys.exit('\nERROR reading particles starfile\n{0} is not a valid star file{1}'.format(sys.argv[1],errmsg))

try:
    overlap = float(sys.argv[3])
    apix = float(sys.argv[4])
    boxsize = float(sys.argv[5])
    print('overlap: {0} px @ {1} a/pix = {2} angstrom'.format(overlap,apix,overlap*apix))
except:
    sys.exit('ERROR: incorrect overlap (px), A/pix, or boxsize\n{0}'.format(errmsg))

labels,header,data,boxdic = read_parts_file(sys.argv[1],sys.argv[2],overlap*apix)

output = open('crYOLO_helix_parts.star','w')
for i in header:
    output.write(i)
output.write('_rlnHelicalTubeID #{0}\n'.format(int(header[-1].split('#')[-1])+1))
output.write('_rlnAngleTiltPrior #{0}\n'.format(int(header[-1].split('#')[-1])+2))
output.write('_rlnAngleRotFlipRatio #{0}\n'.format(int(header[-1].split('#')[-1])+3))
output.write('_rlnAnglePsiFlipRatio #{0}\n'.format(int(header[-1].split('#')[-1])+4))
output.write('_rlnHelicalTrackLengthAngst #{0}\n'.format(int(header[-1].split('#')[-1])+5))
output.write('_rlnAnglePsiPrior #{0}\n'.format(int(header[-1].split('#')[-1])+6))


for i in data:
    mic = i[labels['_rlnMicrographName']]
    partid ='{0}{1}'.format(i[labels['_rlnCoordinateX']].split('.')[0],i[labels['_rlnCoordinateY']].split('.')[0])  
    output.write('\n{0}  {1}'.format('  '.join(i),'  '.join([str(x) for x in boxdic[mic][partid]])))
