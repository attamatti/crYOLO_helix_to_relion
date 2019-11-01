#!/usr/bin/env python

import sys
import glob
import subprocess
import os

errmsg = '\nUSAGE: rln3p1_crYOLO_helix_pick.py <boxfiles directory>'
try:
	yoloboxes = glob.glob('{0}/*.box'.format(sys.argv[1]))
except:
	sys.exit('\nERROR no boxfiles dir specified{0}'.format(errmsg))
if len(yoloboxes) == 0:
	sys.exit('\nERROR: no boxfiles found in {0}{1}'.format(sys.argv[1],errmsg))

def read_yolobox(boxfile):
	if 'YOLObox' not in boxfile:
		boxdata = open(boxfile,'r').readlines()
		parts=[]
		filename = '{0}_YOLObox.box'.format(boxfile.split('.')[0])
		outbox = open(filename,'w')
		for i in boxdata:
			if i.split()[0] == '#helix:':
				boxsize = int(i.split(',')[-1])
			elif '#' not in i:
				outbox.write('{0:0.0f} {1:0.0f} {2} {2}\n'.format(float(i.split()[0])-(0.5*boxsize),float(i.split()[1])-(0.5*boxsize),boxsize))
		outbox.close()
		return(0,1)
	else:
		return(1,0)

skipped,success,fails = 0,0,[]
for i in yoloboxes:
	try:
		s,w = read_yolobox(i)
		skipped +=s
		success+=w
	except:
		fails.append(i)

subprocess.call(['touch','{}/coords_suffix_YOLObox.box'.format('/'.join(yoloboxes[0].split('/')[:-2]))])
print('\nWrote {0} _YOLObox file(s) with {1} errors'.format(success,len(fails)))
if len(fails) > 0:
	print('Errors in processing the following files:')
	for i in fails:
		print(i)
if skipped > 0:
	print('Skipped {0} file(s) because they were already written by this script'.format(skipped))
print('Wrote: {}/coords_suffix_YOLObox.box'.format('/'.join(yoloboxes[0].split('/')[:-2])))
