#!/usr/bin/python
'''Script to create a template rts.in file using get_obs_info.py'''
import jdcal as j
import os
import sys
import subprocess
## sys.path.append("~/Data/Tools")
import tool_lib as t

#input 1 - gps id

obs_id=sys.argv[1]



#subprocess.call("change_db.py curtin",shell=True)
#subprocess.call("make_metafiles.py -g %s -l --rts" %obs_id,shell=True)

obs_info = subprocess.check_output("get_observation_info.py -g %s" %obs_id, stderr=subprocess.STDOUT, shell=True)

for line in obs_info.split('\n'):
	if 'LST' in line:
		lst=float(line[4:11])
	if '(RA,Dec)' in line:
		line=line.split()
#		print line
		RA=float(line[2][1:-1])
		Dec=line[3][:-1]
	if 'Channels'in line:
		line=line.split(',')
		ch = line[0].split()
		base_ch = float(ch[1])
	if 'epoch' in line:
		info=line.split()
		date=info[1].split("/")
		year=date[0][1:]
		month=date[1]
		day=date[2][:-1]
		time=info[2]

		
#print year, month, day, time
		
HA = lst-RA
RA = RA/15.0
HA = HA/15.0

tmp_time = time.split(':')
tmp_min = float(tmp_time[1]*60 + tmp_time[2])/3600.

d1,d2 = j.gcal2jd(year, month, day)
time_frac=t.hour_to_dechr(time)/24.0
#time_frac = (float(tmp_time[0])+tmp_min)/24.
jdate = d1+d2+time_frac
print jdate


template_file = open('./input_cal_temp.in','r').read()
template = template_file.split('\n')
out_file = open('input_cal.in',"w+")
for line in template:
	line_out = line
	if 'ObservationPointCentreHA' in line:
		line_out = line.replace('hh',str(HA))
	if 'ObservationPointCentreDec' in line:
		line_out = line.replace('dd',Dec)
	if 'ObservationImageCentreRA' in line:
		line_out = line.replace('rr',str(RA))
	if 'ObservationImageCentreDec' in line:
		line_out = line.replace('dd',Dec)
	if 'ObservationFrequencyBase' in line:
		line_out = line.replace('ff', '%.3f' %((base_ch*1.28)-0.645))
	if 'ObservationTimeBase' in line:
		line_out = line.replace('tt',str(jdate))
	out_file.write(line_out+"\n")
out_file.close()

template_file = open('./input_img_temp.in','r').read()
template = template_file.split('\n')
out_file = open('input_img.in',"w+")
for line in template:
	line_out = line
	if 'ObservationPointCentreHA' in line:
		line_out = line.replace('hh',str(HA))
	if 'ObservationPointCentreDec' in line:
		line_out = line.replace('dd',Dec)
	if 'ObservationImageCentreRA' in line:
		line_out = line.replace('rr',str(RA))
	if 'ObservationImageCentreDec' in line:
		line_out = line.replace('dd',Dec)
	if 'ObservationFrequencyBase' in line:
		line_out = line.replace('ff', '%.3f' %((base_ch*1.28)-0.645))
	if 'ObservationTimeBase' in line:
		line_out = line.replace('tt',str(jdate))
	out_file.write(line_out+"\n")
out_file.close()
