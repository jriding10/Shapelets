import numpy as np

def hour_to_deg(time):      #converts hh:mm:ss.ss in to degrees, must input as a string
	negtest=time[0]
	time=time.split(':')
	degr=float(time[0])*15.0
	mins=float(time[1])*(15.0/60.0)
	secs=float(time[2])*(15.0/(3600.0))
	if negtest=='-':
		deg=degr-mins-secs
	if negtest!='-':
		deg=degr+mins+secs
	return deg

def dec_to_deg(time):    #converts dd:mm:ss.ss in to degrees, must input as a string
	negtest=time[0]
	time=time.split(':')
	degr=float(time[0])
	mins=float(time[1])*(1.0/60.0)
	secs=float(time[2])*(1.0/(3600.0))
	if negtest=='-':
		deg=degr-mins-secs
	if negtest!='-':
		deg=degr+mins+secs
	return deg

param_infos = open('jennys_params.txt').read().split('--------------------------------------------------------------------------------')
del param_infos[0]

##Now each entry is [name,gauss_info,pos_info]
param_infos = [param_info.split('\n')[1:-1] for param_info in param_infos]

r_to_d = 180.0/np.pi

file_names = ['3C33_model.txt','3C353_model.txt','3C353_model.txt','3C444_model.txt','HerA_model.txt','HydA_model_new.txt','PKS0408-65_model.txt','PKS2153-69_model.txt']

for param_info,file_name in zip(param_infos,file_names):
	#print param_info
	name,gauss_info,pos_info = param_info
	
	#print gauss_info
	major,minor,PA = gauss_info.split()
	##Convert major/minor from radians to arcmins
	major = (float(major)*r_to_d)/60.0
	minor = (float(minor)*r_to_d)/60.0
	##Convert PA from rads to deg
	PA = float(PA)*r_to_d
	
	##hrs,deg,deg-amins,deg
	ra,ra_off,dec,dec_off = pos_info.split()
	
	##Covert to deg, then in to hours decimal as RTS like it like that
	RA = ( hour_to_deg(ra) + float(ra_off) ) / 15.0
	
	DEC = dec_to_deg(dec) + float(dec_off)
	
	out_file = open("%s_RTS.txt" %name,'w+')
	out_file.write("SOURCE %s %.5f %.5f\n" %(name,RA,DEC))
	print "SOURCE %s %.5f %.5f\n" %(name,RA*15.0,DEC)
	
	out_file.write("FREQ 150e+6 flux 0 0 0\n")
	out_file.write("SHAPELET %.5f %.5f %.5f" %(PA,major,minor))
	
	
	shape_coeffs = open(file_name).read().split('\n')
	del shape_coeffs[-1]
	
	for shape_coeff in shape_coeffs:
		c1,c2,c3 = shape_coeff.split()
		out_file.write('\nCOEFF %s %s %s' %(c1,c2,c3))
	
	out_file.close()

