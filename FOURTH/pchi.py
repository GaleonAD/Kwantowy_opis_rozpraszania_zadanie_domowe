import rozp
import scipy.special as scsp
import math

#U=0
U = -100.75
a = 1
h = 1e-2
#l = 0
E = 20

JJ = scsp.sph_jn(4,math.sqrt(2*E)*a*2.)

for l in range (0, 3):
	logfile = open('chi'+str(l), 'w')  	
	
	chii = rozp.CHI_L( E,U,l,h,a )
	n = 0
	maax = 0
	for i in chii :
		if i > maax : 
			maax = i
	
	for i in chii :
		n += 1
		logfile.write( str(n*h) + " " + str(i/maax) + "\n")


logfile.close()
