import numpy as np
import scipy.special
import math


def sigmaI(I,J,N,k):
	delta_I = (J/N)*(J/N)/(1+(J/N)*(J/N))
	return (2*I+1)*delta_I

def sigma_tot(E):
	k = math.sqrt(2*E)
	siT=0.
	l = 0
	lmax =int( 2*k+100 ) #doesn't matter
	JJ = scipy.special.sph_jn(2*lmax,k)
	NN = scipy.special.sph_yn(2*lmax,k)
	while l < lmax:
		siI = sigmaI(l,JJ[0][l],NN[0][l],k)
		siI1 = sigmaI(l,JJ[0][l+1],NN[0][l+1],k)
		siI2 = sigmaI(l,JJ[0][l+2],NN[0][l+2],k)
		siT += siI
		if (siI1+siI2+siI)/siT < 10e-4:
			lmax = l
		
		l += 1
	
	return ( 4*math.pi*siT/(k*k) , lmax )
	

logfile = open('test.log', 'w')                                           #(1)



ENER=0.
while ENER <= 10 :
	ENER += 0.0001		
	k = math.sqrt(2*ENER)
	logfile.write( str(ENER) )
	logfile.write(' ')
	logfile.write( str(sigma_tot(ENER)[0]) )
	logfile.write(' ')
	logfile.write( str(k) )
	logfile.write(' ')
	logfile.write( str(sigma_tot(ENER)[1]) )
	logfile.write('\n')
	

logfile.close()
