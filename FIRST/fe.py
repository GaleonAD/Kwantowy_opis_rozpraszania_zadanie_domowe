import numpy as np
import scipy.special
import math


def sigmaIl0(I,J,N,k,kp,Jp,Np,JKP,JPKP):
	delta_I = math.atan((k/(kp))*math.tan(kp)-k)
#	A = (JKP*k)/(JPKP*kp)
#	delta_I = math.atan( (J-A*Jp)/(N-A*Np) )
	return (2*I+1)*math.sin(delta_I)*math.sin(delta_I)

def sigma_tot(E):
	k = math.sqrt(2*E)
	
	kp = math.sqrt(2*E+1.0)

	siT=0.
	l = 0
	lmax = 1
	JJ = scipy.special.sph_jn(2*lmax,k)
	NN = scipy.special.sph_yn(2*lmax,k)
	JJKP = scipy.special.sph_jn(2*lmax,kp)
	NNKP = scipy.special.sph_yn(2*lmax,kp)
	while l < lmax:
		siI = sigmaIl0(l,JJ[0][l],NN[0][l],k,kp,JJ[1][l],NN[1][l],JJKP[0][l],JJKP[1][l])
		siT += siI
		l += 1
	
	return ( 4*math.pi*siT/(k*k) , lmax )
	

logfile = open('test.log', 'w')                                           #(1)

V=-1.

ENER=0.
while ENER <= 10 :
	ENER += 0.001		
	k = math.sqrt(2*(ENER-V))
	logfile.write( str(ENER) )
	logfile.write(' ')
	logfile.write( str(sigma_tot(ENER)[0]) )
	logfile.write(' ')
	logfile.write( str(k) )
	logfile.write(' ')
	logfile.write( str(sigma_tot(ENER)[1]) )
	logfile.write('\n')
	

logfile.close()
