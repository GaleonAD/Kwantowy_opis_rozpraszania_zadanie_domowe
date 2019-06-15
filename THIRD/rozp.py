#!/usr/bin/env python3

import scipy.special as scsp
import math

# z lenistwa
def gl(E , U, l, r ):
	return 2*E - 2*U - l*(l+1)/(r*r)



def CHI_L( E, U, l, h, a):
	x_0 = h
	N_graniczne = int(a/h) 

	chi =  [ x_0**(l+1) , (x_0+h)**(l+1) ] 

	for n in range(2,N_graniczne+1):
		YN = (2.-5.*h*h*gl(E, U, l, x_0 + (n-1)*h)/6.)*chi[n-1]
		YN_1 = (1.+h*h*gl(E, U, l, x_0 + (n-2)*h)/12.)*chi[n-2]
		chi.append( (YN-YN_1)/(1.+h*h*gl(E, U, l, x_0 + n*h)/12.) )
		
	for n in range(N_graniczne+1,2*N_graniczne+1):
		YN = (2.-5.*h*h*gl(E, 0, l, x_0 + (n-1)*h)/6.)*chi[n-1]
		YN_1 = (1.+h*h*gl(E, 0, l, x_0 + (n-2)*h)/12.)*chi[n-2]
		chi.append( (YN-YN_1)/(1.+h*h*gl(E, 0, l, x_0 + n*h)/12.) )
	return chi



def SIGMA_L( l, E, U, h, a, Jl, Nl, Jlp, Nlp ):
	
	k = math.sqrt( 2*E )
	N_max = 2*int(a/h) - 5
	
	chi_l = CHI_L( E, U, l, h, a)

	y = chi_l[N_max] * h / (chi_l[N_max+1]-chi_l[N_max])
	tgdl = (Jl-k*y*Jlp)/(Nl-k*y*Nlp)

	s2 = tgdl*tgdl/(tgdl*tgdl+1)
	
	return	(2*l+1)*s2*4.*math.pi/(2*E)
	


def SIGMA_TOT( E, U, lmax, h, a ):
	
	k = math.sqrt(2*E)
	N_max = 2*int(a/h) - 5

	JJ = scsp.sph_jn((lmax+2),k*2*a)
	NN = scsp.sph_yn((lmax+2),k*2*a)
	
	sig_l = []
	sig_tot = 0.

	for l in range( 0 , lmax + 1):
		sig_l.append( SIGMA_L(l, E, U, h, a, 
			JJ[0][l], NN[0][l], JJ[1][l], NN[1][l] ))
		sig_tot += sig_l[l]
	
	return [sig_tot , sig_l]

