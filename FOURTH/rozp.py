#!/usr/bin/env python3

import scipy.special as scsp
import math
from numerov import Numerov_solve_f

def CHI_L( E, U, l, h, a):
	
	def glU ( x ):
		if x <= a :
			return 2*E - 2*U - l*(l+1)/(x*x)
		else :
			return 2*E  - l*(l+1)/(x*x)

	def S ( x ):
		return 0.

	x_0 = h
	Nk = 2*int(a/h) + 10

	y0 = (x_0)**(l + 1)
	y1 = y0 + 2*h*(l+1)*(x_0**(l))

	chi = Numerov_solve_f(S, glU, x_0, Nk, y0, y1, h) 
	return chi

def SIGMA_L( l, E, U, h, a, Jl, Nl, Jlp, Nlp ):
	
	k = math.sqrt( 2*E )
	N_max = 2*int(a/h)  	

	
	chi_l = CHI_L( E, U, l, h, a)


	yp = (chi_l[N_max-2] - 8*chi_l[N_max-1] + 
				8*chi_l[N_max+1] - chi_l[N_max+2]) /(12.*h)
	y = chi_l[N_max] 
	
	tgdl = (yp*Jl-k*y*Jlp)/(yp*Nl-k*y*Nlp)

	sin2 = (tgdl**2)/(tgdl**2 + 1)

	return	(2*l+1.)*sin2*4.*math.pi/(2*E)
	


def SIGMA_TOT( E, U, lmax, h, a ):
	
	k = math.sqrt(2*E)

	JJ = scsp.sph_jn((lmax+2),k*2*a)
	NN = scsp.sph_yn((lmax+2),k*2*a)
	
	sig_l = []
	sig_tot = 0.

	for l in range( 0 , lmax + 1):
		sig_l.append( SIGMA_L(l, E, U, h, a, 
			JJ[0][l], NN[0][l], JJ[1][l], NN[1][l] ) )
		sig_tot += sig_l[l]
	
	return [sig_tot , sig_l]

