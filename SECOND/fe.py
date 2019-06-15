import numpy as np
import scipy.special as scsp
import math


# z lenistwa
def gl(E , U, l, r ):
	return 2*E - 2*U - l*(l+1)/(r*r)


# liczenie chi_l(x)
def CHI( E, l_max ) :
	h = 10e-5
	
	x_0 = h
	U = 0.0001
	a=1.

	N_AT_A = int( a/h )
	chi = []

	for l in range(0, l_max+1):
		chi.append( [ x_0**(l+1) , (x_0+h)**(l+1) ] )
	
		
		for n in range(2,N_AT_A+1):
			YN = (2.-5.*h*h*gl(E, U, l, x_0 + (n-1)*h)/6.)*chi[l][n-1]
			YN_1 = (1.+h*h*gl(E, U, l, x_0 + (n-2)*h)/12.)*chi[l][n-2]
			
			chi[l].append( (YN-YN_1)/(1.+h*h*gl(E, U, l, x_0 + n*h)/12.) )
			
		for n in range(N_AT_A+1,2*N_AT_A+1):
			YN = (2.-5.*h*h*gl(E, 0, l, x_0 + (n-1)*h)/6.)*chi[l][n-1]
			YN_1 = (1.+h*h*gl(E, 0, l, x_0 + (n-2)*h)/12.)*chi[l][n-2]
			
			chi[l].append( (YN-YN_1)/(1.+h*h*gl(E, 0, l, x_0 + n*h)/12.) )
					
			 
	return chi

# pojedyncza skladowa przekroju
def sigmal( E, l, chi_l, Jl, Nl, Jlp, Nlp ):
	h = 10e-5
	a = 1.

	k = math.sqrt( 2*E )
	N_max = int(a/h) - 5

	y = chi_l[N_max] * h / (chi_l[N_max+1]-chi_l[N_max])
	tgdl = (Jl-k*y*Jlp)/(Nl-k*y*Nlp)

	s2 = tgdl*tgdl/(tgdl*tgdl+1)
	
	return	(2*l+1)*s2*4.*math.pi/(2*E)


# całkowity
def sigmatot( E ):
	k = math.sqrt(2*E)
	siT = 0.
	lmin = 0
	lmax = 7
	siI = []

	JJ = scsp.sph_jn((lmax+2),k)
	NN = scsp.sph_yn((lmax+2),k)

	chi = CHI(E, lmax)

	
	for l in range(lmin, lmax+1):
		siI.append( sigmal( E, l, chi[l] , 
			JJ[0][l], NN[0][l], JJ[1][l], NN[1][l] ))
		siT += siI[l]
	return (siT,siI,chi)


logfile = open('test.log', 'w')                                           


# tu zalezy co chce liczyc

#EMAX=5.
#div = 10.
#for E in range(1,int(EMAX*div)+1):
#	out = sigmatot(E/div)
#	wypisz = str(E/div) +'\t' + str(out[0])	
#	for i in out[1]:
#		wypisz += '\t' + str(i)
#	
#	wypisz += '\n'
#	logfile.write( wypisz )
#	print (E/div)


chii = CHI(0.0001,5)
n = 0
for i in chii[0] :
	n += 1
	logfile.write( str(n*10e-5) + " " + str(i)+ "\n")
logfile.close()

