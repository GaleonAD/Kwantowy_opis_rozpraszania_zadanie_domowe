from numerov import Numerov_solve_f
import math

def gl ( x ):
	return 0.01*0.01

def S ( x ):
	return 0

h = 0.1
Nk = 10000

wynik =  Numerov_solve_f( S, gl, 0, Nk, 0, h, h   )
for i in range(0,Nk):
	print (str(i*h)+ ' ' +str(wynik[i]) )
