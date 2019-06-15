
import rozp
logfile = open('test.log', 'w')                                           

U = 100
a = 1
h = 10e-3
l = 5
STOP =4000
div = 100.

print (STOP/div)

for it in range(1, STOP+1):
	E = it/div
	print (E)
	out = rozp.SIGMA_TOT( E, U, l, h, a)
	
	wypisz = str(E) +'\t' + str(out[0])	
	for i in out[1]:
		wypisz += '\t' + str(i)
	wypisz += '\n'
	logfile.write( wypisz )

#chii = rozp.CHI_L( E,U,l,h,a )
#n = 0
#for i in chii :
#	n += 1
#	logfile.write( str(n*h) + " " + str(i)+ "\n")


logfile.close()
