import rozp

U = 100
a = 1
h = 10e-4
l = 3
E = 9.7

logfile = open('chi'+str(l), 'w')                                           
chii = rozp.CHI_L( E,U,l,h,a )
n = 0
maax = 0
for i in chii :
	if i > maax : maax = i

for i in chii :
	n += 1
	logfile.write( str(n*h) + " " + str(i/maax)+ "\n")


logfile.close()
