#!/usr/bin/env python3


def Numerov_step (S, gl, y, x, n, h):
	H = h*h/12.
	x_1 = x - h
	x_2 = x - 2*h
	YN_1 = (2.-10.*H*gl(x_1))*y[n - 1]
	YN_2 = - (1. + H*gl(x_2))*y[n - 2]
	SN = H*(S(x) + 10*S(x_1) + S(x_2))
	return (YN_2 + YN_1 + SN)/(1.+ H*gl(x)) 

def Numerov_solve_f (S, gl, x_0, Nk, y_0, y_1, h):
	y = [y_0, y_1]

	for n in range (2, Nk+1):
		y.append(Numerov_step(S, gl, 
			y, x_0 + n*h, n, h))

	return y

def Numerov_solve_b (S, gl, x_k, Nk, y_k, y_k_1, h):
	y = [y_k, y_k_1]

	for n in range (2, Nk+1):
		y.append(Numerov_step(S, gl, 
			y, x_k - n*h, n, -h))

	return y[::-1]

