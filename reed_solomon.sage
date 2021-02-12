#!/usr/bin/env sage

import sys
import copy
import random
from sage.all import *


def gen_error(amnt_service_pos, L):
	seed()
	out = [L[0] for i in range(7-1-(amnt_service_pos//2))]
	for i in range(amnt_service_pos//2):
		out.insert(randint(0, 7-1-(amnt_service_pos//2)+i), L[randint(0,7-1)])
	return out


def encode(vec, L):
	R = PolynomialRing(GF(7), 'x')(vec)	
	return [R(L[5]**i) for i in range(7-1)]


def decode(vec, L):
	R = PolynomialRing(GF(7), 'x')(vec)
	return [R(L[5]**(-i))/L[7-1] for i in range(7-1)]


def berlekamp_messi(S, L):	
	p.<x> = PolynomialRing(GF(7), 'x')
	C = p([1])
	B = p([1])
	r = 0
	l = 0
	while (r < len(S)):
		delta = L[0]
		for i in range(len(C.list())):
			delta += S[r-i]*C.list()[i]

		if (delta != 0):
			T = C - delta * x * B
			if (2*l<=r):
				B = (delta**(-1))*C
				C = T
				l = r-l+1
			else:
				C = T
				B *= x
		else:
			B *= x
		r += 1
	return C


def forni(S, L, poly, amnt_service_pos):
	out = deepcopy(S)
	poly_list = poly.list()
	poly_len = len(poly_list)-1
	for i in range(7-1-amnt_service_pos):
		out.insert(0,L[0])
		for r in range(poly_len):
			out[0] -= out[poly_len-r]*poly_list[r]/poly_list[-1]
	return encode(out, L)


def main():
	L = GF(7).list()
	amnt_service_pos = 7-1-len(sys.argv[1])
	inp = [L[int(i)] for i in sys.argv[1]]
	for i in range(amnt_service_pos):
		inp.append(L[0])
	err = gen_error(amnt_service_pos, L)

	print("input vector: ",inp)
	print("error vector: ",err)
	enc = encode(inp, L)
	print("encode vector:", enc)
	noisy = [enc[i]+err[i] for i in range(7-1)]
	print("noisy vector: ", noisy)
	dec = decode(noisy, L)
	print("decode vector:", dec)
	if (dec[-amnt_service_pos:] != [L[0] for i in range(amnt_service_pos)]):
		poly = berlekamp_messi(dec[-amnt_service_pos:], L)
		print("Berlekamp-Messi poly: ",poly)
		print("list of error positions: "," ".join([(str)(r[0].log(5)) for r in poly.roots()]))
		err = forni(dec[-amnt_service_pos:], L, poly, amnt_service_pos)
		print("forni err vector:", err)
		enc = [noisy[i]-err[i] for i in range(7-1)]
		print("cleaned encode vector:", enc)
		dec = decode(enc, L)
		print("decode vector:", dec)


if __name__ == "__main__":
	main()
