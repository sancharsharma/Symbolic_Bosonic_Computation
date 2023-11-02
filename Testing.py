
import numpy as np
import sympy as sym
import sympy.physics.secondquant as SQ
import Bosons as B
import Simps_Boson as Fun 
import re

from importlib import reload
reload(B)
reload(Fun)

wm = B.SQ_Symbol('omega_m',real=True)
wa = B.SQ_Symbol('omega_a',real=True)
wq = B.SQ_Symbol('omega_q',real=True)
al = B.SQ_Symbol('alpha',real=True)

gam = B.SQ_Symbol('g_am',real=True)
gaq = B.SQ_Symbol('g_aq',real=True)

def defs(name):
	return B.Des(name),B.Cre(name),B.Num(name)
	#sym.Symbol("\hat{n}_{%s}" %name, commutative = False)
	#B.Num(name) 
	#B.Cre(name)*B.Des(name)

m,md,nm = defs('m')
q,qd,nq = defs('q')
a,ad,na = defs('a')

exprs = []
exprs.append(md*m - m*md)
exprs.append(md**2*m*ad*a - m*md**2*a*ad)
exprs.append(md**3*m**9*md**6 - m**9*md**9)
exprs.append(md*(2-3*nm)*(8-9*nq)*m*q*qd+4*a)

expr_numforms = [Fun.NumberForm(e) for e in exprs]
expr_simps = [Fun.NC_Simplify(e) for e in exprs]

## In these expression, the bottleneck of simplification is the polynomial simplification of sympy. Not sure what to do about it!
