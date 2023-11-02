
# Convention: For a bosonic operator \hat{h}, H[x] is h^{\dagger,x} for x>0, and h^{-x} for x<0.

import numpy as np
import sympy as sym
import Bosons as B

from importlib import reload
reload(B)

# Takes in the expression and spits out an equivalent expression of the form: sums over terms like f(numbers)*B1[.]*B2[.]*... 
# The algorithm is explained in lyx file CodeExplained.lyx
# One can output 'expr' which is what you normally want but can also output 'list' which is largely useful for internal computations and is otherwise, not very user-friendly.
def NumberForm(expr,output='expr'):
	
	e = expr.doit() # If there are unevaluated commutators in the expression
	e = e.expand() # We want to deal with only sums of monomials. TODO: This might slow down the process for long expressions. Maybe there is a recursive way to run NumberForm on different brackets and expand only as a possible simplification step.

	# If the input consists of multiple terms, apply the function to all the terms
	if e.is_Add:
		if (output == 'expr'):
			return sym.Add(*[NumberForm(term) for term in e.args])  # Apply to each term
		elif (output == 'list'):
			return [NumberForm(term,'list')[0] for term in e.args]
		else:
			print(output)
			raise ValueError("Unknown output method")
	
	BosNames = _unique_bosons_(e) # Gives a list of the names of all unique bosonic operators appearing in the expression

	comm,ncomm = e.args_cnc() # Separates the commutative and non-commutative part of a monomial term
	coeff = sym.Mul(*comm);  # Multiply the commutative part

	# We want to separate the number part from the annihilation and creation part
	ncomm_nums = []
	ncomm_anncre = {}
	pow_count = {}  # Keeps a record of all bosonic powers.
	# As each bosonic operator commutes with other bosonic operators, we want to separate them from each other
	for bos in BosNames:
		pow_count[bos] = 0
		ncomm_anncre[bos] = sym.sympify(1)  # starting with 1, it will keep track of all terms

	for i,t in enumerate(ncomm): # Each entry in the list ncomm can be of the type B[.] for some boson or a function (defined or otherwise) of number operators.
		
		name,power,flag = _pow_id_(t)  # Dismantle out 't'
		if flag == 'ann': # Annihilation operator
			ncomm_anncre[name] *= t
			pow_count[name] -= power
		elif flag == 'cre': # Creation operator
			ncomm_anncre[name] *= t
			pow_count[name] += power
		else:
			if t.find(B.Cre) == set() and t.find(B.Des) == set():  # TODO: Functions of ann and cre operators are not coded.
				reps = [(B.Num(bos),B.Num(bos) - pow_count[bos]) for bos in BosNames]
				ncomm_nums.append(t.subs(reps))  # The replacement reps ensures replacements like m**2 f(n_m) --> f(n_m + 2) m**2
			else:
				raise ValueError("Functions of annihilation and creation operators not supported yet")

	# End of loop

	try:
		monomial_orderings = [_monomial_num_order_(ncomm_anncre[bos]) for bos in BosNames]
	except:
		print("ncomm_nums:",ncomm_nums)
		print("ncomm_anncre:",ncomm_anncre)
		raise ValueError("Something went wrong")

	poly_mul = sym.Mul(*[t[0] for t in monomial_orderings])
	anncre_mul = sym.Mul(*[t[1] for t in monomial_orderings])

	poly = coeff*sym.Mul(*ncomm_nums)*poly_mul

	if (output == 'expr'):
		return poly*anncre_mul
	elif (output == 'list'):
		return [(poly,anncre_mul)]
	else:
		print(output)
		raise ValueError("Unknown output method")
		

# Simplification function. First, finds the NumberForm. Then, partitions the sum into factors with equal exponents, e.g. terms of the form md*m**2*ad**3*f(nm,na) and md**4*g(nm,na)*ad**4*ad*m**5 have exponents {'m':-1,'a':+3} and will be in the same partition. Finally, does more optimizations on the pre-factor of each partition.
def NC_Simplify(expr):
	
	BosNames = _unique_bosons_(expr) # Gives a list of the names of all unique bosonic operators appearing in the expression
	if len(BosNames) == 0:
		return expr

	e_list = NumberForm(expr,'list')
	if len(e_list) == 0:
		print(expr)
		raise ValueError("Something went wrong")
	
	try:
		separated_factors = _uniques_(e_list,lambda x: x[1])  # Separating the number forms by the anncre-part. 
	except:
		print(e_list)
		raise ValueError("Something went wrong")
	
	simplified_factors = []
	for anncre in separated_factors.keys():
		poly_sum = sym.Add(*[t[0] for t in separated_factors[anncre]])
		poly_sum_simplified = _simplify_(poly_sum,BosNames) # Simplify taking into account that all number operators commute with each other.
		simplified_factors.append(poly_sum_simplified*anncre)

	return sym.Add(*simplified_factors)

# Numerically tests if two exprs are equivalent
def Equivalence(e,f,Hilb=35):
	
	return 'Not Implemented yet'


################### Helper functions


def _pow_id_(t):

	t_atom = t if not t.is_Pow else t.args[0]
	t_pow = 1 if not t.is_Pow else t.args[1]

	if isinstance(t_atom, B.SQ.B):
		flag = 'ann'
	elif isinstance(t_atom, B.SQ.Bd):
		flag = 'cre'
	elif isinstance(t_atom, B.Num):
		flag = 'num'
	else: # TODO: Can certainly add integers, functions, and other stuff in flags
		try:
			return 'unknown',t_pow,t_atom.__class__
		except:
			return 'unknown',t_pow,'unknown'
	
	return t_atom.state.name,t_pow,flag


n_sym = sym.Symbol('n',commutative=False)
# The input should be a monomial of powers of annihilation and creation operators
def _monomial_num_order_(e):
	
	if e.is_Add:
		print(e)
		raise ValueError("The input should be a monomial")
	
	if not e.is_Mul:
		return (1,e)
	
	l = e.args

	if len(l)<2:
		return(1,e)

	Pow_Ids = [_pow_id_(t) for t in l]
	name = Pow_Ids[0][0]
	
	if any([t[0]!=name for t in Pow_Ids]):
		raise ValueError("The input should contain only one type of bosonic operator")
	
	if any([t[2]=='num' for t in Pow_Ids]):
		raise ValueError("The input should not contain number operators")

	if Pow_Ids[-1][2] == 'ann': # If the last term is an annihilation operator, add bd**0 at the end
		Pow_Ids.append((name,0,'cre'))
	if Pow_Ids[0][2] == 'cre': # If the first term is a creation operator, add b**0 at the end
		Pow_Ids.insert(0,(name,0,'ann'))
	
	powers = [t[1] for t in Pow_Ids]
	r = powers[1::2] # All powers of bd
	s = powers[0::2] # All powers of b

	b_ann = B.Des(name)
	b_cre = B.Cre(name)
	n_b = B.Num(name)
	poly = (_poly_(r,s)).subs(n_sym,n_b)
	d = sum(r) - sum(s)
	
	if (d>0):
		LastTerm = b_cre**(+d)
	else:
		LastTerm = b_ann**(-d)
	
	return poly,LastTerm

def _poly_single_(a,b):

	return sym.Mul(*[n_sym+r for r in range(a,b+1)])

def _poly_(r,s):

	rlen = len(r)
	if len(s) != rlen :
		print (r,s)
		raise ValueError("The lengths of the inputs should be same")

	if len(r) == 0:
		raise ValueError("Did you forget to give an input?")

	negs = [n for n in r+s if n<0]
	if len(negs) != 0:
		raise ValueError("There should be no negative numbers")
	
	if rlen == 1:
		r0 = r[0]
		s0 = s[0]

		return _poly_single_(1 + max(0,s0-r0),s0)
	
	# Divide and conquer
	rless = r[0:int(rlen/2)]
	sless = s[0:int(rlen/2)]
	rmore = r[int(rlen/2):]
	smore = s[int(rlen/2):]

	pless = _poly_(rless,sless)
	pmore = _poly_(rmore,smore)

	dless = sum(rless) - sum(sless)
	dmore = sum(rmore) - sum(smore)

	peff = pless * (pmore.subs(n_sym,n_sym-dless))
	if (dless > 0):
		if (dmore > 0):
			return peff
		else :
			return peff*_poly_single_(1-dless,min(0,-dmore-dless))
	else:
		if (dmore > 0):
			return peff*_poly_single_(1 + max(0,-dless-dmore),-dless)
		else:
			return peff


# Outputs the names of all the bosonic operators appearing in an expression
def _unique_bosons_(expr):
	
	s = set()
	for cls in [B.Cre,B.Des,B.Num,B.SQ.B,B.SQ.Bd]: # For some reason find() does not always recognize B.Cre and B.Des
		terms = expr.find(cls)  # Finds all terms of the type cls
		for t in terms:
			s.add(t.state.name)  # Does not repeat terms
	
	return sorted(s)

# The expression should contain only number operators
def _simplify_(expr,bosnames):
	
	reps = []
	rev_reps = []
	for name in bosnames:
		reps.append((B.Num(name),sym.Symbol('n_{%s}'%name)))  # We want to replace number operators by numbers as all the number operators commute with each other
		rev_reps.append((sym.Symbol('n_{%s}'%name),B.Num(name)))

	# Make all number operators commutative
	e = expr.subs(reps)

	# TODO: sym.simplify is very primitive! Would be better to implement our own simplifications based off polynomial simplification functions in sympy. An idea could be to use sym.together() on the whole expr, use sym.cancel() hoping it works well, and then possibly sym.apart(). Seems like a very computationally expensive way to do it. A manual method could be to compute gcds of denominators and figure out somehow which fractions do you actually want to add first.
	
	e_simp = sym.simplify(e)

	#print(expr)
	#print(e)
	#print(e_simp)

	return e_simp.subs(rev_reps)

# If the input to this function is g(nq,nm,na)*ad**3*q**4*md**6, then the output would be (3,-4,6)
def _coll_exp_(e,bosnames):
	if e.is_Add:
		raise ValueError("The input should be a monomial")
	
	if not e.is_Mul:
		return _sing_exp_(e,bosnames)
	
	terms = e.args
	exps = [np.array(_sing_exp_(t),bosnames) for t in terms]

	return (np.sum(exps,axis=0)).tolist()

def _sing_exp_(e,bosnames):

	res = [0]*len(bosnames)

	name,power,flag = _pow_id_(e)

	if flag == 'ann':
		signed_pow = -power
	elif flag == 'cre':
		signed_pow = power
	else:
		return res
	
	pos = bosnames.index(name)
	res[pos] = signed_pow

	return res

def _reverse_args_cnc_(l):
	if len(l) == 0:
		return 0
	mul_facs = [sym.Mul(*(t[0] + t[1])) for t in l] # Remember that '+' here is a string concatenation
	return sym.Add(*mul_facs)

# Partitions the list l s.t. in each partition, cond(*) has the same value
def _uniques_(l,cond):
	
	sep = {}
	for t in l:
		c = cond(t)
		if cond(t) in sep.keys():
			sep[c].append(t)
		else: 
			sep[c] = [t]
	
	return sep


