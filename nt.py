'''
Library of elementary number theory

This library contains some basic number theory functions.

'''
def GCD(a,b):
	'''
		GCD(a,b) returns the gcd of two integers, a and b.

	'''
	try:
		aInt = isinstance(a,int)
		bInt = isinstance(b,int)
		if not (aInt and bInt):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of GCD is not integer.')
	else:
		if a==0:
			return b
		if b==0:
			return a
    
		A = max(abs(a),abs(b))
		B = min(abs(a),abs(b))
		while B!=0:
			tmp = B
			B = A%B
			A = tmp
		return A


def LCM(a,b):
	'''
		LCM(a,b) returns the lcm of two integers, a and b.

	'''
	try:
		aInt = isinstance(a,int)
		bInt = isinstance(b,int)
		if not (aInt and bInt):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of GCD is not integer.')
	else:
		if a==0 or b==0:
			return 0
		lcm = a*b/GCD(a,b)
		if lcm<0:
			return -lcm
		else:
			return lcm


def Factor(n):
	'''
		Factor(n) returns the factor list of n

		Return the factors of n, in the form {a1:b1, a2:b2, ...}, which means
		n has b1 many factor a1's, b2 many factor a2's, ...

	'''
	try:
		if not isinstance(n,int):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of factorization is not integer.') 
	else:
		factors = {}
		if n<0:
			n = -n
			factors[-1]=1
		if n==0 or n==1:
			return {n:1}
		mid = int(n**(1/2.0))
		n_copy = n
		for i in range(2,mid+1):
			k = 0
			if n_copy%i==0:
				while (n_copy%i==0):
					k = k+1
					n_copy = n_copy/i
				factors[i] = k
		if n_copy>1:
			factors[n_copy] = 1			
		return factors
		
		

def FindFactor(n):
	'''
		FindFactor(n) returns a factor of n

		Return a factor of n if n is a composite number.
		When the method fails, return 0.
		The method uses pollard rho factorization algorithm.

	'''
	try:
		if not isinstance(n,int):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of factorization is not integer.') 
	else:
		if n==1:
			return 1
		for i in range(1-n,1):
			x = 2
			y = 2
			d = 1
			while d==1:
				x = (x**2-i)%n
				y = (((y**2-i)%n)**2-i)%n
				d = GCD(y-x,n)
			if d!=n and d!=0:
				return d
		return 0
		
def NumOfFactors(n):
	'''
		NumOfFactors(n) returns the number of distinct factors of n.

		For example, 6 has four factors 1, 2, 3, 6, so NumOfFactors(6)=4.
	'''
	try:
		if not isinstance(n,int):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of factorization is not integer.') 
	else:
		factors = Factor(n)
		num = 1
		if n==1 or n==0:
			return n
		for i in factors.values():
			num = num*(i+1)
		return num

def SumOfFactors(n):
	'''
		SumOfFactors(n) returns sum of distinct factors of n.

		For example, 6 has four factors 1, 2, 3, 6, so 
		SumOfFactors(6) = 1 + 2 + 3 + 6 = 12.
	'''
	try:
		if not isinstance(n,int):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of factorization is not integer.') 
	else:
		factors = Factor(n)
		num = 1
		if n==1 or n==0:
			return n
		for i in factors.keys():
			num = num*(n/i)*(i-1)
		return num

def SquareFree(n):
	''' 
		SquareFree(n) return True if n is square free, and False otherwise.

		Example: 12 has a factor 4, which is a square, so SquareFree(12)
		is False.
	'''
	try:
		if not isinstance(n,int):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of factorization is not integer.') 
	else:
		if n==0 or n==1 or n==-1:
			return False
		factors = Factor(n)
		for i in factors.values():
			if i>=2:
				return False
		return True


def Phi(n):
	'''
		Phi(n) is the Euler totient function, which is the number of
		positive integers not greater than n that are coprime to n.

		For example, in 1, 2, ... , 12, the integers coprime to 12 are
		1, 5, 7, 11, so Phi(12) = 4.
	'''
	try:
		if not isinstance(n,int):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of factorization is not integer.') 
	else:
		factors = Factor(n)
		num = n
		if n==1 or n==0:
			return n
		for i in factors.keys():
			num = num/i*(i-1)
		return num

def LinCongEqn(a,b,m):
	'''
		LinCongEqn(a,b,m) gives all solutions to the congruence equation
			a*x = b (mod m)
		where 0<= x <m
		It returns all solutions if exist, and return [-1] otherwise.

		Example: 2*x = 4 (mod 20) has two solutions, x = 2, x = 12, so
		LinCongEqn(2,4,20) = [2, 12] 
	'''
	try:
		if not (isinstance(a,int) and isinstance(b,int) and\
			 isinstance(m,int)):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of GCD is not integer.')
	else:
		if m<0:
			m=-m
		if m==0:
			return [-1]
		# a * x = b (mod m) has solution if and only if
		# b = 0 (mod d), where d = GCD(a, m)
		d = GCD(a,m)
		if b%d!=0:
			return [-1]
		# otherwise has solution, now solve (a/d) * x = b/d (mod m/d)
		# Euler's theorem says (a/d)^Phi(m/d) = 1 (mod m/d)
		# so x = (a/d)^(Phi(m/d)-1) * (b/d) mod (m/d) is a solution
		# all solutions are x, x+m/d, ... , x+(d-1)*(m/d)  
		soln = []
		a1 = a/d
		b1 = b/d
		m1 = m/d
		x = a1**(Phi(m1)-1)*b1%m1
		for i in range(1,d+1):
			soln.append(x+(i-1)*m1)
		return soln

def PowMod(a,k,n):
	'''
		PowMod(a,k,n) returns a^k mod n, here k could be negative. If
		it doesn't exist, return -1.

		Example: 1 = 4 (mod 3), so 2^(-1) = 1/2 = 4/2 = 2 (mod 3), thus
		PowMod(2,-1,3) = 2
	'''
	try:
		if not (isinstance(a,int) and isinstance(a,int) and\
			 isinstance(a,int)):
			raise TypeError
	except TypeError:
		print('TypeError: Argument of GCD is not integer.')
	else:
		# in case a<0
		a = n + a
		if n<0:
			n=-n
		if n==1:
			return 0
		if k>=0:
			return (a**k)%n
		# if k<0, then let ak = (a^(-k)), now we need to find the inverse
		# of ak (mod n). If gcd(ak,n) is not 1, then no solution; otherwise
		# it has one solution ak^(Phi(n)-1)
		ak = (a**(-k))%n
		if GCD(ak,n)!=1:
			return -1
		return (ak)**(Phi(n)-1)%n
		
		
		

























		
    
    
