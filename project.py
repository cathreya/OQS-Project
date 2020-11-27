from sympy.abc import x,p
from sympy import poly, roots, Matrix, eye
import numpy as np

# Matrix I \otimes \Lamda_T(\rho) - I\lamda
A = Matrix([[-x,0,0,0,0,0,0,0,1/3],
	        [0,-x,0,0,0,1/3,0,0,0],
	        [0,0,1/3-x,0,0,0,0,0,0], 
	        [0,0,0,-x,0,0,0,1/3,0], 
	        [0,0,0,0,1/3-x,0,0,0,0], 
	        [0,1/3,0,0,0,-x,0,0,0], 
	        [0,0,0,0,0,0,1/3-x,0,0], 
	        [0,0,0,1/3,0,0,0,-x,0],
	        [1/3,0,0,0,0,0,0,0,-x]])


print(A.det())

det = poly(A.det())
r = roots(det)

print(r)



B = Matrix([[p/9,0,0,0,0,0,0,0,(1-p)/3],
			   [0,p/9,0,0,0,(1-p)/3,0,0,0],
			   [0,0,p/9+(1-p)/3,0,0,0,0,0,0],
			   [0,0,0,p/9,0,0,0,(1-p)/3,0],
			   [0,0,0,0,p/9+(1-p)/3,0,0,0,0],
			   [0,(1-p)/3,0,0,0,p/9,0,0,0],
			   [0,0,0,0,0,0,p/9+(1-p)/3,0,0],
			   [0,0,0,(1-p)/3,0,0,0,p/9,0],
			   [(1-p)/3,0,0,0,0,0,0,0,p/9]])

lamda = eye(9)*x

B = B - lamda


charEqn = poly(B.det())

pvals = []
results = []

for i in np.linspace(0,1,100):
	pvals.append(i)
	charEqni = charEqn.eval(p,i)
	ri = roots(charEqni)
	if any(t < 0 for t in ri.keys()):
		results.append(0)
	else:
		results.append(1)



l = 0
r = 1

for i in range(1000):
	m = (l+r)/2
	print(m)
	charEqni = charEqn.eval(p,m)
	ri = roots(charEqni)
	if any(t < 0 for t in ri.keys()):
		l = m
	else:
		r = m

print(l)



