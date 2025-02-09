#%%
# Import Dependencies
from sympy import (Expr, Symbol, diff, expand, integrate, sqrt,
                   together, fraction)
from sympy.plotting import plot
from sympy.solvers import solve

#%%
# Create Symbols
s = Symbol('s', real=True, positive=True)
An = Symbol('An', real=True)
cr = Symbol('cr', real=True, positive=True, nonzero=True)
ct = Symbol('ct', real=True, positive=True, nonzero=True)
al0r = Symbol('al0r', real=True)
al0t = Symbol('al0t', real=True)
ali = Symbol('ali', real=True)

a = Symbol('a', real=True)
b = Symbol('b', real=True)
c = Symbol('c', real=True)
d = Symbol('d', real=True)
e = Symbol('e', real=True)
f = Symbol('f', real=True)

#%%
# Elliptical Twist Equation
alg_ell = An*sqrt(1 - s**2)/(cr*(1 - s) + ct*s) + ali + al0r*(1 - s) + al0t*s
alg_ell: Expr = together(expand(alg_ell))
print(f'{alg_ell = }\n')

num_ell, den_ell = fraction(alg_ell)
num_ell: Expr = expand(num_ell)
den_ell: Expr = expand(den_ell)
print(f'{num_ell = }\n')
print(f'{den_ell = }\n')

# num_ell = alg_ell

#%%
# Define the Elliptical Best Fit Equation
syms = [a, b, c, d, e]
num_fit: Expr = a*s**4 + b*s**3 + c*s**2 + d*s + e

eqn1 = num_fit.subs(s, 0) - num_ell.subs(s, 0)
eqn2 = num_fit.subs(s, 1) - num_ell.subs(s, 1)

res: dict[Symbol, Expr] = solve([eqn1, eqn2], syms[-2:], dict=True)[0]

for sym, val in res.items():
    print(f'{sym} = {val}\n')

num_fit = num_fit.subs(res).simplify()
print(f'{num_fit = }\n')

r = (num_fit - num_ell)**2
print(f'{r = }\n')

R = integrate(r, (s, 0, 1)).simplify()
print(f'{R = }\n')

dRds = []

for sym in syms[:-2]:
    dRds.append(diff(R, sym).simplify())

res = solve(dRds, syms[:-2], dict=True)[0]
print(f'{res = }\n')

num_fit = num_fit.subs(res).simplify()
print(f'{num_fit = }\n')

alg_fit = num_fit/den_ell
print(f'{alg_fit = }\n')

#%%
# Plot the Elliptical Best Fit
# alg = (((0.03757357259624004*sin(th))/(0.392 - abs(s)*0.1372)) + (-0.060562925044203235 - abs(s)*-0.060562925044203235) + (0.017087559298996885))
sbs = {
    An: 0.03757357259624004,
    cr: 0.392,
    ct: 0.392 - 0.1372,
    al0r: -0.060562925044203235,
    al0t: 0.0,
    ali: 0.017087559298996885
}

num_fit_eval = num_fit.subs(sbs)
print(f'{num_fit_eval = }\n')

num_ell_eval = num_ell.subs(sbs)
print(f'{num_ell_eval = }\n')

plot(num_fit_eval, num_ell_eval, (s, 0, 1), ylabel='Twist [rad]', xlabel='s')

#%%
# Plot the Elliptical Best Fit
# alg = (((0.03757357259624004*sin(th))/(0.392 - abs(s)*0.1372)) + (-0.060562925044203235 - abs(s)*-0.060562925044203235) + (0.017087559298996885))
sbs = {
    An: 0.03757357259624004,
    cr: 0.392,
    ct: 0.392 - 0.1372,
    al0r: -0.060562925044203235,
    al0t: 0.0,
    ali: 0.017087559298996885
}

alg_fit_eval = alg_fit.subs(sbs)
print(f'{alg_fit_eval = }\n')

alg_ell_eval = alg_ell.subs(sbs)
print(f'{alg_ell_eval = }\n')

plot(alg_fit_eval, alg_ell_eval, (s, 0, 1), ylabel='Twist [rad]', xlabel='s')
