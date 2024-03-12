#%%
# Import Dependencies
from sympy import Symbol, cos, expand, integrate, pi, sin, expand_trig, Rational

#%%
# Create Symbols
b = Symbol('b', real=True, positive=True)
rho = Symbol('rho', real=True, positive=True)
Vinf = Symbol('Vinf', real=True, positive=True)
Sref = Symbol('Sref', real=True, positive=True)
th = Symbol('th', real=True)
n = Symbol('n', real=True, positive=True, nonzero=True, integer=True)
m = Symbol('m', real=True, positive=True, nonzero=True, integer=True)
N = Symbol('N', real=True, positive=True, nonzero=True, integer=True)
AR = Symbol('AR', real=True, positive=True)
thv = Symbol('thv', real=True)
thm = Symbol('thm', real=True)
yv = Symbol('yv', real=True)
A1 = Symbol('A1', real=True)
A2 = Symbol('A2', real=True)
A3 = Symbol('A3', real=True)
A4 = Symbol('A4', real=True)
A5 = Symbol('A5', real=True)
A6 = Symbol('A6', real=True)
A7 = Symbol('A7', real=True)
A8 = Symbol('A8', real=True)

Adct = {1: A1, 2: A2, 3: A3, 4: A4, 5: A5, 6: A6, 7: A7, 8: A8}

#%%
# Create Gamma
s = cos(th)
print(f's = {s}\n')

dsdth = s.diff(th)
print(f'dsdth = {dsdth}\n')

y = b/2*s
print(f'y = {y}\n')

dydth = y.diff(th)
print(f'dydth = {dydth}\n')

Gamma = 2*b*Vinf*sum([Ai*sin(n*th) for n, Ai in Adct.items()])
print(f'Gamma = {Gamma}\n')

l = expand(rho*Vinf*Gamma)
print(f'l = {l}\n')

integrand = expand(l*dydth)
print(f'integrand = {integrand}\n')

L = integrate(integrand, (th, pi, 0))
L = L.simplify().subs(b**2, AR*Sref).simplify()
print(f'L = {L}\n')

CL = L/(rho*Vinf**2/2*Sref)
print(f'CL = {CL}\n')

#%%
# Lift Distribution Coefficient
Cl_factor = Sref/4/b

Cl = {}
for n, An in Adct.items():
    Cln = l.coeff(An)/(rho*Vinf**2/2*Sref)
    Cln = Cln.expand().factor()*Cl_factor
    Cl[n] = Cln
    if n > 0:
        print(f'Cl_A{n} = {Cln}\n')

#%%
# Shear Force Coefficient Integration
Cv_factor = 2/b

Cv = {}
for n, An in Adct.items():
    integrand = Cl[n]*dydth
    Cvn = integrate(integrand, (th, pi, thv))
    Cvn = Cvn.expand().factor()*Cv_factor
    Cv[n] = Cvn
    print(f'Cv_A{n} = {Cvn}\n')

Cv_root = {n: Cv[n].subs(thv, pi/2) for n in Cv}
print(f'Cv_root = {Cv_root}\n')

#%%
# Bending Moment Coefficient Integration
Cm_factor = 2/b

Cm = {}
for n, An in Adct.items():
    integrand = Cv[n]*dydth.subs(th, thv)
    Cmn = integrate(integrand, (thv, pi, thm))
    Cmn = Cmn.expand().factor()*Cm_factor
    Cmn = Cmn.simplify()
    if n == 1:
        Cmn = Cmn.subs(sin(thm)**3, 3*sin(thm)/4 - sin(3*thm)/4).expand()
    Cm[n] = Cmn
    print(f'Cm_A{n} = {Cmn}')

    guess = sin((n+2)*thm)/(n+2)/(n+1)/4
    if n > 1:
        guess -= sin(n*thm)/(n+1)/(n-1)/2
    if n > 2:
        guess += sin((n-2)*thm)/(n-1)/(n-2)/4
    print(f'guess = {guess.expand().simplify()}')
    check = Cmn - guess.expand().simplify()
    check = check.expand().simplify()
    check = expand_trig(check).expand().simplify()
    print(f'check = {check}\n')

Cm_root = {n: Cm[n].subs(thm, pi/2) for n in Cm}
print(f'Cm_root = {Cm_root}\n')

Cm_root_guess = {1: Rational(1, 3), 2: pi/8 - pi/4}

for n in range(3, 9):
    if n % 2 == 0:
        Cm_root_guess[n] = 0
    else:
        Cm_root_guess[n] = Rational(1, (n+2)*(n+1)*4) + Rational(1, (n+1)*(n-1)*2) + Rational(1, (n-1)*(n-2)*4)
    if (n - 1) % 4 == 0:
        Cm_root_guess[n] = -Cm_root_guess[n]

print(f'Cm_root_guess = {Cm_root_guess}\n')

#%%
# Downwash and Drag Calculation
ali = sum([n*Ai*sin(n*th)/sin(th) for n, Ai in Adct.items()])
print(f'ali = {ali}\n')

wi = Vinf*ali
print(f'wi = {wi}\n')

di = Gamma*wi
di = di.expand()
print(f'di = {di}\n')

Di = rho*integrate(di*dydth, (th, pi, 0))
Di = Di.simplify().subs(b**2, AR*Sref).simplify()
print(f'Di = {Di}\n')

CDi = Di/(rho*Vinf**2/2*Sref)
print(f'CDi = {CDi}\n')

#%%
# Induced Drag Coefficient Integration
Cdi_factor = Sref/4/b

Cdi = {}
for n, An in Adct.items():
    Cdin = di.coeff(An)/(rho*Vinf**2/2*Sref)
    Cdin = Cdin.expand().factor()*Cdi_factor
    Cdi[n] = Cdin
    print(f'Cdi_A{n} = {Cdin}\n')

#%%
# Distribution Derivative w.r.t. theta
dGammadth = Gamma.diff(th)
print(f'dGammadth = {dGammadth}\n')

dalidth = ali.diff(th)
print(f'dalidth = {dalidth}\n')

#%%
# Drag Optimisation
dCDidA = {}

for n, An in Adct.items():
    dCDidAn = CDi.diff(An)
    dCDidA[n] = dCDidAn
    print(f'dCDi/dA{n:d} = {dCDidAn}\n')
