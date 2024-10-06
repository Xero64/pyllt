#%%
# Import Dependencies
from typing import List

from IPython.display import display_markdown
from matplotlib.pyplot import figure
from numpy import asarray, linspace, pi, sqrt
from pyllt.liftingline import LiftingLineResult

from pyllt import BellShape, ConstantShape, EllipticalShape, LiftingLine

#%%
# Flow Parameters
lift = 200.0
vel = 25.0
rho = 1.225

#%%
# Lifting Line
cd0 = 0.006
b_e = 4.8
c_e = 0.4

num = 30

b_b_fac = sqrt(1.5)

fac = 3.0

b_fac = linspace(1.0, 1.0 + fac*(b_b_fac - 1.0), num)

b = b_e*b_fac
c = c_e/b_fac

cla_shp = ConstantShape(2*pi)

b_b = b_b_fac*b_e
c_b_shp = EllipticalShape(c_e/b_b_fac)

ll_b = LiftingLine('Bell', b_b, c_b_shp, cd0=cd0)
llr_b = ll_b.return_result_L(lift, vel=vel, rho=rho)
llr_b.set_lift_distribution(lift, BellShape())
llr_b.set_lifting_line_twist()
display_markdown(ll_b)

display_markdown(llr_b)

lls: List[LiftingLine] = []
llrs: List[LiftingLineResult] = []

for bi, ci in zip(b, c):
    c_shp = EllipticalShape(ci)
    ll = LiftingLine(f'LiftingLine {bi:.5f}', bi, c_shp, cd0=cd0)
    lls.append(ll)
    llr = ll.return_result_L(lift, vel=vel, rho=rho)
    llr.minimum_induced_drag_optimum(lift, llr_b.bmr, num=12)
    llr.set_lifting_line_twist()
    llrs.append(llr)


#%%
# Plot Geometry
ax_c = None
for ll in lls:
    ax_c = ll.plot_c(ax_c)

ax_alg = None
for ll in lls:
    ax_alg = ll.plot_alg(ax_alg)

ax_al0 = None
for ll in lls:
    ax_al0 = ll.plot_al0(ax_al0)

ax_cla = None
for ll in lls:
    ax_cla = ll.plot_cla(ax_cla)

#%%
# Plot Results
ax_gamma = None
for llr in llrs:
    ax_gamma = llr.plot_gamma(ax_gamma)

ax_l = None
for llr in llrs:
    ax_l = llr.plot_l(ax_l)

ax_cl = None
for llr in llrs:
    ax_cl = llr.plot_cl(ax_cl)

ax_ali = None
for llr in llrs:
    ax_ali = llr.plot_ali(ax_ali)

ax_ale = None
for llr in llrs:
    ax_ale = llr.plot_ale(ax_ale)

ax_wi = None
for llr in llrs:
    ax_wi = llr.plot_wi(ax_wi)

ax_di = None
for llr in llrs:
    ax_di = llr.plot_di(ax_di)

ax_cdi = None
for llr in llrs:
    ax_cdi = llr.plot_cdi(ax_cdi)

ax_sf = None
for llr in llrs:
    ax_sf = llr.plot_sf(ax_sf)

ax_bm = None
for llr in llrs:
    ax_bm = llr.plot_bm(ax_bm)

#%%
# Plot Span Results
Di0 = llrs[0].Di

Di = asarray([llr.Di/Di0 for llr in llrs])

fig_Di = figure()
ax_Di = fig_Di.gca()
ax_Di.grid(True)
ax_Di.set_xlabel('b')
ax_Di.set_ylabel('Di')
_ = ax_Di.plot(b, Di)

Sref = asarray([ll.area for ll in lls])

fig_Sref = figure()
ax_Sref = fig_Sref.gca()
ax_Sref.grid(True)
ax_Sref.set_xlabel('b')
ax_Sref.set_ylabel('Sref')
_ = ax_Sref.plot(b, Sref)

e = asarray([llr.e for llr in llrs])

fig_e = figure()
ax_e = fig_e.gca()
ax_e.grid(True)
ax_e.set_xlabel('b')
ax_e.set_ylabel('e')
_ = ax_e.plot(b, e)

ar = asarray([ll.ar for ll in lls])

fig_ar = figure()
ax_ar = fig_ar.gca()
ax_ar.grid(True)
ax_ar.set_xlabel('b')
ax_ar.set_ylabel('ar')
_ = ax_ar.plot(b, ar)

are = asarray([ll.ar*llr.e for ll, llr in zip(lls, llrs)])

fig_are = figure()
ax_are = fig_are.gca()
ax_are.grid(True)
ax_are.set_xlabel('b')
ax_are.set_ylabel('ar*e')
_ = ax_are.plot(b, are)

#%%
# 3g Stall Result
llr_bm = ll_b.return_result_L(3*lift, vel=vel, rho=rho)
display_markdown(llr_bm)

#%%
# Plot Results
ax_gamma = llr_bm.plot_gamma()

ax_l = llr_bm.plot_l()

ax_cl = llr_bm.plot_cl()

ax_ali = llr_bm.plot_ali()

ax_ale = llr_bm.plot_ale()

ax_wi = llr_bm.plot_wi()

ax_di = llr_bm.plot_di()

ax_cdi = llr_bm.plot_cdi()

ax_sf = llr_bm.plot_sf()

ax_bm = llr_bm.plot_bm()
