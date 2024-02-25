#%%
# Import Dependencies
from IPython.display import display_markdown
from numpy.linalg import norm
from pyllt import BellShape, LiftingLine, TaperedShape

#%%
# Flow Parameters
vel = 25.0
rho = 1.225

#%%
# Lifting Line
num = 101

b = 3.749 # m
c_r = 0.4 # m
c_t = 0.1 # m
m = 6.577 # kg
W = m*9.0665 # N
vel = 14.13 # m/s
rho = 1.145 # kg/m^3

c_fn = TaperedShape(c_r, c_t)

ll_p = LiftingLine('Prandtl D2', b, num, c_fn)
display_markdown(ll_p)

llr_p = ll_p.return_result_L(W, vel=vel, rho=rho)
llr_p.set_lift_distribution(W, BellShape())
llr_p.set_lifting_line_twist()
display_markdown(llr_p)

display_markdown(ll_p)

#%%
# Plot Geometry
ax_c = None
ax_c = ll_p.plot_c(ax_c)
_ = ax_c.legend()

ax_alg = None
ax_alg = ll_p.plot_alg(ax_alg)
_ = ax_alg.legend()

ax_al0 = None
ax_al0 = ll_p.plot_al0(ax_al0)
_ = ax_al0.legend()

ax_cla = None
ax_cla = ll_p.plot_cla(ax_cla)
_ = ax_cla.legend()

#%%
# Plot Results
ax_gamma = None
ax_gamma = llr_p.plot_gamma(ax_gamma)
_ = ax_gamma.legend()

print(f'norm(gamma_fn - gamma) = {norm(llr_p.gamma_fn(ll_p.s) - llr_p.gamma)}')

ax_l = None
ax_l = llr_p.plot_l(ax_l)
_ = ax_l.legend()

print(f'norm(l_fn - l) = {norm(llr_p.l_fn(ll_p.s) - llr_p.l)}')

ax_cl = None
ax_cl = llr_p.plot_cl(ax_cl)
_ = ax_cl.legend()

print(f'norm(cl_fn - cl) = {norm(llr_p.cl_fn(ll_p.s) - llr_p.cl)}')

ax_ali = None
ax_ali = llr_p.plot_ali(ax_ali)
_ = ax_ali.legend()

print(f'norm(ali_fn - ali) = {norm(llr_p.ali_fn(ll_p.s) - llr_p.ali)}')

ax_ale = None
ax_ale = llr_p.plot_ale(ax_ale)
_ = ax_ale.legend()

print(f'norm(ale_fn - ale) = {norm(llr_p.ale_fn(ll_p.s) - llr_p.ale)}')

ax_wi = None
ax_wi = llr_p.plot_wi(ax_wi)
_ = ax_wi.legend()

print(f'norm(wi_fn - wi) = {norm(llr_p.wi_fn(ll_p.s) - llr_p.wi)}')

ax_di = None
ax_di = llr_p.plot_di(ax_di)
_ = ax_di.legend()

print(f'norm(di_fn - di) = {norm(llr_p.di_fn(ll_p.s) - llr_p.di)}')

ax_cdi = None
ax_cdi = llr_p.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

print(f'norm(cdi_fn - cdi) = {norm(llr_p.cdi_fn(ll_p.s) - llr_p.cdi)}')

ax_sf = None
ax_sf = llr_p.plot_sf(ax_sf)
_ = ax_sf.legend()

# print(f'norm(sf_fn - sf) = {norm(llr_p.sf_fn(ll_p.s) - llr_p.sf)}')

ax_bm = None
ax_bm = llr_p.plot_bm(ax_bm)
_ = ax_bm.legend()

# print(f'norm(bm_fn - bm) = {norm(llr_p.bm_fn(ll_p.s) - llr_p.bm)}')

# #%%
# # Stall Conditions
# al_max = 16.0
# llr_pm = ll_p.return_result_alpha(al_max, name='Elliptical CL Max',
#                                   vel=vel, rho=rho)
# llr_pm.stall()
# display_markdown(llr_pm)

# #%%
# # Plot Geometry
# ax_c = None
# ax_c = ll_p.plot_c(ax_c)
# _ = ax_c.legend()

# ax_alg = None
# ax_alg = ll_p.plot_alg(ax_alg)
# _ = ax_alg.legend()

# ax_al0 = None
# ax_al0 = ll_p.plot_al0(ax_al0)
# _ = ax_al0.legend()

# ax_cla = None
# ax_cla = ll_p.plot_cla(ax_cla)
# ax_cla = llr_pm.plot_cla(ax_cla)
# _ = ax_cla.legend()

# #%%
# # Plot Results
# ax_gamma = None
# ax_gamma = llr_p.plot_gamma(ax_gamma)
# _ = ax_gamma.legend()

# ax_l = None
# ax_l = llr_p.plot_l(ax_l)
# _ = ax_l.legend()

# ax_cl = None
# ax_cl = llr_p.plot_cl(ax_cl)
# _ = ax_cl.legend()

# ax_ali = None
# ax_ali = llr_p.plot_ali(ax_ali)
# _ = ax_ali.legend()

# ax_ale = None
# ax_ale = llr_p.plot_ale(ax_ale)
# _ = ax_ale.legend()

# ax_wi = None
# ax_wi = llr_p.plot_wi(ax_wi)
# _ = ax_wi.legend()

# ax_di = None
# ax_di = llr_p.plot_di(ax_di)
# _ = ax_di.legend()

# ax_cdi = None
# ax_cdi = llr_p.plot_cdi(ax_cdi)
# _ = ax_cdi.legend()

# ax_sf = None
# ax_sf = llr_p.plot_sf(ax_sf)
# _ = ax_sf.legend()

# ax_bm = None
# ax_bm = llr_p.plot_bm(ax_bm)
# _ = ax_bm.legend()

# #%%
# # Create Polars
# al_degs = linspace(-25.0, 25.0, 101)

# llp_em = ll_p.return_polar(al_degs).stall()

# axp_CL = None
# axp_CL = llp_em.plot_CL(axp_CL)
# _ = axp_CL.legend()

# axp_CDi = None
# axp_CDi = llp_em.plot_CDi(axp_CDi)
# _ = axp_CDi.legend()

# axp_pol = None
# axp_pol = llp_em.plot_polar(axp_pol)
# _ = axp_pol.legend()
