#%%
# Import Dependencies
from IPython.display import display_markdown
from numpy import linspace, pi, radians
from pyllt import BellShape, ConstantShape, EllipticalShape, LiftingLine, TaperedShape

#%%
# Flow Parameters
vel = 25.0
rho = 1.225

#%%
# Lifting Line
b = 4.8
c_r = 0.9
c_t = 0.6

b_b = 1.2247*b

cla_shp = ConstantShape(2*pi)

clmax_shp = ConstantShape(1.4)
clmin_shp = ConstantShape(-1.4)

c_e_shp = EllipticalShape(0.4)
alg_shp_0 = ConstantShape(0.0)
al0_shp_0 = ConstantShape(0.0)

c_b_shp = EllipticalShape(0.4)

ll_e = LiftingLine('Elliptical', b, c_e_shp)
display_markdown(ll_e)

clmax_shp = ConstantShape(1.4)
clmin_shp = ConstantShape(-1.4)

ll_e.clmax_shp = clmax_shp
ll_e.clmin_shp = clmin_shp

fac = ll_e.area/((c_r + c_t)/2*b)

c_r = c_r*fac
c_t = c_t*fac

c_t_shp = TaperedShape(c_r, c_t)
alg_shp = TaperedShape(radians(4.0), radians(2.0))
al0_shp = TaperedShape(radians(-2.0), radians(0.0))

ll_t = LiftingLine('Tapered', b, c_t_shp, alg_shp, al0_shp)
display_markdown(ll_t)

ll_t.clmax_shp = clmax_shp
ll_t.clmin_shp = clmin_shp

c_c_shp = c_e_shp.to_constant()

ll_c = LiftingLine('Constant', b, c_c_shp, alg_shp, al0_shp, cla_shp)
display_markdown(ll_c)

ll_b = LiftingLine('Bell', b_b, c_b_shp)
llr_b = ll_b.return_result_L(200.0, vel=vel, rho=rho)
llr_b.set_lift_distribution(200.0, BellShape())
llr_b.set_lifting_line_twist()
display_markdown(ll_b)

print(f'll_b.b/ll_e.b = {ll_b.b/ll_e.b:6f}')
print(f'll_b.area/ll_e.area = {ll_b.area/ll_e.area:6f}')

#%%
# Results
llr_e = ll_e.return_result_L(200.0, vel=vel, rho=rho)
display_markdown(llr_e)

display_markdown(llr_b)

print(f'llr_b.L/llr_e.L = {llr_b.L/llr_e.L:6f}')
print(f'llr_b.CL/llr_e.CL = {llr_b.CL/llr_e.CL:6f}')
print(f'llr_b.Di/llr_e.Di = {llr_b.Di/llr_e.Di:6f}')

llr_t = ll_t.return_result_L(200.0, vel=vel, rho=rho)
display_markdown(llr_t)

llr_c = ll_c.return_result_L(200.0, vel=vel, rho=rho)
display_markdown(llr_c)

#%%
# Plot Geometry
ax_c = None
ax_c = ll_e.plot_c(ax_c)
ax_c = ll_t.plot_c(ax_c)
ax_c = ll_c.plot_c(ax_c)
ax_c = ll_b.plot_c(ax_c)

_ = ax_c.legend()

ax_alg = None
ax_alg = ll_e.plot_alg(ax_alg)
ax_alg = ll_t.plot_alg(ax_alg)
ax_alg = ll_c.plot_alg(ax_alg)
ax_alg = ll_b.plot_alg(ax_alg)
_ = ax_alg.legend()

ax_al0 = None
ax_al0 = ll_e.plot_al0(ax_al0)
ax_al0 = ll_t.plot_al0(ax_al0)
ax_al0 = ll_c.plot_al0(ax_al0)
ax_al0 = ll_b.plot_al0(ax_al0)
_ = ax_al0.legend()

ax_cla = None
ax_cla = ll_e.plot_cla(ax_cla)
ax_cla = ll_t.plot_cla(ax_cla)
ax_cla = ll_c.plot_cla(ax_cla)
ax_cla = ll_b.plot_cla(ax_cla)
_ = ax_cla.legend()

#%%
# Plot Results
ax_gamma = None
ax_gamma = llr_e.plot_gamma(ax_gamma)
ax_gamma = llr_t.plot_gamma(ax_gamma)
ax_gamma = llr_c.plot_gamma(ax_gamma)
ax_gamma = llr_b.plot_gamma(ax_gamma)
_ = ax_gamma.legend()

ax_l = None
ax_l = llr_e.plot_l(ax_l)
ax_l = llr_t.plot_l(ax_l)
ax_l = llr_c.plot_l(ax_l)
ax_l = llr_b.plot_l(ax_l)
_ = ax_l.legend()

ax_cl = None
ax_cl = llr_e.plot_cl(ax_cl)
ax_cl = llr_t.plot_cl(ax_cl)
ax_cl = llr_c.plot_cl(ax_cl)
ax_cl = llr_b.plot_cl(ax_cl)
_ = ax_cl.legend()

ax_ali = None
ax_ali = llr_e.plot_ali(ax_ali)
ax_ali = llr_t.plot_ali(ax_ali)
ax_ali = llr_c.plot_ali(ax_ali)
ax_ali = llr_b.plot_ali(ax_ali)
_ = ax_ali.legend()

ax_ale = None
ax_ale = llr_e.plot_ale(ax_ale)
ax_ale = llr_t.plot_ale(ax_ale)
ax_ale = llr_c.plot_ale(ax_ale)
ax_ale = llr_b.plot_ale(ax_ale)
_ = ax_ale.legend()

ax_wi = None
ax_wi = llr_e.plot_wi(ax_wi)
ax_wi = llr_t.plot_wi(ax_wi)
ax_wi = llr_c.plot_wi(ax_wi)
ax_wi = llr_b.plot_wi(ax_wi)
_ = ax_wi.legend()

ax_di = None
ax_di = llr_e.plot_di(ax_di)
ax_di = llr_t.plot_di(ax_di)
ax_di = llr_c.plot_di(ax_di)
ax_di = llr_b.plot_di(ax_di)
_ = ax_di.legend()

ax_cdi = None
ax_cdi = llr_e.plot_cdi(ax_cdi)
ax_cdi = llr_t.plot_cdi(ax_cdi)
ax_cdi = llr_c.plot_cdi(ax_cdi)
ax_cdi = llr_b.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

ax_sf = None
ax_sf = llr_e.plot_sf(ax_sf)
ax_sf = llr_t.plot_sf(ax_sf)
ax_sf = llr_c.plot_sf(ax_sf)
ax_sf = llr_b.plot_sf(ax_sf)
_ = ax_sf.legend()

ax_bm = None
ax_bm = llr_e.plot_bm(ax_bm)
ax_bm = llr_t.plot_bm(ax_bm)
ax_bm = llr_c.plot_bm(ax_bm)
ax_bm = llr_b.plot_bm(ax_bm)
_ = ax_bm.legend()

#%%
# Stall Conditions
al_max = 16.0
llr_em = ll_e.return_result_alpha(al_max, name='Elliptical CL Max',
                                  vel=vel, rho=rho)
llr_em.stall()
display_markdown(llr_em)

llr_tm = ll_t.return_result_alpha(al_max, name='TaperedShape CL Max',
                                  vel=vel, rho=rho)
llr_tm.stall()
display_markdown(llr_tm)

#%%
# Plot Geometry
ax_c = None
ax_c = ll_e.plot_c(ax_c)
ax_c = ll_t.plot_c(ax_c)

_ = ax_c.legend()

ax_alg = None
ax_alg = ll_e.plot_alg(ax_alg)
ax_alg = ll_t.plot_alg(ax_alg)
_ = ax_alg.legend()

ax_al0 = None
ax_al0 = ll_e.plot_al0(ax_al0)
ax_al0 = ll_t.plot_al0(ax_al0)
_ = ax_al0.legend()

ax_cla = None
ax_cla = ll_e.plot_cla(ax_cla)
ax_cla = llr_em.plot_cla(ax_cla)
ax_cla = ll_t.plot_cla(ax_cla)
ax_cla = llr_tm.plot_cla(ax_cla)
_ = ax_cla.legend()

#%%
# Plot Results
ax_gamma = None
ax_gamma = llr_e.plot_gamma(ax_gamma)
ax_gamma = llr_em.plot_gamma(ax_gamma)
ax_gamma = llr_t.plot_gamma(ax_gamma)
ax_gamma = llr_tm.plot_gamma(ax_gamma)
_ = ax_gamma.legend()

ax_l = None
ax_l = llr_e.plot_l(ax_l)
ax_l = llr_em.plot_l(ax_l)
ax_l = llr_t.plot_l(ax_l)
ax_l = llr_tm.plot_l(ax_l)
_ = ax_l.legend()

ax_cl = None
ax_cl = llr_e.plot_cl(ax_cl)
ax_cl = llr_em.plot_cl(ax_cl)
ax_cl = llr_t.plot_cl(ax_cl)
ax_cl = llr_tm.plot_cl(ax_cl)
_ = ax_cl.legend()

ax_ali = None
ax_ali = llr_e.plot_ali(ax_ali)
ax_ali = llr_em.plot_ali(ax_ali)
ax_ali = llr_t.plot_ali(ax_ali)
ax_ali = llr_tm.plot_ali(ax_ali)
_ = ax_ali.legend()

ax_ale = None
ax_ale = llr_e.plot_ale(ax_ale)
ax_ale = llr_em.plot_ale(ax_ale)
ax_ale = llr_t.plot_ale(ax_ale)
ax_ale = llr_tm.plot_ale(ax_ale)
_ = ax_ale.legend()

ax_wi = None
ax_wi = llr_e.plot_wi(ax_wi)
ax_wi = llr_em.plot_wi(ax_wi)
ax_wi = llr_t.plot_wi(ax_wi)
ax_wi = llr_tm.plot_wi(ax_wi)
_ = ax_wi.legend()

ax_di = None
ax_di = llr_e.plot_di(ax_di)
ax_di = llr_em.plot_di(ax_di)
ax_di = llr_t.plot_di(ax_di)
ax_di = llr_tm.plot_di(ax_di)
_ = ax_di.legend()

ax_cdi = None
ax_cdi = llr_e.plot_cdi(ax_cdi)
ax_cdi = llr_em.plot_cdi(ax_cdi)
ax_cdi = llr_t.plot_cdi(ax_cdi)
ax_cdi = llr_tm.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

ax_sf = None
ax_sf = llr_e.plot_sf(ax_sf)
ax_sf = llr_em.plot_sf(ax_sf)
ax_sf = llr_t.plot_sf(ax_sf)
ax_sf = llr_tm.plot_sf(ax_sf)
_ = ax_sf.legend()

ax_bm = None
ax_bm = llr_e.plot_bm(ax_bm)
ax_bm = llr_em.plot_bm(ax_bm)
ax_bm = llr_t.plot_bm(ax_bm)
ax_bm = llr_tm.plot_bm(ax_bm)
_ = ax_bm.legend()

#%%
# Create Polars
al_degs = linspace(-25.0, 25.0, 101)

llp_em = ll_e.return_polar(al_degs).stall()

llp_tm = ll_t.return_polar(al_degs).stall()

axp_CL = None
axp_CL = llp_em.plot_CL(axp_CL)
axp_CL = llp_tm.plot_CL(axp_CL)
_ = axp_CL.legend()

axp_CDi = None
axp_CDi = llp_em.plot_CDi(axp_CDi)
axp_CDi = llp_tm.plot_CDi(axp_CDi)
_ = axp_CDi.legend()

axp_pol = None
axp_pol = llp_em.plot_polar(axp_pol)
axp_pol = llp_tm.plot_polar(axp_pol)
_ = axp_pol.legend()

#%%
# Plot CL
ax_cl = None
ax_cl = llp_em.results[al_degs[-1]].plot_cl(ax_cl)
ax_cl = llp_tm.results[al_degs[-1]].plot_cl(ax_cl)
_ = ax_cl.legend()
