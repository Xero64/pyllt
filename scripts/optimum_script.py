#%%
# Import Dependencies
from IPython.display import display_markdown
from numpy import pi

from pyllt import BellShape, ConstantShape, EllipticalShape, LiftingLine

#%%
# Flow Parameters
lift = 200.0
vel = 25.0
rho = 1.225

#%%
# Lifting Line
cd0 = 0.02
b = 4.8

b_fac = 1.2247

b_b = b_fac*b

cla_shp = ConstantShape(2*pi)

clmax_shp = ConstantShape(1.4)
clmin_shp = ConstantShape(-1.4)

c_e_shp = EllipticalShape(0.4)

c_b_shp = EllipticalShape(0.4)

c_o_shp = EllipticalShape(0.4/b_fac)

ll_e = LiftingLine('Elliptical', b, c_e_shp, cd0=cd0)
display_markdown(ll_e)

ll_e.clmax_shp = clmax_shp
ll_e.clmin_shp = clmin_shp

llr_e = ll_e.return_result_L(lift, vel=vel, rho=rho)

ll_b = LiftingLine('Bell', b_b, c_b_shp, cd0=cd0)

ll_b.clmax_shp = clmax_shp
ll_b.clmin_shp = clmin_shp

llr_b = ll_b.return_result_L(lift, vel=vel, rho=rho)
llr_b.set_lift_distribution(lift, BellShape())
llr_b.set_lifting_line_twist()
display_markdown(ll_b)

print(f'll_b.b/ll_e.b = {ll_b.b/ll_e.b:6f}')
print(f'll_b.area/ll_e.area = {ll_b.area/ll_e.area:6f}')

ll_o = LiftingLine('Optimum', b_b, c_o_shp, cd0=cd0)

ll_o.clmax_shp = clmax_shp
ll_o.clmin_shp = clmin_shp

llr_o = ll_o.return_result_L(lift, vel=vel, rho=rho)
llr_o.minimum_induced_drag_optimum(lift, llr_b.bmr, num=12)
llr_o.set_lifting_line_twist()
display_markdown(ll_o)

print(f'llr_b.L/llr_e.L = {llr_b.L/llr_e.L:6f}')
print(f'llr_b.CL/llr_e.CL = {llr_b.CL/llr_e.CL:6f}')
print(f'llr_b.Di/llr_e.Di = {llr_b.Di/llr_e.Di:6f}')
print(f'llr_b.bmr/llr_e.bmr = {llr_b.Di/llr_e.Di:6f}')

#%%
# Results
display_markdown(llr_e)

display_markdown(llr_b)

display_markdown(llr_o)

#%%
# Plot Geometry
ax_c = None
ax_c = ll_e.plot_c(ax_c)
ax_c = ll_b.plot_c(ax_c)
ax_c = ll_o.plot_c(ax_c)
_ = ax_c.legend()

ax_alg = None
ax_alg = ll_e.plot_alg(ax_alg)
ax_alg = ll_b.plot_alg(ax_alg)
ax_alg = ll_o.plot_alg(ax_alg)
_ = ax_alg.legend()

ax_al0 = None
ax_al0 = ll_e.plot_al0(ax_al0)
ax_al0 = ll_b.plot_al0(ax_al0)
ax_al0 = ll_o.plot_al0(ax_al0)
_ = ax_al0.legend()

ax_cla = None
ax_cla = ll_e.plot_cla(ax_cla)
ax_cla = ll_b.plot_cla(ax_cla)
ax_cla = ll_o.plot_cla(ax_cla)
_ = ax_cla.legend()

#%%
# Plot Results
ax_gamma = None
ax_gamma = llr_e.plot_gamma(ax_gamma)
ax_gamma = llr_b.plot_gamma(ax_gamma)
ax_gamma = llr_o.plot_gamma(ax_gamma)
_ = ax_gamma.legend()

ax_l = None
ax_l = llr_e.plot_l(ax_l)
ax_l = llr_b.plot_l(ax_l)
ax_l = llr_o.plot_l(ax_l)
_ = ax_l.legend()

ax_cl = None
ax_cl = llr_e.plot_cl(ax_cl)
ax_cl = llr_b.plot_cl(ax_cl)
ax_cl = llr_o.plot_cl(ax_cl)
_ = ax_cl.legend()

ax_ali = None
ax_ali = llr_e.plot_ali(ax_ali)
ax_ali = llr_b.plot_ali(ax_ali)
ax_ali = llr_o.plot_ali(ax_ali)
_ = ax_ali.legend()

ax_ale = None
ax_ale = llr_e.plot_ale(ax_ale)
ax_ale = llr_b.plot_ale(ax_ale)
ax_ale = llr_o.plot_ale(ax_ale)
_ = ax_ale.legend()

ax_wi = None
ax_wi = llr_e.plot_wi(ax_wi)
ax_wi = llr_b.plot_wi(ax_wi)
ax_wi = llr_o.plot_wi(ax_wi)
_ = ax_wi.legend()

ax_di = None
ax_di = llr_e.plot_di(ax_di)
ax_di = llr_b.plot_di(ax_di)
ax_di = llr_o.plot_di(ax_di)
_ = ax_di.legend()

ax_cdi = None
ax_cdi = llr_e.plot_cdi(ax_cdi)
ax_cdi = llr_b.plot_cdi(ax_cdi)
ax_cdi = llr_o.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

ax_sf = None
ax_sf = llr_e.plot_sf(ax_sf)
ax_sf = llr_b.plot_sf(ax_sf)
ax_sf = llr_o.plot_sf(ax_sf)
_ = ax_sf.legend()

ax_bm = None
ax_bm = llr_e.plot_bm(ax_bm)
ax_bm = llr_b.plot_bm(ax_bm)
ax_bm = llr_o.plot_bm(ax_bm)
_ = ax_bm.legend()
