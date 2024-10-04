#%%
# Import Dependencies
from IPython.display import display_markdown
from numpy import linspace, pi, radians
from pyllt import (BellShape, ConstantShape, EllipticalShape, LiftingLine,
                   TaperedShape)

#%%
# Flow Parameters
n = 3.0
lift = 300.0
vel = 25.0
vels = linspace(10.0, 35.0, 50)
rho = 1.225

#%%
# Lifting Line
cd0 = 0.006
# b_e = 2.617170337 # Minimum Power at 25 m/s
# c_c = 0.582630605 # Minimum Power at 25 m/s
b_e = 4.533 # Maximum L/D at 25 m/s
c_c = 0.336 # Maximum L/D at 25 m/s
c_r = c_c*1.1
c_t = c_c*0.9
c_e = (c_r + c_t)*2/pi

b_b = 5.234
c_b = 0.291*4/pi

cla_shp = ConstantShape(2*pi)

clmax_shp = ConstantShape(1.4)
clmin_shp = ConstantShape(-1.4)

c_e_shp = EllipticalShape(c_e)
alg_shp_0 = ConstantShape(0.0)
al0_shp_0 = ConstantShape(0.0)

c_b_shp = EllipticalShape(c_b)

ll_e = LiftingLine('Elliptical', b_e, c_e_shp, cd0=cd0)
display_markdown(ll_e)

ll_e.clmax_shp = clmax_shp
ll_e.clmin_shp = clmin_shp

c_t_shp = TaperedShape(c_r, c_t)
alg_shp = TaperedShape(radians(4.0), radians(2.0))
al0_shp = TaperedShape(radians(-2.0), radians(0.0))

ll_t = LiftingLine('Tapered', b_e, c_t_shp, alg_shp, al0_shp, cd0=cd0)
display_markdown(ll_t)

ll_t.clmax_shp = clmax_shp
ll_t.clmin_shp = clmin_shp

c_c_shp = c_e_shp.to_constant()

ll_c = LiftingLine('Constant', b_e, c_c_shp, alg_shp, al0_shp, cd0=cd0)
display_markdown(ll_c)

ll_c.clmax_shp = clmax_shp
ll_c.clmin_shp = clmin_shp

ll_b = LiftingLine('Bell', b_b, c_b_shp, cd0=cd0)
llr_b = ll_b.return_result_L(lift, vel=vel, rho=rho)
llr_b.set_lift_distribution(lift, BellShape())
llr_b.set_lifting_line_twist()
display_markdown(ll_b)

ll_b.clmax_shp = clmax_shp
ll_b.clmin_shp = clmin_shp

print(f'll_b.b/ll_e.b = {ll_b.b/ll_e.b:6f}')
print(f'll_b.area/ll_e.area = {ll_b.area/ll_e.area:6f}')

#%%
# Results
llr_e = ll_e.return_result_L(lift, vel=vel, rho=rho)
display_markdown(llr_e)

display_markdown(llr_b)

print(f'llr_b.L/llr_e.L = {llr_b.L/llr_e.L:6f}')
print(f'llr_b.CL/llr_e.CL = {llr_b.CL/llr_e.CL:6f}')
print(f'llr_b.Di/llr_e.Di = {llr_b.Di/llr_e.Di:6f}')

llr_t = ll_t.return_result_L(lift, vel=vel, rho=rho)
display_markdown(llr_t)

llr_c = ll_c.return_result_L(lift, vel=vel, rho=rho)
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
# Create Level Flight Plots
llf_e = ll_e.return_level_flight(lift, rho, vels)
llf_t = ll_t.return_level_flight(lift, rho, vels)
llf_c = ll_c.return_level_flight(lift, rho, vels)
llf_b = ll_b.return_level_flight(lift, rho, vels)

#%%
# Level Flight Plots
axf_CL = None
axf_CL = llf_e.plot_CL(axf_CL)
axf_CL = llf_t.plot_CL(axf_CL)
axf_CL = llf_c.plot_CL(axf_CL)
axf_CL = llf_b.plot_CL(axf_CL)
_ = axf_CL.legend()

axf_CL = None
axf_CL = llf_e.plot_CD(axf_CL)
axf_CL = llf_t.plot_CD(axf_CL)
axf_CL = llf_c.plot_CD(axf_CL)
axf_CL = llf_b.plot_CD(axf_CL)
_ = axf_CL.legend()

axf_LoD = None
axf_LoD = llf_e.plot_LoD(axf_LoD)
axf_LoD = llf_t.plot_LoD(axf_LoD)
axf_LoD = llf_c.plot_LoD(axf_LoD)
axf_LoD = llf_b.plot_LoD(axf_LoD)
_ = axf_LoD.legend()

axf_pwr = None
axf_pwr = llf_e.plot_pwr(axf_pwr)
axf_pwr = llf_t.plot_pwr(axf_pwr)
axf_pwr = llf_c.plot_pwr(axf_pwr)
axf_pwr = llf_b.plot_pwr(axf_pwr)
_ = axf_pwr.legend()

#%%
# Create Level Flight Plots
lln_e = ll_e.return_level_flight(n*lift, rho, vels).stall()
lln_t = ll_t.return_level_flight(n*lift, rho, vels).stall()
lln_c = ll_c.return_level_flight(n*lift, rho, vels).stall()
lln_b = ll_b.return_level_flight(n*lift, rho, vels).stall()

#%%
# Level Flight Plots
axf_CL = None
axf_CL = lln_e.plot_CL(axf_CL)
axf_CL = lln_t.plot_CL(axf_CL)
axf_CL = lln_c.plot_CL(axf_CL)
axf_CL = lln_b.plot_CL(axf_CL)
_ = axf_CL.legend()

axf_CL = None
axf_CL = lln_e.plot_CD(axf_CL)
axf_CL = lln_t.plot_CD(axf_CL)
axf_CL = lln_c.plot_CD(axf_CL)
axf_CL = lln_b.plot_CD(axf_CL)
_ = axf_CL.legend()

axf_LoD = None
axf_LoD = lln_e.plot_LoD(axf_LoD)
axf_LoD = lln_t.plot_LoD(axf_LoD)
axf_LoD = lln_c.plot_LoD(axf_LoD)
axf_LoD = lln_b.plot_LoD(axf_LoD)
_ = axf_LoD.legend()

axf_pwr = None
axf_pwr = lln_e.plot_pwr(axf_pwr)
axf_pwr = lln_t.plot_pwr(axf_pwr)
axf_pwr = lln_c.plot_pwr(axf_pwr)
axf_pwr = lln_b.plot_pwr(axf_pwr)
_ = axf_pwr.legend()

axf_L = None
axf_L = lln_e.plot_custom('L', axf_L)
axf_L = lln_t.plot_custom('L', axf_L)
axf_L = lln_c.plot_custom('L', axf_L)
axf_L = lln_b.plot_custom('L', axf_L)
_ = axf_L.legend()
