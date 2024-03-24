#%%
# Import Dependencies
from IPython.display import display_markdown
from numpy import pi, sqrt, divide
from pyllt import BellShape, ConstantShape, EllipticalShape, LiftingLine
from matplotlib.pyplot import figure

#%%
# Flow Parameters
vel = 25.0
rho = 1.225

#%%
# Lifting Line
cd0 = 0.02
b = 4.8
c_r = 0.9
c_t = 0.6

b_b = sqrt(1.5)*b

cla_shp = ConstantShape(2*pi)

clmax_shp = ConstantShape(1.4)
clmin_shp = ConstantShape(-1.4)

c_b_shp = EllipticalShape(0.4)

ll_b = LiftingLine('Bell', b_b, c_b_shp, cd0=cd0)

llr_b = ll_b.return_result_L(200.0, vel=vel, rho=rho)
llr_b.set_lift_distribution(200.0, BellShape())
llr_b.set_lifting_line_twist()
display_markdown(ll_b)

ll_b.clmax_shp = clmax_shp
ll_b.clmin_shp = clmin_shp

#%%
# Stall Conditions
llr_bm1 = ll_b.return_result_alpha(16.0, name=f'Bell AoA 16.0 deg',
                                   vel=vel, rho=rho)
llr_bm1.stall(display=True)
display_markdown(llr_bm1)

llr_bm2 = ll_b.return_result_alpha(12.0, name='Bell AoA 12.0 deg',
                                   vel=vel, rho=rho)

llr_bm2.cla = llr_bm1.cla

# llr_bm2.stall(display=True)
display_markdown(llr_bm2)

#%%
# Plot Geometry
ax_c = None
ax_c = ll_b.plot_c(ax_c)
_ = ax_c.legend()

ax_alg = None
ax_alg = ll_b.plot_alg(ax_alg)
_ = ax_alg.legend()

ax_al0 = None
ax_al0 = ll_b.plot_al0(ax_al0)
_ = ax_al0.legend()

ax_cla = None
ax_cla = llr_b.plot_cla(ax_cla)
ax_cla = llr_bm1.plot_cla(ax_cla)
ax_cla = llr_bm2.plot_cla(ax_cla)
_ = ax_cla.legend()

#%%
# Plot Results
ax_gamma = None
ax_gamma = llr_b.plot_gamma(ax_gamma)
ax_gamma = llr_bm1.plot_gamma(ax_gamma)
ax_gamma = llr_bm2.plot_gamma(ax_gamma)
_ = ax_gamma.legend()

ax_l = None
ax_l = llr_b.plot_l(ax_l)
ax_l = llr_bm1.plot_l(ax_l)
ax_l = llr_bm2.plot_l(ax_l)
_ = ax_l.legend()

ax_cl = None
ax_cl = llr_b.plot_cl(ax_cl)
ax_cl = llr_bm1.plot_cl(ax_cl)
ax_cl = llr_bm2.plot_cl(ax_cl)
_ = ax_cl.legend()

ax_ali = None
ax_ali = llr_b.plot_ali(ax_ali)
ax_ali = llr_bm1.plot_ali(ax_ali)
ax_ali = llr_bm2.plot_ali(ax_ali)
_ = ax_ali.legend()

ax_ale = None
ax_ale = llr_b.plot_ale(ax_ale)
ax_ale = llr_bm1.plot_ale(ax_ale)
ax_ale = llr_bm2.plot_ale(ax_ale)
_ = ax_ale.legend()

ax_wi = None
ax_wi = llr_b.plot_wi(ax_wi)
ax_wi = llr_bm1.plot_wi(ax_wi)
ax_wi = llr_bm2.plot_wi(ax_wi)
_ = ax_wi.legend()

ax_di = None
ax_di = llr_bm1.plot_di(ax_di)
ax_di = llr_bm2.plot_di(ax_di)
_ = ax_di.legend()

ax_cdi = None
ax_cdi = llr_bm1.plot_cdi(ax_cdi)
ax_cdi = llr_bm1.plot_cdi(ax_cdi)
ax_cdi = llr_bm2.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

ax_sf = None
ax_sf = llr_b.plot_sf(ax_sf)
ax_sf = llr_bm1.plot_sf(ax_sf)
ax_sf = llr_bm2.plot_sf(ax_sf)
_ = ax_sf.legend()

ax_bm = None
ax_bm = llr_b.plot_bm(ax_bm)
ax_bm = llr_bm1.plot_bm(ax_bm)
ax_bm = llr_bm2.plot_bm(ax_bm)
_ = ax_bm.legend()

#%%
# Stall at AoA of 12 deg
cla_def = ll_b.cla_shp(ll_b.spacing)
cl_max = ll_b.clmax_shp(ll_b.spacing)
al0 = ll_b.al0_shp(ll_b.spacing)
ale2 = llr_bm2.ale_shp(ll_b.spacing)

cla1 = llr_bm1.cla
cla2 = llr_bm2.cla

cl2 = cla2*ale2

cl_def = cla_def*ale2
cl_def[cl_def > cl_max] = cl_max[cl_def > cl_max]

dcl_1 = cl_max - cl2
dcl_1[dcl_1 > 0.0] = 0.0

dcl_2 = cl_def - cl2

ax_ale = None
ax_ale = llr_bm2.plot_ale(ax_ale)
_ = ax_ale.legend()

fig = figure()
ax = fig.gca()
ax.grid(True)
ax.plot(ll_b.y, cla_def*ale2, label='cla_def*ale2')
ax.plot(ll_b.y, cl_def, label='cl_def')
ax.plot(ll_b.y, cl2, label='cl2')
_ = ax.legend()

fig = figure()
ax = fig.gca()
ax.grid(True)
ax.plot(ll_b.y, dcl_1, label='dcl_1')
ax.plot(ll_b.y, dcl_2, label='dcl_2')
_ = ax.legend()

dcla_2 = divide(dcl_2, ale2)

fig = figure()
ax = fig.gca()
ax.grid(True)
ax.plot(ll_b.y, dcla_2, label='dcla')
_ = ax.legend()

cla2 = cla2 + dcla_2
cla2[cla2 > cla_def] = cla_def[cla2 > cla_def]

fig = figure()
ax = fig.gca()
ax.grid(True)
ax.plot(ll_b.y, cla1, label='llr_bm1.cla')
ax.plot(ll_b.y, cla2, label='cla2')
_ = ax.legend()

# dcl2a = cl2 - cl_max
# dcl2a[dcl2a < 0.0] = 0.0

# dcl2b = cl_max - cla_def*ale2
# # dcl2b[dcl2b < 0.0] = 0.0

# dcla = cla_def - llr_bm2.cla
# dcla[dcl2b < 0.0] = 0.0

# fig = figure()
# ax = fig.gca()
# ax.plot(ll_b.y, dcl2a, label='dcl2a')
# ax.plot(ll_b.y, dcl2b, label='dcl2b')
# ax.plot(ll_b.y, dcla, label='dcla')
# _ = ax.legend()
