#%%
# Import Dependencies
from IPython.display import display_markdown
from numpy import linspace
from pyllt import BellShape, LiftingLine, TaperedShape, ConstantShape

#%%
# Lifting Line
b = 3.749 # m
c_r = 0.4 # m
c_t = 0.1 # m
m = 6.577 # kg
W = m*9.0665 # N
vel = 12.3 # m/s
rho = 1.225 # kg/m^3

c_shp = TaperedShape(c_r, c_t)

ll_p = LiftingLine('Prandtl D2', b, c_shp)
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

alg = [8.3274, 8.5524, 8.7259, 8.8441, 8.9030, 8.8984, 8.8257, 8.6801,
       8.4565, 8.1492, 7.7522, 7.2592, 6.6634, 5.9579, 5.1362, 4.1927,
       3.1253, 1.9394, 0.6589, -0.6417, -1.6726]

num = len(alg)

y = [float(i)*b/2/(num-1) for i in range(num)]

ax_alg = None
ax_alg = ll_p.plot_alg(ax_alg)
ax_alg.plot(y, alg, label='Actual')
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

ax_l = None
ax_l = llr_p.plot_l(ax_l)
_ = ax_l.legend()

ax_cl = None
ax_cl = llr_p.plot_cl(ax_cl)
_ = ax_cl.legend()

ax_ali = None
ax_ali = llr_p.plot_ali(ax_ali)
_ = ax_ali.legend()

ax_ale = None
ax_ale = llr_p.plot_ale(ax_ale)
_ = ax_ale.legend()

ax_wi = None
ax_wi = llr_p.plot_wi(ax_wi)
_ = ax_wi.legend()

ax_di = None
ax_di = llr_p.plot_di(ax_di)
_ = ax_di.legend()

ax_cdi = None
ax_cdi = llr_p.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

ax_sf = None
ax_sf = llr_p.plot_sf(ax_sf)
_ = ax_sf.legend()

ax_bm = None
ax_bm = llr_p.plot_bm(ax_bm)
_ = ax_bm.legend()

#%%
# Stall Conditions
al_max = 8.0

clmax_shp = ConstantShape(1.4)
clmin_shp = ConstantShape(-1.4)

ll_p.clmax_shp = clmax_shp
ll_p.clmin_shp = clmin_shp

llr_pm = ll_p.return_result_alpha(al_max, name='Prandtl D2 CL Max',
                                  vel=vel, rho=rho)
llr_pm.stall()
display_markdown(llr_pm)

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
ax_cla = llr_pm.plot_cla(ax_cla)
_ = ax_cla.legend()

#%%
# Plot Results
ax_gamma = None
ax_gamma = llr_p.plot_gamma(ax_gamma)
_ = ax_gamma.legend()

ax_l = None
ax_l = llr_p.plot_l(ax_l)
_ = ax_l.legend()

ax_cl = None
ax_cl = llr_p.plot_cl(ax_cl)
_ = ax_cl.legend()

ax_ali = None
ax_ali = llr_p.plot_ali(ax_ali)
_ = ax_ali.legend()

ax_ale = None
ax_ale = llr_p.plot_ale(ax_ale)
_ = ax_ale.legend()

ax_wi = None
ax_wi = llr_p.plot_wi(ax_wi)
_ = ax_wi.legend()

ax_di = None
ax_di = llr_p.plot_di(ax_di)
_ = ax_di.legend()

ax_cdi = None
ax_cdi = llr_p.plot_cdi(ax_cdi)
_ = ax_cdi.legend()

ax_sf = None
ax_sf = llr_p.plot_sf(ax_sf)
_ = ax_sf.legend()

ax_bm = None
ax_bm = llr_p.plot_bm(ax_bm)
_ = ax_bm.legend()

#%%
# Create Polars
al_degs = linspace(-32.0, 18.0, 101)

llp_em = ll_p.return_polar(al_degs).stall()

axp_CL = None
axp_CL = llp_em.plot_CL(axp_CL)
_ = axp_CL.legend()

axp_CDi = None
axp_CDi = llp_em.plot_CDi(axp_CDi)
_ = axp_CDi.legend()

axp_pol = None
axp_pol = llp_em.plot_polar(axp_pol)
_ = axp_pol.legend()
