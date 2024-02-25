#%%
# Import Dependencies
from pyllt import BellShape, EllipticalShape, GeneralShape

#%%
# Compare Shapes
esd = EllipticalShape().normalise_area()
bsd = BellShape().normalise_area()
gsd = GeneralShape([1.0, 0.0, 1/3, 0.0, 1/5, 0.0, 1/7, 0.0, 1/9]).normalise_area()

ax = None
ax = esd.plot(ax=ax, label='Elliptical')
ax = bsd.plot(ax=ax, label='Bell')
ax = gsd.plot(ax=ax, label='General')
_ = ax.legend()

print(f'esd.area = {esd.area}\n')
print(f'bsd.area = {bsd.area}\n')
print(f'gsd.area = {gsd.area}\n')
