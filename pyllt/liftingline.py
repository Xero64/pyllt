from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from matplotlib.pyplot import figure
from numpy import (absolute, asarray, degrees, divide, fill_diagonal, full, pi,
                   radians, zeros)
from numpy.linalg import norm, solve
from py2md.classes import MDHeading, MDReport, MDTable

from .shape import (ConstantShape, EllipticalShape, GeneralShape,
                    InducedAngleShape)
from .spacing import CosineSpacing

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy.typing import NDArray

    from .shape import Shape
    from .spacing import Spacing


class LiftingLineResult():

    liftingline: 'LiftingLine' = None
    al_deg: float = None
    vel: float = None
    rho: float = None
    name: str = None
    _num: int = None
    _al_rad: float = None
    _al_shp: 'Shape' = None
    _spacing: 'Spacing' = None
    _s: 'NDArray' = None
    _th: 'NDArray' = None
    _cla: 'NDArray' = None
    _cos_n_th: dict[int, 'NDArray'] = None
    _sin_n_th: dict[int, 'NDArray'] = None
    _y: 'NDArray' = None
    _An: GeneralShape = None
    _n: 'NDArray' = None
    _q: float = None
    _gamma_shp: 'Shape' = None
    _gamma: 'NDArray' = None
    _ali_shp: 'Shape' = None
    _ali: 'NDArray' = None
    _wi_shp: 'Shape' = None
    _wi: 'NDArray' = None
    _ale_shp: 'Shape' = None
    _ale: 'NDArray' = None
    _l_shp: 'Shape' = None
    _l: 'NDArray' = None
    _cl_shp: 'Shape' = None
    _cl: 'NDArray' = None
    _di_shp: 'Shape' = None
    _di: 'NDArray' = None
    _cdi_shp: 'Shape' = None
    _cdi: 'NDArray' = None
    _d: 'NDArray' = None
    _CL: float = None
    _L: float = None
    _CDi: float = None
    _CD0: float = None
    _CD: float = None
    _Di: float = None
    _D0: float = None
    _D: float = None
    _delta: float = None
    _e: float = None
    _dcl: 'NDArray' = None
    _norm_dcl: float = None
    _sf: 'NDArray' = None
    _bm: 'NDArray' = None
    _bmr: float = None
    _pwr: float = None
    _LoD: float = None

    def __init__(self, liftingline: 'LiftingLine',
                 al_deg: float, **kwargs: dict[str, Any]) -> None:

        self.liftingline = liftingline
        self.al_deg = al_deg

        self.vel = kwargs.get('vel', 1.0)
        self.rho = kwargs.get('rho', 1.0)
        self.name = kwargs.get('name', self.liftingline.name)
        self._num = kwargs.get('num', self.liftingline.num)

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def num(self) -> int:
        if self._num is None:
            self._num = self.liftingline.num
        return self._num

    @num.setter
    def num(self, num: int) -> None:
        self.reset()
        self._num = num

    @property
    def al_rad(self) -> float:
        if self._al_rad is None:
            self._al_rad = radians(self.al_deg)
        return self._al_rad

    @al_rad.setter
    def al_rad(self, al_rad: float) -> None:
        self.reset()
        self._al_rad = al_rad

    @property
    def al_shp(self) -> 'Shape':
        if self._al_shp is None:
            self._al_shp = ConstantShape(self.al_rad)
        return self._al_shp

    @property
    def spacing(self) -> 'Spacing':
        if self._spacing is None:
            self._spacing = CosineSpacing(self.num)
        return self._spacing

    @spacing.setter
    def spacing(self, spacing: 'Spacing') -> None:
        self.reset()
        self.num = spacing.num
        self._spacing = spacing

    @property
    def y(self) -> 'NDArray':
        if self._y is None:
            self._y = self.liftingline.b/2*self.spacing.s
        return self._y

    @property
    def cla(self) -> 'NDArray':
        if self._cla is None:
            self._cla = self.liftingline.cla_shp(self.liftingline.spacing)
        return self._cla

    @cla.setter
    def cla(self, cla: 'NDArray') -> None:
        if cla.size != self.liftingline.num:
            raise ValueError('cla is not correct size.')
        self._cla = cla

    @property
    def An(self) -> GeneralShape:
        if self._An is None:
            if self.stall:
                Ana, An0 = self.liftingline.solve(self.cla)
            else:
                Ana = self.liftingline.Ana
                An0 = self.liftingline.An0
            self._An = Ana*self.al_rad + An0
        return self._An

    @An.setter
    def An(self, An: GeneralShape) -> None:
        self.reset()
        self._An = An

    @property
    def n(self) -> 'NDArray':
        if self._n is None:
            self._n = asarray(self.An.n)
        return self._n

    @property
    def q(self) -> float:
        if self._q is None:
            self._q = 0.5*self.rho*self.vel**2
        return self._q

    @q.setter
    def q(self, q: float) -> None:
        self._q = q

    @property
    def gamma_shp(self) -> 'Shape':
        if self._gamma_shp is None:
            self._gamma_shp = 2*self.liftingline.b*self.vel*self.An
        return self._gamma_shp

    @property
    def gamma(self) -> 'NDArray':
        if self._gamma is None:
            self._gamma = self.gamma_shp(self.spacing)
        return self._gamma

    @property
    def ali_shp(self) -> InducedAngleShape:
        if self._ali_shp is None:
            self._ali_shp = InducedAngleShape(self.An.An)
        return self._ali_shp

    @property
    def ali(self) -> 'NDArray':
        if self._ali is None:
            self._ali = self.ali_shp(self.spacing)
        return self._ali

    @property
    def wi_shp(self) -> 'Shape':
        if self._wi_shp is None:
            self._wi_shp = self.vel*self.ali_shp
        return self._wi_shp

    @property
    def wi(self) -> 'NDArray':
        if self._wi is None:
            self._wi = self.wi_shp(self.spacing)
        return self._wi

    @property
    def ale_shp(self) -> 'Shape':
        if self._ale_shp is None:
            al_rad_shp = ConstantShape(self.al_rad)
            alg_shp = self.liftingline.alg_shp
            al0_shp = self.liftingline.al0_shp
            self._ale_shp = al_rad_shp + alg_shp - al0_shp - self.ali_shp
        return self._ale_shp

    @property
    def ale(self) -> 'NDArray':
        if self._ale is None:
            self._ale = self.ale_shp(self.spacing)
        return self._ale

    @property
    def l_shp(self) -> 'Shape':
        if self._l_shp is None:
            self._l_shp = self.rho*self.vel*self.gamma_shp
        return self._l_shp

    @property
    def l(self) -> 'NDArray':
        if self._l is None:
            self._l = self.l_shp(self.spacing)
        return self._l

    @property
    def cl_shp(self) -> 'Shape':
        if self._cl_shp is None:
            self._cl_shp = self.l_shp/self.liftingline.c_shp/self.q
        return self._cl_shp

    @property
    def cl(self) -> float:
        if self._cl is None:
            self._cl = self.cl_shp(self.spacing)
        return self._cl

    @property
    def di_shp(self) -> 'Shape':
        if self._di_shp is None:
            self._di_shp = self.l_shp*self.ali_shp
        return self._di_shp

    @property
    def di(self) -> 'NDArray':
        if self._di is None:
            self._di = self.di_shp(self.spacing)
        return self._di

    @property
    def cdi_shp(self) -> 'Shape':
        if self._cdi_shp is None:
            self._cdi_shp = self.di_shp/self.liftingline.c_shp/self.q
        return self._cdi_shp

    @property
    def cdi(self) -> float:
        if self._cdi is None:
            self._cdi = self.cdi_shp(self.spacing)
        return self._cdi

    @property
    def sf(self) -> 'NDArray':
        if self._sf is None:
            sfoqSar = zeros(self.spacing.th.size)
            for n, An in self.An.items():
                if n == 1:
                    sfoqSar += An*(self.spacing.sinth*self.spacing.costh + pi - self.spacing.th)
                else:
                    n2m1 = n**2 - 1
                    cos_n_th = self.spacing.cos_n_th(n)
                    sin_n_th = self.spacing.sin_n_th(n)
                    term = n*self.spacing.sinth*cos_n_th - sin_n_th*self.spacing.costh
                    sfoqSar += 2*An/n2m1*term
            q = self.q
            Sref = self.liftingline.area
            ar = self.liftingline.ar
            self._sf = sfoqSar*q*Sref*ar
            check = self.spacing.s > 0.0
            self._sf[check] -= self.L
        return self._sf

    @property
    def bm(self) -> 'NDArray':
        if self._bm is None:
            bmoqSarb = zeros(self.spacing.th.size)
            pimth = pi - self.spacing.th
            for n, An in self.An.items():
                if n == 1:
                    term = pimth*self.spacing.costh/2 + 3*self.spacing.sinth/8 + self.spacing.sin_n_th(3)/24
                    bmoqSarb += An*term
                elif n == 2:
                    term = pimth/4 - self.spacing.sin_n_th(2)/6 + self.spacing.sin_n_th(4)/48
                    bmoqSarb += An*term
                else:
                    nm1 = n - 1
                    nm2 = n - 2
                    np1 = n + 1
                    np2 = n + 2
                    term = self.spacing.sin_n_th(nm2)/nm2/nm1/4
                    term -= self.spacing.sin_n_th(n)/nm1/np1/2
                    term += self.spacing.sin_n_th(np2)/np2/np1/4
                    bmoqSarb += An*term
            q = self.q
            Sref = self.liftingline.area
            ar = self.liftingline.ar
            b = self.liftingline.b
            self._bm = bmoqSarb*q*Sref*ar*b
            check = self.spacing.s > 0.0
            self._bm[check] -= self.L*self.y[check]
        return self._bm

    @property
    def CL(self) -> float:
        if self._CL is None:
            self._CL = pi*self.liftingline.ar*self.An[0]
        return self._CL

    @CL.setter
    def CL(self, CL: float) -> None:
        self._CL = CL

    @property
    def L(self) -> float:
        if self._L is None:
            self._L = self.q*self.liftingline.area*self.CL
        return self._L

    @L.setter
    def L(self, L: float) -> None:
        self._L = L

    @property
    def delta(self) -> float:
        if self._delta is None:
            if self.CL == 0.0:
                self._delta = 0.0
            else:
                self._delta = (self.An[1:]**2*self.n[1:]).sum()/self.An[0]**2
        return self._delta

    @property
    def e(self) -> float:
        if self._e is None:
            self._e = 1.0/(1.0 + self.delta)
        return self._e

    @property
    def CDi(self) -> float:
        if self._CDi is None:
            self._CDi = self.CL**2/(pi*self.liftingline.ar*self.e)
        return self._CDi

    @property
    def CD0(self) -> float:
        if self._CD0 is None:
            self._CD0 = self.liftingline.cd0
        return self._CD0

    @property
    def CD(self) -> float:
        if self._CD is None:
            self._CD = self.CD0 + self.CDi
        return self._CD

    @property
    def Di(self) -> float:
        if self._Di is None:
            self._Di = self.q*self.liftingline.area*self.CDi
        return self._Di

    @property
    def D0(self) -> float:
        if self._D0 is None:
            self._D0 = self.q*self.liftingline.area*self.CD0
        return self._D0

    @property
    def D(self) -> float:
        if self._D is None:
            self._D = self.D0 + self.Di
        return self._D

    @property
    def LoD(self) -> float:
        if self._LoD is None:
            if self.D != 0.0:
                self._LoD = self.L/self.D
            else:
                if self.L == 0.0:
                    self._LoD = 0.0
                else:
                    self._LoD = float('inf')
        return self._LoD

    @property
    def pwr(self) -> float:
        if self._pwr is None:
            self._pwr = self.vel*self.D
        return self._pwr

    @property
    def bmr(self) -> float:
        if self._bmr is None:
            ind = self.liftingline.num//2 + 1
            self._bmr = self.bm[ind]
        return self._bmr

    def stall(self, norm_tol: float = 1e-6, ale_tol: float = 1e-6,
              num_iter: int = 100, display: bool = False) -> 'LiftingLineResult':

        ll = self.liftingline

        cla = self.cla

        clmax = ll.clmax_shp(ll.spacing)
        clmin = ll.clmin_shp(ll.spacing)

        count = 0

        norm_dcl = float('inf')

        while norm_dcl >= norm_tol:

            Ana, An0 = ll.solve(cla)
            An = Ana*self.al_rad + An0

            gamma_shp = 2*ll.b*self.vel*An
            l_shp = self.rho*self.vel*gamma_shp
            cl_shp = l_shp/ll.c_shp/self.q

            ali_shp = An*self.n/EllipticalShape()
            ale_shp = self.al_shp + ll.alg_shp - ll.al0_shp - ali_shp

            cl = cl_shp(ll.spacing)
            ale = ale_shp(ll.spacing)

            ale_chk = absolute(ale) > ale_tol

            dclp = cl - clmax
            dclp_chk = dclp <= 0.0
            dclp[dclp_chk] = 0.0

            dcln = cl - clmin
            dcln_chk = dcln >= 0.0
            dcln[dcln_chk] = 0.0

            dcl = dclp + dcln

            norm_dcl = norm(dcl)

            dcla = zeros(cla.shape)

            divide(dcl, ale, where=ale_chk, out=dcla)

            cla -= dcla

            count = count + 1

            if display:
                print(f'Iteration {count:d}: Norm dcl = {norm_dcl}')

            if count > num_iter:
                err = f'Failed to converge in {num_iter:d} iterations.'
                raise ValueError(err)

        self.reset()
        self.An = An
        self.cla = cla

        if display:
            print(f'Converged in {count:d} iterations.')

        return self

    def set_lift_distribution(self, L: float, An: GeneralShape) -> None:
        Sref = self.liftingline.area
        ar = self.liftingline.ar
        factor = L/self.q/Sref/ar/An.area/2
        self.An = An*factor

    def set_lifting_line_twist(self) -> None:
        b = self.liftingline.b
        c_shp = self.liftingline.c_shp
        al0_shp = self.liftingline.al0_shp
        ali_shp = InducedAngleShape(self.An.An)
        temp = self.An*2/pi*b/c_shp
        alg_shp = temp + al0_shp + ali_shp
        self.liftingline.alg_shp = alg_shp
        self.liftingline.reset()
        al_rad = (self.CL - self.liftingline.CL0)/self.liftingline.CLa
        self.al_deg = degrees(al_rad)
        self.reset(excl = ['_An'])

    def minimum_induced_drag_optimum(self, Lspec: float, Mspec: float = None,
                                     num: int = 3, display: bool = False) -> 'LiftingLineResult':

        dCDdA = asarray([n for n in range(1, num+1, 2)])
        numsym = dCDdA.size

        qS = self.q*self.liftingline.area
        ar = self.liftingline.ar
        b = self.liftingline.b
        Amat = zeros((numsym, numsym))
        fill_diagonal(Amat, dCDdA)
        Bmat = zeros(numsym)
        Amat[0, 0] = 1.0
        Bmat[0] = Lspec/qS/pi/ar
        if Mspec is not None:
            CMopiar = Mspec/qS/b/ar
            dCMdA = [1/3]
            for n in range(3, num+1, 2):
                dCMdAn = 1/(n+2)/(n+1)/4 + 1/(n+1)/(n-1)/2 + 1/(n-1)/(n-2)/4
                if (n - 1) % 4 == 0:
                    dCMdAn = -dCMdAn
                dCMdA.append(dCMdAn)
            Amat[1, :] = asarray(dCMdA)
            Bmat[1] = CMopiar
        if display:
            print(f'Amat = \n{Amat}\n')
            print(f'Bmat = {Bmat}\n')
        Cmat = solve(Amat, Bmat)
        if display:
            print(f'Cmat = {Cmat}\n')
        An = zeros(num)
        An[::2] = Cmat
        self.reset()
        self.An = GeneralShape(An)
        if display:
            print(f'An = {An}\n')
        return self

    def plot_gamma(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Circulation Distribution - $\Gamma$ (m/s)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.gamma, **kwargs)
        return ax

    def plot_ali(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Angle Distribution - $\alpha_i$ (deg)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, degrees(self.ali), **kwargs)
        return ax

    def plot_wi(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Wash Distribution $w_i$ (m/s)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.wi, **kwargs)
        return ax

    def plot_ale(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Effective Angle Distribution $\alpha_{eff}$ (deg)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, degrees(self.ale), **kwargs)
        return ax

    def plot_cla(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Slope - $c_{la}$ (1/rad)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.liftingline.y, self.cla, **kwargs)
        return ax

    def plot_l(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Distribution - l (N/m)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.l, **kwargs)
        return ax

    def plot_cl(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Distribution - $c_l$')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.cl, **kwargs)
        return ax

    def plot_di(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Drag Distribution - $d_i$ (N/m)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.di, **kwargs)
        return ax

    def plot_cdi(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Drag Coefficient - $c_{di}$')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.cdi, **kwargs)
        return ax

    def plot_sf(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Shear Force - $SF$ (N)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.sf, **kwargs)
        return ax

    def plot_bm(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Bending Moment - $BM$ (N.m)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.bm, **kwargs)
        return ax

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        heading = MDHeading(f'Lifting Line Result for {self.name:s}', 2)
        report.add_object(heading)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('Alpha (deg)', '.3f', data=[self.al_deg])
        table.add_column('V (m/s)', '.3f', data=[self.vel])
        table.add_column('&rho; (kg/m<sup>3</sup>)', '.3f', data=[self.rho])
        table.add_column('q (Pa)', '.3f', data=[self.q])
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('C<sub>L</sub>', '.6f', data=[self.CL])
        table.add_column('C<sub>Di</sub>', '.6f', data=[self.CDi])
        table.add_column('C<sub>D0</sub>', '.6f', data=[self.CD0])
        table.add_column('C<sub>D</sub>', '.6f', data=[self.CD])
        table.add_column('&delta;', '.6f', data=[self.delta])
        table.add_column('e', '.6f', data=[self.e])
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('L (N)', '.6f', data=[self.L])
        table.add_column('D<sub>i</sub> (N)', '.6f', data=[self.Di])
        table.add_column('D<sub>0</sub> (N)', '.6f', data=[self.D0])
        table.add_column('D (N)', '.6f', data=[self.D])
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('BM Root (N.m)', '.6f', data=[self.bmr])
        table.add_column('Power (W)', '.3f', data=[self.pwr])
        table.add_column('Lift to Drag Ratio', '.3f', data=[self.LoD])
        report.add_object(table)
        return report

    def _repr_markdown_(self) -> MDReport:
        return self.to_mdobj()._repr_markdown_()


class LiftingLinePolar():

    name: str = None
    al_degs: 'NDArray' = None
    liftingline: 'LiftingLine' = None
    _num: int = None
    _spacing: 'Spacing' = None
    _results: dict[float, LiftingLineResult] = None
    _CL: 'NDArray' = None
    _CDi: 'NDArray' = None
    _CD: 'NDArray' = None

    def __init__(self, liftingline: 'LiftingLine', al_degs: 'NDArray',
                 **kwargs: dict[str, Any]) -> None:

        self.liftingline = liftingline
        self.al_degs = al_degs
        self.name = kwargs.get('name', self.liftingline.name)
        self._num = kwargs.get('num', self.liftingline.num)

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def num(self) -> int:
        if self._num is None:
            self._num = self.liftingline.num
        return self._num

    @num.setter
    def num(self, num: int) -> None:
        self.reset()
        self.num = num

    @property
    def spacing(self) -> 'Spacing':
        if self._spacing is None:
            self._spacing = CosineSpacing(self.num)
        return self._spacing

    @property
    def results(self) -> dict[float, LiftingLineResult]:
        if self._results is None:
            self._results = {}
            for al_deg in self.al_degs:
                name = f'{self.name:s} {al_deg:.3f}'
                result = self.liftingline.return_result_alpha(al_deg, name=name)
                result.spacing = self.spacing
                self._results[al_deg] = result
        return self._results

    @property
    def CL(self) -> 'NDArray':
        if self._CL is None:
            self._CL = asarray([result.CL for result in self.results.values()])
        return self._CL

    @property
    def CDi(self) -> 'NDArray':
        if self._CDi is None:
            self._CDi = asarray([result.CDi for result in self.results.values()])
        return self._CDi

    @property
    def CD(self) -> 'NDArray':
        if self._CD is None:
            self._CD = asarray([result.CD for result in self.results.values()])
        return self._CD

    def stall(self, norm_tol: float = 1e-6, ale_tol: float = 1e-6,
              num_iter: int = 100, display: bool = False) -> 'LiftingLinePolar':

        for al_deg, result in self.results.items():
            result.stall(norm_tol, ale_tol, num_iter, display)
            self._results[al_deg] = result

        self.reset(excl=['_results'])

        return self

    def plot_CL(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Angle of Attack - $\alpha$ (deg)')
            ax.set_ylabel(r'Lift Coefficient - $C_L$')
        ax.plot(self.al_degs, self.CL, label=self.name)
        return ax

    def plot_CDi(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Angle of Attack - $\alpha$ (deg)')
            ax.set_ylabel(r'Induced Drag Coefficient - $C_{Di}$')
        ax.plot(self.al_degs, self.CDi, label=self.name)
        return ax

    def plot_polar(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Drag Coefficient - $C_D$')
            ax.set_ylabel(r'Lift Coefficient - $C_L$')
            ax.set_title(r'Drag Polar')
        ax.plot(self.CD, self.CL, label=self.name)
        return ax


class LiftingLineLevelFlight():

    name: str = None
    lift: float = None
    rho: float = None
    vels: 'NDArray' = None
    liftingline: 'LiftingLine' = None
    _num: int = None
    _spacing: 'Spacing' = None
    _results: dict[float, LiftingLineResult] = None
    _CL: 'NDArray' = None
    _CD: 'NDArray' = None
    _LoD: 'NDArray' = None
    _pwr: 'NDArray' = None

    def __init__(self, liftingline: 'LiftingLine', lift: float, rho: float,
                 vels: 'NDArray', **kwargs: dict[str, Any]) -> None:

        self.liftingline = liftingline
        self.lift = lift
        self.rho = rho
        self.vels = vels
        self.name = kwargs.get('name', self.liftingline.name)
        self._num = kwargs.get('num', self.liftingline.num)

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def num(self) -> int:
        if self._num is None:
            self._num = self.liftingline.num
        return self._num

    @num.setter
    def num(self, num: int) -> None:
        self.reset()
        self.num = num

    @property
    def spacing(self) -> 'Spacing':
        if self._spacing is None:
            self._spacing = CosineSpacing(self.num)
        return self._spacing

    @property
    def results(self) -> dict[float, LiftingLineResult]:
        if self._results is None:
            self._results = {}
            for vel in self.vels:
                name = f'{self.name:s} {vel:.3f}'
                result = self.liftingline.return_result_L(self.lift, vel=vel,
                                                          rho=self.rho, name=name)
                result.spacing = self.spacing
                self._results[vel] = result
        return self._results

    @property
    def CL(self) -> 'NDArray':
        if self._CL is None:
            self._CL = asarray([result.CL for result in self.results.values()])
        return self._CL

    @property
    def CD(self) -> 'NDArray':
        if self._CD is None:
            self._CD = asarray([result.CD for result in self.results.values()])
        return self._CD

    @property
    def LoD(self) -> 'NDArray':
        if self._LoD is None:
            self._LoD = asarray([result.LoD for result in self.results.values()])
        return self._LoD

    @property
    def pwr(self) -> 'NDArray':
        if self._pwr is None:
            self._pwr = asarray([result.pwr for result in self.results.values()])
        return self._pwr

    def stall(self, norm_tol: float = 1e-6, ale_tol: float = 1e-6,
              num_iter: int = 100, display: bool = False) -> 'LiftingLineLevelFlight':

        for vel, result in self.results.items():
            result.stall(norm_tol, ale_tol, num_iter, display)
            self._results[vel] = result

        self.reset(excl=['_results'])

        return self

    def plot_LoD(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Velocity - V (m/s)')
            ax.set_ylabel(r'Lift to Drag Ratio - $\frac{L}{D}$')
        ax.plot(self.vels, self.LoD, label=self.name)
        return ax

    def plot_pwr(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Velocity - V (m/s)')
            ax.set_ylabel(r'Power - P (W)')
        ax.plot(self.vels, self.pwr, label=self.name)
        return ax

    def plot_CL(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Velocity - V (m/s)')
            ax.set_ylabel(r'Lift Coefficient - $C_L$')
        ax.plot(self.vels, self.CL, label=self.name)
        return ax

    def plot_CD(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Velocity - V (m/s)')
            ax.set_ylabel(r'Drag Coefficient - $C_D$')
        ax.plot(self.vels, self.CD, label=self.name)
        return ax

    def plot_custom(self, attr: str, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Velocity - V (m/s)')
            ax.set_ylabel(attr)
        vals = asarray([getattr(result, attr) for result in self.results.values()])
        ax.plot(self.vels, vals, label=self.name)
        return ax


class LiftingLine():

    name: str = None
    b: float = None
    c_shp: 'Shape' = None
    alg_shp: 'Shape' = None
    al0_shp: 'Shape' = None
    cla_shp: 'Shape' = None
    num: int = None
    clmax_shp: 'Shape' = None
    clmin_shp: 'Shape' = None
    cd0: float = None
    _clmax: 'NDArray' = None
    _clmin: 'NDArray' = None
    _area: float = None
    _mac: float = None
    _ar: float = None
    _spacing: 'Spacing' = None
    _y: 'NDArray' = None
    _c: 'NDArray' = None
    _alg: 'NDArray' = None
    _al0: 'NDArray' = None
    _cla: 'NDArray' = None
    _Ana: GeneralShape = None
    _An0: GeneralShape = None
    _CLa: float = None
    _CL0: float = None
    _aL0: float = None

    def __init__(self, name: str, b: float, c_shp: 'Shape', /,
                 alg_shp: 'Shape' = ConstantShape(0.0),
                 al0_shp: 'Shape' = ConstantShape(0.0),
                 cla_shp: 'Shape' = ConstantShape(2*pi),
                 cd0: float = 0.0, num: int = 101) -> None:
        self.name = name
        self.b = b
        self.c_shp = c_shp
        self.alg_shp = alg_shp
        self.al0_shp = al0_shp
        self.cla_shp = cla_shp
        self.cd0 = cd0
        self.num = 2*(num//2) + 1

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def area(self) -> float:
        if self._area is None:
            self._area = self.c_shp.area*self.b/2
        return self._area

    @area.setter
    def area(self, area: float) -> None:
        self._area = area

    @property
    def mac(self) -> float:
        if self._mac is None:
            self._mac = self.area/self.b
        return self._mac

    @property
    def ar(self) -> float:
        if self._ar is None:
            self._ar = self.b**2/self.area
        return self._ar

    @ar.setter
    def ar(self, ar: float) -> None:
        self._ar = ar

    @property
    def spacing(self) -> 'Spacing':
        if self._spacing is None:
            self._spacing = CosineSpacing(self.num, sol=True)
        return self._spacing

    @property
    def y(self) -> 'NDArray':
        if self._y is None:
            self._y = self.b/2*self.spacing.s
        return self._y

    @property
    def c(self) -> 'NDArray':
        if self._c is None:
            self._c = self.c_shp(self.spacing)
        return self._c

    @property
    def alg(self) -> 'NDArray':
        if self._alg is None:
            self._alg = self.alg_shp(self.spacing)
        return self._alg

    @property
    def al0(self) -> 'NDArray':
        if self._al0 is None:
            self._al0 = self.al0_shp(self.spacing)
        return self._al0

    @property
    def cla(self) -> 'NDArray':
        if self._cla is None:
            self._cla = self.cla_shp(self.spacing)
        return self._cla

    @cla.setter
    def cla(self, cla: 'NDArray') -> None:
        if cla.size != self.num:
            raise ValueError('cla is not correct size.')
        self._cla = cla

    @property
    def clmax(self) -> 'NDArray':
        if self.clmax_shp is None:
            self._clmax = full(self.num, float('inf'))
        else:
            self._clmax = self.clmax_shp(self.spacing)
        return self._clmax

    @property
    def clmin(self) -> 'NDArray':
        if self.clmin_shp is None:
            self._clmin = full(self.num, -float('inf'))
        else:
            self._clmin = self.clmin_shp(self.spacing)
        return self._clmin

    def solve(self, cla: 'NDArray | None' = None) -> tuple[GeneralShape,
                                                           GeneralShape] | None:

        if cla is None:
            claco4b = self.cla*self.c/4/self.b
        else:
            claco4b = cla*self.c/4/self.b

        lhs = zeros((self.num, self.num))
        rhs = zeros((self.num, 2))

        sinth = self.spacing.sinth

        for n in range(1, self.num+1):
            sin_nth = self.spacing.sin_n_th(n)
            n_claco8s = n*claco4b
            lhs[:, n-1] = sin_nth*(sinth + n_claco8s)

        rhs[:, 0] = claco4b*sinth
        rhs[:, 1] = rhs[:, 0]*(self.alg - self.al0)

        res = solve(lhs, rhs)

        if cla is None:

            self.Ana = GeneralShape(res[:, 0])
            self.An0 = GeneralShape(res[:, 1])

        else:

            Ana = GeneralShape(res[:, 0])
            An0 = GeneralShape(res[:, 1])

            return Ana, An0

    @property
    def Ana(self) -> GeneralShape:
        if self._Ana is None:
            self.solve()
        return self._Ana

    @Ana.setter
    def Ana(self, Ana: GeneralShape) -> None:
        if Ana.size != self.num:
            raise ValueError('Ana is not correct size.')
        self._Ana = Ana

    @property
    def An0(self) -> GeneralShape:
        if self._An0 is None:
            self.solve()
        return self._An0

    @An0.setter
    def An0(self, An0: GeneralShape) -> None:
        if An0.size != self.num:
            raise ValueError('An0 is not correct size.')
        self._An0 = An0

    @property
    def CLa(self) -> float:
        if self._CLa is None:
            self._CLa = pi*self.ar*self.Ana[0]
        return self._CLa

    @property
    def CL0(self) -> float:
        if self._CL0 is None:
            self._CL0 = pi*self.ar*self.An0[0]
        return self._CL0

    @property
    def aL0(self) -> float:
        if self._aL0 is None:
            self._aL0 = -self.CL0/self.CLa
        return self._aL0

    def copy(self, incl_attr: bool = False) -> 'LiftingLine':
        ll = LiftingLine(self.name, self.c_shp, self.alg_shp, self.al0_shp,
                         self.cla_shp, self.b, self.num)
        ll.clmax_shp = self.clmax_shp
        ll.clmin_shp = self.clmin_shp
        if incl_attr:
            for attr in self.__dict__:
                ll.__dict__[attr] = self.__dict__[attr]
        return ll

    def plot_c(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Chord - c (m)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.c, **kwargs)
        return ax

    def plot_al0(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Zero Lift Angle of Attack -w $\alpha_{l0}$ (deg)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, degrees(self.al0), **kwargs)
        return ax

    def plot_alg(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Spanwise Twist - $\alpha_g$ (deg)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, degrees(self.alg), **kwargs)
        return ax

    def plot_al(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'$\alpha_g$ & $\alpha_{l0}$ (deg)')
        ax.plot(self.y, degrees(self.alg), label='Geometric Twist')
        ax.plot(self.y, degrees(self.al0), label='Aerodynamic Twist')
        return ax

    def plot_cla(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Slope - $c_{la}$ (1/rad)')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.cla, **kwargs)
        return ax

    def plot_An0(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Zero Lift Distribution Shape - $A_{n0}$')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.An0(self.spacing), **kwargs)
        return ax

    def plot_Ana(self, ax: 'Axes' = None, **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Distribution Shape Variation - $A_{na}$')
        kwargs.setdefault('label', self.name)
        ax.plot(self.y, self.Ana(self.spacing), **kwargs)
        return ax

    def return_result_alpha(self, al_deg: float, **kwargs: dict[str, Any]) -> 'LiftingLineResult':
        res = LiftingLineResult(self, al_deg, **kwargs)
        return res

    def return_result_CL(self, CL: float, **kwargs: dict[str, Any]) -> 'LiftingLineResult':
        al_rad = (CL - self.CL0)/self.CLa
        al_deg = degrees(al_rad)
        res = self.return_result_alpha(al_deg, **kwargs)
        res.al_rad = al_rad
        res.CL = CL
        return res

    def return_result_L(self, L: float, **kwargs: dict[str, Any]) -> 'LiftingLineResult':
        rho = kwargs.get('rho', 1.0)
        vel = kwargs.get('vel', 1.0)
        q = rho*vel**2/2
        res = self.return_result_CL(L/self.area/q, **kwargs)
        res.q = q
        res.L = L
        return res

    def return_polar(self, al_degs: Iterable[float],
                     **kwargs: dict[str, Any]) -> LiftingLinePolar:
        al_degs = asarray(al_degs)
        ll_pol = LiftingLinePolar(self, al_degs, **kwargs)
        return ll_pol

    def return_level_flight(self, lift: float, rho: float, vels: Iterable[float],
                            **kwargs: dict[str, Any]) -> LiftingLineLevelFlight:
        vels = asarray(vels)
        ll_lev = LiftingLineLevelFlight(self, lift, rho, vels, **kwargs)
        return ll_lev

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        heading = MDHeading(f'Lifting Line for {self.name:s}', 2)
        report.add_object(heading)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('Span (m)', '.3f', data=[self.b])
        table.add_column('MAC (m)', '.3f', data=[self.mac])
        table.add_column('Area (m<sup>2</sup>)', '.3f', data=[self.area])
        table.add_column('Aspect Ratio', '.3f', data=[self.ar])
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('C<sub>L&alpha;</sub>', '.3f', data=[self.CLa])
        table.add_column('C<sub>L0</sub>', '.3f', data=[self.CL0])
        table.add_column('&alpha;<sub>L0</sub> (deg)', '.3f', data=[degrees(self.aL0)])
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('Number of Terms', 'd', data=[self.num])
        report.add_object(table)
        return report

    def _repr_markdown_(self) -> MDReport:
        return self.to_mdobj()._repr_markdown_()
