from typing import TYPE_CHECKING, Any, Dict, Iterable, Optional, Tuple

from matplotlib.pyplot import figure
from numpy import (absolute, asarray, cos, degrees, divide, full, linspace, pi,
                   radians, sin, zeros)
from numpy.linalg import norm, solve
from py2md.classes import MDHeading, MDReport, MDTable

from .shape import ConstantShape, GeneralShape, EllipticalShape

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy import ndarray

    from .shape import Shape


def CL_fn(ar: float , An: GeneralShape) -> float:
    return pi*ar*An[0]


class LiftingLineResult():

    liftingline: 'LiftingLine' = None
    al_deg: float = None
    vel: float = None
    rho: float = None
    name: str = None
    _al_rad: float = None
    _al_fn: 'Shape' = None
    _An: GeneralShape = None
    _n: 'ndarray' = None
    _q: float = None
    _gamma_fn: 'Shape' = None
    _gamma: 'ndarray' = None
    _ali_fn: 'Shape' = None
    _ali: 'ndarray' = None
    _wi_fn: 'Shape' = None
    _wi: 'ndarray' = None
    _ale_fn: 'Shape' = None
    _ale: 'ndarray' = None
    _l_fn: 'Shape' = None
    _l: 'ndarray' = None
    _cl_fn: 'Shape' = None
    _cl: 'ndarray' = None
    _di_fn: 'Shape' = None
    _di: 'ndarray' = None
    _cdi_fn: 'Shape' = None
    _cdi: float = None
    _CL: float = None
    _L: float = None
    _CDi: float = None
    _Di: float = None
    _delta: float = None
    _e: float = None
    _dcl: 'ndarray' = None
    _norm_dcl: float = None
    _th: 'ndarray' = None
    _cla: 'ndarray' = None
    _cos_n_th: Dict[int, 'ndarray'] = None
    _sin_n_th: Dict[int, 'ndarray'] = None
    _y: 'ndarray' = None
    _sf: 'ndarray' = None
    _bm: 'ndarray' = None
    _bmr: float = None

    def __init__(self, liftingline: 'LiftingLine',
                 al_deg: float, **kwargs: Dict[str, Any]) -> None:

        self.liftingline = liftingline
        self.al_deg = al_deg

        self.vel = kwargs.get('vel', 1.0)
        self.rho = kwargs.get('rho', 1.0)
        self.name = kwargs.get('name', self.liftingline.name)

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def al_rad(self) -> float:
        if self._al_rad is None:
            self._al_rad = radians(self.al_deg)
        return self._al_rad

    @al_rad.setter
    def al_rad(self, al_rad: float) -> None:
        self._al_rad = al_rad

    @property
    def al_fn(self) -> 'Shape':
        if self._al_fn is None:
            self._al_fn = ConstantShape(self.al_rad)
        return self._al_fn

    @property
    def cla(self) -> 'ndarray':
        if self._cla is None:
            self._cla = self.liftingline.cla_fn(self.liftingline.s)
        return self._cla

    @cla.setter
    def cla(self, cla: 'ndarray') -> None:
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
    def n(self) -> 'ndarray':
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
    def gamma_fn(self) -> 'Shape':
        if self._gamma_fn is None:
            self._gamma_fn = 2*self.liftingline.b*self.vel*self.An
        return self._gamma_fn

    @property
    def gamma(self) -> 'ndarray':
        if self._gamma is None:
            self._gamma = 2*self.vel*self.liftingline.b*self.An(self.liftingline.s)
        return self._gamma

    @property
    def ali_fn(self) -> 'Shape':
        if self._ali_fn is None:
            nAn = self.An*self.n
            self._ali_fn = nAn/EllipticalShape()
        return self._ali_fn

    @property
    def ali(self) -> 'ndarray':
        if self._ali is None:
            sinth = sin(self.liftingline.th)
            nAn = self.An*self.n
            self._ali = nAn(self.liftingline.s)/sinth
        return self._ali

    @property
    def wi_fn(self) -> 'Shape':
        if self._wi_fn is None:
            self._wi_fn = self.vel*self.ali_fn
        return self._wi_fn

    @property
    def wi(self) -> 'ndarray':
        if self._wi is None:
            self._wi = self.vel*self.ali
        return self._wi

    @property
    def ale_fn(self) -> 'Shape':
        if self._ale_fn is None:
            al_rad_fn = ConstantShape(self.al_rad)
            alg_fn = self.liftingline.alg_fn
            al0_fn = self.liftingline.al0_fn
            self._ale_fn = al_rad_fn + alg_fn - al0_fn - self.ali_fn
        return self._ale_fn

    @property
    def ale(self) -> 'ndarray':
        if self._ale is None:
            self._ale = self.al_rad + self.liftingline.alg - self.liftingline.al0 - self.ali
        return self._ale

    @property
    def l_fn(self) -> 'Shape':
        if self._l_fn is None:
            self._l_fn = self.rho*self.vel*self.gamma_fn
        return self._l_fn

    @property
    def l(self) -> 'ndarray':
        if self._l is None:
            self._l = self.rho*self.vel*self.gamma
        return self._l

    @property
    def cl_fn(self) -> 'Shape':
        if self._cl_fn is None:
            self._cl_fn = self.l_fn/self.liftingline.c_fn/self.q
        return self._cl_fn

    @property
    def cl(self) -> float:
        if self._cl is None:
            self._cl = self.l/self.liftingline.c/self.q
        return self._cl

    @property
    def di_fn(self) -> 'Shape':
        if self._di_fn is None:
            self._di_fn = self.l_fn*self.ali_fn
        return self._di_fn

    @property
    def di(self) -> 'ndarray':
        if self._di is None:
            self._di = self.l*self.ali
        return self._di

    @property
    def cdi_fn(self) -> 'Shape':
        if self._cdi_fn is None:
            self._cdi_fn = self.di_fn/self.liftingline.c_fn/self.q
        return self._cdi_fn

    @property
    def cdi(self) -> float:
        if self._cdi is None:
            self._cdi = self.di/self.liftingline.c/self.q
        return self._cdi

    @property
    def CL(self) -> float:
        if self._CL is None:
            self._CL = CL_fn(self.liftingline.ar, self.An)
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
    def CDi(self) -> float:
        if self._CDi is None:
            self._CDi = (pi*self.liftingline.ar*self.An**2*self.n).sum()
        return self._CDi

    @property
    def Di(self) -> float:
        if self._Di is None:
            self._Di = self.q*self.liftingline.area*self.CDi
        return self._Di

    @property
    def delta(self) -> float:
        if self._delta is None:
            self._delta = (self.An[1:]**2*self.n[1:]).sum()/self.An[0]**2
        return self._delta

    @property
    def e(self) -> float:
        if self._e is None:
            self._e = 1.0/(1.0 + self.delta)
        return self._e

    @property
    def dcl(self) -> 'ndarray':
        if self._dcl is None:
            dclp = self.cl - self.liftingline.clmax
            dcln = self.cl - self.liftingline.clmin
            dclp[dclp < 0.0] = 0.0
            dcln[dcln > 0.0] = 0.0
            self._dcl = dclp + dcln
        return self._dcl

    @property
    def norm_dcl(self) -> float:
        if self._norm_dcl is None:
            self._norm_dcl = norm(self.dcl)
        return self._norm_dcl

    @property
    def th(self) -> 'ndarray':
        if self._th is None:
            if self.liftingline.num % 2 == 0:
                self._th = self.liftingline.th
            else:
                self._th = zeros(self.liftingline.num + 3)
                ind = self.liftingline.num//2 + 1
                self._th[0] = pi
                self._th[1:ind+1] = self.liftingline.th[:ind]
                self._th[ind+1:-1] = self.liftingline.th[ind-1:]
        return self._th

    def cos_n_th(self, n: int) -> 'ndarray':
        if self._cos_n_th is None:
            self._cos_n_th = {}
        if n not in self._cos_n_th:
            self._cos_n_th[n] = cos(n*self.th)
        return self._cos_n_th[n]

    def sin_n_th(self, n: int) -> 'ndarray':
        if self._sin_n_th is None:
            self._sin_n_th = {}
        if n not in self._sin_n_th:
            self._sin_n_th[n] = sin(n*self.th)
        return self._sin_n_th[n]

    @property
    def costh(self) -> 'ndarray':
        return self.cos_n_th(1)

    @property
    def sinth(self) -> 'ndarray':
        return self.sin_n_th(1)

    @property
    def y(self) -> 'ndarray':
        if self._y is None:
            self._y = self.liftingline.b/2*self.costh
        return self._y

    @property
    def s(self) -> 'ndarray':
        return self.costh

    @property
    def sf(self) -> 'ndarray':
        if self._sf is None:
            sin_2_th = self.sin_n_th(2)
            sfoqSar = self.An[0]*(sin_2_th/2 + pi - self.th)
            for nm1 in range(1, self.An.size):
                n = nm1 + 1
                n2m1 = n**2 - 1
                cos_n_th = self.cos_n_th(n)
                sin_n_th = self.sin_n_th(n)
                term = n*self.sinth*cos_n_th - sin_n_th*self.costh
                sfoqSar += 2*self.An[nm1]/n2m1*term
            self._sf = sfoqSar*self.liftingline.ar*self.q*self.liftingline.area
            ind = self.liftingline.num//2 + 1
            self._sf[ind+1:] -= self.L
        return self._sf

    @property
    def bm(self) -> 'ndarray':
        if self._bm is None:
            pimth = pi - self.th
            term = pimth*self.costh/2 + 3*self.sinth/8 + self.sin_n_th(3)/24
            bmoqSarb = self.An[0]*term
            term = pimth/4 - self.sin_n_th(2)/6 + self.sin_n_th(4)/48
            bmoqSarb += self.An[1]*term
            for nm1 in range(2, self.An.size):
                n = nm1 + 1
                nm2 = n - 2
                np1 = n + 1
                np2 = n + 2
                term = self.sin_n_th(nm2)/nm2/nm1/4
                term -= self.sin_n_th(n)/nm1/np1/2
                term += self.sin_n_th(np2)/np2/np1/4
                bmoqSarb += self.An[nm1]*term
            qSarb = self.q*self.liftingline.area*self.liftingline.ar*self.liftingline.b
            self._bm = bmoqSarb*qSarb
            ind = self.liftingline.num//2 + 1
            self._bm[ind+1:] -= self.L*self.y[ind+1:]
        return self._bm

    @property
    def bmr(self) -> float:
        if self._bmr is None:
            ind = self.liftingline.num//2 + 1
            self._bmr = self.bm[ind]
        return self._bmr

    def stall(self, norm_tol: float = 1e-6, ale_tol: float = 1e-6,
              num_iter: int = 100, display: bool = False) -> 'LiftingLineResult':

        cla = self.cla

        count = 0

        if display:
            print(f'Iteration {count:d}: Norm dcl = {self.norm_dcl}')

        while self.norm_dcl > norm_tol:

            al_chk = absolute(self.ale) > ale_tol

            dcla = zeros(cla.shape)

            divide(self.dcl, self.ale, where=al_chk, out=dcla)

            cla -= dcla

            self.reset()

            self.cla = cla

            count = count + 1

            if display:
                print(f'Iteration {count:d}: Norm dcl = {self.norm_dcl}')

            if count > num_iter:
                err = f'Failed to converge in {num_iter:d} iterations.'
                raise ValueError(err)

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
        c_fn = self.liftingline.c_fn
        al0_fn = self.liftingline.al0_fn
        ali_fn = GeneralShape(self.An.An*self.An.n)/GeneralShape([1.0])
        temp = self.An/c_fn
        alg_fn = temp*2*b/pi + al0_fn + ali_fn
        self.liftingline.alg_fn = alg_fn
        self.liftingline.reset(excl=['_area', '_ar', '_th', '_s', '_y', '_c', '_al0'])

    def plot_gamma(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Circulation Distribution - $\Gamma$ (m/s)')
        ax.plot(self.liftingline.y, self.gamma, label=self.name)
        return ax

    def plot_ali(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Angle Distribution - $\alpha_i$ (deg)')
        ax.plot(self.liftingline.y, degrees(self.ali), label=self.name)
        return ax

    def plot_wi(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Wash Distribution $w_i$ (m/s)')
        ax.plot(self.liftingline.y, self.wi, label=self.name)
        return ax

    def plot_ale(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Effective Angle Distribution $\alpha_{eff}$ (deg)')
        ax.plot(self.liftingline.y, degrees(self.ale), label=self.name)
        return ax

    def plot_cla(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Slope - $c_{la}$ (1/rad)')
        ax.plot(self.liftingline.y, self.cla, label=self.name)
        return ax

    def plot_l(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Distribution - l (N/m)')
        ax.plot(self.liftingline.y, self.l, label=self.name)
        return ax

    def plot_cl(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Distribution - $c_l$')
        ax.plot(self.liftingline.y, self.cl, label=self.name)
        return ax

    def plot_di(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Drag Distribution - $d_i$ (N/m)')
        ax.plot(self.liftingline.y, self.di, label=self.name)
        return ax

    def plot_cdi(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Drag Coefficient - $c_{di}$')
        ax.plot(self.liftingline.y, self.cdi, label=self.name)
        return ax

    def plot_sf(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Shear Force - $SF$ (N)')
        ax.plot(self.y, self.sf, label=self.name)
        return ax

    def plot_bm(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Bending Moment - $BM$ (N.m)')
        ax.plot(self.y, self.bm, label=self.name)
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
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('CL', '.6f', data=[self.CL])
        table.add_column('CDi', '.6f', data=[self.CDi])
        table.add_column('&delta;', '.6f', data=[self.delta])
        table.add_column('e', '.6f', data=[self.e])
        report.add_object(table)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('L (N)', '.6f', data=[self.L])
        table.add_column('Di (N)', '.6f', data=[self.Di])
        table.add_column('BM Root (N.m)', '.6f', data=[self.bmr])
        report.add_object(table)
        return report

    def _repr_markdown_(self) -> MDReport:
        return self.to_mdobj()._repr_markdown_()


class LiftingLinePolar():
    name: str = None
    al_degs: 'ndarray' = None
    liftingline: 'LiftingLine' = None
    _results: Dict[float, LiftingLineResult] = None
    _CL: 'ndarray' = None
    _CDi: 'ndarray' = None

    def __init__(self, liftingline: 'LiftingLine', al_degs: 'ndarray',
                 **kwargs: Dict[str, Any]) -> None:

        self.liftingline = liftingline
        self.al_degs = al_degs

        self.name = kwargs.get('name', self.liftingline.name)

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def results(self) -> Dict[float, LiftingLineResult]:
        if self._results is None:
            self._results = {}
            for al_deg in self.al_degs:
                name = f'{self.name:s} {al_deg:.3f}'
                result = self.liftingline.return_result_alpha(al_deg, name=name)
                self._results[al_deg] = result
        return self._results

    @property
    def CL(self) -> 'ndarray':
        if self._CL is None:
            self._CL = asarray([result.CL for result in self.results.values()])
        return self._CL

    @property
    def CDi(self) -> 'ndarray':
        if self._CDi is None:
            self._CDi = asarray([result.CDi for result in self.results.values()])
        return self._CDi

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
            ax.set_xlabel(r'Induced Drag Coefficient - $C_{Di}$')
            ax.set_ylabel(r'Lift Coefficient - $C_L$')
            ax.set_title(r'Drag Polar')
        ax.plot(self.CDi, self.CL, label=self.name)
        return ax


class LiftingLine():

    name: str = None
    c_fn: 'Shape' = None
    alg_fn: 'Shape' = None
    al0_fn: 'Shape' = None
    cla_fn: 'Shape' = None
    b: float = None
    num: int = None
    _area: float = None
    _ar: float = None
    _th: 'ndarray' = None
    _s: 'ndarray' = None
    _y: 'ndarray' = None
    _c: 'ndarray' = None
    _alg: 'ndarray' = None
    _al0: 'ndarray' = None
    _cla: 'ndarray' = None
    _Ana: GeneralShape = None
    _An0: GeneralShape = None
    _CLa: float = None
    _CL0: float = None
    _aL0: float = None
    clmax_fn: 'Shape' = None
    _clmax: 'ndarray' = None
    clmin_fn: 'Shape' = None
    _clmin: 'ndarray' = None

    def __init__(self, name: str, b: float, num: int, c_fn: 'Shape',
                 alg_fn: 'Shape' = ConstantShape(0.0),
                 al0_fn: 'Shape' = ConstantShape(0.0),
                 cla_fn: 'Shape' = ConstantShape(2*pi)) -> None:
        self.name = name
        self.c_fn = c_fn
        self.alg_fn = alg_fn
        self.al0_fn = al0_fn
        self.cla_fn = cla_fn
        self.b = b
        self.num = num

    def reset(self, excl: Iterable[str] = []) -> None:
        for attr in self.__dict__:
            if attr.startswith('_') and attr not in excl:
                self.__dict__[attr] = None

    @property
    def area(self) -> float:
        if self._area is None:
            self._area = self.c_fn.area*self.b/2
        return self._area

    @area.setter
    def area(self, area: float) -> None:
        self._area = area

    @property
    def ar(self) -> float:
        if self._ar is None:
            self._ar = self.b**2/self.area
        return self._ar

    @ar.setter
    def ar(self, ar: float) -> None:
        self._ar = ar

    @property
    def th(self) -> 'ndarray':
        if self._th is None:
            self._th = linspace(pi, 0.0, self.num+2)[1:-1]
        return self._th

    @th.setter
    def th(self, th: 'ndarray') -> None:
        if th.size != self.num:
            raise ValueError('th is not correct size.')
        self._th = th

    @property
    def s(self) -> 'ndarray':
        if self._s is None:
            self._s = cos(self.th)
        return self._s

    @s.setter
    def s(self, s: 'ndarray') -> None:
        if s.size != self.num:
            raise ValueError('s is not correct size.')
        self._s = s

    @property
    def y(self) -> 'ndarray':
        if self._y is None:
            self._y = self.b/2*self.s
        return self._y

    @y.setter
    def y(self, y: 'ndarray') -> None:
        if y.size != self.num:
            raise ValueError('y is not correct size.')
        self._y = y

    @property
    def c(self) -> 'ndarray':
        if self._c is None:
            self._c = self.c_fn(self.s)
        return self._c

    @c.setter
    def c(self, c: 'ndarray') -> None:
        if c.size != self.num:
            raise ValueError('c is not correct size.')
        self._c = c

    @property
    def alg(self) -> 'ndarray':
        if self._alg is None:
            self._alg = self.alg_fn(self.s)
        return self._alg

    @alg.setter
    def alg(self, alg: 'ndarray') -> None:
        if alg.size != self.num:
            raise ValueError('alg is not correct size.')
        self._alg = alg

    @property
    def al0(self) -> 'ndarray':
        if self._al0 is None:
            self._al0 = self.al0_fn(self.s)
        return self._al0

    @al0.setter
    def al0(self, al0: 'ndarray') -> None:
        if al0.size != self.num:
            raise ValueError('al0 is not correct size.')
        self._al0 = al0

    @property
    def cla(self) -> 'ndarray':
        if self._cla is None:
            self._cla = self.cla_fn(self.s)
        return self._cla

    @cla.setter
    def cla(self, cla: 'ndarray') -> None:
        if cla.size != self.num:
            raise ValueError('cla is not correct size.')
        self._cla = cla

    @property
    def clmax(self) -> 'ndarray':
        if self.clmax_fn is None:
            self._clmax = full(self.num, float('inf'))
        else:
            self._clmax = self.clmax_fn(self.s)
        return self._clmax

    @clmax.setter
    def clmax(self, clmax: 'ndarray') -> None:
        if clmax.size != self.num:
            raise ValueError('clmax is not correct size.')
        self._clmax = clmax

    @property
    def clmin(self) -> 'ndarray':
        if self.clmin_fn is None:
            self._clmin = full(self.num, -float('inf'))
        else:
            self._clmin = self.clmin_fn(self.s)
        return self._clmin

    @clmin.setter
    def clmin(self, clmin: 'ndarray') -> None:
        if clmin.size != self.num:
            raise ValueError('clmin is not correct size.')
        self._clmin = clmin

    def solve(self, cla: Optional['ndarray'] = None) -> Optional[Tuple[GeneralShape, GeneralShape]]:

        if cla is None:
            claco4b = self.cla*self.c/4/self.b
        else:
            claco4b = cla*self.c/4/self.b

        lhs = zeros((self.num, self.num))
        rhs = zeros((self.num, 2))

        sinth = sin(self.th)

        for n in range(1, self.num+1):
            sin_nth = sin(n*self.th)
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
    def Ana(self) -> 'ndarray':
        if self._Ana is None:
            self.solve()
        return self._Ana

    @Ana.setter
    def Ana(self, Ana: 'ndarray') -> None:
        if Ana.size != self.num:
            raise ValueError('Ana is not correct size.')
        self._Ana = Ana

    @property
    def An0(self) -> 'ndarray':
        if self._An0 is None:
            self.solve()
        return self._An0

    @An0.setter
    def An0(self, An0: 'ndarray') -> None:
        if An0.size != self.num:
            raise ValueError('An0 is not correct size.')
        self._An0 = An0

    @property
    def CLa(self) -> float:
        if self._CLa is None:
            self._CLa = CL_fn(self.ar, self.Ana)
        return self._CLa

    @property
    def CL0(self) -> float:
        if self._CL0 is None:
            self._CL0 = CL_fn(self.ar, self.An0)
        return self._CL0

    @property
    def aL0(self) -> float:
        if self._aL0 is None:
            self._aL0 = -self.CL0/self.CLa
        return self._aL0

    def copy(self, incl_attr: bool = False) -> 'LiftingLine':
        ll = LiftingLine(self.name, self.c_fn, self.alg_fn, self.al0_fn,
                         self.cla_fn, self.b, self.num)
        ll.clmax_fn = self.clmax_fn
        ll.clmin_fn = self.clmin_fn
        if incl_attr:
            for attr in self.__dict__:
                ll.__dict__[attr] = self.__dict__[attr]
        return ll

    def plot_c(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Chord - c (m)')
        ax.plot(self.y, self.c, label=self.name)
        return ax

    def plot_al0(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Zero Lift Angle of Attack -w $\alpha_{l0}$ (deg)')
        ax.plot(self.y, degrees(self.al0), label=self.name)
        return ax

    def plot_alg(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Spanwise Twist - $\alpha_g$ (deg)')
        ax.plot(self.y, degrees(self.alg), label=self.name)
        return ax

    def plot_cla(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Slope - $c_{la}$ (1/rad)')
        ax.plot(self.y, self.cla, label=self.name)
        return ax

    def return_result_alpha(self, al_deg: float, **kwargs: Dict[str, Any]) -> 'LiftingLineResult':
        res = LiftingLineResult(self, al_deg, **kwargs)
        return res

    def return_result_CL(self, CL: float, **kwargs: Dict[str, Any]) -> 'LiftingLineResult':
        al_rad = (CL - self.CL0)/self.CLa
        al_deg = degrees(al_rad)
        res = self.return_result_alpha(al_deg, **kwargs)
        res.al_rad = al_rad
        res.CL = CL
        return res

    def return_result_L(self, L: float, **kwargs: Dict[str, Any]) -> 'LiftingLineResult':
        rho = kwargs.get('rho', 1.0)
        vel = kwargs.get('vel', 1.0)
        q = rho*vel**2/2
        res = self.return_result_CL(L/self.area/q, **kwargs)
        res.q = q
        res.L = L
        return res

    def return_polar(self, al_degs: Iterable[float],
                     **kwargs: Dict[str, Any]) -> LiftingLinePolar:
        al_degs = asarray(al_degs)
        ll_pol = LiftingLinePolar(self, al_degs, **kwargs)
        return ll_pol

    def to_mdobj(self) -> MDReport:
        report = MDReport()
        heading = MDHeading(f'Lifting Line for {self.name:s}', 2)
        report.add_object(heading)
        table = MDTable()
        table.add_column('Name', 's', data=[self.name])
        table.add_column('Span (m)', '.1f', data=[self.b])
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
