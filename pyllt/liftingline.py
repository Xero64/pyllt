from typing import TYPE_CHECKING, Any, Dict, Iterable, Optional, Tuple

from matplotlib.pyplot import figure
from numpy import (absolute, arccos, asarray, cos, degrees, divide, flip, full,
                   hstack, linspace, pi, radians, sin, zeros)
from numpy.linalg import norm, solve
from py2md.classes import MDHeading, MDReport, MDTable

from .shape import ConstantShape, EllipticalShape, GeneralShape, InducedAngleShape

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy import ndarray

    from .shape import Shape


class LiftingLineResult():

    liftingline: 'LiftingLine' = None
    al_deg: float = None
    vel: float = None
    rho: float = None
    name: str = None
    _num: int = None
    _al_rad: float = None
    _al_shp: 'Shape' = None
    _s: 'ndarray' = None
    _th: 'ndarray' = None
    _cla: 'ndarray' = None
    _cos_n_th: Dict[int, 'ndarray'] = None
    _sin_n_th: Dict[int, 'ndarray'] = None
    _y: 'ndarray' = None
    _An: GeneralShape = None
    _n: 'ndarray' = None
    _q: float = None
    _gamma_shp: 'Shape' = None
    _gamma: 'ndarray' = None
    _ali_shp: 'Shape' = None
    _ali: 'ndarray' = None
    _wi_shp: 'Shape' = None
    _wi: 'ndarray' = None
    _ale_shp: 'Shape' = None
    _ale: 'ndarray' = None
    _l_shp: 'Shape' = None
    _l: 'ndarray' = None
    _cl_shp: 'Shape' = None
    _cl: 'ndarray' = None
    _di_shp: 'Shape' = None
    _di: 'ndarray' = None
    _cdi_shp: 'Shape' = None
    _cdi: 'ndarray' = None
    _d: 'ndarray' = None
    _CL: float = None
    _L: float = None
    _CDi: float = None
    _Di: float = None
    _D0: float = None
    _D: float = None
    _delta: float = None
    _e: float = None
    _dcl: 'ndarray' = None
    _norm_dcl: float = None
    _sf: 'ndarray' = None
    _bm: 'ndarray' = None
    _BM_root: float = None

    def __init__(self, liftingline: 'LiftingLine',
                 al_deg: float, **kwargs: Dict[str, Any]) -> None:

        self.liftingline = liftingline
        self.al_deg = al_deg

        self.vel = kwargs.get('vel', 1.0)
        self.rho = kwargs.get('rho', 1.0)
        self.name = kwargs.get('name', self.liftingline.name)
        self.num = kwargs.get('num', self.liftingline.num)

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
        self.reset()

    @property
    def al_shp(self) -> 'Shape':
        if self._al_shp is None:
            self._al_shp = ConstantShape(self.al_rad)
        return self._al_shp

    @property
    def num(self) -> int:
        if self._num is None:
            self._num = self.liftingline.num // 2
        return self._num

    @num.setter
    def num(self, num: int) -> None:
        self._num = num
        self.reset()

    @property
    def s(self) -> 'ndarray':
        if self._s is None:
            num = self.liftingline.num//2 + 1
            th = linspace(pi/2, 0.0, num)
            s = cos(th)
            self._s = hstack((-flip(s), s))
        return self._s

    @property
    def th(self) -> 'ndarray':
        if self._th is None:
            self._th = arccos(self.s)
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
            self._y = self.liftingline.b/2*self.s
        return self._y

    @property # Needs rework
    def cla(self) -> 'ndarray':
        if self._cla is None:
            self._cla = self.liftingline.cla_shp(self.liftingline.s)
        return self._cla

    @cla.setter # Needs rework
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
    def gamma_shp(self) -> 'Shape':
        if self._gamma_shp is None:
            self._gamma_shp = 2*self.liftingline.b*self.vel*self.An
        return self._gamma_shp

    @property
    def gamma(self) -> 'ndarray':
        if self._gamma is None:
            self._gamma = self.gamma_shp(self.s)
        return self._gamma

    @property
    def ali_shp(self) -> 'Shape': # Needs its own shape
        if self._ali_shp is None:
            self._ali_shp = InducedAngleShape(self.An)
        return self._ali_shp

    @property
    def ali(self) -> 'ndarray':
        if self._ali is None:
            self._ali = self.ali_shp(self.s)
        return self._ali

    @property
    def wi_shp(self) -> 'Shape':
        if self._wi_shp is None:
            self._wi_shp = self.vel*self.ali_shp
        return self._wi_shp

    @property
    def wi(self) -> 'ndarray':
        if self._wi is None:
            self._wi = self.wi_shp(self.s)
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
    def ale(self) -> 'ndarray':
        if self._ale is None:
            self._ale = self.ale_shp(self.s)
        return self._ale

    @property
    def l_shp(self) -> 'Shape':
        if self._l_shp is None:
            self._l_shp = self.rho*self.vel*self.gamma_shp
        return self._l_shp

    @property
    def l(self) -> 'ndarray':
        if self._l is None:
            self._l = self.l_shp(self.s)
        return self._l

    @property
    def cl_shp(self) -> 'Shape':
        if self._cl_shp is None:
            self._cl_shp = self.l_shp/self.liftingline.c_shp/self.q
        return self._cl_shp

    @property
    def cl(self) -> float:
        if self._cl is None:
            self._cl = self.cl_shp(self.s)
        return self._cl

    @property
    def di_shp(self) -> 'Shape':
        if self._di_shp is None:
            self._di_shp = self.l_shp*self.ali_shp
        return self._di_shp

    @property
    def di(self) -> 'ndarray':
        if self._di is None:
            self._di = self.di_shp(self.s)
        return self._di

    @property
    def cdi_shp(self) -> 'Shape':
        if self._cdi_shp is None:
            self._cdi_shp = self.di_shp/self.liftingline.c_shp/self.q
        return self._cdi_shp

    @property
    def cdi(self) -> float:
        if self._cdi is None:
            self._cdi = self.cdi_shp(self.s)
        return self._cdi

    @property
    def sf(self) -> 'ndarray':
        if self._sf is None:
            sfoqSar = zeros(self.th.size)
            for n, An in self.An.items():
                if n == 1:
                    sfoqSar += An*(self.sinth*self.costh + pi - self.th)
                else:
                    n2m1 = n**2 - 1
                    cos_n_th = self.cos_n_th(n)
                    sin_n_th = self.sin_n_th(n)
                    term = n*self.sinth*cos_n_th - sin_n_th*self.costh
                    sfoqSar += 2*An/n2m1*term
            q = self.q
            Sref = self.liftingline.area
            ar = self.liftingline.ar
            self._sf = sfoqSar*q*Sref*ar
            self._sf[self.s > 0.0] -= self.L
        return self._sf

    @property
    def bm(self) -> 'ndarray':
        if self._bm is None:
            bmoqSarb = zeros(self.th.size)
            pimth = pi - self.th
            for n, An in self.An.items():
                if n == 1:
                    term = pimth*self.costh/2 + 3*self.sinth/8 + self.sin_n_th(3)/24
                    bmoqSarb += An*term
                elif n == 2:
                    term = pimth/4 - self.sin_n_th(2)/6 + self.sin_n_th(4)/48
                    bmoqSarb += An*term
                else:
                    nm1 = n - 1
                    nm2 = n - 2
                    np1 = n + 1
                    np2 = n + 2
                    term = self.sin_n_th(nm2)/nm2/nm1/4
                    term -= self.sin_n_th(n)/nm1/np1/2
                    term += self.sin_n_th(np2)/np2/np1/4
                    bmoqSarb += An*term
            q = self.q
            Sref = self.liftingline.area
            ar = self.liftingline.ar
            b = self.liftingline.b
            self._bm = bmoqSarb*q*Sref*ar*b
            check = self.s > 0.0
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
    def Di(self) -> float:
        if self._Di is None:
            self._Di = self.q*self.liftingline.area*self.CDi
        return self._Di

    @property
    def BM_root(self) -> float:
        if self._BM_root is None:
            ind = self.liftingline.num//2 + 1
            self._BM_root = self.bm[ind]
        return self._BM_root

    def stall(self, norm_tol: float = 1e-6, ale_tol: float = 1e-6,
              num_iter: int = 100, display: bool = False) -> 'LiftingLineResult':

        cla = self.cla

        clmax = self.liftingline.clmax_shp(self.liftingline.s)
        clmin = self.liftingline.clmin_shp(self.liftingline.s)

        count = 0

        norm_dcl = float('inf')

        while norm_dcl >= norm_tol:

            Ana, An0 = self.liftingline.solve(cla)
            An = Ana*self.al_rad + An0

            gamma_shp = 2*self.liftingline.b*self.vel*An
            l_shp = self.rho*self.vel*gamma_shp
            cl_shp = l_shp/self.liftingline.c_shp/self.q

            ali_shp = An*self.n/EllipticalShape()
            ale_shp = self.al_shp + self.liftingline.alg_shp - self.liftingline.al0_shp - ali_shp

            cl = cl_shp(self.liftingline.s)
            ale = ale_shp(self.liftingline.s)

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

            # cla_chk = logical_and(dclp_chk, dcln_chk)
            # cla[cla_chk] = self.liftingline.cla_shp(self.liftingline.s[cla_chk])
            # self.reset()

            # self.cla = cla

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
        ali_shp = GeneralShape(self.An.An*self.An.n)/GeneralShape([1.0])
        temp = self.An/c_shp
        alg_shp = temp*2*b/pi + al0_shp + ali_shp
        self.liftingline.alg_shp = alg_shp
        self.liftingline.reset(excl=['_area', '_ar', '_th', '_s', '_y', '_c', '_al0'])

    def plot_gamma(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Circulation Distribution - $\Gamma$ (m/s)')
        ax.plot(self.y, self.gamma, label=self.name)
        return ax

    def plot_ali(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Angle Distribution - $\alpha_i$ (deg)')
        ax.plot(self.y, degrees(self.ali), label=self.name)
        return ax

    def plot_wi(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Wash Distribution $w_i$ (m/s)')
        ax.plot(self.y, self.wi, label=self.name)
        return ax

    def plot_ale(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Effective Angle Distribution $\alpha_{eff}$ (deg)')
        ax.plot(self.y, degrees(self.ale), label=self.name)
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
        ax.plot(self.y, self.l, label=self.name)
        return ax

    def plot_cl(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Lift Coefficient Distribution - $c_l$')
        ax.plot(self.y, self.cl, label=self.name)
        return ax

    def plot_di(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Drag Distribution - $d_i$ (N/m)')
        ax.plot(self.y, self.di, label=self.name)
        return ax

    def plot_cdi(self, ax: 'Axes' = None) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel(r'Spanwise Coordinate - y (m)')
            ax.set_ylabel(r'Induced Drag Coefficient - $c_{di}$')
        ax.plot(self.y, self.cdi, label=self.name)
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
        table.add_column('q (Pa)', '.3f', data=[self.q])
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
        table.add_column('BM Root (N.m)', '.6f', data=[self.BM_root])
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
    b: float = None
    c_shp: 'Shape' = None
    alg_shp: 'Shape' = None
    al0_shp: 'Shape' = None
    cla_shp: 'Shape' = None
    num: int = None
    clmax_shp: 'Shape' = None
    _clmax: 'ndarray' = None
    clmin_shp: 'Shape' = None
    _clmin: 'ndarray' = None
    cd0_shp: 'Shape' = None
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

    def __init__(self, name: str, b: float, c_shp: 'Shape',
                 alg_shp: 'Shape' = ConstantShape(0.0),
                 al0_shp: 'Shape' = ConstantShape(0.0),
                 cla_shp: 'Shape' = ConstantShape(2*pi),
                 cd0_shp: 'Shape' = ConstantShape(0.0),
                 num: int = 101) -> None:
        self.name = name
        self.b = b
        self.c_shp = c_shp
        self.alg_shp = alg_shp
        self.al0_shp = al0_shp
        self.cla_shp = cla_shp
        self.cd0_shp = cd0_shp
        self.num = num

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
            self._c = self.c_shp(self.s)
        return self._c

    @c.setter
    def c(self, c: 'ndarray') -> None:
        if c.size != self.num:
            raise ValueError('c is not correct size.')
        self._c = c

    @property
    def alg(self) -> 'ndarray':
        if self._alg is None:
            self._alg = self.alg_shp(self.s)
        return self._alg

    @alg.setter
    def alg(self, alg: 'ndarray') -> None:
        if alg.size != self.num:
            raise ValueError('alg is not correct size.')
        self._alg = alg

    @property
    def al0(self) -> 'ndarray':
        if self._al0 is None:
            self._al0 = self.al0_shp(self.s)
        return self._al0

    @al0.setter
    def al0(self, al0: 'ndarray') -> None:
        if al0.size != self.num:
            raise ValueError('al0 is not correct size.')
        self._al0 = al0

    @property
    def cla(self) -> 'ndarray':
        if self._cla is None:
            self._cla = self.cla_shp(self.s)
        return self._cla

    @cla.setter
    def cla(self, cla: 'ndarray') -> None:
        if cla.size != self.num:
            raise ValueError('cla is not correct size.')
        self._cla = cla

    @property
    def clmax(self) -> 'ndarray':
        if self.clmax_shp is None:
            self._clmax = full(self.num, float('inf'))
        else:
            self._clmax = self.clmax_shp(self.s)
        return self._clmax

    @clmax.setter
    def clmax(self, clmax: 'ndarray') -> None:
        if clmax.size != self.num:
            raise ValueError('clmax is not correct size.')
        self._clmax = clmax

    @property
    def clmin(self) -> 'ndarray':
        if self.clmin_shp is None:
            self._clmin = full(self.num, -float('inf'))
        else:
            self._clmin = self.clmin_shp(self.s)
        return self._clmin

    @clmin.setter
    def clmin(self, clmin: 'ndarray') -> None:
        if clmin.size != self.num:
            raise ValueError('clmin is not correct size.')
        self._clmin = clmin

    @property
    def cd0(self) -> 'ndarray':
        if self._cd0 is None:
            self._cd0 = self.cd0_shp(self.s)
        return self._cd0

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
