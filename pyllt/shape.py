from collections.abc import Iterable
from typing import TYPE_CHECKING, Any

from matplotlib.pyplot import figure
from numpy import (absolute, asarray, divide, full, logical_and, logical_not,
                   ndarray, pi, radians, zeros)

from .spacing import CosineSpacing

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from matplotlib.axes import Axes

    from .spacing import Spacing


class Shape():

    _area: float = None

    def __init__(self) -> None:
        pass

    @property
    def area(self) -> None:
        return self._area

    def __call__(self) -> None:
        return None

    def __add__(self, other: 'Shape') -> 'AddShape | Shape':
        if isinstance(other, Shape):
            return AddShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(other)}')

    def __radd__(self, other: 'Shape') -> 'AddShape | Shape':
        if isinstance(other, Shape):
            return AddShape(other, self)
        else:
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(other)}')

    def __sub__(self, other: 'Shape') -> 'SubShape | Shape':
        if isinstance(other, Shape):
            return SubShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(other)}')

    def __rsub__(self, other: 'Shape') -> 'SubShape | Shape':
        if isinstance(other, Shape):
            return SubShape(other, self)
        else:
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(other)}')

    def __mul__(self, other: 'Shape | float | int') -> 'MulShape | Shape':
        if isinstance(other, Shape):
            return MulShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(other)}')

    def __rmul__(self, other: 'Shape | float | int') -> 'MulShape | Shape':
        if isinstance(other, Shape):
            return MulShape(other, self)
        else:
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(other)}')

    def __div__(self, other: 'Shape | float | int') -> 'DivShape | Shape':
        if isinstance(other, Shape):
            return DivShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(other)}')

    def __truediv__(self, other: 'Shape | float | int') -> 'DivShape | Shape':
        if isinstance(other, Shape):
            return DivShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(other)}')

    def plot(self, ax: 'Axes | None' = None, num: int = 100,
             **kwargs: dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel('s')
            ax.set_ylabel('d')
        spc = CosineSpacing(num)
        d = self(spc)
        ax.plot(spc.s, d, **kwargs)
        return ax

    def __repr__(self) -> str:
        return f'<{self.__str__():s}>'

    def __str__(self) -> str:
        return 'Shape()'


class AddShape(Shape):
    arg1: Shape = None
    arg2: Shape = None

    def __init__(self, arg1: Shape, arg2: Shape) -> None:
        self.arg1 = arg1
        self.arg2 = arg2

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        arg1 = self.arg1(spc)
        arg2 = self.arg2(spc)
        return arg1 + arg2

    def __mul__(self, other: float | int) -> 'MulShape | AddShape':
        if isinstance(other, (float, int)):
            return AddShape(self.arg1*other, self.arg2*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: float | int) -> 'MulShape | AddShape':
        if isinstance(other, (float, int)):
            return AddShape(other*self.arg1, other*self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: float | int) -> 'DivShape | AddShape':
        if isinstance(other, (float, int)):
            return AddShape(self.arg1/other, self.arg2/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: float | int) -> 'DivShape | AddShape':
        if isinstance(other, (float, int)):
            return AddShape(self.arg1/other, self.arg2/other)
        else:
            return super().__truediv__(other)

    def __str__(self) -> str:
        return f'AddShape({self.arg1}, {self.arg2})'


class SubShape(Shape):
    arg1: Shape = None
    arg2: Shape = None

    def __init__(self, arg1: Shape, arg2: Shape) -> None:
        self.arg1 = arg1
        self.arg2 = arg2

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        arg1 = self.arg1(spc)
        arg2 = self.arg2(spc)
        return arg1 - arg2

    def __mul__(self, other: float | int) -> 'MulShape | SubShape':
        if isinstance(other, (float, int)):
            return SubShape(self.arg1*other, self.arg2*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: float | int) -> 'MulShape | SubShape':
        if isinstance(other, (float, int)):
            return SubShape(other*self.arg1, other*self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: float | int) -> 'DivShape | SubShape':
        if isinstance(other, (float, int)):
            return SubShape(self.arg1/other, self.arg2/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: float | int) -> 'DivShape | SubShape':
        if isinstance(other, (float, int)):
            return SubShape(self.arg1/other, self.arg2/other)
        else:
            return super().__truediv__(other)

    def __str__(self) -> str:
        return f'SubShape({self.arg1}, {self.arg2})'


class MulShape(Shape):
    arg1: Shape = None
    arg2: Shape = None

    def __init__(self, arg1: Shape, arg2: Shape) -> None:
        self.arg1 = arg1
        self.arg2 = arg2

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        arg1 = self.arg1(spc)
        arg2 = self.arg2(spc)
        return arg1*arg2

    def __mul__(self, other: float | int) -> 'MulShape':
        if isinstance(other, (float, int)):
            return MulShape(self.arg1, self.arg2*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Shape | float | int) -> 'MulShape':
        if isinstance(other, (float, int)):
            return MulShape(other*self.arg1, self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Shape | float | int) -> 'DivShape | MulShape':
        if isinstance(other, (float, int)):
            return MulShape(self.arg1, self.arg2/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Shape | float | int) -> 'DivShape | MulShape':
        if isinstance(other, (float, int)):
            return MulShape(self.arg1, self.arg2/other)
        else:
            return super().__truediv__(other)

    def __str__(self) -> str:
        return f'MulShape({self.arg1}, {self.arg2})'


class DivShape(Shape):
    arg1: Shape = None
    arg2: Shape = None

    def __init__(self, arg1: Shape, arg2: Shape) -> None:
        self.arg1 = arg1
        self.arg2 = arg2

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        result = zeros(spc.shape)
        res1 = self.arg1(spc)
        res2 = self.arg2(spc)
        res1_eql0 = absolute(res1) < 1e-12
        res2_eql0 = absolute(res2) < 1e-12
        both_eql0 = logical_and(res1_eql0, res2_eql0)
        both_not0 = logical_not(both_eql0)
        divide(res1, res2, out=result, where=both_not0)
        if isinstance(res2_eql0, bool):
            anycheck = res2_eql0
        else:
            res2_eql0 = asarray(res2_eql0)
            anycheck = res2_eql0.any()
        if anycheck:
            if isinstance(self.arg1, GeneralShape) and isinstance(self.arg2, GeneralShape):
                lhopital = zeros(spc.shape)
                der1 = self.arg1.derivative(spc)
                der2 = self.arg2.derivative(spc)
                check_der = logical_and(der1 != 0.0, der2 != 0.0)
                divide(der1, der2, out=lhopital, where=check_der)
                result[both_eql0] = lhopital[both_eql0]
        return result

    def __mul__(self, other: float | int) -> 'MulShape | DivShape':
        if isinstance(other, (float, int)):
            return DivShape(self.arg1*other, self.arg2)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: float | int) -> 'MulShape | DivShape':
        if isinstance(other, (float, int)):
            return DivShape(other*self.arg1, self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: float | int) -> 'DivShape':
        if isinstance(other, (float, int)):
            return DivShape(self.arg1/other, self.arg2)
        else:
            return super().__div__(other)

    def __truediv__(self, other: float | int) -> 'DivShape':
        if isinstance(other, (float, int)):
            return DivShape(self.arg1/other, self.arg2)
        else:
            return super().__truediv__(other)

    def __str__(self) -> str:
        return f'DivShape({self.arg1}, {self.arg2})'


class ConstantShape(Shape):

    value: float = None

    def __init__(self, value: float) -> None:
        self.value = value

    @property
    def area(self) -> None:
        if self._area is None:
            self._area = 2*self.value
        return self._area

    def as_radians(self) -> 'TaperedShape':
        self.value = radians(self.value)
        return self

    def __add__(self, other: 'ConstantShape | Shape') -> 'ConstantShape | AddShape':
        if isinstance(other, ConstantShape):
            return ConstantShape(self.value + other.value)
        else:
            return super().__add__(other)

    def __radd__(self, other: 'ConstantShape | Shape') -> 'ConstantShape | AddShape':
        if isinstance(other, ConstantShape):
            return ConstantShape(self.value + other.value)
        else:
            return super().__radd__(other)

    def __sub__(self, other: 'ConstantShape | Shape') -> 'ConstantShape | SubShape':
        if isinstance(other, ConstantShape):
            return ConstantShape(self.value - other.value)
        else:
            return super().__sub__(other)

    def __rsub__(self, other: 'ConstantShape | Shape') -> 'ConstantShape | SubShape':
        if isinstance(other, ConstantShape):
            return ConstantShape(other.value - self.value)
        else:
            return super().__rsub__(other)

    def __mul__(self, other: Shape | float | int) -> 'ConstantShape | MulShape':
        if isinstance(other, (float, int)):
            return ConstantShape(self.value*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Shape | float | int) -> 'ConstantShape | MulShape':
        if isinstance(other, (float, int)):
            return ConstantShape(other*self.value)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Shape | float | int) -> 'ConstantShape | DivShape':
        if isinstance(other, (float, int)):
            return ConstantShape(self.value/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Shape | float | int) -> 'ConstantShape | DivShape':
        if isinstance(other, (float, int)):
            return ConstantShape(self.value/other)
        else:
            return super().__truediv__(other)

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        return full(spc.shape, self.value)

    def __str__(self) -> str:
        return f'ConstantShape({self.value})'


class TaperedShape(Shape):

    value_root: float = None
    value_tip: float = None

    def __init__(self, value_root: float, value_tip: float) -> None:
        self.value_root = value_root
        self.value_tip = value_tip

    @property
    def area(self) -> None:
        if self._area is None:
            self._area = self.value_root + self.value_tip
        return self._area

    def as_radians(self) -> 'TaperedShape':
        self.value_root = radians(self.value_root)
        self.value_tip = radians(self.value_tip)
        return self

    def __add__(self, other: 'TaperedShape | Shape') -> 'TaperedShape | AddShape':
        if isinstance(other, TaperedShape):
            return TaperedShape(self.value_root + other.value_root, self.value_tip + other.value_tip)
        else:
            return super().__add__(other)

    def __radd__(self, other: 'TaperedShape | Shape') -> 'TaperedShape | AddShape':
        if isinstance(other, TaperedShape):
            return TaperedShape(self.value_root + other.value_root, self.value_tip + other.value_tip)
        else:
            return super().__radd__(other)

    def __sub__(self, other: 'TaperedShape | Shape') -> 'TaperedShape | SubShape':
        if isinstance(other, TaperedShape):
            return TaperedShape(self.value_root - other.value_root, self.value_tip - other.value_tip)
        else:
            return super().__sub__(other)

    def __rsub__(self, other: 'TaperedShape | Shape') -> 'TaperedShape | SubShape':
        if isinstance(other, TaperedShape):
            return TaperedShape(other.value_root - self.value_root, other.value_tip - self.value_tip)
        else:
            return super().__rsub__(other)

    def __mul__(self, other: 'Shape | NDArray | float') -> 'TaperedShape | MulShape':
        if isinstance(other, (float, int)):
            return TaperedShape(self.value_root*other, self.value_tip*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: 'Shape | NDArray | float') -> 'TaperedShape | MulShape':
        if isinstance(other, (float, int)):
            return TaperedShape(other*self.value_root, other*self.value_tip)
        else:
            return super().__rmul__(other)

    def __div__(self, other: 'Shape | NDArray | float') -> 'TaperedShape | DivShape':
        if isinstance(other, (float, int)):
            return TaperedShape(self.value_root/other, self.value_tip/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: 'Shape | NDArray | float') -> 'TaperedShape | DivShape':
        if isinstance(other, (float, int)):
            return TaperedShape(self.value_root/other, self.value_tip/other)
        else:
            return super().__truediv__(other)

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        return self.value_root - absolute(spc.s)*(self.value_root - self.value_tip)

    def __str__(self) -> str:
        return f'TaperedShape({self.value_root}, {self.value_tip})'


class GeneralShape(Shape):

    An: 'NDArray' = None
    _n: Iterable[int] = None

    def __init__(self, An: 'NDArray') -> None:
        self.An = asarray(An)

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr.startswith('_'):
                setattr(self, attr, None)

    def resize(self, size: int) -> None:
        An = zeros(size)
        An[:self.An.size] = self.An
        self.An = An
        self.reset()

    @property
    def n(self) -> Iterable[int]:
        if self._n is None:
            self._n = range(1, self.An.size + 1)
        return self._n

    @property
    def area(self) -> float:
        if self._area is None:
            self._area = self.An[0]*pi/2
        return self._area

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        result = zeros(spc.shape)
        for n, An in self.items():
            if An != 0.0:
                result += An*spc.sin_n_th(n)
        return result

    def derivative(self, spc: 'Spacing') -> 'NDArray':
        result = zeros(spc.shape)
        for n, An in self.items():
            if An != 0.0:
                result += An*spc.n_cos_n_th(n)
        return result

    @property
    def size(self) -> int:
        return self.An.size

    def __getitem__(self, n: int | slice) -> 'float | NDArray':
        return self.An[n]

    def __setitem__(self, n: int | slice, An: 'float | NDArray') -> None:
        self.An[n] = An

    def __mul__(self, other: 'Shape | NDArray | float') -> 'GeneralShape | MulShape':
        if isinstance(other, (int, float, ndarray)):
            return GeneralShape(self.An*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: 'Shape | NDArray | float') -> 'GeneralShape | MulShape':
        if isinstance(other, (int, float, ndarray)):
            return GeneralShape(other*self.An)
        else:
            return super().__rmul__(other)

    def __div__(self, other: 'Shape | NDArray | float') -> 'GeneralShape | DivShape':
        if isinstance(other, (float, ndarray)):
            return GeneralShape(self.An/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: 'Shape | NDArray | float') -> 'GeneralShape | DivShape':
        if isinstance(other, (float, ndarray)):
            return GeneralShape(self.An/other)
        else:
            return super().__truediv__(other)

    def __add__(self, other: 'GeneralShape | Shape') -> 'GeneralShape | AddShape':
        if isinstance(other, GeneralShape):
            if self.size == other.size:
                result = GeneralShape(self.An + other.An)
            elif self.size > other.size:
                result = GeneralShape(self.An[:other.size] + other.An)
            elif self.size < other.size:
                result = GeneralShape(self.An + other.An[:self.size])
            if result.size == 1:
                result.__class__ = EllipticalShape
            elif result.size == 3:
                if result.An[1] == 0.0:
                    result.__class__ = BellShape
            return result
        else:
            return super().__add__(other)

    def __radd__(self, other: 'GeneralShape | Shape') -> 'GeneralShape | AddShape':
        if isinstance(other, GeneralShape):
            if self.size == other.size:
                result = GeneralShape(other.An + self.An)
            elif self.size > other.size:
                result = GeneralShape(other.An[:self.size] + self.An)
            elif self.size < other.size:
                result = GeneralShape(other.An + self.An[:other.size])
            if result.size == 1:
                result.__class__ = EllipticalShape
            elif result.size == 3:
                if result.An[1] == 0.0:
                    result.__class__ = BellShape
            return result
        else:
            return super().__radd__(other)

    def __sub__(self, other: 'GeneralShape | Shape') -> 'GeneralShape | SubShape':
        if isinstance(other, GeneralShape):
            if self.size == other.size:
                result = GeneralShape(self.An - other.An)
            elif self.size > other.size:
                result = GeneralShape(self.An[:other.size] - other.An)
            elif self.size < other.size:
                result = GeneralShape(self.An - other.An[:self.size])
            if result.size == 1:
                result.__class__ = EllipticalShape
            elif result.size == 3:
                if result.An[1] == 0.0:
                    result.__class__ = BellShape
            return result
        else:
            return super().__sub__(other)

    def __rsub__(self, other: 'GeneralShape | Shape') -> 'GeneralShape | SubShape':
        if isinstance(other, GeneralShape):
            if self.size == other.size:
                result = GeneralShape(other.An - self.An)
            elif self.size > other.size:
                result = GeneralShape(other.An[:self.size] - self.An)
            elif self.size < other.size:
                result = GeneralShape(other.An - self.An[:other.size])
            if result.size == 1:
                result.__class__ = EllipticalShape
            elif result.size == 3:
                if result.An[1] == 0.0:
                    result.__class__ = BellShape
            return result
        else:
            return super().__rsub__(other)

    def __pow__(self, other: float) -> 'GeneralShape':
        if isinstance(other, (float, int)):
            return GeneralShape(self.An**other)
        else:
            raise TypeError(f'unsupported operand type(s) for ** or pow(): {type(self)} and {type(other)}')

    def sum(self) -> float:
        return self.An.sum()

    def normalise_area(self) -> 'GeneralShape':
        result = GeneralShape(self.An/self.area)
        result.__class__ = self.__class__
        return result

    def to_constant(self) -> ConstantShape:
        return ConstantShape(self.area/2)

    def keys(self) -> Iterable[int]:
        return self.n

    def values(self) -> 'NDArray':
        return self.An

    def items(self) -> Iterable[tuple[int, float]]:
        return zip(self.n, self.An)

    def __str__(self) -> str:
        return f'GeneralShape({self.An})'


class EllipticalShape(GeneralShape):

    def __init__(self, scale: float = 1.0) -> None:
        An = zeros(1)
        An[0] = scale
        super().__init__(An)

    def __str__(self) -> str:
        return f'EllipticalShape({self.An})'


class BellShape(GeneralShape):

    def __init__(self, scale: float = 1.0) -> None:
        An = zeros(3)
        An[0] = scale
        An[2] = -scale/3
        super().__init__(An)

    def __str__(self) -> str:
        return f'BellShape({self.An})'


class InducedAngleShape(Shape):

    An: 'NDArray' = None
    _n: Iterable[int] = None

    def __init__(self, An: 'NDArray') -> None:
        self.An = asarray(An)

    def reset(self) -> None:
        for attr in self.__dict__:
            if attr.startswith('_'):
                setattr(self, attr, None)

    def resize(self, size: int) -> None:
        An = zeros(size)
        An[:self.An.size] = self.An
        self.An = An
        self.reset()

    @property
    def n(self) -> Iterable[int]:
        if self._n is None:
            self._n = range(1, self.An.size + 1)
        return self._n

    def __call__(self, spc: 'Spacing') -> 'NDArray':
        result = zeros(spc.shape)
        sinth = spc.sinth
        sinth_not0 = absolute(sinth) > 1e-12
        sinth_eql0 = logical_not(sinth_not0)
        if isinstance(sinth_eql0, bool):
            anycheck = sinth_eql0
        else:
            sinth_eql0 = asarray(sinth_eql0)
            anycheck = sinth_eql0.any()
        if anycheck:
            costh = spc.costh
            costh_not0 = absolute(costh) > 1e-12
        for n, An in self.items():
            if n == 1:
                result += An
            else:
                if An != 0.0:
                    temp1 = zeros(spc.shape)
                    divide(spc.n_sin_n_th(n), sinth, out=temp1, where=sinth_not0)
                    if anycheck:
                        temp2 = zeros(spc.shape)
                        divide(spc.n2_cos_n_th(n), costh, out=temp2, where=costh_not0)
                        temp1[sinth_eql0] = temp2[sinth_eql0]
                    result += An*temp1
        return result

    def derivative(self, spc: 'Spacing') -> 'NDArray':
        result = zeros(spc.shape)
        sinth = spc.sinth
        sinth2 = sinth**2
        costh = spc.costh
        check = absolute(sinth) > 1e-12
        for n, An in self.items():
            if n > 1:
                if An != 0.0:
                    temp1 = zeros(spc.shape)
                    divide(spc.n2_cos_n_th(n), sinth, out=temp1, where=check)
                    temp2 = zeros(spc.shape)
                    divide(spc.n_sin_n_th(n)*costh, sinth2, out=temp2, where=check)
                    result += An*(temp1 - temp2)
        return result

    @property
    def size(self) -> int:
        return self.An.size

    def __getitem__(self, n: int | slice) -> 'float | NDArray':
        return self.An[n]

    def __setitem__(self, n: int | slice, An: 'float | NDArray') -> None:
        self.An[n] = An

    def __mul__(self, other: 'Shape | NDArray | float') -> 'InducedAngleShape | MulShape':
        if isinstance(other, (int, float, ndarray)):
            return InducedAngleShape(self.An*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: 'Shape | NDArray | float') -> 'InducedAngleShape | MulShape':
        if isinstance(other, (int, float, ndarray)):
            return InducedAngleShape(other*self.An)
        else:
            return super().__rmul__(other)

    def __div__(self, other: 'Shape | NDArray | float') -> 'InducedAngleShape | DivShape':
        if isinstance(other, (float, ndarray)):
            return InducedAngleShape(self.An/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: 'Shape | NDArray | float') -> 'InducedAngleShape | DivShape':
        if isinstance(other, (float, ndarray)):
            return InducedAngleShape(self.An/other)
        else:
            return super().__truediv__(other)

    def __add__(self, other: 'InducedAngleShape | Shape') -> 'InducedAngleShape | AddShape':
        if isinstance(other, InducedAngleShape):
            if self.size == other.size:
                result = InducedAngleShape(self.An + other.An)
            elif self.size > other.size:
                result = InducedAngleShape(self.An[:other.size] + other.An)
            elif self.size < other.size:
                result = InducedAngleShape(self.An + other.An[:self.size])
            return result
        else:
            return super().__add__(other)

    def __radd__(self, other: 'InducedAngleShape | Shape') -> 'InducedAngleShape | AddShape':
        if isinstance(other, InducedAngleShape):
            if self.size == other.size:
                result = InducedAngleShape(other.An + self.An)
            elif self.size > other.size:
                result = InducedAngleShape(other.An[:self.size] + self.An)
            elif self.size < other.size:
                result = InducedAngleShape(other.An + self.An[:other.size])
            return result
        else:
            return super().__radd__(other)

    def __sub__(self, other: 'InducedAngleShape | Shape') -> 'InducedAngleShape | SubShape':
        if isinstance(other, InducedAngleShape):
            if self.size == other.size:
                result = InducedAngleShape(self.An - other.An)
            elif self.size > other.size:
                result = InducedAngleShape(self.An[:other.size] - other.An)
            elif self.size < other.size:
                result = InducedAngleShape(self.An - other.An[:self.size])
            return result
        else:
            return super().__sub__(other)

    def __rsub__(self, other: 'InducedAngleShape | Shape') -> 'InducedAngleShape | SubShape':
        if isinstance(other, InducedAngleShape):
            if self.size == other.size:
                result = InducedAngleShape(other.An - self.An)
            elif self.size > other.size:
                result = InducedAngleShape(other.An[:self.size] - self.An)
            elif self.size < other.size:
                result = InducedAngleShape(other.An - self.An[:other.size])
            return result
        else:
            return super().__rsub__(other)

    def __pow__(self, other: float) -> 'InducedAngleShape':
        if isinstance(other, (float, int)):
            return InducedAngleShape(self.An**other)
        else:
            raise TypeError(f'unsupported operand type(s) for ** or pow(): {type(self)} and {type(other)}')

    def sum(self) -> float:
        return self.An.sum()

    def keys(self) -> Iterable[int]:
        return self.n

    def values(self) -> 'NDArray':
        return self.An

    def items(self) -> Iterable[tuple[int, float]]:
        return zip(self.n, self.An)

    def __str__(self) -> str:
        return f'InducedAngleShape({self.An})'
