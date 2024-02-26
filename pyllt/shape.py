#%%
# Import Dependencies
from typing import TYPE_CHECKING, Any, Dict, Iterable, Optional, Union

from matplotlib.pyplot import figure
from numpy import (absolute, arccos, asarray, cos, divide, full, linspace,
                   logical_and, ndarray, pi, sin, zeros)

if TYPE_CHECKING:
    from matplotlib.axes import Axes


class Shape():

    _area: float = None

    def __init__(self) -> None:
        pass

    @property
    def area(self) -> None:
        return self._area

    def __call__(self) -> None:
        return None

    def __add__(self, other: 'Shape') -> Union['AddShape', 'Shape']:
        if isinstance(other, Shape):
            return AddShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(other)}')

    def __radd__(self, other: 'Shape') -> Union['AddShape', 'Shape']:
        if isinstance(other, Shape):
            return AddShape(other, self)
        else:
            raise TypeError(f'unsupported operand type(s) for +: {type(self)} and {type(other)}')

    def __sub__(self, other: 'Shape') -> Union['SubShape', 'Shape']:
        if isinstance(other, Shape):
            return SubShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(other)}')

    def __rsub__(self, other: 'Shape') -> Union['SubShape', 'Shape']:
        if isinstance(other, Shape):
            return SubShape(other, self)
        else:
            raise TypeError(f'unsupported operand type(s) for -: {type(self)} and {type(other)}')

    def __mul__(self, other: Union['Shape', float, int]) -> Union['MulShape', 'Shape']:
        if isinstance(other, Shape):
            return MulShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(other)}')

    def __rmul__(self, other: Union['Shape', float, int]) -> Union['MulShape', 'Shape']:
        if isinstance(other, Shape):
            return MulShape(other, self)
        else:
            raise TypeError(f'unsupported operand type(s) for *: {type(self)} and {type(other)}')

    def __div__(self, other: Union['Shape', float, int]) -> Union['DivShape', 'Shape']:
        if isinstance(other, Shape):
            return DivShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(other)}')

    def __truediv__(self, other: Union['Shape', float, int]) -> Union['DivShape', 'Shape']:
        if isinstance(other, Shape):
            return DivShape(self, other)
        else:
            raise TypeError(f'unsupported operand type(s) for /: {type(self)} and {type(other)}')

    def plot(self, ax: Optional['Axes'] = None, num: int = 100,
             **kwargs: Dict[str, Any]) -> 'Axes':
        if ax is None:
            fig = figure()
            ax = fig.gca()
            ax.grid(True)
            ax.set_xlabel('s')
            ax.set_ylabel('d')
        th = linspace(pi, 0, num)
        s = cos(th)
        d = self(s)
        ax.plot(s, d, **kwargs)
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

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        arg1 = self.arg1(s)
        arg2 = self.arg2(s)
        return arg1 + arg2

    def __mul__(self, other: Union[float, int]) -> Union['MulShape', 'AddShape']:
        if isinstance(other, (float, int)):
            return AddShape(self.arg1*other, self.arg2*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[float, int]) -> Union['MulShape', 'AddShape']:
        if isinstance(other, (float, int)):
            return AddShape(other*self.arg1, other*self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[float, int]) -> Union['DivShape', 'AddShape']:
        if isinstance(other, (float, int)):
            return AddShape(self.arg1/other, self.arg2/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[float, int]) -> Union['DivShape', 'AddShape']:
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

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        arg1 = self.arg1(s)
        arg2 = self.arg2(s)
        return arg1 - arg2

    def __mul__(self, other: Union[float, int]) -> Union['MulShape', 'SubShape']:
        if isinstance(other, (float, int)):
            return SubShape(self.arg1*other, self.arg2*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[float, int]) -> Union['MulShape', 'SubShape']:
        if isinstance(other, (float, int)):
            return SubShape(other*self.arg1, other*self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[float, int]) -> Union['DivShape', 'SubShape']:
        if isinstance(other, (float, int)):
            return SubShape(self.arg1/other, self.arg2/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[float, int]) -> Union['DivShape', 'SubShape']:
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

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        arg1 = self.arg1(s)
        arg2 = self.arg2(s)
        return arg1*arg2

    def __mul__(self, other: Union[float, int]) -> 'MulShape':
        if isinstance(other, (float, int)):
            return MulShape(self.arg1, self.arg2*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[Shape, float, int]) -> 'MulShape':
        if isinstance(other, (float, int)):
            return MulShape(other*self.arg1, self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[Shape, float, int]) -> Union['DivShape', 'MulShape']:
        if isinstance(other, (float, int)):
            return MulShape(self.arg1, self.arg2/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[Shape, float, int]) -> Union['DivShape', 'MulShape']:
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

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        result = zeros(s.shape)
        arg1 = self.arg1(s)
        arg2 = self.arg2(s)
        arg1_not_zero = arg1 != 0.0
        arg2_not_zero = arg2 != 0.0
        check = logical_and(arg1_not_zero, arg2_not_zero)
        divide(arg1, arg2, out=result, where=check)
        return result

    def __mul__(self, other: Union[float, int]) -> Union['MulShape', 'DivShape']:
        if isinstance(other, (float, int)):
            return DivShape(self.arg1*other, self.arg2)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[float, int]) -> Union['MulShape', 'DivShape']:
        if isinstance(other, (float, int)):
            return DivShape(other*self.arg1, self.arg2)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[float, int]) -> 'DivShape':
        if isinstance(other, (float, int)):
            return DivShape(self.arg1/other, self.arg2)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[float, int]) -> 'DivShape':
        if isinstance(other, (float, int)):
            return DivShape(self.arg1/other, self.arg2)
        else:
            return super().__truediv__(other)

    def __str__(self) -> str:
        return f'DivShape({self.arg1}, {self.arg2})'


class ConstantShape(Shape):

    c: float = None

    def __init__(self, c: float) -> None:
        self.c = c

    @property
    def area(self) -> None:
        if self._area is None:
            self._area = 2*self.c
        return self._area

    def __add__(self, other: Union['ConstantShape', Shape]) -> Union['ConstantShape', AddShape]:
        if isinstance(other, ConstantShape):
            return ConstantShape(self.c + other.c)
        else:
            return super().__add__(other)

    def __radd__(self, other: Union['ConstantShape', Shape]) -> Union['ConstantShape', AddShape]:
        if isinstance(other, ConstantShape):
            return ConstantShape(self.c + other.c)
        else:
            return super().__radd__(other)

    def __sub__(self, other: Union['ConstantShape', Shape]) -> Union['ConstantShape', SubShape]:
        if isinstance(other, ConstantShape):
            return ConstantShape(self.c - other.c)
        else:
            return super().__sub__(other)

    def __rsub__(self, other: Union['ConstantShape', Shape]) -> Union['ConstantShape', SubShape]:
        if isinstance(other, ConstantShape):
            return ConstantShape(other.c - self.c)
        else:
            return super().__rsub__(other)

    def __mul__(self, other: Union[float, int, Shape]) -> Union['ConstantShape', MulShape]:
        if isinstance(other, (float, int)):
            return ConstantShape(self.c*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[float, int, Shape]) -> Union['ConstantShape', MulShape]:
        if isinstance(other, (float, int)):
            return ConstantShape(other*self.c)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[float, int, Shape]) -> Union['ConstantShape', DivShape]:
        if isinstance(other, (float, int)):
            return ConstantShape(self.c/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[float, int, Shape]) -> Union['ConstantShape', DivShape]:
        if isinstance(other, (float, int)):
            return ConstantShape(self.c/other)
        else:
            return super().__truediv__(other)

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        if isinstance(s, float):
            result = self.c
        else:
            result = full(s.shape, self.c)
        return result

    def __str__(self) -> str:
        return f'ConstantShape({self.c})'


class TaperedShape(Shape):

    c_root: float = None
    c_tip: float = None

    def __init__(self, c_root: float, c_tip: float) -> None:
        self.c_root = c_root
        self.c_tip = c_tip

    @property
    def area(self) -> None:
        if self._area is None:
            self._area = self.c_root + self.c_tip
        return self._area

    def __add__(self, other: Union['TaperedShape', Shape]) -> Union['TaperedShape', AddShape]:
        if isinstance(other, TaperedShape):
            return TaperedShape(self.c_root + other.c_root, self.c_tip + other.c_tip)
        else:
            return super().__add__(other)

    def __radd__(self, other: Union['TaperedShape', Shape]) -> Union['TaperedShape', AddShape]:
        if isinstance(other, TaperedShape):
            return TaperedShape(self.c_root + other.c_root, self.c_tip + other.c_tip)
        else:
            return super().__radd__(other)

    def __sub__(self, other: Union['TaperedShape', Shape]) -> Union['TaperedShape', SubShape]:
        if isinstance(other, TaperedShape):
            return TaperedShape(self.c_root - other.c_root, self.c_tip - other.c_tip)
        else:
            return super().__sub__(other)

    def __rsub__(self, other: Union['TaperedShape', Shape]) -> Union['TaperedShape', SubShape]:
        if isinstance(other, TaperedShape):
            return TaperedShape(other.c_root - self.c_root, other.c_tip - self.c_tip)
        else:
            return super().__rsub__(other)

    def __mul__(self, other: Union[float, ndarray, Shape]) -> Union['TaperedShape', MulShape]:
        if isinstance(other, (float, int)):
            return TaperedShape(self.c_root*other, self.c_tip*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[float, ndarray, Shape]) -> Union['TaperedShape', MulShape]:
        if isinstance(other, (float, int)):
            return TaperedShape(other*self.c_root, other*self.c_tip)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[float, ndarray, Shape]) -> Union['TaperedShape', DivShape]:
        if isinstance(other, (float, int)):
            return TaperedShape(self.c_root/other, self.c_tip/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[float, ndarray, Shape]) -> Union['TaperedShape', DivShape]:
        if isinstance(other, (float, int)):
            return TaperedShape(self.c_root/other, self.c_tip/other)
        else:
            return super().__truediv__(other)

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        return self.c_root - absolute(s)*(self.c_root - self.c_tip)

    def __str__(self) -> str:
        return f'TaperedShape({self.c_root}, {self.c_tip})'


class GeneralShape(Shape):

    An: ndarray = None
    _n: Iterable[int] = None

    def __init__(self, An: ndarray) -> None:
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

    def __call__(self, s: Union[ndarray, float]) -> Union[ndarray, float]:
        th = arccos(s)
        if isinstance(s, float):
            result = 0.0
        else:
            result = zeros(s.shape)
        for n, An in self.items():
            result += An*sin(n*th)
        return result

    @property
    def size(self) -> int:
        return self.An.size

    def __getitem__(self, n: Union[int, slice]) -> Union[float, ndarray]:
        return self.An[n]

    def __setitem__(self, n: Union[int, slice], An: Union[float, ndarray]) -> None:
        self.An[n] = An

    def __mul__(self, other: Union[float, ndarray, Shape]) -> Union['GeneralShape', MulShape]:
        if isinstance(other, (int, float, ndarray)):
            return GeneralShape(self.An*other)
        else:
            return super().__mul__(other)

    def __rmul__(self, other: Union[float, ndarray, Shape]) -> Union['GeneralShape', MulShape]:
        if isinstance(other, (int, float, ndarray)):
            return GeneralShape(other*self.An)
        else:
            return super().__rmul__(other)

    def __div__(self, other: Union[float, ndarray, Shape]) -> Union['GeneralShape', DivShape]:
        if isinstance(other, (float, ndarray)):
            return GeneralShape(self.An/other)
        else:
            return super().__div__(other)

    def __truediv__(self, other: Union[float, ndarray, Shape]) -> Union['GeneralShape', DivShape]:
        if isinstance(other, (float, ndarray)):
            return GeneralShape(self.An/other)
        else:
            return super().__truediv__(other)

    def __add__(self, other: Union['GeneralShape', Shape]) -> Union['GeneralShape', AddShape]:
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

    def __radd__(self, other: Union['GeneralShape', Shape]) -> Union['GeneralShape', AddShape]:
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

    def __sub__(self, other: Union['GeneralShape', Shape]) -> Union['GeneralShape', SubShape]:
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

    def __rsub__(self, other: Union['GeneralShape', Shape]) -> Union['GeneralShape', SubShape]:
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

    def values(self) -> ndarray:
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
