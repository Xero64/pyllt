#%%
# Import Dependencies
from typing import TYPE_CHECKING, Iterable, Union

from numpy import absolute, arccos, full, pi, sin, zeros

if TYPE_CHECKING:
    from numpy import ndarray
    from numpy.typing import DTypeLike


class Distribution():

    _area: float = None

    def __init__(self) -> None:
        pass

    @property
    def area(self) -> None:
        return self._area

    def __call__(self) -> None:
        return None



class ConstantSpan(Distribution):

    c: float = None

    def __init__(self, c: float) -> None:
        self.c = c

    @property
    def area(self) -> None:
        if self._area is None:
            self._area = 2*self.c
        return self._area

    def __call__(self, s: Union['ndarray', float]) -> Union['ndarray', float]:
        if isinstance(s, float):
            result = self.c
        else:
            result = full(s.shape, self.c)
        return result


class TaperedSpan(Distribution):

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

    def __call__(self, s: Union['ndarray', float]) -> Union['ndarray', float]:
        return self.c_root - absolute(s)*(self.c_root - self.c_tip)


class GeneralSpan(Distribution):

    An: 'ndarray' = None
    _n: Iterable[int] = None

    def __init__(self, An: 'ndarray') -> None:
        self.An = An

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

    def __call__(self, s: Union['ndarray', float]) -> Union['ndarray', float]:
        th = arccos(s)
        if isinstance(s, float):
            result = 0.0
        else:
            result = zeros(s.shape)
        for i, Ai in enumerate(self.An):
            result += Ai*sin((i+1)*th)
        return result

    @property
    def size(self) -> int:
        return self.An.size

    @property
    def shape(self) -> 'ndarray':
        return self.An.shape

    @property
    def ndim(self) -> int:
        return self.An.ndim

    @property
    def dtype(self) -> 'DTypeLike':
        return self.An.dtype

    def __getitem__(self, n: Union[int, slice]) -> Union[float, 'ndarray']:
        return self.An[n]

    def __setitem__(self, n: Union[int, slice], An: Union[float, 'ndarray']) -> None:
        self.An[n] = An

    def __mul__(self, other: Union[float, 'ndarray']) -> 'GeneralSpan':
        return GeneralSpan(self.An*other)

    def __rmul__(self, other: Union[float, 'ndarray']) -> 'GeneralSpan':
        return GeneralSpan(self.An*other)

    def __add__(self, other: 'GeneralSpan') -> 'GeneralSpan':
        return GeneralSpan(self.An + other.An)

    def __radd__(self, other: 'GeneralSpan') -> 'GeneralSpan':
        return GeneralSpan(other.An + self.An)

    def __pow__(self, other: float) -> 'GeneralSpan':
        return GeneralSpan(self.An**other)

    def sum(self) -> float:
        return self.An.sum()

    def to_constant(self) -> ConstantSpan:
        return ConstantSpan(self.area/2)

    def keys(self) -> Iterable[int]:
        return self.n

    def values(self) -> 'ndarray':
        return self.An

    def items(self) -> Iterable[tuple[int, float]]:
        return zip(self.n, self.An)

    def __repr__(self) -> str:
        return f'<GeneralSpan({self.An})>'

    def __str__(self) -> str:
        return f'GeneralSpan: An={self.An}'


class EllipticalSpan(GeneralSpan):

    def __init__(self, scale: float = 1.0) -> None:
        An = zeros(1)
        An[0] = scale
        super().__init__(An)

    def __add__(self, other: GeneralSpan) -> GeneralSpan:
        An = other.An
        An[0] += self.An[0]
        return GeneralSpan(An)

    def __radd__(self, other: GeneralSpan) -> GeneralSpan:
        An = other.An
        An[0] += self.An[0]
        return GeneralSpan(An)

    def __repr__(self) -> str:
        return f'<EllipticalSpan({self.An})>'

    def __str__(self) -> str:
        return f'EllipticalSpan: An={self.An}'


class BellSpan(GeneralSpan):

    def __init__(self, scale: float = 1.0) -> None:
        An = zeros(3)
        An[0] = scale
        An[2] = -scale/3
        super().__init__(An)

    def __repr__(self) -> str:
        return f'<BellSpan({self.An})>'

    def __str__(self) -> str:
        return f'BellSpan: An={self.An}'
