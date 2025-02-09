from typing import TYPE_CHECKING

from numpy import absolute, arccos, cos, flip, hstack, linspace, pi, sin

if TYPE_CHECKING:
    from numpy.typing import NDArray

class Spacing():

    num: int = None
    sol: bool = None
    _s: 'NDArray' = None
    _th: 'NDArray' = None
    _cos_n_th: dict[int, 'NDArray'] = None
    _sin_n_th: dict[int, 'NDArray'] = None
    _n_sin_n_th: dict[int, 'NDArray'] = None
    _n_cos_n_th: dict[int, 'NDArray'] = None
    _n2_sin_n_th: dict[int, 'NDArray'] = None
    _n2_cos_n_th: dict[int, 'NDArray'] = None
    _abs_s: 'NDArray' = None

    def __init__(self, num: int = 101, sol: bool = False) -> None:
        self.num = num
        self.sol = sol

    @property
    def s(self) -> 'NDArray':
        if self._s is None:
            if self.sol:
                th = linspace(pi, 0.0, self.num+2)
                s = cos(th)
                self._s = s[1:-1]
            else:
                num = self.num//2 + 1
                th = linspace(pi/2, 0.0, num)
                s = cos(th)
                self._s = hstack((-flip(s), s))
        return self._s

    @property
    def abs_s(self) -> 'NDArray':
        if self._abs_s is None:
            self._abs_s = absolute(self.s)
        return self._abs_s

    @property
    def th(self) -> 'NDArray':
        if self._th is None:
            self._th = arccos(self.s)
        return self._th

    @property
    def shape(self) -> int | tuple[int, ...]:
        return self.s.shape

    def cos_n_th(self, n: int) -> 'NDArray':
        if self._cos_n_th is None:
            self._cos_n_th = {}
        if n not in self._cos_n_th:
            self._cos_n_th[n] = cos(n*self.th)
        return self._cos_n_th[n]

    def sin_n_th(self, n: int) -> 'NDArray':
        if self._sin_n_th is None:
            self._sin_n_th = {}
        if n not in self._sin_n_th:
            self._sin_n_th[n] = sin(n*self.th)
        return self._sin_n_th[n]

    @property
    def costh(self) -> 'NDArray':
        return self.cos_n_th(1)

    @property
    def sinth(self) -> 'NDArray':
        return self.sin_n_th(1)

    def n_cos_n_th(self, n: int) -> 'NDArray':
        if self._n_cos_n_th is None:
            self._n_cos_n_th = {}
        if n not in self._n_cos_n_th:
            self._n_cos_n_th[n] = n*self.cos_n_th(n)
        return self._n_cos_n_th[n]

    def n_sin_n_th(self, n: int) -> 'NDArray':
        if self._n_sin_n_th is None:
            self._n_sin_n_th = {}
        if n not in self._n_sin_n_th:
            self._n_sin_n_th[n] = n*self.sin_n_th(n)
        return self._n_sin_n_th[n]

    def n2_cos_n_th(self, n: int) -> 'NDArray':
        if self._n2_cos_n_th is None:
            self._n2_cos_n_th = {}
        if n not in self._n2_cos_n_th:
            self._n2_cos_n_th[n] = n*self.n_cos_n_th(n)
        return self._n2_cos_n_th[n]

    def n2_sin_n_th(self, n: int) -> 'NDArray':
        if self._n2_sin_n_th is None:
            self._n2_sin_n_th = {}
        if n not in self._n2_sin_n_th:
            self._n2_sin_n_th[n] = n*self.n_sin_n_th(n)
        return self._n2_sin_n_th[n]

    def stringify(self) -> str:
        return f'arccos(s)'


class CosineSpacing(Spacing):

    pass
