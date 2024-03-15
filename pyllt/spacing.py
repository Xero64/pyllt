from typing import TYPE_CHECKING, Dict

from numpy import arccos, cos, flip, hstack, linspace, pi, sin, asarray

if TYPE_CHECKING:
    from numpy import ndarray

class Spacing():

    num: int = None
    sol: bool = None
    _s: 'ndarray' = None
    _th: 'ndarray' = None
    _cos_n_th: Dict[int, 'ndarray'] = None
    _sin_n_th: Dict[int, 'ndarray'] = None

    def __init__(self, num: int = 101, sol: bool = False) -> None:
        self.num = num
        self.sol = sol

    @property
    def s(self) -> 'ndarray':
        if self._s is None:
            if self.sol:
                th = linspace(pi, 0.0, self.num)
                s = cos(th)
                self._s = s[1:-1]
            else:
                num = self.num//2 + 1
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


class CosineSpacing(Spacing):

    pass
