"""
Calcium includes a simple Python interface implemented using ``ctypes``.

Introduction
------------------------------------------------------------------------

Setup and usage
..................

Make sure the Calcium library and its dependencies are built
and in the path of the system's dynamic library loader.
Then make sure that ``<calcium_source_dir>/pycalcium`` is in the
Python path, for example by adding it to ``PYTHONPATH``,
adding it to ``sys.path``,
or simply starting Python inside the ``pycalcium`` directory.

Import the module and run perform a calculation:

    >>> import calcium
    >>> calcium.ca(1) / 3
    0.333333 {1/3}
    >>> calcium.exp(calcium.pi * calcium.i / 2)
    1.00000*I {a where a = I [a^2+1=0]}

If you don't mind polluting the global namespace, import everything:

    >>> from calcium import *
    >>> exp(pi*i/2) + ca(1)/3
    0.333333 + 1.00000*I {(3*a+1)/3 where a = I [a^2+1=0]}

Current limitations
.....................

* Does not support creating new context objects or modifying
  the settings of a context object.
* Leaks memory (for example, when printing).
* Because ``ctypes`` is used, this is not as efficient as a
  Cython wrapper. This interface should be used for testing
  and not for absolute performance.
* Does not wrap various functions.

API documentation
------------------------------------------------------------------------

Functions are available both as regular functions and as methods
on the ``ca`` class.

"""

import ctypes
import ctypes.util
import sys

libcalcium_path = ctypes.util.find_library('calcium')
libcalcium = ctypes.CDLL(libcalcium_path)

libarb_path = ctypes.util.find_library('arb')
libarb = ctypes.CDLL(libarb_path)

libflint_path = ctypes.util.find_library('flint')
libflint = ctypes.CDLL(libflint_path)

T_TRUE = 0
T_FALSE = 1
T_UNKNOWN = 2

class _fmpz_struct(ctypes.Structure):
    _fields_ = [('val', ctypes.c_long)]


class ca_struct(ctypes.Structure):
    """Low-level wrapper for ca_struct, for internal use by ctypes."""
    _fields_ = [('data', ctypes.c_long * 5)]

class ca_ctx_struct(ctypes.Structure):
    """Low-level wrapper for ca_ctx_struct, for internal use by ctypes."""
    # todo: use the real size
    _fields_ = [('content', ctypes.c_ulong * 64)]

class ca_mat_struct(ctypes.Structure):
    """Low-level wrapper for ca_mat_struct, for internal use by ctypes."""
    _fields_ = [('entries', ctypes.c_void_p),
                ('r', ctypes.c_long),
                ('c', ctypes.c_long),
                ('rows', ctypes.c_void_p)]


class ca_ctx:
    """
    Python class wrapping the ca_ctx_t context object.
    Currently only supports a global instance.
    """

    def __init__(self):
        self._data = ca_ctx_struct()
        self._ref = ctypes.byref(self._data)
        libcalcium.ca_ctx_init(self._ref)

    def __del__(self):
        # Python calls ctx_default.__del__ before all ca instances are
        # destroyed, despite the fact that the ca instances contain
        # references to ctx_default.
        # WHY??????????????
        # libcalcium.ca_ctx_clear(self._ref)
        pass

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

ctx_default = ca_ctx()


class ca:
    """
    Python class wrapping the ca_t type for numbers.

    Examples::

        >>> ca(1)
        1
        >>> ca()
        0
        >>> ca(0)
        0
        >>> ca(-5)
        -5
        >>> ca(2.25)
        2.25000 {9/4}
        >>> ca(1) / 3
        0.333333 {1/3}
        >>> (-1) ** ca(0.5)
        1.00000*I {a where a = I [a^2+1=0]}
        >>> ca(1-2j)
        1.00000 - 2.00000*I {-2*a+1 where a = I [a^2+1=0]}
        >>> ca(0.1)           # be careful with float input!
        0.100000 {3602879701896397/36028797018963968}
        >>> ca(float("inf"))
        +Infinity
        >>> ca(float("nan"))
        Unknown
        >>> 3 < pi < ca(22)/7
        True

    """

    def __init__(self, val=0):
        self._ctx_python = ctx_default
        self._ctx = self._ctx_python._ref
        self._data = ca_struct()
        self._ref = ctypes.byref(self._data)
        libcalcium.ca_init(self, self._ctx)
        if val is not 0:
            typ = type(val)
            if typ is int:
                b = sys.maxsize
                if -b <= val <= b:
                    libcalcium.ca_set_si(self, val, self._ctx)
                else:
                    n = _fmpz_struct()
                    nref = ctypes.byref(n)
                    libflint.fmpz_init(nref)
                    libflint.fmpz_set_str(nref, ctypes.c_char_p(str(val).encode('ascii')), 10)
                    libcalcium.ca_set_fmpz(self, nref, self._ctx)
                    libflint.fmpz_clear(nref)
            elif typ is ca:
                libcalcium.ca_set(self, val, self._ctx)
            elif typ is float:
                libcalcium.ca_set_d(self, val, self._ctx)
            elif typ is complex:
                libcalcium.ca_set_d_d(self, val.real, val.imag, self._ctx)
            else:
                raise TypeError

    def __del__(self):
        libcalcium.ca_clear(self, self._ctx)

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

    @staticmethod
    def inf():
        """
        The special value positive infinity.

        Examples::

            >>> inf == ca.inf()    # alias for the method
            True
            >>> inf
            +Infinity
            >>> -inf
            -Infinity
            >>> abs(-inf)
            +Infinity
            >>> inf + inf
            +Infinity
            >>> (-inf) + (-inf)
            -Infinity
            >>> -inf * inf
            -Infinity
            >>> inf / inf
            Undefined
            >>> 1 / inf
            0
            >>> re(inf)
            Undefined
            >>> im(inf)
            Undefined
            >>> sign(inf)
            1
            >>> sign((1+i)*inf) == (1+i)/sqrt(2)
            True

        """
        res = ca()
        libcalcium.ca_pos_inf(res, res._ctx)
        return res

    @staticmethod
    def uinf():
        """
        The special value unsigned infinity.

        Examples::

            >>> uinf == ca.uinf()    # alias for the method
            True
            >>> uinf
            UnsignedInfinity
            >>> abs(uinf)
            +Infinity
            >>> -3 * uinf
            UnsignedInfinity
            >>> 1/uinf
            0
            >>> sign(uinf)
            Undefined
            >>> uinf * uinf
            UnsignedInfinity
            >>> uinf / uinf
            Undefined
            >>> uinf + uinf
            Undefined
            >>> re(uinf)
            Undefined
            >>> im(uinf)
            Undefined

        """
        res = ca()
        libcalcium.ca_uinf(res, res._ctx)
        return res

    @staticmethod
    def undefined():
        """
        The special value undefined.

        Examples::

            >>> undefined == ca.undefined()    # alias for the method
            True
            >>> undefined
            Undefined
            >>> undefined == undefined
            True
            >>> undefined * 0
            Undefined

        """
        res = ca()
        libcalcium.ca_undefined(res, res._ctx)
        return res

    @staticmethod
    def unknown():
        """
        The special meta-value unknown.
        This meta-value is not comparable.

        Examples::

            >>> unknown == unknown
            Traceback (most recent call last):
              ...
            ValueError: unable to decide predicate: equal
            >>> unknown == 0
            Traceback (most recent call last):
              ...
            ValueError: unable to decide predicate: equal
            >>> unknown == undefined
            Traceback (most recent call last):
              ...
            ValueError: unable to decide predicate: equal
            >>> unknown + unknown
            Unknown
            >>> unknown + 3
            Unknown
            >>> unknown + undefined
            Undefined

        """
        res = ca()
        libcalcium.ca_unknown(res, res._ctx)
        return res

    @staticmethod
    def pi():
        """
        The constant pi.

        Examples::

            >>> pi == ca.pi()    # alias for the method
            True
            >>> pi
            3.14159 {a where a = 3.14159 [Pi]}
            >>> sin(pi/6)
            0.500000 {1/2}
            >>> (pi - 3) ** 3
            0.00283872 {a^3-9*a^2+27*a-27 where a = 3.14159 [Pi]}

        """
        res = ca()
        libcalcium.ca_pi(res, res._ctx)
        return res

    @staticmethod
    def euler():
        """
        Euler's constant (gamma).

        Examples::

            >>> euler == ca.euler()    # alias for the method
            True
            >>> euler
            0.577216 {a where a = 0.577216 [Euler]}
            >>> exp(euler)
            1.78107 {a where a = 1.78107 [Exp(0.577216 {b})], b = 0.577216 [Euler]}

        """
        res = ca()
        libcalcium.ca_euler(res, res._ctx)
        return res

    @staticmethod
    def i():
        """
        The imaginary unit.

        Examples::

            >>> i == ca.i()    # alias for the method
            True
            >>> i == I == j    # extra aliases for convenience
            True
            >>> i
            1.00000*I {a where a = I [a^2+1=0]}
            >>> i**2
            -1
            >>> i**3
            -1.00000*I {-a where a = I [a^2+1=0]}
            >>> abs(i)
            1
            >>> sign(i)
            1.00000*I {a where a = I [a^2+1=0]}
            >>> abs(sqrt(1+i) / sqrt(1-i))
            1
            >>> re(i), im(i)
            (0, 1)

        """
        res = ca()
        libcalcium.ca_i(res, res._ctx)
        return res

    def __repr__(self):
        s = libcalcium.ca_get_str(self, self._ctx)
        res = str(s, 'ascii')
        #libflint.flint_free(s)
        return res

    def __str__(self):
        return self.__repr__()

    def __bool__(self):
        """

        Does it work?

        """
        t = libcalcium.ca_check_is_zero(self, self._ctx)
        if t == T_TRUE:
            return False
        if t == T_FALSE:
            return True
        raise ValueError("unable to decide predicate: is_zero")

    def __eq__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        t = libcalcium.ca_check_equal(self, other, self._ctx)
        if t == T_TRUE:
            return True
        if t == T_FALSE:
            return False
        raise ValueError("unable to decide predicate: equal")

    def __ne__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        t = libcalcium.ca_check_equal(self, other, self._ctx)
        if t == T_TRUE:
            return False
        if t == T_FALSE:
            return True
        raise ValueError("unable to decide predicate: equal")

    def __le__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        t = libcalcium.ca_check_le(self, other, self._ctx)
        if t == T_TRUE:
            return True
        if t == T_FALSE:
            return False
        raise ValueError("unable to decide predicate: le")

    def __lt__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        t = libcalcium.ca_check_lt(self, other, self._ctx)
        if t == T_TRUE:
            return True
        if t == T_FALSE:
            return False
        raise ValueError("unable to decide predicate: lt")

    def __ge__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        t = libcalcium.ca_check_ge(self, other, self._ctx)
        if t == T_TRUE:
            return True
        if t == T_FALSE:
            return False
        raise ValueError("unable to decide predicate: ge")

    def __gt__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        t = libcalcium.ca_check_gt(self, other, self._ctx)
        if t == T_TRUE:
            return True
        if t == T_FALSE:
            return False
        raise ValueError("unable to decide predicate: gt")

    def __add__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_add(res, self, other, self._ctx)
        return res

    __radd__ = __add__

    def __sub__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_sub(res, self, other, self._ctx)
        return res

    def __rsub__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_sub(res, other, self, self._ctx)
        return res

    def __mul__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_mul(res, self, other, self._ctx)
        return res

    __rmul__ = __mul__

    def __truediv__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_div(res, self, other, self._ctx)
        return res

    def __rtruediv__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_div(res, other, self, self._ctx)
        return res

    def __floordiv__(self, other):
        return (self / other).floor()

    def __rfloordiv__(self, other):
        return (other / self).floor()

    def __pow__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_pow(res, self, other, self._ctx)
        return res

    def __rpow__(self, other):
        if type(self) is not type(other):
            try:
                other = ca(other)
            except TypeError:
                return NotImplemented
        if self._ctx_python is not self._ctx_python:
            raise ValueError("different context objects!")
        res = ca()
        libcalcium.ca_pow(res, other, self, self._ctx)
        return res

    def __abs__(self):
        res = ca()
        libcalcium.ca_abs(res, self, self._ctx)
        return res

    def __neg__(self):
        res = ca()
        libcalcium.ca_neg(res, self, self._ctx)
        return res

    def __pos__(self):
        res = ca()
        libcalcium.ca_set(res, self, self._ctx)
        return res

    def re(self):
        """
        Real part.

        Examples::

            >>> re(2+3j) == ca(2+3j).re()
            True
            >>> re(2+3*i)
            2

        """
        res = ca()
        libcalcium.ca_re(res, self, self._ctx)
        return res

    def im(self):
        """
        Imaginary part.

        Examples::

            >>> im(2+3j) == ca(2+3j).im()    # alias for the method
            True
            >>> im(2+3*i)
            3

        """
        res = ca()
        libcalcium.ca_im(res, self, self._ctx)
        return res

    def conj(self):
        """
        Complex conjugate.

        Examples::

            >>> conj(1j) == conjugate(1j) == ca(1j).conj() == ca(1j).conjugate()    # alias for the method
            True
            >>> conj(2+i)
            2.00000 - 1.00000*I {-a+2 where a = I [a^2+1=0]}
            >>> conj(pi)
            3.14159 {a where a = 3.14159 [Pi]}

        """
        res = ca()
        libcalcium.ca_conjugate(res, self, self._ctx)
        return res

    conjugate = conj

    def floor(self):
        """
        Floor function.

        Examples::

            >>> floor(3) == ca(3).floor()    # alias for the method
            True
            >>> floor(pi)
            3
            >>> floor(-pi)
            -4

        """
        res = ca()
        libcalcium.ca_floor(res, self, self._ctx)
        return res

    def ceil(self):
        """
        Ceiling function.

        Examples::

            >>> ceil(3) == ca(3).ceil()    # alias for the method
            True
            >>> ceil(pi)
            4
            >>> ceil(-pi)
            -3

        """
        res = ca()
        libcalcium.ca_ceil(res, self, self._ctx)
        return res

    def sgn(self):
        """
        Sign function.

        Examples::

            >>> sgn(2) == sign(2) == ca(2).sgn()    # aliases for the method
            True
            >>> sign(0)
            0
            >>> sign(sqrt(2))
            1
            >>> sign(-sqrt(2))
            -1
            >>> sign(-sqrt(-2))
            -1.00000*I {-a where a = I [a^2+1=0]}

        """
        res = ca()
        libcalcium.ca_sgn(res, self, self._ctx)
        return res

    sign = sgn

    def sqrt(self):
        """
        Principal square root.

        Examples::

            >>> sqrt(2) == ca(2).sqrt()    # alias for the method
            True
            >>> sqrt(0)
            0
            >>> sqrt(1)
            1
            >>> sqrt(2)
            1.41421 {a where a = 1.41421 [a^2-2=0]}
            >>> sqrt(-1)
            1.00000*I {a where a = I [a^2+1=0]}
            >>> sqrt(inf)
            +Infinity
            >>> sqrt(-inf)
            +I * Infinity
            >>> sqrt(uinf)
            UnsignedInfinity
            >>> sqrt(undefined)
            Undefined
            >>> sqrt(unknown)
            Unknown

        """
        res = ca()
        libcalcium.ca_sqrt(res, self, self._ctx)
        return res

    def log(self):
        """
        Natural logarithm.

        Examples::

            >>> log(2) == ca(2).log()    # alias for the method
            True
            >>> log(1)
            0
            >>> log(exp(2))
            2
            >>> log(2)
            0.693147 {a where a = 0.693147 [Log(2)]}
            >>> log(4) == 2*log(2)
            True
            >>> log(1+sqrt(2)) / log(3+2*sqrt(2))
            0.500000 {1/2}
            >>> log(ca(10)**(10**30)) / log(10)
            1.00000e+30 {1000000000000000000000000000000}
            >>> log(-1)
            3.14159*I {a where a = 3.14159*I [Log(-1)]}
            >>> log(i)
            1.57080*I {(a*b)/2 where a = 3.14159 [Pi], b = I [b^2+1=0]}
            >>> log(0)
            -Infinity
            >>> log(inf)
            +Infinity
            >>> log(-inf)
            +Infinity
            >>> log(uinf)
            +Infinity
            >>> log(undefined)
            Undefined
            >>> log(unknown)
            Unknown

        """
        res = ca()
        libcalcium.ca_log(res, self, self._ctx)
        return res

    def exp(self):
        """
        Exponential function.

        Examples::

            >>> exp(0)
            1
            >>> exp(1)
            2.71828 {a where a = 2.71828 [Exp(1)]}
            >>> exp(-1)
            0.367879 {a where a = 0.367879 [Exp(-1)]}
            >>> exp(-1) * exp(1)
            1
            >>> exp(7*pi*i/2)
            -1.00000*I {-a where a = I [a^2+1=0]}
            >>> exp(inf)
            +Infinity
            >>> exp(-inf)
            0
            >>> exp(uinf)
            Undefined
            >>> exp(undefined)
            Undefined
            >>> exp(unknown)
            Unknown

        """
        res = ca()
        libcalcium.ca_exp(res, self, self._ctx)
        return res

    def erf(self):
        """
        Error function.

        Examples::

            >>> erf(0)
            0
            >>> erf(1)
            0.842701 {a where a = 0.842701 [Erf(1)]}
            >>> erf(inf)
            1
            >>> erf(-inf)
            -1
            >>> erf(i*inf)
            +I * Infinity
            >>> erf(-i*inf)
            -I * Infinity
            >>> erf(uinf)
            Undefined

        """
        res = ca()
        libcalcium.ca_erf(res, self, self._ctx)
        return res

    def erfc(self):
        """
        Complementary error function.

        Examples::

            >>> erfc(inf)
            0
            >>> erfc(-inf)
            2
            >>> erfc(1000)
            1.86004e-434298 {a where a = 1.86004e-434298 [Erfc(1000)]}
            >>> erfc(i*inf)
            -I * Infinity
            >>> erfc(-i*inf)
            +I * Infinity
            >>> erfc(sqrt(2)) + erf(sqrt(2))
            1
            >>> erfc(uinf)
            Undefined

        """
        res = ca()
        libcalcium.ca_erfc(res, self, self._ctx)
        return res

    def erfi(self):
        """
        Imaginary error function.

        Examples::

            >>> erfi(0)
            0
            >>> erfi(1)
            1.65043 {a where a = 1.65043 [Erfi(1)]}
            >>> erfi(inf)
            +Infinity
            >>> erfi(-inf)
            -Infinity
            >>> erfi(i*inf)
            1.00000*I {a where a = I [a^2+1=0]}
            >>> erfi(-i*inf)
            -1.00000*I {-a where a = I [a^2+1=0]}
            >>> erf(2)**2 + erfi(i*2)**2
            0

        """
        res = ca()
        libcalcium.ca_erfi(res, self, self._ctx)
        return res

    def gamma(self):
        """
        Gamma function.

        Examples::

            >>> [gamma(n) for n in range(1,11)]
            [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880]
            >>> gamma(0)
            UnsignedInfinity
            >>> 1 / gamma(0)
            0
            >>> gamma(0.5)
            1.77245 {a where a = 1.77245 [Sqrt(3.14159 {b})], b = 3.14159 [Pi]}
            >>> gamma(0.5)**2 == pi
            True
            >>> pi * gamma(pi) / gamma(pi+1)
            1

        """
        res = ca()
        libcalcium.ca_gamma(res, self, self._ctx)
        return res

class ca_mat:
    """
    Python class wrapping the ca_mat_t type for matrices.

        >>> ca_mat(2,3)
        ca_mat of size 2 x 3
        [0, 0, 0]
        [0, 0, 0]
        >>> ca_mat([[1,2],[3,4],[5,6]])
        ca_mat of size 3 x 2
        [1, 2]
        [3, 4]
        [5, 6]
        >>> ca_mat(2, 5, range(10))
        ca_mat of size 2 x 5
        [0, 1, 2, 3, 4]
        [5, 6, 7, 8, 9]

    """

    def __init__(self, *args):
        self._ctx_python = ctx_default
        self._ctx = self._ctx_python._ref
        self._data = ca_mat_struct()
        self._ref = ctypes.byref(self._data)

        if len(args) == 1:
            val = args[0]
            if isinstance(val, ca_mat):
                m = val.nrows()
                n = val.ncols()
                libcalcium.ca_mat_init(self, m, n, self._ctx)
                libcalcium.ca_mat_set(self, val, self._ctx)
            elif isinstance(val, (list, tuple)):
                m = len(val)
                n = 0
                if m != 0:
                    if not isinstance(val[0], (list, tuple)):
                        raise TypeError("single input to ca_mat must be a list of lists")
                    n = len(val[0])
                    for i in range(1, m):
                        if len(val[i]) != n:
                            raise ValueError("input rows have different lengths")
                libcalcium.ca_mat_init(self, m, n, self._ctx)
                for i in range(m):
                    row = val[i]
                    for j in range(n):
                        x = ca(row[j])
                        libcalcium.ca_set(libcalcium.ca_mat_entry_ptr(self, i, j, self._ctx), x, self._ctx)
            else:
                raise TypeError("cannot create ca_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            assert m >= 0
            assert n >= 0
            libcalcium.ca_mat_init(self, m, n, self._ctx)
        elif len(args) == 3:
            m, n, entries = args
            assert m >= 0
            assert n >= 0
            libcalcium.ca_mat_init(self, m, n, self._ctx)
            entries = list(entries)
            if len(entries) != m*n:
                raise ValueError("list of entries has the wrong length")
            for i in range(m):
                for j in range(n):
                    x = ca(entries[i*n + j])
                    libcalcium.ca_set(libcalcium.ca_mat_entry_ptr(self, i, j, self._ctx), x, self._ctx)
        else:
            raise ValueError("ca_mat: expected 1-3 arguments")

    def __del__(self):
        libcalcium.ca_mat_clear(self, self._ctx)

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

    def nrows(self):
        return self._data.r

    def ncols(self):
        return self._data.c

    def entries(self):
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i in range(m):
            for j in range(n):
                L[i*n + j] = self[i, j]
        return L

    def __iter__(self):
        m = self.nrows()
        n = self.ncols()
        for i in range(m):
            for j in range(n):
                yield self[i, j]

    def table(self):
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n : (i+1)*n] for i in range(m)]

    # supports mpmath conversions
    tolist = table

    def __repr__(self):
        s = "ca_mat of size %i x %i\n" % (self.nrows(), self.ncols())
        nrows = self.nrows()
        ncols = self.ncols()

        def matrix_to_str(tab):
            if len(tab) == 0 or len(tab[0]) == 0:
                return "[]"
            tab = [[str(c) for c in row] for row in tab]
            widths = []
            for i in range(len(tab[0])):
                w = max([len(row[i]) for row in tab])
                widths.append(w)
            for i in range(len(tab)):
                tab[i] = [s.rjust(widths[j]) for j, s in enumerate(tab[i])]
                tab[i] = "[" + (", ".join(tab[i])) + "]"
            return "\n".join(tab)

        return s + matrix_to_str(self.table())

    __str__ = __repr__

    def __getitem__(self, ij):
        i, j = ij
        nrows = self.nrows()
        ncols = self.ncols()
        assert 0 <= i < nrows
        assert 0 <= j < ncols
        x = ca()
        libcalcium.ca_set(x, libcalcium.ca_mat_entry_ptr(self, i, j, self._ctx), self._ctx)
        return x

    def __setitem__(self, ij, x):
        i, j = ij
        nrows = self.nrows()
        ncols = self.ncols()
        assert 0 <= i < nrows
        assert 0 <= j < ncols
        x = ca(x)
        libcalcium.ca_set(libcalcium.ca_mat_entry_ptr(self, i, j, self._ctx), x, self._ctx)


def re(x):
    return ca(x).re()

def im(x):
    return ca(x).im()

def sgn(x):
    return ca(x).sgn()

def sign(x):
    return ca(x).sign()

def floor(x):
    return ca(x).floor()

def ceil(x):
    return ca(x).ceil()

def conj(x):
    return ca(x).conj()

def conjugate(x):
    return ca(x).conjugate()

def sqrt(x):
    return ca(x).sqrt()

def log(x):
    return ca(x).log()

def exp(x):
    return ca(x).exp()

def erf(x):
    return ca(x).erf()

def erfc(x):
    return ca(x).erfc()

def erfi(x):
    return ca(x).erfi()

def gamma(x):
    return ca(x).gamma()

def fac(x):
    """
    Alias for gamma(x+1).

    Examples::

        >>> fac(10)
        3.62880e+6 {3628800}

    """
    return (ca(x)+1).gamma()

def cos(x):
    """
    The cosine function is not yet implemented in Calcium.
    This placeholder function evaluates the cosine function
    using complex exponentials.

    Examples::

        >>> cos(0)
        1
        >>> cos(pi)
        -1
        >>> cos(pi/2)
        0
        >>> cos(pi/3)
        0.500000 {1/2}
        >>> cos(pi/6)**2
        0.750000 {3/4}
        >>> cos(1)**2 + sin(1)**2
        1

    """
    ix = ca(x)*i
    y = exp(ix)
    return (y + 1/y)/2

def sin(x):
    """
    The sine function is not yet implemented in Calcium.
    This placeholder function evaluates the sine function
    using complex exponentials.

    Examples::

        >>> sin(0)
        0
        >>> sin(pi)
        0
        >>> sin(pi/2)
        1
        >>> sin(pi/6)
        0.500000 {1/2}
        >>> sin(sqrt(2))**2 + cos(sqrt(2))**2
        1
        >>> sin(3 + pi) + sin(3)
        0

    """
    ix = ca(x)*i
    y = exp(ix)
    return (y - 1/y)/(2*i)

def tan(x):
    """
    The tangent function is not yet implemented in Calcium.
    This placeholder function evaluates the tangent function
    using complex exponentials.

    """
    return sin(x)/cos(x)

def atan(x):
    """
    The inverse tangent function is not yet implemented in Calcium.
    This placeholder function evaluates the inverse tangent function
    using complex logarithms.

    Examples::

        >>> atan(0)
        0
        >>> 4 * atan(1) == pi
        True

    """
    return (-i/2)*log((i-x)/(i+x))

def cosh(x):
    """
    The hyperbolic cosine function is not yet implemented in Calcium.
    This placeholder function evaluates the hyperbolic cosine function
    using exponentials.

    Examples::

        >>> cosh(1)
        1.54308 {(a^2+1)/(2*a) where a = 2.71828 [Exp(1)]}

    """
    y = exp(x)
    return (y + 1/y)/2

def sinh(x):
    """
    The hyperbolic sine function is not yet implemented in Calcium.
    This placeholder function evaluates the hyperbolic sine function
    using exponentials.

    Examples::

        >>> sinh(1)
        1.17520 {(a^2-1)/(2*a) where a = 2.71828 [Exp(1)]}

    """
    y = exp(x)
    return (y - 1/y)/2

def tanh(x):
    """
    The hyperbolic tangent function is not yet implemented in Calcium.
    This placeholder function evaluates the hyperbolic tangent function
    using exponentials.

    Examples::

        >>> tanh(1)
        0.761594 {(a^2-1)/(a^2+1) where a = 2.71828 [Exp(1)]}

    """
    return sinh(x)/cosh(x)


#class allocated_c_char_p(ctypes.c_char_p):
#    def __del__(self):
#        libflint.flint_free(self)

libcalcium.ca_set_si.argtypes = ca, ctypes.c_long, ca_ctx
libcalcium.ca_set_d.argtypes = ca, ctypes.c_double, ca_ctx
libcalcium.ca_set_d_d.argtypes = ca, ctypes.c_double, ctypes.c_double, ca_ctx

libcalcium.ca_get_str.argtypes = ca, ca_ctx
libcalcium.ca_get_str.restype = ctypes.c_char_p

libflint.fmpz_set_str.argtypes = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int


i = j = I = ca.i()
pi = ca.pi()
euler = ca.euler()
e = E = ca(1).exp()

inf = ca.inf()
uinf = ca.uinf()
undefined = ca.undefined()
unknown = ca.unknown()


def test_floor_ceil():
    assert floor(sqrt(2)) == 1
    assert ceil(sqrt(2)) == 2

def test_power_identities():
    assert (sqrt(2)**sqrt(2))**sqrt(2) == 2
    assert (sqrt(-2)**sqrt(2))**sqrt(2) == -2
    assert (sqrt(3)**sqrt(3))**sqrt(3) == 3*sqrt(3)
    assert sqrt(-pi)**2 == -pi

def test_log():
    assert log(1+pi) - log(pi) - log(1+1/pi) == 0
    assert log(log(-log(log(exp(exp(-exp(exp(3)))))))) == 3

def test_exp():
    assert exp(pi*i) + 1 == 0
    assert exp(pi*i) == -1
    assert exp(log(2)*log(3)) > 2
    assert e**2 == exp(2)

def test_erf():
    assert erf(2*log(sqrt(ca(1)/2-sqrt(2)/4))+log(4)) - erf(log(2-sqrt(2))) == 0
    assert 1-erf(pi)-erfc(pi) == 0
    assert erf(sqrt(2))**2 + erfi(sqrt(-2))**2 == 0

def test_gudermannian():
    def gd(x):
        return 2*atan(exp(x))-pi/2
    assert sin(gd(1)) == tanh(1)
    assert tan(gd(1)) == sinh(1)
    assert sin(gd(sqrt(2))) == tanh(sqrt(2))
    assert tan(gd(1)/2) - tanh(ca(1)/2) == 0

def test_gamma():
    assert gamma(1) == 1
    assert gamma(0.5) == sqrt(pi)
    assert 1/gamma(0) == 0
    assert gamma(sqrt(2)*sqrt(3)) == gamma(sqrt(6))
    #assert gamma(pi+1)/gamma(pi) == pi
    assert gamma(pi)/gamma(pi-1) == pi-1


def test_xfail():
    # Test some simplifications that are known not to work yet.
    # When a case starts working, we will get a test failure so we can
    # catch it and add it to the working tests
    def gd(x):
        return 2*atan(exp(x))-pi/2

    def expect_error(f):
        try:
            f()
        except ValueError:
            return
        raise AssertionError

    # expect_error(lambda: tan(gd(1)/2) - tanh(ca(1)/2) == 0)

if __name__ == "__main__":
    from time import time
    print("Testing pycalcium_ctypes")
    print("----------------------------------------------------------")
    for fname in dir():
        if fname.startswith("test_"):
            print(fname + "...", end="")
            t1 = time()
            globals()[fname]()
            t2 = time()
            print("PASS", end="     ")
            print("%.2f" % (t2-t1))
    print("----------------------------------------------------------")
    import doctest
    doctest.testmod(verbose=True)
    print("----------------------------------------------------------")
