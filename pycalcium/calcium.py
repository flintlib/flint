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
    Python class wrapping the ca_t type as a number object.

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

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

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
        res = ca()
        libcalcium.ca_re(res, self, self._ctx)
        return res

    def im(self):
        res = ca()
        libcalcium.ca_im(res, self, self._ctx)
        return res

    def conj(self):
        res = ca()
        libcalcium.ca_conjugate(res, self, self._ctx)
        return res

    conjugate = conj

    def floor(self):
        res = ca()
        libcalcium.ca_floor(res, self, self._ctx)
        return res

    def ceil(self):
        res = ca()
        libcalcium.ca_ceil(res, self, self._ctx)
        return res

    def sgn(self):
        res = ca()
        libcalcium.ca_sgn(res, self, self._ctx)
        return res

    sign = sgn

    def sqrt(self):
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
        res = ca()
        libcalcium.ca_exp(res, self, self._ctx)
        return res

    def erf(self):
        res = ca()
        libcalcium.ca_erf(res, self, self._ctx)
        return res

    def erfc(self):
        res = ca()
        libcalcium.ca_erfc(res, self, self._ctx)
        return res

    def erfi(self):
        res = ca()
        libcalcium.ca_erfi(res, self, self._ctx)
        return res

    def gamma(self):
        res = ca()
        libcalcium.ca_gamma(res, self, self._ctx)
        return res


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
    return (ca(x)+1).gamma()

def cos(x):
    ix = ca(x)*i
    y = exp(ix)
    return (y + 1/y)/2

def sin(x):
    ix = ca(x)*i
    y = exp(ix)
    return (y - 1/y)/(2*i)

def tan(x):
    return sin(x)/cos(x)

def atan(x):
    return (-i/2)*log((i-x)/(i+x))

def cosh(x):
    y = exp(x)
    return (y + 1/y)/2

def sinh(x):
    y = exp(x)
    return (y - 1/y)/2

def tanh(x):
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
