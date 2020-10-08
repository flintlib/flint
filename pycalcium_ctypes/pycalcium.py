"""
Simple ctypes wrapper around Calcium mainly intended for test code.

Currently leaks memory, doesn't support multiple context objects...

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
    _fields_ = [('head', ctypes.c_ulong),
                ('data0', ctypes.c_ulong),
                ('data1', ctypes.c_ulong),
                ('data2', ctypes.c_ulong),
                ('data3', ctypes.c_ulong)]


class ca_ctx_struct(ctypes.Structure):
    # todo: use the real size
    _fields_ = [('content', ctypes.c_ulong * 64)]


class ca_ctx:

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
    def pi():
        res = ca()
        libcalcium.ca_pi(res, res._ctx)
        return res

    @staticmethod
    def euler():
        res = ca()
        libcalcium.ca_euler(res, res._ctx)
        return res

    @staticmethod
    def i():
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
