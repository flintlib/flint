import ctypes
import ctypes.util
import sys

libflint_path = ctypes.util.find_library('flint')
libflint = ctypes.CDLL(libflint_path)

libcalcium_path = ctypes.util.find_library('calcium')
libcalcium = ctypes.CDLL(libcalcium_path)

libarb_path = ctypes.util.find_library('arb')
libarb = ctypes.CDLL(libarb_path)

libgr_path = ctypes.util.find_library('genericrings')
libgr = ctypes.CDLL(libgr_path)

T_TRUE = 0
T_FALSE = 1
T_UNKNOWN = 2

GR_SUCCESS = 0
GR_DOMAIN = 1
GR_UNABLE = 2

c_slong = ctypes.c_long
c_ulong = ctypes.c_ulong

class flint_rand_struct(ctypes.Structure):
    # todo: use the real size
    _fields_ = [('data', c_slong * 16)]

_flint_rand = flint_rand_struct()
libflint.flint_randinit(ctypes.byref(_flint_rand))

class fmpz_struct(ctypes.Structure):
    _fields_ = [('val', c_slong)]

class fmpq_struct(ctypes.Structure):
    _fields_ = [('num', c_slong),
                ('den', c_slong)]

class fmpz_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', ctypes.c_long),
                ('length', ctypes.c_long)]

class fmpq_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', ctypes.c_long),
                ('length', ctypes.c_long),
                ('den', ctypes.c_long)]

class arb_struct(ctypes.Structure):
    _fields_ = [('data', ctypes.c_long * 6)]

class acb_struct(ctypes.Structure):
    _fields_ = [('real', arb_struct),
                ('imag', arb_struct)]

class fexpr_struct(ctypes.Structure):
    _fields_ = [('data', ctypes.c_void_p),
                ('alloc', ctypes.c_long)]

class qqbar_struct(ctypes.Structure):
    _fields_ = [('poly', fmpz_poly_struct),
                ('enclosure', acb_struct)]

class ca_struct(ctypes.Structure):
    _fields_ = [('data', ctypes.c_long * 5)]

class gr_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', ctypes.c_long),
                ('length', ctypes.c_long)]

class psl2z_struct(ctypes.Structure):
    _fields_ = [('a', ctypes.c_long), ('b', ctypes.c_long),
                ('c', ctypes.c_long), ('d', ctypes.c_long)]


# todo: efficiently
def fmpz_to_python_int(xref):
    ptr = libflint.fmpz_get_str(None, 10, xref)
    try:
        return int(ctypes.cast(ptr, ctypes.c_char_p).value.decode())
    finally:
        libflint.flint_free(ptr)

# todo
def fmpq_set_python(cref, x):
    assert isinstance(x, int) and -sys.maxsize <= x <= sys.maxsize
    libflint.fmpq_set_si(cref, x, 1)


class Undecidable(NotImplementedError):
    pass

class gr_ctx_struct(ctypes.Structure):
    # todo: use the real size
    _fields_ = [('content', ctypes.c_ulong * 64)]


libflint.flint_malloc.restype = ctypes.c_void_p
libflint.flint_free.argtypes = (ctypes.c_void_p,)
libflint.fmpz_set_str.argtypes = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int
libflint.fmpz_get_str.argtypes = ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(fmpz_struct)
libflint.fmpz_get_str.restype = ctypes.c_void_p

libgr.gr_heap_init.argtypes = (ctypes.POINTER(gr_ctx_struct),)
libgr.gr_heap_init.restype = ctypes.c_void_p

libgr.gr_set_si.argtypes = (ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(gr_ctx_struct))
libgr.gr_add_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(gr_ctx_struct))
libgr.gr_sub_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(gr_ctx_struct))
libgr.gr_mul_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(gr_ctx_struct))
libgr.gr_div_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(gr_ctx_struct))
libgr.gr_pow_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.POINTER(gr_ctx_struct))


libgr.gr_set_str.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmp.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmpabs.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_heap_clear.argtypes = (ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

_add_methods = [libgr.gr_add, libgr.gr_add_ui, libgr.gr_add_si, libgr.gr_add_fmpz, libgr.gr_add_fmpq]
_sub_methods = [libgr.gr_sub, libgr.gr_sub_ui, libgr.gr_sub_si, libgr.gr_sub_fmpz, libgr.gr_sub_fmpq]
_mul_methods = [libgr.gr_mul, libgr.gr_mul_ui, libgr.gr_mul_si, libgr.gr_mul_fmpz, libgr.gr_mul_fmpq]
_div_methods = [libgr.gr_div, libgr.gr_div_ui, libgr.gr_div_si, libgr.gr_div_fmpz, libgr.gr_div_fmpq]
_pow_methods = [libgr.gr_pow, libgr.gr_pow_ui, libgr.gr_pow_si, libgr.gr_pow_fmpz, libgr.gr_pow_fmpq]

_gr_logic = 0

class LogicContext(object):
    """
    Handle the result of predicates (experimental):

        >>> a = (RR_arb(1) / 3) * 3
        >>> a
        [1.00000000000000 +/- 3.89e-16]
        >>> with strict_logic:
        ...     a == 1
        ...
        Traceback (most recent call last):
          ...
        Undecidable: x == y cannot be decided for x = [1.00000000000000 +/- 3.89e-16], y = 1.000000000000000 over Real numbers (arb, prec = 53)
        >>> with pessimistic_logic:
        ...     a == 1
        ...
        False
        >>> with optimistic_logic:
        ...     a == 1
        ...
        True
    """

    def __init__(self, value):
        self.logic = value
    def __enter__(self):
        global _gr_logic
        self.original = _gr_logic
        _gr_logic = self.logic
    def __exit__(self, type, value, traceback):
        global _gr_logic
        _gr_logic = self.original

strict_logic = LogicContext(0)
pessimistic_logic = LogicContext(-1)
optimistic_logic = LogicContext(1)


class gr_ctx:

    def __init__(self):
        self._data = gr_ctx_struct()
        self._ref = ctypes.byref(self._data)
        self._str = None

    def _repr(self):
        if self._str is None:
            arr = ctypes.c_char_p()
            if libgr.gr_ctx_get_str(ctypes.byref(arr), self._ref) != GR_SUCCESS:
                raise NotImplementedError
            try:
                self._str = ctypes.cast(arr, ctypes.c_char_p).value.decode("ascii")
            finally:
                libflint.flint_free(arr)
        return self._str

    def __call__(self, value=None, random=False):
        return self._elem_type(value, self, random)

    def __repr__(self):
        return self._repr()

    def __del__(self):
        # todo: refcounting
        # libgr.gr_ctx_clear(self._ref)
        pass

class gr_elem:

    def __init__(self, val=None, context=None, random=False):
        assert context is not None
        if context is None:
            raise ValueError("a context object is needed")
        else:
            self._ctx_python = context
        self._ctx = self._ctx_python._ref
        self._data = self._struct_type()
        self._ref = ctypes.byref(self._data)
        libgr.gr_init(self._ref, self._ctx)
        if val is not None:
            typ = type(val)
            if typ is int:
                b = sys.maxsize
                if -b <= val <= b:
                    status = libgr.gr_set_si(self._ref, val, self._ctx)
                    if status:
                        if status & GR_UNABLE: raise NotImplementedError()
                        if status & GR_DOMAIN: raise ValueError()
                else:
                    n = fmpz_struct()
                    nref = ctypes.byref(n)
                    libflint.fmpz_init(nref)
                    libflint.fmpz_set_str(nref, ctypes.c_char_p(str(val).encode('ascii')), 10)
                    status = libgr.gr_set_fmpz(self._ref, nref, self._ctx)
                    if status:
                        if status & GR_UNABLE: raise NotImplementedError()
                        if status & GR_DOMAIN: raise ValueError()
                    libflint.fmpz_clear(nref)
            elif isinstance(val, gr_elem):
                status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
                if status:
                    if status & GR_UNABLE: raise NotImplementedError()
                    if status & GR_DOMAIN: raise ValueError()
            elif typ is str:
                status = libgr.gr_set_str(self._ref, ctypes.c_char_p(str(val).encode('ascii')), self._ctx)
                if status:
                    if status & GR_UNABLE: raise NotImplementedError()
                    if status & GR_DOMAIN: raise ValueError()
            else:
                raise NotImplementedError(f"unable to create {type(self)} from {type(val)}")
        elif random:
            libgr.gr_randtest(self._ref, ctypes.byref(_flint_rand), self._ctx)

    def __del__(self):
        libgr.gr_clear(self._ref, self._ctx)

    def parent(self):
        return self._ctx_python

    def __repr__(self):
        arr = ctypes.c_char_p()
        if libgr.gr_get_str(ctypes.byref(arr), self._ref, self._ctx) != GR_SUCCESS:
            raise NotImplementedError
        try:
            return ctypes.cast(arr, ctypes.c_char_p).value.decode("ascii")
        finally:
            libflint.flint_free(arr)

    @staticmethod
    def _binary_coercion(self, other):
        elem_type = type(self)
        other_type = type(other)
        if elem_type is not other_type:
            if not isinstance(other, gr_elem):
                other = self.parent()(other)
            elif not isinstance(self, gr_elem):
                self = other.parent()(self)
        if self._ctx_python is not other._ctx_python:
            c = libgr.gr_ctx_cmp_coercion(self._ctx, other._ctx)
            if c >= 0:
                other = self.parent()(other)
            else:
                self = other.parent()(self)
        return self, other

    @staticmethod
    def _binary_op(self, other, op, rstr):
        self, other = gr_elem._binary_coercion(self, other)
        res = type(self)(None, self._ctx_python)
        status = op(res._ref, self._ref, other._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"{rstr} is not implemented for x = {self}, y = {other} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self}, y = {other} over {self.parent()}")
        return res

    @staticmethod
    def _binary_op2(self, other, ops, rstr):
        if type(other) is int and -sys.maxsize <= other <= sys.maxsize:
            res = type(self)(None, self._ctx_python)
            status = ops[2](res._ref, self._ref, other, self._ctx)
        else:
            self, other = gr_elem._binary_coercion(self, other)
            res = type(self)(None, self._ctx_python)
            status = ops[0](res._ref, self._ref, other._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"{rstr} is not implemented for x = {self}, y = {other} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self}, y = {other} over {self.parent()}")
        return res

    @staticmethod
    def _binary_predicate(self, other, op, rstr):
        self, other = gr_elem._binary_coercion(self, other)
        truth = op(self._ref, other._ref, self._ctx)
        if truth == T_TRUE: return True
        if truth == T_FALSE: return False
        if _gr_logic == 1: return True
        if _gr_logic == -1: return False
        raise Undecidable(f"{rstr} cannot be decided for x = {self}, y = {other} over {self.parent()}")

    @staticmethod
    def _unary_op(self, op, rstr):
        elem_type = type(self)
        res = elem_type(None, self._ctx_python)
        status = op(res._ref, self._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"{rstr} cannot be computed for x = {self} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self} over {self.parent()}")
        return res

    def __eq__(self, other):
        return self._binary_predicate(self, other, libgr.gr_equal, "x == y")

    def __ne__(self, other):
        return self._binary_predicate(self, other, libgr.gr_not_equal, "x != y")

    def _cmp(self, other):
        self, other = gr_elem._binary_coercion(self, other)
        c = (ctypes.c_int * 1)()
        status = libgr.gr_cmp(c, self._ref, other._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise Undecidable(f"unable to compare x = {self} and y = {other} in {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"ordering not defined for x = {self} and y = {other} in {self.parent()}")
        return c[0]

    def __lt__(self, other):
        return gr_elem._cmp(self, other) < 0

    def __le__(self, other):
        return gr_elem._cmp(self, other) <= 0

    def __gt__(self, other):
        return gr_elem._cmp(self, other) > 0

    def __ge__(self, other):
        return gr_elem._cmp(self, other) >= 0

    def __neg__(self):
        return self._unary_op(self, libgr.gr_neg, "-x")

    def __pos__(self):
        return self

    def __abs__(self):
        return self._unary_op(self, libgr.gr_abs, "abs(x)")

    def __add__(self, other):
        return self._binary_op(self, other, libgr.gr_add, "x + y")

    def __radd__(self, other):
        return self._binary_op(other, self, libgr.gr_add, "x + y")

    def __sub__(self, other):
        return self._binary_op(self, other, libgr.gr_sub, "x - y")

    def __rsub__(self, other):
        return self._binary_op(other, self, libgr.gr_sub, "x - y")

    def __mul__(self, other):
        return self._binary_op(self, other, libgr.gr_mul, "x * y")

    def __rmul__(self, other):
        return self._binary_op(other, self, libgr.gr_mul, "x * y")

    def __truediv__(self, other):
        return self._binary_op(self, other, libgr.gr_div, "x / y")

    def __rtruediv__(self, other):
        return self._binary_op(other, self, libgr.gr_div, "x / y")

    def __pow__(self, other):
        return self._binary_op2(self, other, _pow_methods, "x ** y")

    def __rpow__(self, other):
        return self._binary_op2(other, self, _pow_methods, "x ** y")

    def sqrt(self):
        return self._unary_op(self, libgr.gr_sqrt, "sqrt(x)")

    def abs(self):
        return self._unary_op(self, libgr.gr_abs, "abs(x)")

    def conj(self):
        return self._unary_op(self, libgr.gr_conj, "conj(x)")

    def re(self):
        return self._unary_op(self, libgr.gr_re, "re(x)")

    def im(self):
        return self._unary_op(self, libgr.gr_im, "im(x)")

    def sgn(self):
        return self._unary_op(self, libgr.gr_sgn, "sgn(x)")

    def csgn(self):
        return self._unary_op(self, libgr.gr_csgn, "csgn(x)")

    def exp(self):
        return self._unary_op(self, libgr.gr_exp, "exp(x)")

    def log(self):
        return self._unary_op(self, libgr.gr_log, "log(x)")

    def sin(self):
        return self._unary_op(self, libgr.gr_sin, "sin(x)")

    def cos(self):
        return self._unary_op(self, libgr.gr_cos, "cos(x)")

    def atan(self):
        return self._unary_op(self, libgr.gr_atan, "atan(x)")


class IntegerRing_fmpz(gr_ctx):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fmpz(self._ref)
        self._elem_type = fmpz

class RationalField_fmpq(gr_ctx):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fmpq(self._ref)
        self._elem_type = fmpq

class ComplexAlgebraicField_qqbar(gr_ctx):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_qqbar(self._ref)
        self._elem_type = qqbar

class RealAlgebraicField_qqbar(gr_ctx):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_real_qqbar(self._ref)
        self._elem_type = qqbar

class gr_arb_ctx(gr_ctx):

    @property
    def prec(self):
        return libgr.gr_ctx_arb_get_prec(self._ref)

    @prec.setter
    def prec(self, prec):
        libgr.gr_ctx_arb_set_prec(self._ref, prec)


class RealField_arb(gr_arb_ctx):
    def __init__(self, prec=53):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_real_arb(self._ref, prec)
        self._elem_type = arb

class ComplexField_acb(gr_arb_ctx):
    def __init__(self, prec=53):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_acb(self._ref, prec)
        self._elem_type = acb

_ca_options = [
    "verbose",
    "print_flags",
    "mpoly_ord",
    "prec_limit",
    "qqbar_deg_limit",
    "low_prec",
    "smooth_limit",
    "lll_prec",
    "pow_limit",
    "use_gb",
    "gb_length_limit",
    "gb_poly_length_limit",
    "gb_poly_bits_limit",
    "vieta_limit",
    "trig_form"]

class gr_ctx_ca(gr_ctx):

    def _set_options(self, kwargs):
        for w in kwargs:
            i = _ca_options.index(w)
            if i == -1:
                raise ValueError(f"unknown option {w}")
            libgr.gr_ctx_ca_set_option(self._ref, i, kwargs[w])

    def options(self):
        opts = {_ca_options[i] : libgr.gr_ctx_ca_get_option(self._ref, i) for i in range(len(_ca_options))}
        return opts

class RealAlgebraicField_ca(gr_ctx_ca):
    def __init__(self, **kwargs):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_real_algebraic_ca(self._ref)
        self._elem_type = ca
        self._set_options(kwargs)

class ComplexAlgebraicField_ca(gr_ctx_ca):
    def __init__(self, **kwargs):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_algebraic_ca(self._ref)
        self._elem_type = ca
        self._set_options(kwargs)

class RealField_ca(gr_ctx_ca):
    def __init__(self, **kwargs):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_real_ca(self._ref)
        self._elem_type = ca
        self._set_options(kwargs)

class ComplexField_ca(gr_ctx_ca):
    def __init__(self, **kwargs):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_ca(self._ref)
        self._elem_type = ca
        self._set_options(kwargs)


class PolynomialRing_gr_poly(gr_ctx):
    def __init__(self, coefficient_ring):
        assert isinstance(coefficient_ring, gr_ctx)
        if libgr.gr_ctx_is_ring(coefficient_ring._ref) != T_TRUE:
            raise ValueError("coefficient structure must be a ring")
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_polynomial(self._ref, coefficient_ring._ref)
        self._coefficient_ring = coefficient_ring
        self._elem_type = gr_poly


class fmpz(gr_elem):
    _struct_type = fmpz_struct

    def __int__(self):
        return fmpz_to_python_int(self._ref)


class fmpq(gr_elem):
    _struct_type = fmpq_struct

class qqbar(gr_elem):
    _struct_type = qqbar_struct

class ca(gr_elem):
    _struct_type = ca_struct

class arb(gr_elem):
    _struct_type = arb_struct

class acb(gr_elem):
    _struct_type = acb_struct

class gr_poly(gr_elem):
    _struct_type = gr_poly_struct

    def __init__(self, val=None, context=None, random=False):
        # todo: also iterables
        if isinstance(val, (list, tuple)):
            gr_elem.__init__(self, None, context)
            coefficient_ring = self.parent()._coefficient_ring
            val = [coefficient_ring(c) for c in val]
            for i in range(len(val)):
                status = libgr.gr_poly_set_coeff_scalar(self._ref, i, val[i]._ref, coefficient_ring._ref)
                if status:
                    raise NotImplementedError
        else:
            gr_elem.__init__(self, val, context)
            # todo: refactor
            if random:
                libgr.gr_randtest(self._ref, ctypes.byref(_flint_rand), self._ctx)

    def __len__(self):
        return self._data.length

    def __getitem__(self, i):
        n = len(self)
        R = self.parent()._coefficient_ring
        c = R()
        status = libgr.gr_poly_get_coeff_scalar(c._ref, self._ref, i, R._ref)
        if status:
            raise NotImplementedError
        return c


class ModularGroup_psl2z(gr_ctx_ca):
    def __init__(self, **kwargs):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_psl2z(self._ref)
        self._elem_type = psl2z

    # todo: C function
    def generators(self):
        S = self()
        T = self()
        S._data.a = 0
        S._data.b = -1
        S._data.c = 1
        S._data.d = 0
        T._data.b = 1
        return (S, T)

class psl2z(gr_elem):
    _struct_type = psl2z_struct



ZZ = IntegerRing_fmpz()
QQ = RationalField_fmpq()
AA = RealAlgebraicField_qqbar()
QQbar = ComplexAlgebraicField_qqbar()
RR = RR_arb = RealField_arb()
CC = CC_acb = ComplexField_acb()
RR_ca = RealField_ca()
CC_ca = ComplexField_ca()

ZZx = PolynomialRing_gr_poly(ZZ)
QQx = PolynomialRing_gr_poly(QQ)
RRx_ca = PolynomialRing_gr_poly(RR_ca)
CCx_ca = PolynomialRing_gr_poly(CC_ca)
RRx = RRx_arb = PolynomialRing_gr_poly(RR_arb)
CCx = CCx_acb = PolynomialRing_gr_poly(CC_acb)

PSL2Z = ModularGroup_psl2z()


def test_all():

    def raises(f, exception):
        try:
            f()
        except exception:
            return True
        raise False

    x = ZZ(23)
    y = ZZ(-1)
    assert str(x) == "23"
    assert x.parent() is ZZ
    assert int(x) == 23
    assert x + y == ZZ(22)
    assert x - y == ZZ(24)
    assert x * y == ZZ(-23)
    assert -x == ZZ(-23)

    assert ZZ(3) != 4
    assert ZZ(3) <= 5
    assert ZZ(3) > 2

    x = QQ(-10000000000000000000075) / QQ(3)
    assert str(x) == "-10000000000000000000075/3"
    assert x.parent() is QQ

    x = QQbar(-2)
    y = QQbar(1) / QQbar(3)
    assert x.parent() is QQbar
    xy = x ** y
    assert (xy ** QQbar(3)) == QQbar(-2)
    assert str(xy) == "Root a = 0.629961 + 1.09112*I of a^3+2"
    i = QQbar(-1) ** (QQ(1)/2)
    assert str(i) == 'Root a = 1.00000*I of a^2+1'
    assert str(-i) == 'Root a = -1.00000*I of a^2+1'
    assert str(1-i) == 'Root a = 1.00000 - 1.00000*I of a^2-2*a+2'
    assert raises(lambda: i > 0, ValueError)
    assert QQ(-3)/2 < i**2 < QQ(1)/2

    assert abs(QQ(-5)) == QQ(5)
    assert QQ(8) ** (QQ(1) / QQ(3)) == QQ(2)
    assert raises(lambda: QQ(2) ** (QQ(1) / QQ(3)), ValueError)

    assert QQ(1) + 2 == QQ(3)
    assert 2 + QQ(1) == QQ(3)
    assert QQ(1) + ZZ(5) == QQ(6)
    assert (QQ(1) + ZZ(5)).parent() is QQ
    assert raises(lambda: ZZ(1) / 2, ValueError)
    assert raises(lambda: (-1) ** (QQ(1) / 2), ValueError)
    assert ((-1) ** (QQbar(1) / 2)) ** 2 == QQbar(-1)

    f = ZZx([1,2,3]) + QQx([1,2])
    assert f == ZZx([2,4,3])
    assert f.parent() is QQx
    assert RRx([1,QQ(2),AA(3)]) != ZZx([1,2,3,4])
    assert RRx([1,QQ(2),AA(3),4]) == ZZx([1,2,3,4])
    assert ZZx(3) + ZZx(2) == ZZx([5])
    assert ZZx(3) + 2 == ZZx([5])

if __name__ == "__main__":
    from time import time
    print("Testing flint_ctypes")
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
    doctest.testmod(optionflags=doctest.FAIL_FAST, verbose=True)
    print("----------------------------------------------------------")

