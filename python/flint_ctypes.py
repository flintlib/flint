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

if sys.maxsize < 2**32:
    FLINT_BITS = 32
else:
    FLINT_BITS = 64
    if ctypes.sizeof(c_slong) == 4:
        c_slong = ctypes.c_longlong
        c_ulong = ctypes.c_ulonglong
        assert ctypes.sizeof(c_slong) == 8
        assert ctypes.sizeof(c_ulong) == 8

UWORD_MAX = (1<<FLINT_BITS)-1
WORD_MAX = (1<<(FLINT_BITS-1))-1
WORD_MIN = -(1<<(FLINT_BITS-1))

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

class fmpzi_struct(ctypes.Structure):
    _fields_ = [('real', c_slong),
                ('imag', c_slong)]

class fmpz_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong)]

class fmpq_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong),
                ('den', c_slong)]

class arf_struct(ctypes.Structure):
    _fields_ = [('data', c_slong * 4)]

class acf_struct(ctypes.Structure):
    _fields_ = [('real', arf_struct),
                ('imag', arf_struct)]

class arb_struct(ctypes.Structure):
    _fields_ = [('data', c_slong * 6)]

class acb_struct(ctypes.Structure):
    _fields_ = [('real', arb_struct),
                ('imag', arb_struct)]

class fexpr_struct(ctypes.Structure):
    _fields_ = [('data', ctypes.c_void_p),
                ('alloc', c_slong)]

class qqbar_struct(ctypes.Structure):
    _fields_ = [('poly', fmpz_poly_struct),
                ('enclosure', acb_struct)]

class ca_struct(ctypes.Structure):
    _fields_ = [('data', c_slong * 5)]

class nmod_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong),
                ('n', c_ulong),
                ('ninv', c_ulong),
                ('nnorm', c_slong)]

class fq_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong)]

class fq_nmod_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong),
                ('n', c_ulong),
                ('ninv', c_ulong),
                ('nnorm', c_slong)]

class fq_zech_struct(ctypes.Structure):
    _fields_ = [('n', ctypes.c_ulong)]

class gr_vec_struct(ctypes.Structure):
    _fields_ = [('entries', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong)]

class gr_poly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong)]

class psl2z_struct(ctypes.Structure):
    _fields_ = [('a', c_slong), ('b', c_slong),
                ('c', c_slong), ('d', c_slong)]

class dirichlet_char_struct(ctypes.Structure):
    _fields_ = [('n', c_ulong),
                ('log', ctypes.POINTER(c_ulong))]

class perm_struct(ctypes.Structure):
    _fields_ = [('entries', ctypes.POINTER(c_slong))]

class gr_mat_struct(ctypes.Structure):
    _fields_ = [('entries', ctypes.c_void_p),
                ('r', c_slong),
                ('c', c_slong),
                ('rows', ctypes.c_void_p)]


# todo: efficiently
def fmpz_to_python_int(xref):
    ptr = libflint.fmpz_get_str(None, 10, xref)
    try:
        return int(ctypes.cast(ptr, ctypes.c_char_p).value.decode())
    finally:
        libflint.flint_free(ptr)

# todo
def fmpq_set_python(cref, x):
    assert isinstance(x, int) and WORD_MIN <= x <= WORD_MAX
    libflint.fmpq_set_si(cref, x, 1)


class Undecidable(NotImplementedError):
    pass

class gr_ctx_struct(ctypes.Structure):
    _fields_ = [('content', ctypes.c_char * libgr.gr_ctx_sizeof_ctx())]


libflint.flint_malloc.restype = ctypes.c_void_p
libflint.flint_free.argtypes = (ctypes.c_void_p,)
libflint.fmpz_set_str.argtypes = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int
libflint.fmpz_get_str.argtypes = ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(fmpz_struct)
libflint.fmpz_get_str.restype = ctypes.c_void_p

libgr.gr_heap_init.argtypes = (ctypes.POINTER(gr_ctx_struct),)
libgr.gr_heap_init.restype = ctypes.c_void_p

libgr.gr_set_si.argtypes = (ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_add_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_sub_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_mul_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_div_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_pow_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))

libgr.gr_set_d.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(gr_ctx_struct))

libgr.gr_set_str.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmp.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmpabs.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))


libgr.gr_heap_clear.argtypes = (ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_ctx_init_dirichlet_group.argtypes = (ctypes.POINTER(gr_ctx_struct), c_ulong)

_add_methods = [libgr.gr_add, libgr.gr_add_si, libgr.gr_add_fmpz, libgr.gr_add_other, libgr.gr_other_add]
_sub_methods = [libgr.gr_sub, libgr.gr_sub_si, libgr.gr_sub_fmpz, libgr.gr_sub_other, libgr.gr_other_sub]
_mul_methods = [libgr.gr_mul, libgr.gr_mul_si, libgr.gr_mul_fmpz, libgr.gr_mul_other, libgr.gr_other_mul]
_div_methods = [libgr.gr_div, libgr.gr_div_si, libgr.gr_div_fmpz, libgr.gr_div_other, libgr.gr_other_div]
_pow_methods = [libgr.gr_pow, libgr.gr_pow_si, libgr.gr_pow_fmpz, libgr.gr_pow_other, libgr.gr_other_pow]

_gr_logic = 0

class Truth:

    def __init__(self, value):
        self.value = value

    def __repr__(self):
        if self.value == T_TRUE:
            return "TRUE"
        if self.value == T_FALSE:
            return "FALSE"
        return "UNKNOWN"

    def __bool__(self):
        if self.value == T_UNKNOWN:
            raise ValueError("unknown truth value")
        return self.value == T_TRUE



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
        Undecidable: unable to decide x == y for x = [1.00000000000000 +/- 3.89e-16], y = 1.000000000000000 over Real numbers (arb, prec = 53)
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
none_logic = LogicContext(2)
triple_logic = LogicContext(3)

def set_logic(which_logic):
    global _gr_logic
    _gr_logic = which_logic.logic


class gr_ctx:

    def __init__(self):
        self._data = gr_ctx_struct()
        self._ref = ctypes.byref(self._data)
        self._str = None
        self._refcount = 1

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

    def __call__(self, *args, **kwargs):
        kwargs['context'] = self
        return self._elem_type(*args, **kwargs)

    def __repr__(self):
        return self._repr()

    def __del__(self):
        self._decrement_refcount()

    def _decrement_refcount(self):
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)

    # todo
    def i(self):
        return self().i()

    def pi(self):
        return self().pi()


class gr_elem:
    """
    Base class for elements.
    """

    @staticmethod
    def _default_context():
        return None

    def __init__(self, val=None, context=None, random=False):
        if context is None:
            context = self._default_context()
            if context is None:
                raise ValueError("a context object is needed")
        self._ctx_python = context
        self._ctx = self._ctx_python._ref
        self._data = self._struct_type()
        self._ref = ctypes.byref(self._data)
        libgr.gr_init(self._ref, self._ctx)
        self._ctx_python._refcount += 1
        if val is not None:
            typ = type(val)
            status = GR_UNABLE
            if typ is int:
                if WORD_MIN <= val <= WORD_MAX:
                    status = libgr.gr_set_si(self._ref, val, self._ctx)
                else:
                    n = fmpz_struct()
                    nref = ctypes.byref(n)
                    libflint.fmpz_init(nref)
                    libflint.fmpz_set_str(nref, ctypes.c_char_p(str(val).encode('ascii')), 10)
                    status = libgr.gr_set_fmpz(self._ref, nref, self._ctx)
                    libflint.fmpz_clear(nref)
            elif isinstance(val, gr_elem):
                status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
            elif typ is str:
                status = libgr.gr_set_str(self._ref, ctypes.c_char_p(str(val).encode('ascii')), self._ctx)
            elif typ is float:
                status = libgr.gr_set_d(self._ref, val, self._ctx)
            elif typ is complex:
                # todo
                x = context(val.real) + context(val.imag) * context.i()
                status = libgr.gr_set(self._ref, x._ref, self._ctx)
            elif hasattr(val, "_gr_elem_"):
                val = val._gr_elem_(context)
                assert val.parent() is context
                status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
            else:
                status = GR_UNABLE
            if status:
                if status & GR_UNABLE: raise NotImplementedError(f"unable to create element of {self.parent()} from {val} of type {type(val)}")
                if status & GR_DOMAIN: raise ValueError(f"{val} is not defined in {self.parent()}")
        elif random:
            libgr.gr_randtest(self._ref, ctypes.byref(_flint_rand), self._ctx)

    def __del__(self):
        libgr.gr_clear(self._ref, self._ctx)
        self._ctx_python._decrement_refcount()

    def parent(self):
        """
        Return the parent object of this element.

            >>> ZZ(0).parent()
            Integer ring (fmpz)
            >>> ZZ(0).parent() is ZZ
            True
        """
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
        res = type(self)(context=self._ctx_python)
        status = op(res._ref, self._ref, other._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"unable to compute {rstr} for x = {self}, y = {other} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self}, y = {other} over {self.parent()}")
        return res

    @staticmethod
    def _binary_op2(self, other, ops, rstr):
        self_type = type(self)
        other_type = type(other)
        if self_type is other_type and self._ctx_python is other._ctx_python:
            res = type(self)(context=self._ctx_python)
            status = ops[0](res._ref, self._ref, other._ref, self._ctx)
        elif isinstance(self, gr_elem) and isinstance(other, gr_elem):
            c = libgr.gr_ctx_cmp_coercion(self._ctx, other._ctx)
            if c >= 0:
                # other -> self
                # print("trying", other, "into", self)
                res = type(self)(context=self._ctx_python)
                status = ops[3](res._ref, self._ref, other._ref, other._ctx, self._ctx)
            else:
                # self -> other
                # print("trying", self, "into", other)
                res = type(other)(context=other._ctx_python)
                status = ops[4](res._ref, self._ref, self._ctx, other._ref, other._ctx)
            # needed?
            if status:
                if c >= 0:
                    other = self.parent()(other)
                else:
                    self = other.parent()(self)
                res = type(self)(context=self._ctx_python)
                status = ops[0](res._ref, self._ref, other._ref, self._ctx)
        elif other_type is int:
            if WORD_MIN <= other <= WORD_MAX:   # todo: efficient code from left also
                res = type(self)(context=self._ctx_python)
                status = ops[1](res._ref, self._ref, other, self._ctx)
            else:
                other = ZZ(other)
                res = type(self)(context=self._ctx_python)
                status = ops[2](res._ref, self._ref, other._ref, self._ctx)
        elif self_type is int:
            return other._binary_op2(ZZ(self), other, ops, rstr)
        else:
            if not isinstance(other, gr_elem):
                other = self.parent()(other)
            elif not isinstance(self, gr_elem):
                self = other.parent()(self)
            return self._binary_op2(self, other, ops, rstr)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"unable to compute {rstr} for x = {self}, y = {other} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self}, y = {other} over {self.parent()}")
        return res

    @staticmethod
    def _unary_predicate(self, op, rstr):
        truth = op(self._ref, self._ctx)
        if _gr_logic == 3:
            return Truth(truth)
        if truth == T_TRUE: return True
        if truth == T_FALSE: return False
        if _gr_logic == 1: return True
        if _gr_logic == -1: return False
        if _gr_logic == 2: return None
        raise Undecidable(f"unable to decide {rstr} for x = {self} over {self.parent()}")

    @staticmethod
    def _binary_predicate(self, other, op, rstr):
        self, other = gr_elem._binary_coercion(self, other)
        truth = op(self._ref, other._ref, self._ctx)
        if _gr_logic == 3:
            return Truth(truth)
        if truth == T_TRUE: return True
        if truth == T_FALSE: return False
        if _gr_logic == 1: return True
        if _gr_logic == -1: return False
        if _gr_logic == 2: return None
        raise Undecidable(f"unable to decide {rstr} for x = {self}, y = {other} over {self.parent()}")

    @staticmethod
    def _unary_op(self, op, rstr):
        elem_type = type(self)
        res = elem_type(context=self._ctx_python)
        status = op(res._ref, self._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"unable to compute {rstr} for x = {self} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self} over {self.parent()}")
        return res

    @staticmethod
    def _unary_op_get_fmpz(self, op, rstr):
        res = ZZ()
        status = op(res._ref, self._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"unable to compute {rstr} for x = {self} over {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self} over {self.parent()}")
        return res

    @staticmethod
    def _constant(self, op, rstr):
        elem_type = type(self)
        res = elem_type(context=self._ctx_python)
        status = op(res._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"unable to compute {rstr} in {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined in {self.parent()}")
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
        return self._binary_op2(self, other, _add_methods, "x + y")

    def __radd__(self, other):
        return self._binary_op2(other, self, _add_methods, "x + y")

    def __sub__(self, other):
        return self._binary_op2(self, other, _sub_methods, "x - y")

    def __rsub__(self, other):
        return self._binary_op2(other, self, _sub_methods, "x - y")

    def __mul__(self, other):
        return self._binary_op2(self, other, _mul_methods, "x * y")

    def __rmul__(self, other):
        return self._binary_op2(other, self, _mul_methods, "x * y")

    def __truediv__(self, other):
        return self._binary_op2(self, other, _div_methods, "x / y")

    def __rtruediv__(self, other):
        return self._binary_op2(other, self, _div_methods, "x / y")

    def __pow__(self, other):
        return self._binary_op2(self, other, _pow_methods, "x ** y")

    def __rpow__(self, other):
        return self._binary_op2(other, self, _pow_methods, "x ** y")

    def __floordiv__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_div, "x // y")

    def __rfloordiv__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_div, "x // y")

    def __mod__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_rem, "x % y")

    def __rmod__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_rem, "x % y")

    def is_invertible(self):
        """
        Return whether self has a multiplicative inverse in its domain.

            >>> 
            >>> ZZ(3).is_invertible()
            False
            >>> ZZ(-1).is_invertible()
            True
        """
        return self._unary_predicate(self, libgr.gr_is_invertible, "is_invertible")

    def divides(self, other):
        """
        Return whether self divides other.

            >>> ZZ(5).divides(10)
            True
            >>> ZZ(5).divides(12)
            False
        """
        return self._binary_predicate(self, other, libgr.gr_divides, "divides")

    def gcd(self, other):
        """
        Greatest common divisor.

            >>> ZZ(24).gcd(30)
            6
        """
        return self._binary_op(self, other, libgr.gr_gcd, "gcd")

    def lcm(self, other):
        """
        Least common multiple.

            >>> ZZ(24).lcm(30)
            120
        """
        return self._binary_op(self, other, libgr.gr_lcm, "lcm")

    def factor(self):
        """
        Returns a factorization of self as a tuple (prefactor, factors, exponents).

            >>> ZZ(-120).factor()
            (-1, [2, 3, 5], [3, 1, 1])

        """
        elem_type = type(self)
        c = elem_type(context=self._ctx_python)
        factors = Vec(self._ctx_python)()
        exponents = ZZ_vec()
        # print("c", c)
        # print("factors", factors)
        # print("c", exponents)
        status = libgr.gr_factor(c._ref, factors._ref, exponents._ref, self._ref, 0, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (c, factors, exponents)

    def is_square(self):
        """
        Return whether self is a perfect square in its domain.

            >>> ZZ(3).is_square()
            False
            >>> ZZ(4).is_square()
            True
        """
        return self._unary_predicate(self, libgr.gr_is_square, "is_square")

    def __index__(self):
        n = fmpz_struct()
        nref = ctypes.byref(n)
        libflint.fmpz_init(nref)
        status = libgr.gr_get_fmpz(nref, self._ref, self._ctx)
        v = fmpz_to_python_int(nref)
        libflint.fmpz_clear(nref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"unable to convert x = {self} to integer in {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"x = {self} is not an integer in {self.parent()}")
        return v

    def __int__(self):
        return self.trunc().__index__()

    def __float__(self):
        c = (ctypes.c_double * 1)()
        status = libgr.gr_get_d(c, self._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"x = {self} is not an integer in {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"x = {self} is not a float in {self.parent()}")
        return c[0]

    def inv(self):
        """
        Multiplicative inverse of this element.

            >>> QQ(3).inv()
            1/3
            >>> QQ(0).inv()
            Traceback (most recent call last):
              ...
            ValueError: inv(x) is not defined for x = 0 over Rational field (fmpq)

        """
        return self._unary_op(self, libgr.gr_inv, "inv(x)")

    def sqrt(self):
        """
        Square root of this element.

            >>> ZZ(4).sqrt()
            2
            >>> ZZ(2).sqrt()
            Traceback (most recent call last):
              ...
            ValueError: sqrt(x) is not defined for x = 2 over Integer ring (fmpz)
            >>> QQbar(2).sqrt()
            Root a = 1.41421 of a^2-2
            >>> (QQ(25)/16).sqrt()
            5/4
            >>> QQbar(-1).sqrt()
            Root a = 1.00000*I of a^2+1
            >>> RR(-1).sqrt()
            Traceback (most recent call last):
              ...
            ValueError: sqrt(x) is not defined for x = -1.000000000000000 over Real numbers (arb, prec = 53)
            >>> RF(-1).sqrt()
            nan

        """
        return self._unary_op(self, libgr.gr_sqrt, "sqrt(x)")

    def rsqrt(self):
        """
        Reciprocal square root of this element.

            >>> QQ(25).rsqrt()
            1/5
        """
        return self._unary_op(self, libgr.gr_rsqrt, "rsqrt(x)")

    def floor(self):
        r"""
        Floor function: closest integer in the direction of `-\infty`.

            >>> (QQ(3) / 2).floor()
            1
            >>> (QQ(3) / 2).ceil()
            2
            >>> (QQ(3) / 2).nint()
            2
            >>> (QQ(3) / 2).trunc()
            1
        """
        return self._unary_op(self, libgr.gr_floor, "floor(x)")

    def ceil(self):
        r"""
        Ceiling function: closest integer in the direction of `+\infty`.

            >>> (QQ(3) / 2).ceil()
            2
        """
        return self._unary_op(self, libgr.gr_ceil, "ceil(x)")

    def trunc(self):
        r"""
        Truncate to integer: closest integer in the direction of zero.

            >>> (QQ(3) / 2).trunc()
            1
        """
        return self._unary_op(self, libgr.gr_trunc, "trunc(x)")

    def nint(self):
        r"""
        Nearest integer function: nearest integer, rounding to
        even on a tie.

            >>> (QQ(3) / 2).nint()
            2
        """
        return self._unary_op(self, libgr.gr_nint, "nint(x)")

    def abs(self):
        return self._unary_op(self, libgr.gr_abs, "abs(x)")

    def i(self):
        """
        Imaginary unit as an element of the same domain as self.

            >>> QQbar().i()
            Root a = 1.00000*I of a^2+1
            >>> QQ().i()
            Traceback (most recent call last):
              ...
            ValueError: i is not defined in Rational field (fmpq)

        """
        return self._constant(self, libgr.gr_i, "i")

    def conj(self):
        """
        Complex conjugate.

            >>> QQbar().i().conj()
            Root a = -1.00000*I of a^2+1
            >>> CC(-2).log().conj()
            ([0.693147180559945 +/- 4.12e-16] + [-3.141592653589793 +/- 3.39e-16]*I)
            >>> QQ(3).conj()
            3
        """
        return self._unary_op(self, libgr.gr_conj, "conj(x)")

    def re(self):
        """
        Real part.

            >>> QQ(1).re()
            1
            >>> (QQbar(-1) ** (QQ(1) / 3)).re()
            1/2
        """
        return self._unary_op(self, libgr.gr_re, "re(x)")

    def im(self):
        """
        Imaginary part.

            >>> QQ(1).im()
            0
            >>> (QQbar(-1) ** (QQ(1) / 3)).im()
            Root a = 0.866025 of 4*a^2-3
        """
        return self._unary_op(self, libgr.gr_im, "im(x)")

    def sgn(self):
        """
        Sign function.

            >>> QQ(-5).sgn()
            -1
            >>> CC(-10).sqrt().sgn()
            1.000000000000000*I
        """
        return self._unary_op(self, libgr.gr_sgn, "sgn(x)")

    def csgn(self):
        """
        Real-valued extension of the sign function: gives
        the sign of the real part when nonzero, and the sign of the
        imaginary part when on the imaginary axis.

            >>> QQbar(-10).sqrt().csgn()
            1
            >>> (-QQbar(-10).sqrt()).csgn()
            -1
        """
        return self._unary_op(self, libgr.gr_csgn, "csgn(x)")

    def pi(self):
        """
        The number pi as an element of the same domain as self.

            >>> RR().pi()
            [3.141592653589793 +/- 3.39e-16]
            >>> QQbar().pi()
            Traceback (most recent call last):
              ...
            ValueError: pi is not defined in Complex algebraic numbers (qqbar)
        """
        return self._constant(self, libgr.gr_pi, "pi")

    def exp(self):
        """
        Exponential function.

            >>> RR(1).exp()
            [2.718281828459045 +/- 5.41e-16]
            >>> RR_ca(1).exp()
            2.71828 {a where a = 2.71828 [Exp(1)]}
            >>> QQ(0).exp()
            1
            >>> QQ(1).exp()
            Traceback (most recent call last):
              ...
            NotImplementedError: unable to compute exp(x) for x = 1 over Rational field (fmpq)

        """
        return self._unary_op(self, libgr.gr_exp, "exp(x)")

    def log(self):
        return self._unary_op(self, libgr.gr_log, "log(x)")

    def sin(self):
        return self._unary_op(self, libgr.gr_sin, "sin(x)")

    def cos(self):
        return self._unary_op(self, libgr.gr_cos, "cos(x)")

    def tan(self):
        return self._unary_op(self, libgr.gr_tan, "tan(x)")

    def sinh(self):
        return self._unary_op(self, libgr.gr_sinh, "sinh(x)")

    def cosh(self):
        return self._unary_op(self, libgr.gr_cosh, "cosh(x)")

    def tanh(self):
        return self._unary_op(self, libgr.gr_tanh, "tanh(x)")

    def atan(self):
        return self._unary_op(self, libgr.gr_atan, "atan(x)")

    def exp_pi_i(self):
        r"""
        `\exp(\pi i x)` evaluated at self.

            >>> (QQbar(1) / 3).exp_pi_i()
            Root a = 0.500000 + 0.866025*I of a^2-a+1
            >>> (QQbar(2).sqrt()).exp_pi_i()
            Traceback (most recent call last):
              ...
            ValueError: exp_pi_i(x) is not defined for x = Root a = 1.41421 of a^2-2 over Complex algebraic numbers (qqbar)
        """
        return self._unary_op(self, libgr.gr_exp_pi_i, "exp_pi_i(x)")

    def log_pi_i(self):
        r"""
        `\log(x) / (\pi i)` evaluated at self.

            >>> (QQbar(-1) ** (QQbar(7) / 5)).log_pi_i()
            -3/5
            >>> (QQbar(1) / 2).log_pi_i()
            Traceback (most recent call last):
              ...
            ValueError: log_pi_i(x) is not defined for x = 1/2 over Complex algebraic numbers (qqbar)
        """
        return self._unary_op(self, libgr.gr_log_pi_i, "log_pi_i(x)")

    def sin_pi(self):
        r"""
        `\sin(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).sin_pi()
            Root a = 0.866025 of 4*a^2-3
        """
        return self._unary_op(self, libgr.gr_sin_pi, "sin_pi(x)")

    def cos_pi(self):
        r"""
        `\cos(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).cos_pi()
            1/2
        """
        return self._unary_op(self, libgr.gr_cos_pi, "cos_pi(x)")

    def tan_pi(self):
        r"""
        `\tan(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).tan_pi()
            Root a = 1.73205 of a^2-3
        """
        return self._unary_op(self, libgr.gr_tan_pi, "tan_pi(x)")

    def cot_pi(self):
        r"""
        `\cot(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).cot_pi()
            Root a = 0.577350 of 3*a^2-1
        """
        return self._unary_op(self, libgr.gr_cot_pi, "cot_pi(x)")

    def sec_pi(self):
        r"""
        `\sec(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).sec_pi()
            2
        """
        return self._unary_op(self, libgr.gr_sec_pi, "sec_pi(x)")

    def csc_pi(self):
        r"""
        `\csc(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).csc_pi()
            Root a = 1.15470 of 3*a^2-4
        """
        return self._unary_op(self, libgr.gr_csc_pi, "csc_pi(x)")

    def asin_pi(self):
        return self._unary_op(self, libgr.gr_asin_pi, "asin_pi(x)")

    def acos_pi(self):
        return self._unary_op(self, libgr.gr_acos_pi, "acos_pi(x)")

    def atan_pi(self):
        return self._unary_op(self, libgr.gr_atan_pi, "atan_pi(x)")

    def acot_pi(self):
        return self._unary_op(self, libgr.gr_acot_pi, "acot_pi(x)")

    def asec_pi(self):
        return self._unary_op(self, libgr.gr_asec_pi, "asec_pi(x)")

    def acsc_pi(self):
        return self._unary_op(self, libgr.gr_acsc_pi, "acsc_pi(x)")

    def erf(self):
        return self._unary_op(self, libgr.gr_erf, "erf(x)")

    def erfi(self):
        return self._unary_op(self, libgr.gr_erfi, "erfi(x)")

    def erfc(self):
        return self._unary_op(self, libgr.gr_erfc, "erfc(x)")

    def gamma(self):
        return self._unary_op(self, libgr.gr_gamma, "gamma(x)")

    def lgamma(self):
        return self._unary_op(self, libgr.gr_lgamma, "lgamma(x)")

    def rgamma(self):
        return self._unary_op(self, libgr.gr_rgamma, "lgamma(x)")

    def digamma(self):
        return self._unary_op(self, libgr.gr_digamma, "digamma(x)")

    def zeta(self):
        return self._unary_op(self, libgr.gr_zeta, "zeta(x)")


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

class GaussianIntegerRing_fmpzi(gr_ctx):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fmpzi(self._ref)
        self._elem_type = fmpzi

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
        gr_ctx.__init__(self)
        #if libgr.gr_ctx_is_ring(coefficient_ring._ref) != T_TRUE:
        #    raise ValueError("coefficient structure must be a ring")
        libgr.gr_ctx_init_polynomial(self._ref, coefficient_ring._ref)
        coefficient_ring._refcount += 1
        self._coefficient_ring = coefficient_ring
        self._elem_type = gr_poly

    def __del__(self):
        self._coefficient_ring._decrement_refcount()

class fmpz(gr_elem):

    _struct_type = fmpz_struct

    @staticmethod
    def _default_context():
        return ZZ

    def __index__(self):
        return fmpz_to_python_int(self._ref)

    def __int__(self):
        return fmpz_to_python_int(self._ref)

    def is_prime(self):
        return bool(libflint.fmpz_is_prime(self._ref))

class fmpq(gr_elem):
    _struct_type = fmpq_struct

    @staticmethod
    def _default_context():
        return QQ

class fmpzi(gr_elem):
    _struct_type = fmpzi_struct

    @staticmethod
    def _default_context():
        return ZZi

class qqbar(gr_elem):
    _struct_type = qqbar_struct

    @staticmethod
    def _default_context():
        return QQbar

class ca(gr_elem):
    _struct_type = ca_struct

    @staticmethod
    def _default_context():
        return CC_ca

class arb(gr_elem):
    _struct_type = arb_struct

    @staticmethod
    def _default_context():
        return RR_arb

class acb(gr_elem):
    _struct_type = acb_struct

    @staticmethod
    def _default_context():
        return CC_acb

class gr_arf_ctx(gr_ctx):

    @property
    def prec(self):
        return libgr.gr_ctx_arf_get_prec(self._ref)

    @prec.setter
    def prec(self, prec):
        libgr.gr_ctx_arf_set_prec(self._ref, prec)


class RealFloat_arf(gr_arf_ctx):
    def __init__(self, prec=53):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_real_float_arf(self._ref, prec)
        self._elem_type = arf

class ComplexFloat_acf(gr_arf_ctx):
    def __init__(self, prec=53):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_float_acf(self._ref, prec)
        self._elem_type = acf

class arf(gr_elem):
    _struct_type = arf_struct

    @staticmethod
    def _default_context():
        return RF

    def __hash__(self):
        # todo
        return hash(float(str(self)))

class acf(gr_elem):
    _struct_type = acf_struct

    @staticmethod
    def _default_context():
        return CF



"""
.. function:: int gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)
.. function:: int gr_ctx_fq_degree(slong * deg, gr_ctx_t ctx)
.. function:: int gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)
"""

class FiniteField_base(gr_ctx):

    def prime(self):
        res = ZZ()
        status = libgr.gr_ctx_fq_prime(res._ref, self._ref, self._ref)
        assert not status
        return res

    def degree(self):
        res = ZZ()
        c = c_slong()
        status = libgr.gr_ctx_fq_degree(ctypes.byref(c), self._ref, self._ref)
        assert not status
        libflint.fmpz_set_si(res._ref, c)
        return res

    def order(self):
        res = ZZ()
        status = libgr.gr_ctx_fq_order(res._ref, self._ref, self._ref)
        assert not status
        return res


class FiniteField_fq(FiniteField_base):
    def __init__(self, p, n):
        gr_ctx.__init__(self)
        p = ZZ(p)
        n = int(n)
        assert p.is_prime()
        assert n >= 1
        libgr.gr_ctx_init_fq(self._ref, p._ref, n, None)
        self._elem_type = fq

class FiniteField_fq_nmod(FiniteField_base):
    def __init__(self, p, n):
        gr_ctx.__init__(self)
        p = ZZ(p)
        n = int(n)
        assert p.is_prime()
        assert n >= 1
        libgr.gr_ctx_init_fq_nmod(self._ref, p._ref, n, None)
        self._elem_type = fq_nmod

class FiniteField_fq_zech(FiniteField_base):
    def __init__(self, p, n):
        gr_ctx.__init__(self)
        p = ZZ(p)
        n = int(n)
        assert p.is_prime()
        assert n >= 1
        libgr.gr_ctx_init_fq_zech(self._ref, p._ref, n, None)
        self._elem_type = fq_zech


class fq_elem(gr_elem):

    def gen(self):
        return self._constant(self, libgr.gr_fq_gen, "gen")

    #def frobenius(self):
    #    return self._binary_op_si(self, libgr.gr_fq_frobenius, "frobenius")

    def multiplicative_order(self):
        return self._unary_op_get_fmpz(self, libgr.gr_fq_multiplicative_order, "multiplicative_order")

    def norm(self):
        return self._unary_op_get_fmpz(self, libgr.gr_fq_norm, "norm")

    def trace(self):
        return self._unary_op_get_fmpz(self, libgr.gr_fq_trace, "trace")

    def is_primitive(self):
        return self._unary_predicate(self, libgr.gr_fq_is_primitive, "is_primitive")

    def pth_root(self):
        return self._unary_op(self, libgr.gr_fq_pth_root, "pth_root")


class fq(fq_elem):
    _struct_type = fq_struct

class fq_nmod(fq_elem):
    _struct_type = fq_nmod_struct

class fq_zech(fq_elem):
    _struct_type = fq_zech_struct



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

    def __call__(self, x, algorithm=None):
        f_R = self.parent()._coefficient_ring
        x_R = x.parent()
        res = x_R()
        if f_R is x_R:
            if algorithm is None:
                status = libgr.gr_poly_evaluate(res._ref, self._ref, x._ref, x_R._ref, f_R._ref)
            elif algorithm == "rectangular":
                status = libgr.gr_poly_evaluate_rectangular(res._ref, self._ref, x._ref, x_R._ref, f_R._ref)
            else:
                raise ValueError
        else:
            if algorithm is None:
                status = libgr.gr_poly_evaluate_other_horner(res._ref, self._ref, x._ref, x_R._ref, f_R._ref)
            elif algorithm == "rectangular":
                status = libgr.gr_poly_evaluate_other_rectangular(res._ref, self._ref, x._ref, x_R._ref, f_R._ref)
            else:
                raise ValueError
        if status:
            raise NotImplementedError
        return res


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


class DirichletGroup_dirichlet_char(gr_ctx_ca):
    def __init__(self, q, **kwargs):
        # todo: automatic range checking with ctypes int -> c_ulong cast?
        if q <= 0:
            raise ValueError(f"modulus must not be zero")
        if q > UWORD_MAX:
            raise NotImplementedError(f"only word-size moduli are supported")
        gr_ctx.__init__(self)
        status = libgr.gr_ctx_init_dirichlet_group(self._ref, q)
        if status & GR_UNABLE: raise NotImplementedError(f"modulus with prime factor p > 10^12 is not currently supported")
        if status & GR_DOMAIN: raise ValueError(f"modulus must not be zero")
        self._elem_type = dirichlet_char

class dirichlet_char(gr_elem):
    _struct_type = dirichlet_char_struct


class SymmetricGroup_perm(gr_ctx_ca):
    def __init__(self, n, **kwargs):
        # todo: automatic range checking with ctypes int -> c_ulong cast?
        if n < 0:
            raise ValueError(f"n must be positive")
        if n > WORD_MAX:
            raise NotImplementedError(f"only word-size moduli n are supported")
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_perm(self._ref, n)
        self._elem_type = perm

class perm(gr_elem):
    _struct_type = perm_struct



class Mat(gr_ctx):
    """
    Parent class for matrix domains.

    There are two kinds of matrix domains:

    - Mat(R), the set of matrices of any size over the domain R.
    - Mat(R, n, m), the set of n x m matrices over the domain R.
      If R is a ring and n = m, then this is also a ring.

    While Mat(R) may be more convenient, e.g. for representing linear
    transformations of arbitrary dimension under a single parent,
    fixed-shape matrix domains have advantages such as allowing
    automatic conversion from scalars to scalar matrices of the
    right size.

        >>> Mat(ZZ)
        Matrices (any shape) over Integer ring (fmpz)
        >>> Mat(ZZ, 2)
        Ring of 2 x 2 matrices over Integer ring (fmpz)
        >>> Mat(ZZ, 2, 3)
        Space of 2 x 3 matrices over Integer ring (fmpz)
        >>> Mat(ZZ)([[1, 2, 3], [4, 5, 6]])
        [[1, 2, 3],
        [4, 5, 6]]
        >>> Mat(ZZ, 2, 2)(5)
        [[5, 0],
        [0, 5]]

    """

    def __init__(self, element_domain, nrows=None, ncols=None):
        assert isinstance(element_domain, gr_ctx)
        assert (nrows is None) or (0 <= nrows <= WORD_MAX)
        assert (ncols is None) or (0 <= ncols <= WORD_MAX)
        gr_ctx.__init__(self)
        if nrows is None and ncols is None:
            libgr.gr_ctx_init_matrix_domain(self._ref, element_domain._ref)
        else:
            if ncols is None:
                ncols = nrows
            libgr.gr_ctx_init_matrix_space(self._ref, element_domain._ref, nrows, ncols)
        self._element_ring = element_domain
        self._elem_type = gr_mat
        self._element_ring._refcount += 1

    def __del__(self):
        self._element_ring._decrement_refcount()

def MatrixRing(element_ring, n):
    assert isinstance(element_ring, gr_ctx)
    assert 0 <= n <= WORD_MAX
    if libgr.gr_ctx_is_ring(element_ring._ref) != T_TRUE:
        raise ValueError("element structure must be a ring")
    return Mat(element_ring, n)


class gr_mat(gr_elem):

    _struct_type = gr_mat_struct

    def __init__(self, *args, **kwargs):
        context = kwargs['context']
        gr_elem.__init__(self, None, context)
        element_ring = context._element_ring
        if kwargs.get('random'):
            libgr.gr_randtest(self._ref, ctypes.byref(_flint_rand), self._ctx)
            return

        if len(args) == 1:
            val = args[0]
            if val is not None:
                status = GR_UNABLE
                if isinstance(val, (list, tuple)):
                    m = len(val)
                    n = 0
                    if m != 0:
                        if not isinstance(val[0], (list, tuple)):
                            raise TypeError("single input to gr_mat must be a list of lists")
                        n = len(val[0])
                        for i in range(1, m):
                            if len(val[i]) != n:
                                raise ValueError("input rows have different lengths")
                    status = libgr._gr_mat_check_resize(self._ref, m, n, self._ctx)
                    if not status:
                        for i in range(m):
                            row = val[i]
                            for j in range(n):
                                x = element_ring(row[j])
                                ijptr = libgr.gr_mat_entry_ptr(self._ref, i, j, x._ctx)
                                status |= libgr.gr_set(ijptr, x._ref, x._ctx)
                elif libgr.gr_ctx_matrix_is_fixed_size(self._ctx) == T_TRUE:
                    if not isinstance(val, gr_elem):
                        val = element_ring(val)
                    status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
                elif isinstance(val, gr_mat):
                    status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
                if status:
                    if status & GR_UNABLE: raise NotImplementedError
                    if status & GR_DOMAIN: raise ValueError
        elif len(args) in (2, 3):
            if len(args) == 2:
                m, n = args
                entries = None
            else:
                m, n, entries = args
                entries = list(entries)
                if len(entries) != m*n:
                    raise ValueError("list of entries has the wrong length")
            status = libgr._gr_mat_check_resize(self._ref, m, n, self._ctx)
            if status:
                if status & GR_UNABLE: raise NotImplementedError
                if status & GR_DOMAIN: raise ValueError("wrong matrix shape for this domain")
            if entries is None:
                status = libgr.gr_mat_zero(self._ref, element_ring._ref)
                if status:
                    if status & GR_UNABLE: raise NotImplementedError
                    if status & GR_DOMAIN: raise ValueError
            else:
                for i in range(m):
                    for j in range(n):
                        x = element_ring(entries[i*n + j])
                        ijptr = libgr.gr_mat_entry_ptr(self._ref, i, j, x._ctx)
                        status = libgr.gr_set(ijptr, x._ref, x._ctx)
                        if status:
                            if status & GR_UNABLE: raise NotImplementedError
                            if status & GR_DOMAIN: raise ValueError

    def nrows(self):
        return self._data.r

    def ncols(self):
        return self._data.c

    def shape(self):
        return (self._data.r, self._data.c)

    def __getitem__(self, ij):
        i, j = ij
        i = int(i)
        j = int(j)
        assert 0 <= i < self.nrows()
        assert 0 <= j < self.ncols()
        element_ring = self.parent()._element_ring
        res = element_ring()
        ijptr = libgr.gr_mat_entry_ptr(self._ref, i, j, res._ctx)
        status = libgr.gr_set(res._ref, ijptr, res._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def __setitem__(self, ij, v):
        i, j = ij
        i = int(i)
        j = int(j)
        assert 0 <= i < self.nrows()
        assert 0 <= j < self.ncols()
        element_ring = self.parent()._element_ring
        # todo: avoid copy
        x = element_ring(v)
        ijptr = libgr.gr_mat_entry_ptr(self._ref, i, j, x._ctx)
        status = libgr.gr_set(ijptr, x._ref, x._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return x

    def det(self, algorithm=None):
        element_ring = self.parent()._element_ring
        res = element_ring()
        if algorithm is None:
            status = libgr.gr_mat_det(res._ref, self._ref, element_ring._ref)
        elif algorithm == "lu":
            status = libgr.gr_mat_det_lu(res._ref, self._ref, element_ring._ref)
        elif algorithm == "fflu":
            status = libgr.gr_mat_det_fflu(res._ref, self._ref, element_ring._ref)
        elif algorithm == "berkowitz":
            status = libgr.gr_mat_det_berkowitz(res._ref, self._ref, element_ring._ref)
        elif algorithm == "cofactor":
            status = libgr.gr_mat_det_cofactor(res._ref, self._ref, element_ring._ref)
        else:
            raise ValueError("unknown algorithm")
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def pascal(self, triangular=0):
        element_ring = self.parent()._element_ring
        res = self.parent()()
        status = libgr.gr_mat_pascal(res._ref, triangular, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def stirling(self, kind=0):
        element_ring = self.parent()._element_ring
        res = self.parent()()
        status = libgr.gr_mat_stirling(res._ref, kind, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def hilbert(self):
        element_ring = self.parent()._element_ring
        res = self.parent()()
        status = libgr.gr_mat_hilbert(res._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def hadamard(self):
        element_ring = self.parent()._element_ring
        res = self.parent()()
        status = libgr.gr_mat_hadamard(res._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def charpoly(self, R, algorithm=None):
        mat_ring = self.parent()
        poly_ring = R
        element_ring = mat_ring._element_ring
        poly_element_ring = poly_ring._coefficient_ring
        assert element_ring is poly_element_ring
        res = poly_ring()
        if algorithm is None:
            status = libgr.gr_mat_charpoly(res._ref, self._ref, element_ring._ref)
        elif algorithm == "berkowitz":
            status = libgr.gr_mat_charpoly_berkowitz(res._ref, self._ref, element_ring._ref)
        elif algorithm == "gauss":
            status = libgr.gr_mat_charpoly_gauss(res._ref, self._ref, element_ring._ref)
        elif algorithm == "householder":
            status = libgr.gr_mat_charpoly_householder(res._ref, self._ref, element_ring._ref)
        elif algorithm == "danilevsky":
            status = libgr.gr_mat_charpoly_danilevsky(res._ref, self._ref, element_ring._ref)
        elif algorithm == "faddeev":
            status = libgr.gr_mat_charpoly_faddeev(res._ref, None, self._ref, element_ring._ref)
        elif algorithm == "faddeev_bsgs":
            status = libgr.gr_mat_charpoly_faddeev_bsgs(res._ref, None, self._ref, element_ring._ref)
        else:
            raise ValueError("unknown algorithm")
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def hessenberg(self, algorithm=None):
        """
        Return this matrix reduced to upper Hessenberg form::

            >>> B = Mat(QQ, 3, 3)([[4, 2, 3], [-1, 5, -3], [-4, 1, 2]]);
            >>> B.hessenberg()
            [[4, 14, 3],
            [-1, -7, -3],
            [0, 37, 14]]

        Options:
        - algorithm: ``None`` (default), ``"gauss"`` or ``"householder"``

        """
        element_ring = self.parent()._element_ring
        res = self.parent()()
        if algorithm is None:
            status = libgr.gr_mat_hessenberg(res._ref, self._ref, element_ring._ref)
        elif algorithm == "gauss":
            status = libgr.gr_mat_hessenberg_gauss(res._ref, self._ref, element_ring._ref)
        elif algorithm == "householder":
            status = libgr.gr_mat_hessenberg_householder(res._ref, self._ref, element_ring._ref)
        else:
            raise ValueError("unknown algorithm")
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res


    #def __getitem__(self, i):
    #    pass



libgr.gr_mat_entry_ptr.argtypes = (ctypes.c_void_p, c_slong, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_mat_entry_ptr.restype = ctypes.POINTER(ctypes.c_char)

libgr.gr_vec_entry_ptr.restype = ctypes.POINTER(ctypes.c_char)


class Vec(gr_ctx):
    """
    Parent class for vector domains.
    """

    def __init__(self, element_domain, n=None):
        assert isinstance(element_domain, gr_ctx)
        assert (n is None) or (0 <= n <= WORD_MAX)
        gr_ctx.__init__(self)
        if n is None:
            libgr.gr_ctx_init_vector_gr_vec(self._ref, element_domain._ref)
        else:
            libgr.gr_ctx_init_vector_space_gr_vec(self._ref, element_domain._ref, n)
        self._element_ring = element_domain
        self._elem_type = gr_vec
        self._element_ring._refcount += 1

    def __del__(self):
        self._element_ring._decrement_refcount()



class gr_vec(gr_elem):

    _struct_type = gr_vec_struct

    def __init__(self, *args, **kwargs):
        context = kwargs['context']
        gr_elem.__init__(self, None, context)
        element_ring = context._element_ring
        if kwargs.get('random'):
            libgr.gr_randtest(self._ref, ctypes.byref(_flint_rand), self._ctx)
            return

        if len(args) == 1:
            val = args[0]
            if val is not None:
                status = GR_UNABLE
                if isinstance(val, (list, tuple)):
                    n = len(val)
                    status = libgr._gr_vec_check_resize(self._ref, n, self._ctx)
                    if not status:
                        for i in range(n):
                            x = element_ring(val[i])
                            iptr = libgr.gr_vec_entry_ptr(self._ref, i, x._ctx)
                            status |= libgr.gr_set(iptr, x._ref, x._ctx)
                elif isinstance(val, gr_elem):
                    status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
                if status:
                    if status & GR_UNABLE: raise NotImplementedError
                    if status & GR_DOMAIN: raise ValueError

    def __len__(self):
        return self._data.length

    '''
    def __getitem__(self, ij):
        i, j = ij
        i = int(i)
        j = int(j)
        assert 0 <= i < self.nrows()
        assert 0 <= j < self.ncols()
        element_ring = self.parent()._element_ring
        res = element_ring()
        ijptr = libgr.gr_mat_entry_ptr(self._ref, i, j, res._ctx)
        status = libgr.gr_set(res._ref, ijptr, res._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def __setitem__(self, ij, v):
        i, j = ij
        i = int(i)
        j = int(j)
        assert 0 <= i < self.nrows()
        assert 0 <= j < self.ncols()
        element_ring = self.parent()._element_ring
        # todo: avoid copy
        x = element_ring(v)
        ijptr = libgr.gr_mat_entry_ptr(self._ref, i, j, x._ctx)
        status = libgr.gr_set(ijptr, x._ref, x._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return x
    '''



PolynomialRing = PolynomialRing_gr_poly

ZZ = IntegerRing_fmpz()
QQ = RationalField_fmpq()
ZZi = GaussianIntegerRing_fmpzi()
AA = RealAlgebraicField_qqbar()
AA_ca = RealAlgebraicField_ca()
QQbar = ComplexAlgebraicField_qqbar()
QQbar_ca = ComplexAlgebraicField_ca()
RR = RR_arb = RealField_arb()
CC = CC_acb = ComplexField_acb()
RR_ca = RealField_ca()
CC_ca = ComplexField_ca()

RF = RealFloat_arf()
CF = ComplexFloat_acf()

ZZ_vec = Vec(ZZ)

ZZx = PolynomialRing_gr_poly(ZZ)
QQx = PolynomialRing_gr_poly(QQ)
RRx_ca = PolynomialRing_gr_poly(RR_ca)
CCx_ca = PolynomialRing_gr_poly(CC_ca)
RRx = RRx_arb = PolynomialRing_gr_poly(RR_arb)
CCx = CCx_acb = PolynomialRing_gr_poly(CC_acb)

ModularGroup = ModularGroup_psl2z
DirichletGroup = DirichletGroup_dirichlet_char

PSL2Z = ModularGroup()
SymmetricGroup = SymmetricGroup_perm


def raises(f, exception):
    try:
        f()
    except exception:
        return True
    return False

def test_perm():
    S = SymmetricGroup(3)
    M = Mat(ZZ)
    A = M([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    assert S(A).parent() is S
    assert S(A).inv() == S(A.inv())
    assert raises(lambda: S(-A), ValueError)
    assert raises(lambda: S(M([[0, 1, 0], [1, 0, 0], [0, 0, 0]])), ValueError)
    assert raises(lambda: S(M([[0, 1, 0], [1, 0, 0], [1, 0, 0]])), ValueError)
    assert raises(lambda: S(M([[0, 1, 0], [1, 0, 0], [0, 1, 1]])), ValueError)

def test_psl2z():
    M = Mat(ZZ)
    A = M([[2, 1], [5, 3]])
    a = PSL2Z(A)
    assert a.parent() is PSL2Z
    assert a == PSL2Z(-A)
    assert a.inv() == PSL2Z(A.inv())
    assert raises(lambda: PSL2Z(M([[1], [2]])), ValueError)
    assert raises(lambda: PSL2Z(M([[1, 3, 4], [4, 5, 6]])), ValueError)
    assert raises(lambda: PSL2Z(M([[1, 2], [3, 4]])), ValueError)

def test_matrix():
    M = Mat(ZZ, 2)
    I = M([[1, 0], [0, 1]])
    assert M(1) == M(ZZ(1)) == I == M(2, 2, [1, 0, 0, 1])
    assert raises(lambda: M(3, 1, [1, 2, 3]), ValueError)
    assert 2 * I == M(2 * I) == I + I == 1 + I == I + 1

    assert Mat(ZZ)([[1],[3]]) * ZZ(5) == Mat(ZZ)([[5],[15]])

    M = Mat(ZZ)
    A = M([[1,2,3],[4,5,6]])
    assert A == M(2, 3, [1,2,3,4,5,6])
    assert A == M(2, 3, [1,2,QQ(3),4,5,6])
    assert Mat(ZZ, 2, 3)(A) == A
    assert Mat(QQ, 2, 3)(A) == A
    assert M(2, 1) == M([[0], [0]])
    assert raises(lambda: M(2, 1, [1,2,3]), ValueError)
    assert raises(lambda: M([[QQ(1)/3]]), ValueError)
    assert raises(lambda: Mat(ZZ, 3, 1)(A), ValueError)
    assert Mat(QQ, 2)(M([[1, 2], [3, 4]])) ** 2 == M([[7,10],[15,22]])

    A[1, 2] = 10
    assert A == M([[1,2,3],[4,5,10]])
    assert A[0,1] == 2
    assert raises(lambda: A[3,4], Exception)
    assert raises(lambda: A.__setitem__((3, 4), 1), Exception)

def test_fq():
    Fq = FiniteField_fq(3, 5)
    x = Fq(random=True)
    y = Fq(random=True)
    assert 3*(x+y) == 4*x+3*y-x
    assert Fq.prime() == 3
    assert Fq.degree() == 5
    assert Fq.order() == 243
    assert x.pth_root() ** 3 == x
    assert (x**2).sqrt() in (x, -x)

def test_floor_ceil_trunc_nint():
    assert ZZ(3).floor() == 3
    assert ZZ(3).ceil() == 3
    assert ZZ(3).trunc() == 3
    assert ZZ(3).nint() == 3

    assert QQ(3).floor() == 3
    assert QQ(3).ceil() == 3
    assert QQ(3).trunc() == 3
    assert QQ(3).nint() == 3

    for R in [QQ, QQbar, QQbar_ca, AA, AA_ca, RR, RR_ca, CC, CC_ca, RF]:
        x = R(3) / 2
        assert x.floor() == 1
        assert x.ceil() == 2
        assert x.trunc() == 1
        assert (-x).floor() == -2
        assert (-x).ceil() == -1
        assert (-x).trunc() == -1
        assert x.nint() == 2
        assert (-x).nint() == -2
        assert (x+1).nint() == 2
        assert (x+2).nint() == 4

    for R in [QQbar, QQbar_ca, CC, CC_ca]:
        x = R(3) / 2 + R().i()
        assert x.floor() == 1
        assert x.ceil() == 2
        assert x.trunc() == 1
        assert (-x).floor() == -2
        assert (-x).ceil() == -1
        assert (-x).trunc() == -1
        assert x.nint() == 2
        assert (-x).nint() == -2
        assert (x+1).nint() == 2
        assert (x+2).nint() == 4

def test_zz():
    assert ZZ(1).factor() == (1, [], [])
    assert ZZ(0).factor() == (0, [], [])
    assert (-ZZ(12)).factor() == (-1, [2, 3], [2, 1])

def test_qq():
    assert QQ(1).factor() == (1, [], [])
    assert QQ(0).factor() == (0, [], [])
    assert (-QQ(12)/175).factor() == (-1, [2, 3, 5, 7], [2, 1, -2, -1])

def test_qqbar():
    a = (-23 + 5*ZZi.i())
    assert ZZi(QQbar(a**2).sqrt()) == -a

def test_arb():
    a = arb(2.5)
    assert a  == arb("2.5")
    b = acb(2.5)
    assert a == b
    c = acb(2.5+1j)
    assert c == b + 1j
    assert raises(lambda: arb(2.5+1j), ValueError)

def test_all():

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

    assert ZZx(QQ(5)) == 5

    v = f(ZZ(3))
    assert v == 41
    assert v.parent() is ZZ

    QM2 = Mat(QQ,2,2)
    A = QM2([[1,2],[3,4]])
    v = f(A)
    assert v == QM2([[27,38],[57,84]])
    assert v.parent() is QM2

    A = Mat(RR,2,2)([[1,2],[3,4]])
    B = ZZx(list(range(10)))(A, algorithm="rectangular")
    assert B == Mat(QQ,2,2)([[9596853, 13986714], [20980071, 30576924]])
    assert B.parent() is A.parent()

    assert CF(2+3j) * (1+1j) == CF((2+3j) * (1+1j))

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

