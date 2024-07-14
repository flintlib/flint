import ctypes
import ctypes.util
import sys
import functools

libflint_path = ctypes.util.find_library('flint')
if libflint_path == None:
    raise ValueError('Could not find libflint.')
libflint = ctypes.CDLL(libflint_path)
libcalcium = libarb = libgr = libflint

T_TRUE = 0
T_FALSE = 1
T_UNKNOWN = 2

GR_SUCCESS = 0
GR_DOMAIN = 1
GR_UNABLE = 2

HUGE_LENGTH = 2**40

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

def set_num_threads(n):
    assert n >= 1
    assert n <= 65536
    libflint.flint_set_num_threads(n)

def get_num_threads():
    return libflint.flint_get_num_threads()

class FlintException(Exception):

    __module__ = Exception.__module__

    def __str__(self):
        if isinstance(self.args[0], str):
            return self.args[0]
        ctx, status, rstr, args = self.args[0]
        rstr2 = rstr.replace("$", "")
        argnames = []
        for i, c in enumerate(rstr):
            if c == "$":
                argname = ""
                for j in range(i + 1, len(rstr)):
                    if rstr[j].isalnum():
                        argname += rstr[j]
                    else:
                        break
                argnames.append(argname)
        if status & GR_UNABLE:
            s = "failed to compute " + rstr2 + " in " + "{" + str(ctx) + "}"
        else:
            s = rstr2 + " is not an element of " + "{" + str(ctx) + "}"
        if args:
            s += " for "
            for i, arg in enumerate(args):
                s += "{"
                if argnames:
                    s += argnames[i] + " = "
                else:
                    s += "input "
                argstr = str(arg)
                if len(argstr) > 200:
                    argstr = argstr[:80] + ("{{{...}}}") + argstr[-80:]
                s += argstr
                s += "}"
                if i < len(args) - 1:
                    s += ", "
        return s

class FlintDomainError(ValueError, FlintException):
    """
    Raised when an operation does not have a well-defined result in the target domain.
    """
    __module__ = Exception.__module__

class FlintUnableError(NotImplementedError, FlintException):
    """
    Raised when an operation cannot be performed because the algorithm is not implemented
    or there is insufficient precision, memory, etc.
    """
    __module__ = Exception.__module__


def _handle_error(ctx, status, rstr, *args):
    if status & GR_UNABLE:
        raise FlintUnableError((ctx, status, rstr, args))
    else:
        raise FlintDomainError((ctx, status, rstr, args))



class flint_rand_struct(ctypes.Structure):
    _fields_ = [('__gmp_state', ctypes.c_void_p),
                ('__randval', c_ulong),
                ('__randval2', c_ulong)]

_flint_rand = flint_rand_struct()
libflint.flint_rand_init(ctypes.byref(_flint_rand))

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

class fmpz_mpoly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('exp', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong),
                ('bits', c_slong)]

class fmpz_mpoly_q_struct(ctypes.Structure):
    _fields_ = [('num', fmpz_mpoly_struct),
                ('den', fmpz_mpoly_struct)]


# todo: actually a union
class nf_elem_struct(ctypes.Structure):
    _fields_ = [('poly', fmpq_poly_struct)]

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

class nmod_struct(ctypes.Structure):
    _fields_ = [('val', c_ulong)]

# todo: want different structure for each size
class mpn_mod_struct(ctypes.Structure):
    _fields_ = [('val', c_ulong * 16)]

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

class gr_series_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('alloc', c_slong),
                ('length', c_slong),
                ('error', c_slong)]

class gr_mpoly_struct(ctypes.Structure):
    _fields_ = [('coeffs', ctypes.c_void_p),
                ('exps', ctypes.c_void_p),
                ('length', c_slong),
                ('bits', c_slong),
                ('coeffs_alloc', c_slong),
                ('exps_alloc', c_slong)]

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

libgr.gr_ctx_data_as_ptr.argtypes = (ctypes.c_void_p,)
libgr.gr_ctx_data_as_ptr.restype = ctypes.c_void_p

libgr.gr_set_si.argtypes = (ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_add_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_sub_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_mul_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_div_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_pow_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))

libgr.gr_set_d.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(gr_ctx_struct))

libgr.gr_set_str.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str_n.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmp.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmpabs.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_heap_clear.argtypes = (ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_ctx_init_nmod.argtypes = (ctypes.POINTER(gr_ctx_struct), c_ulong)
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
        Undecidable: unable to decide x == y for x = [1.00000000000000 +/- 3.89e-16], y = 1 over Real numbers (arb, prec = 53)
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

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

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

    def _ctx_predicate(self, op, rstr):
        truth = op(self._ref)
        if _gr_logic == 3:
            return Truth(truth)
        if truth == T_TRUE: return True
        if truth == T_FALSE: return False
        if _gr_logic == 1: return True
        if _gr_logic == -1: return False
        if _gr_logic == 2: return None
        raise Undecidable(f"unable to decide {rstr} for ctx = {self}")

    def is_ring(self):
        """
        Return whether this structure is a ring.

            >>> RR.is_ring()
            True
            >>> PolynomialRing(QQbar).is_ring()
            True
            >>> RF.is_ring()      # floats do not satisfy the ring laws
            False
            >>> Mat(ZZ, 2).is_ring()
            True
            >>> Mat(ZZ, 2, 3).is_ring()   # nonrectangular matrices
            False
            >>> Mat(ZZ).is_ring()       # matrices of mixed shape do not form a ring
            False
            >>> Mat(RF, 2).is_ring()
            False
            >>> Vec(RR, 3).is_ring()
            True
            >>> Vec(RR).is_ring()       # vectors of mixed length do not form a ring
            False
            >>> Vec(RF, 3).is_ring()
            False
            >>> Vec(RF, 0).is_ring()    # empty vectors form a ring
            True
            >>> DirichletGroup(3).is_ring()
            False
            >>> PSL2Z.is_ring()
            False
            >>> SymmetricGroup(5).is_ring()
            False
            >>> PolynomialRing(RF).is_ring()
            False
            >>> PowerSeriesRing(RF).is_ring()
            False
            >>> PowerSeriesModRing(RF, 1).is_ring()
            False
            >>> PowerSeriesModRing(RF, 0).is_ring()    # is the zero ring
            True
        """
        return self._ctx_predicate(libflint.gr_ctx_is_ring, "is_ring")

    def is_commutative_ring(self):
        """
        Return whether this structure is a commutative ring.

            >>> QQbar.is_commutative_ring()
            True
            >>> CC.is_commutative_ring()
            True
            >>> PolynomialRing(QQ).is_commutative_ring()
            True
            >>> Mat(ZZ, 2).is_commutative_ring()
            False
            >>> PolynomialRing(Mat(ZZ, 2)).is_commutative_ring()
            False
            >>> PowerSeriesRing(QQ).is_commutative_ring()
            True
            >>> PowerSeriesRing(Mat(ZZ, 2)).is_commutative_ring()
            False
            >>> Mat(ZZ, 0).is_commutative_ring()
            True
            >>> Mat(ZZ, 1).is_commutative_ring()
            True
            >>> Mat(ZZ).is_commutative_ring()
            False
            >>> Mat(ZZ, 0, 1).is_commutative_ring()
            False
            >>> Mat(ZZmod(2), 2, 2).is_commutative_ring()
            False
            >>> Mat(ZZmod(1), 2, 2).is_commutative_ring()
            True
            >>> Mat(PowerSeriesModRing(ZZ, 0), 2, 2).is_commutative_ring()
            True
            >>> Mat(RF, 1).is_commutative_ring()
            False
            >>> Mat(RF, 0).is_commutative_ring()
            True
            >>> Vec(ZZ, 3).is_commutative_ring()
            True
            >>> Vec(Mat(ZZ, 2), 3).is_commutative_ring()
            False
            >>> FractionField_fmpz_mpoly_q(3).is_commutative_ring()
            True

        """
        return self._ctx_predicate(libflint.gr_ctx_is_commutative_ring, "is_commutative_ring")

    def is_zero_ring(self):
        """
        Return whether this structure is the zero ring.

            >>> ZZ.is_zero_ring()
            False
            >>> ZZmod(1).is_zero_ring()
            True
            >>> Mat(ZZ, 0).is_zero_ring()
            True
            >>> PowerSeriesModRing(ZZ, 0).is_zero_ring()
            True
            >>> Vec(ZZ, 0).is_zero_ring()
            True
        """
        return self._ctx_predicate(libflint.gr_ctx_is_zero_ring, "is_zero_ring")

    def is_integral_domain(self):
        """
        Return whether this structure is an integral domain.

            >>> ZZ.is_integral_domain()
            True
            >>> ZZx.is_integral_domain()
            True
            >>> PowerSeriesModRing(ZZ, 3).is_integral_domain()
            False

        """
        return self._ctx_predicate(libflint.gr_ctx_is_integral_domain, "is_integral_domain")


    def is_field(self):
        """
        Return whether this structure is a field.

            >>> ZZ.is_field()
            False
            >>> QQ.is_field()
            True

        This check is intended to be fast, and some residue rings may
        not perform a primality test automatically since this would be
        expensive. Rather, the user should set a flag manually in the
        constructor for such rings:

            >>> IntegersMod_fmpz_mod(2**257+1).is_field()
            Traceback (most recent call last):
              ...
            Undecidable: unable to decide is_field for ctx = Integers mod 231584178474632390847141970017375815706539969331281128078915168015826259279873 (fmpz)
            >>> IntegersMod_fmpz_mod(2**257+1, n_is_prime=True).is_field()
            True

        """
        return self._ctx_predicate(libflint.gr_ctx_is_field, "is_field")

    def _set_gen_name(self, s):
        status = libflint.gr_ctx_set_gen_name(self._ref, ctypes.c_char_p(str(s).encode('ascii')))
        self._str = None
        assert not status

    def _set_gen_names(self, s):
        arr = (ctypes.c_char_p * len(s))()
        for i in range(len(s)):
            arr[i] = ctypes.c_char_p(str(s[i]).encode('ascii'))
        status = libflint.gr_ctx_set_gen_names(self._ref, arr)
        self._str = None
        assert not status

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

    @property
    def prec(self):
        p = c_slong()
        status = libgr.gr_ctx_get_real_prec(ctypes.byref(p), self._ref)
        assert not status
        return p.value

    @prec.setter
    def prec(self, prec):
        status = libgr.gr_ctx_set_real_prec(self._ref, prec)
        assert not status

    # constants, sequences etc. with elements in this parent
    # todo: element shortcuts to allow both RR.pi() and RR().pi()

    @staticmethod
    def _constant(ctx, op, rstr):
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr)
        return res

    @staticmethod
    def _as_ui(x):
        type_x = type(x)
        if type_x is not int:
            if type_x is not fmpz:
                x = ZZ(x)
            x = int(x)
        assert 0 <= x <= UWORD_MAX
        return x

    @staticmethod
    def _as_si(x):
        type_x = type(x)
        if type_x is not int:
            if type_x is not fmpz:
                x = ZZ(x)
            x = int(x)
        assert WORD_MIN <= x <= WORD_MAX
        return x

    @staticmethod
    def _as_fmpz(x):
        if type(x) is not fmpz:
            x = ZZ(x)
        return x

    def _as_elem(ctx, x):
        if type(x) is not ctx._elem_type or x._ctx_python is not ctx:
            x = ctx(x)
        return x

    def _as_vec(ctx, x):
        if type(x) is not gr_vec or x._ctx_python._element_ring is not ctx:
            x = Vec(ctx)(x)
        return x

    def _unary_op(ctx, x, op, rstr):
        x = ctx._as_elem(x)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _unary_unary_op(ctx, x, op, rstr):
        x = ctx._as_elem(x)
        res1 = ctx._elem_type(context=ctx)
        res2 = ctx._elem_type(context=ctx)
        status = op(res1._ref, res2._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res1, res2

    def _unary_op_with_flag(ctx, x, flag, op, rstr):
        x = ctx._as_elem(x)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, flag, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _unary_unary_op_with_flag(ctx, x, flag, op, rstr):
        x = ctx._as_elem(x)
        res1 = ctx._elem_type(context=ctx)
        res2 = ctx._elem_type(context=ctx)
        status = op(res1._ref, res2._ref, x._ref, flag, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res1, res2

    def _unary_op_fmpz(ctx, x, op, rstr):
        x = ctx._as_fmpz(x)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _unary_op_with_fmpz_fmpq_overloads(ctx, x, op, op_ui=None, op_fmpz=None, op_fmpq=None, rstr=None):
        type_x = type(x)
        res = ctx._elem_type(context=ctx)
        if type_x is not ctx._elem_type or x._ctx_python is not ctx:
            if type_x is fmpq and op_fmpq is not None:
                status = op_fmpq(res._ref, x._ref, ctx._ref)
            elif type_x is fmpz and op_fmpz is not None:
                status = op_fmpz(res._ref, x._ref, ctx._ref)
            elif type_x is int and op_fmpz is not None:
                x = ZZ(x)
                status = op_fmpz(res._ref, x._ref, ctx._ref)
            elif type_x in (fmpz, int) and op_ui is not None:
                try:
                    x = ctx._as_ui(x)
                    op_ui.argtypes = (ctypes.c_void_p, c_ulong, ctypes.c_void_p)
                    status = op_ui(res._ref, x, ctx._ref)
                except:
                    x = ctx(x)
                    status = op(res._ref, x._ref, ctx._ref)
            else:
                x = ctx(x)
                status = op(res._ref, x._ref, ctx._ref)
        else:
            status = op(res._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _binary_op_with_overloads(ctx, x, y, op, op_ui=None, op_fmpz=None, op_fmpq=None, fmpz_op=None, rstr=None):
        type_x = type(x)
        if fmpz_op is not None and (type_x is fmpz or type_x is int):
            x = ZZ(x)
            type_y = type(y)
            if type_y is not ctx._elem_type or y._ctx_python is not ctx:
                y = ctx(y)
            res = ctx._elem_type(context=ctx)
            status = fmpz_op(res._ref, x._ref, y._ref, ctx._ref)
        else:
            if type_x is not ctx._elem_type or x._ctx_python is not ctx:
                x = ctx(x)
            type_y = type(y)
            res = ctx._elem_type(context=ctx)
            if type_y is not ctx._elem_type or y._ctx_python is not ctx:
                if type_y is fmpq and op_fmpq is not None:
                    status = op_fmpq(res._ref, x._ref, y._ref, ctx._ref)
                elif type_y is fmpz and op_fmpz is not None:
                    status = op_fmpz(res._ref, x._ref, y._ref, ctx._ref)
                elif type_y is int and op_fmpz is not None:
                    y = ZZ(y)
                    status = op_fmpz(res._ref, x._ref, y._ref, ctx._ref)
                elif type_y in (fmpz, int) and op_ui is not None:
                    y = ctx._as_ui(y)
                    op_ui.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_ulong, ctypes.c_void_p)
                    status = op_ui(res._ref, x._ref, y, ctx._ref)
                else:
                    y = ctx(y)
                    status = op(res._ref, x._ref, y._ref, ctx._ref)
            else:
                status = op(res._ref, x._ref, y._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _binary_op(ctx, x, y, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y)
        return res

    def _binary_binary_op(ctx, x, y, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        res1 = ctx._elem_type(context=ctx)
        res2 = ctx._elem_type(context=ctx)
        status = op(res1._ref, res2._ref, x._ref, y._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y)
        return res1, res2

    def _binary_op_with_flag(ctx, x, y, flag, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, flag, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y)
        return res

    def _ternary_op(ctx, x, y, z, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        z = ctx._as_elem(z)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, z._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y, z)
        return res

    def _ternary_op_with_flag(ctx, x, y, z, flag, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        z = ctx._as_elem(z)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, z._ref, flag, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y, z)
        return res

    def _quaternary_op(ctx, x, y, z, w, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        z = ctx._as_elem(z)
        w = ctx._as_elem(w)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, z._ref, w._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y, z, w)
        return res

    def _quaternary_op_with_flag(ctx, x, y, z, w, flag, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        z = ctx._as_elem(z)
        w = ctx._as_elem(w)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, z._ref, w._ref, flag, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y, z, w)
        return res

    def _ternary_unary_op(ctx, x, op, rstr):
        x = ctx._as_elem(x)
        res1 = ctx._elem_type(context=ctx)
        res2 = ctx._elem_type(context=ctx)
        res3 = ctx._elem_type(context=ctx)
        status = op(res1._ref, res2._ref, res3._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return (res1, res2, res3)

    def _quaternary_unary_op(ctx, x, op, rstr):
        x = ctx._as_elem(x)
        res1 = ctx._elem_type(context=ctx)
        res2 = ctx._elem_type(context=ctx)
        res3 = ctx._elem_type(context=ctx)
        res4 = ctx._elem_type(context=ctx)
        status = op(res1._ref, res2._ref, res3._ref, res4._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return (res1, res2, res3, res4)

    def _quaternary_binary_op(ctx, x, y, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_elem(y)
        res1 = ctx._elem_type(context=ctx)
        res2 = ctx._elem_type(context=ctx)
        res3 = ctx._elem_type(context=ctx)
        res4 = ctx._elem_type(context=ctx)
        status = op(res1._ref, res2._ref, res3._ref, res4._ref, x._ref, y._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y)
        return (res1, res2, res3, res4)

    def _binary_op_fmpz(ctx, x, y, op, rstr):
        x = ctx._as_elem(x)
        y = ctx._as_fmpz(y)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, y._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y)
        return res

    def _ui_binary_op(ctx, n, x, op, rstr):
        n = ctx._as_ui(n)
        x = ctx._as_elem(x)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, n, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, n, x)
        return res

    def _op_fmpz(ctx, x, op, rstr):
        x = ctx._as_fmpz(x)
        res = ctx._elem_type(context=ctx)
        status = op(res._ref, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _op_ui(ctx, x, op, rstr):
        x = ctx._as_ui(x)
        res = ctx._elem_type(context=ctx)
        op.argtypes = (ctypes.c_void_p, c_ulong, ctypes.c_void_p)
        status = op(res._ref, x, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x)
        return res

    def _op_uiui(ctx, x, y, op, rstr):
        x = ctx._as_ui(x)
        y = ctx._as_ui(y)
        res = ctx._elem_type(context=ctx)
        op.argtypes = (ctypes.c_void_p, c_ulong, c_ulong, ctypes.c_void_p)
        status = op(res._ref, x, y, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, y)
        return res

    def _op_vec_ctx(ctx, op, rstr):
        op.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
        res = Vec(ctx)()
        status = op(res._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr)
        return res

    def _op_vec_len(ctx, n, op, rstr):
        n = ctx._as_si(n)
        assert n >= 0
        assert n <= HUGE_LENGTH
        op.argtypes = (ctypes.c_void_p, c_slong, ctypes.c_void_p)
        res = Vec(ctx)()
        assert not libgr.gr_vec_set_length(res._ref, n, ctx._ref)
        status = op(libgr.gr_vec_entry_ptr(res._ref, 0, ctx._ref), n, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, n)
        return res

    def _op_vec_arg_len(ctx, x, n, op, rstr):
        x = ctx._as_elem(x)
        n = ctx._as_si(n)
        assert n >= 0
        assert n <= HUGE_LENGTH
        op.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.c_void_p)
        res = Vec(ctx)()
        assert not libgr.gr_vec_set_length(res._ref, n, ctx._ref)
        status = op(libgr.gr_vec_entry_ptr(res._ref, 0, ctx._ref), x._ref, n, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, n)
        return res

    def _op_vec_ui_len(ctx, x, n, op, rstr):
        x = ctx._as_ui(x)
        n = ctx._as_si(n)
        assert n >= 0
        assert n <= HUGE_LENGTH
        op.argtypes = (ctypes.c_void_p, c_ulong, c_slong, ctypes.c_void_p)
        res = Vec(ctx)()
        assert not libgr.gr_vec_set_length(res._ref, n, ctx._ref)
        status = op(libgr.gr_vec_entry_ptr(res._ref, 0, ctx._ref), x, n, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, n)
        return res

    def _op_vec_fmpz_len(ctx, x, n, op, rstr):
        x = ctx._as_fmpz(x)
        n = ctx._as_si(n)
        assert n >= 0
        assert n <= HUGE_LENGTH
        op.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.c_void_p)
        res = Vec(ctx)()
        assert not libgr.gr_vec_set_length(res._ref, n, ctx._ref)
        status = op(libgr.gr_vec_entry_ptr(res._ref, 0, ctx._ref), x._ref, n, ctx._ref)
        if status:
            _handle_error(ctx, status, rstr, x, n)
        return res

    def gen(ctx):
        """
        Gives a generator of this ring.

            >>> ZZi.gen()
            I
            >>> ZZx.gen()
            x
            >>> FiniteField_fq(3, 5).gen() ** 9
            2*a^4+a+2
            >>> (1 + PowerSeriesModRing(QQ, 5).gen())**10
            1 + 10*x + 45*x^2 + 120*x^3 + 210*x^4 (mod x^5)
            >>> (1 + PowerSeriesRing(QQ, 5).gen())**10
            1 + 10*x + 45*x^2 + 120*x^3 + 210*x^4 + O(x^5)

        """
        return ctx._constant(ctx, libgr.gr_gen, "gen")

    def gens(ctx, recursive=False):
        """
        Gives a vector of generators of this ring. If recursive=True,
        includes the generators of all base rings.

            >>> ZZx.gens()
            [x]
            >>> PolynomialRing(QQ, "v").gens()
            [v]
            >>> FiniteField_fq(3, 2).gens()
            [a]
            >>> NumberField(ZZx.gen()**2+1).gens()
            [a]
            >>> ZZ.gens()
            []
            >>> PolynomialRing(PolynomialRing(ZZi, "x"), "y").gens(recursive=True)
            [I, x, y]
            >>> PolynomialRing_fmpz_mpoly(3).gens()
            [x1, x2, x3]
            >>> PolynomialRing(PowerSeriesRing(ZZ, var="b"), "t").gens()
            [t]
            >>> PolynomialRing(PowerSeriesRing(ZZ, var="b"), "t").gens(recursive=True)
            [b, t]
            >>> PowerSeriesRing(ZZx, 3, var="y").gens()
            [y]
            >>> PowerSeriesRing(ZZx, 3, var="y").gens(recursive=True)
            [x, y]

        """
        if recursive:
            return ctx._op_vec_ctx(libgr.gr_gens_recursive, "gens")
        else:
            return ctx._op_vec_ctx(libgr.gr_gens, "gens")

    def zero(ctx):
        """
        The zero element of this domain.

            >>> ZZ.zero()
            0
            >>> Vec(ZZ,3).zero()
            [0, 0, 0]
            >>> Vec(ZZ).zero()
            Traceback (most recent call last):
              ...
            FlintDomainError: zero is not an element of {Vectors (any length) over Integer ring (fmpz)}
        """
        return ctx._constant(ctx, libgr.gr_zero, "zero")

    def one(ctx):
        """
        The one element of this domain.

            >>> ZZ.one()
            1
            >>> Vec(ZZ,3).one()
            [1, 1, 1]
            >>> Vec(ZZ).one()
            Traceback (most recent call last):
              ...
            FlintDomainError: one is not an element of {Vectors (any length) over Integer ring (fmpz)}
        """
        return ctx._constant(ctx, libgr.gr_one, "one")

    def neg_one(ctx):
        """
        The negative one element of this domain.

            >>> ZZ.neg_one()
            -1
            >>> Vec(ZZ,3).neg_one()
            [-1, -1, -1]
            >>> Vec(ZZ).neg_one()
            Traceback (most recent call last):
              ...
            FlintDomainError: neg_one is not an element of {Vectors (any length) over Integer ring (fmpz)}
        """
        return ctx._constant(ctx, libgr.gr_neg_one, "neg_one")

    def i(ctx):
        """
        Imaginary unit as an element of this domain.

            >>> QQbar.i()
            Root a = 1.00000*I of a^2+1
            >>> QQ.i()
            Traceback (most recent call last):
              ...
            FlintDomainError: i is not an element of {Rational field (fmpq)}
        """
        return ctx._constant(ctx, libgr.gr_i, "i")

    def pi(ctx):
        """
        The number pi as an element of this domain.

            >>> RR.pi()
            [3.141592653589793 +/- 3.39e-16]
            >>> QQbar.pi()
            Traceback (most recent call last):
              ...
            FlintDomainError: pi is not an element of {Complex algebraic numbers (qqbar)}
            >>> Vec(RR, 2).pi()
            [[3.141592653589793 +/- 3.39e-16], [3.141592653589793 +/- 3.39e-16]]

        """
        return ctx._constant(ctx, libgr.gr_pi, "pi")

    def euler(ctx):
        """
        Euler's constant as an element of this domain.

            >>> RR.euler()
            [0.5772156649015329 +/- 9.00e-17]

        We do not know whether Euler's constant is rational:

            >>> QQ.euler()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute euler in {Rational field (fmpq)}
        """
        return ctx._constant(ctx, libgr.gr_euler, "euler")

    def catalan(ctx):
        """
        Catalan's constant as an element of this domain.

            >>> RR.catalan()
            [0.915965594177219 +/- 1.23e-16]
        """
        return ctx._constant(ctx, libgr.gr_catalan, "catalan")

    def khinchin(ctx):
        """
        Khinchin's constant as an element of this domain.

            >>> RR.khinchin()
            [2.685452001065306 +/- 6.82e-16]
        """
        return ctx._constant(ctx, libgr.gr_khinchin, "khinchin")

    def glaisher(ctx):
        """
        Khinchin's constant as an element of this domain.

            >>> RR.glaisher()
            [1.282427129100623 +/- 6.02e-16]
        """
        return ctx._constant(ctx, libgr.gr_glaisher, "glaisher")

    def inv(ctx, x):
        """
        Multiplicative inverse.

            >>> ZZ.inv(2)
            Traceback (most recent call last):
              ...
            FlintDomainError: inv(x) is not an element of {Integer ring (fmpz)} for {x = 2}
            >>> QQ.inv(2)
            1/2
            >>> RR.inv(2)
            0.5000000000000000
        """
        return ctx._unary_op(x, libgr.gr_inv, "inv($x)")

    def sqrt(ctx, x):
        """
        Square root.

            >>> ZZ(25).sqrt()
            5
            >>> RR(10).sqrt()
            [3.162277660168379 +/- 5.23e-16]
        """
        return ctx._unary_op(x, libgr.gr_sqrt, "sqrt($x)")

    def rsqrt(ctx, x):
        return ctx._unary_op(x, libgr.gr_rsqrt, "rsqrt($x)")

    def numerator(ctx, x):
        return ctx._unary_op(x, libgr.gr_numerator, "numerator($x)")

    def denominator(ctx, x):
        return ctx._unary_op(x, libgr.gr_denominator, "denominator($x)")

    def floor(ctx, x):
        return ctx._unary_op(x, libgr.gr_floor, "floor($x)")

    def ceil(ctx, x):
        return ctx._unary_op(x, libgr.gr_ceil, "ceil($x)")

    def trunc(ctx, x):
        return ctx._unary_op(x, libgr.gr_trunc, "trunc($x)")

    def nint(ctx, x):
        return ctx._unary_op(x, libgr.gr_nint, "nint($x)")

    def abs(ctx, x):
        return ctx._unary_op(x, libgr.gr_abs, "abs($x)")

    def conj(ctx, x):
        return ctx._unary_op(x, libgr.gr_conj, "conj($x)")

    def re(ctx, x):
        return ctx._unary_op(x, libgr.gr_re, "re($x)")

    def im(ctx, x):
        return ctx._unary_op(x, libgr.gr_im, "im($x)")

    def sgn(ctx, x):
        return ctx._unary_op(x, libgr.gr_sgn, "sgn($x)")

    def csgn(ctx, x):
        return ctx._unary_op(x, libgr.gr_csgn, "csgn($x)")

    def arg(ctx, x):
        return ctx._unary_op(x, libgr.gr_arg, "arg($x)")

    def min(ctx, x, y):
        """
            >>> QQ.min(QQ(1)/3, QQ(1)/4)
            1/4
            >>> RR.min(3, RR.pi())
            3.000000000000000
            >>> RR.min(RR("11 +/- 1"), RR("12 +/- 3"))
            [1e+1 +/- 2.01]
            >>> CC.min(2, 3)
            2.000000000000000
            >>> CC.min(2, CC.i())
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute min(x, y) in {Complex numbers (acb, prec = 53)} for {x = 2.000000000000000}, {y = 1.000000000000000*I}

        """
        return ctx._binary_op(x, y, libgr.gr_min, "min($x, $y)")

    def max(ctx, x, y):
        """
            >>> QQ.max(QQ(1)/3, QQ(1)/4)
            1/3
            >>> RR.max(3, RR.pi())
            [3.141592653589793 +/- 5.61e-16]
            >>> RR.max(RR("10 +/- 1"), RR("9 +/- 3"))
            [1e+1 +/- 2.01]
            >>> CC.max(2, 3)
            3.000000000000000
            >>> CC.max(2, CC.i())
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute max(x, y) in {Complex numbers (acb, prec = 53)} for {x = 2.000000000000000}, {y = 1.000000000000000*I}
        """
        return ctx._binary_op(x, y, libgr.gr_max, "max($x, $y)")

    def inf(ctx):
        """
        Positive infinity (for extended number sets which support it).

            >>> CF.inf()
            +inf
            >>> ComplexExtended_ca().inf()
            +Infinity
            >>> RR.inf()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute inf in {Real numbers (arb, prec = 53)}
        """
        return ctx._constant(ctx, libgr.gr_pos_inf, "inf")

    def neg_inf(ctx):
        """
        Negative infinity (for extended number sets which support it).

            >>> CF.neg_inf()
            -inf
            >>> ComplexExtended_ca().neg_inf()
            -Infinity
            >>> RR.neg_inf()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute neg_inf in {Real numbers (arb, prec = 53)}
        """
        return ctx._constant(ctx, libgr.gr_neg_inf, "neg_inf")

    def uinf(ctx):
        """
        Unsigned infinity (for extended number sets which support it).

            >>> CF.uinf()
            Traceback (most recent call last):
              ...
            FlintDomainError: uinf is not an element of {Complex floating-point numbers (acf, prec = 53)}
            >>> ComplexExtended_ca().uinf()
            UnsignedInfinity
        """
        return ctx._constant(ctx, libgr.gr_uinf, "uinf")

    def undefined(ctx):
        """
        Undefined value (for extended number sets which support it).

            >>> CF.undefined()
            (nan + nan*I)
            >>> ComplexExtended_ca().undefined()
            Undefined
        """
        return ctx._constant(ctx, libgr.gr_undefined, "undefined")

    def unknown(ctx):
        """
        Unknown value (for enclosures and other types which support it).

            >>> ComplexExtended_ca().unknown()
            Unknown
        """
        return ctx._constant(ctx, libgr.gr_unknown, "unknown")

    def mul_2exp(ctx, x, y):
        return ctx._binary_op_fmpz(x, y, libgr.gr_mul_2exp_fmpz, "mul_2exp($x, $y)")

    def exp(ctx, x):
        """
        Exponential function.

            >>> RR.exp(1)
            [2.718281828459045 +/- 5.41e-16]
            >>> RR(1).exp()
            [2.718281828459045 +/- 5.41e-16]

        Matrix exponentials:

            >>> MatRR.exp([[1,2],[3,4]])
            [[[51.96895619870500 +/- 8.39e-15], [74.7365645670032 +/- 1.48e-14]],
            [[112.1048468505048 +/- 2.77e-14], [164.0738030492098 +/- 2.90e-14]]]
            >>> MatCC.exp([[1,2+1j],[3,4]])
            [[([44.75490138773069 +/- 7.60e-15] + [36.14044247515163 +/- 7.33e-15]*I), ([58.81526453925295 +/- 7.84e-15] + [61.26937858302805 +/- 7.52e-15]*I)],
            [([107.3399445969204 +/- 3.89e-14] + [38.23409557608189 +/- 8.58e-15]*I), ([152.0948459846510 +/- 7.30e-14] + [74.3745380512335 +/- 2.88e-14]*I)]]
            >>> Mat(CC_ca)([[1,2],[0,3]]).exp()
            [[2.71828 {a where a = 2.71828 [Exp(1)]}, 17.3673 {b^3-b where a = 20.0855 [Exp(3)], b = 2.71828 [Exp(1)]}],
            [0, 20.0855 {a where a = 20.0855 [Exp(3)]}]]
            >>> Mat(CC_ca)([[1,2],[3,4]]).exp()[0,0]
            51.9690 {(-a*c+11*a+b*c+11*b)/22 where a = 215.354 [Exp(5.37228 {(c+5)/2})], b = 0.689160 [Exp(-0.372281 {(-c+5)/2})], c = 5.74456 [c^2-33=0]}
            >>> Mat(CC_ca)([[0,0,1],[1,0,0],[0,1,0]]).exp().det()
            1
            >>> Mat(CC_ca)([[0,1,0,0,0],[0,0,2,0,0],[0,0,0,3,0],[0,0,0,0,4],[0,0,0,0,0]]).exp()
            [[1, 1, 1, 1, 1],
            [0, 1, 2, 3, 4],
            [0, 0, 1, 3, 6],
            [0, 0, 0, 1, 4],
            [0, 0, 0, 0, 1]]
            >>> MatQQ([[0,1,0,0,0],[0,0,2,0,0],[0,0,0,3,0],[0,0,0,0,4],[0,0,0,0,0]]).exp()
            [[1, 1, 1, 1, 1],
            [0, 1, 2, 3, 4],
            [0, 0, 1, 3, 6],
            [0, 0, 0, 1, 4],
            [0, 0, 0, 0, 1]]


        """
        return ctx._unary_op(x, libgr.gr_exp, "exp($x)")

    def exp2(ctx, x):
        return ctx._unary_op(x, libgr.gr_exp2, "exp2($x)")

    def exp10(ctx, x):
        return ctx._unary_op(x, libgr.gr_exp10, "exp10($x)")

    def expm1(ctx, x):
        """
            >>> RR.expm1(1)
            [1.718281828459045 +/- 3.19e-16]
        """
        return ctx._unary_op(x, libgr.gr_expm1, "expm1($x)")

    def exp_pi_i(ctx, x):
        """
            >>> QQbar.exp_pi_i(QQ(1) / 3)
            Root a = 0.500000 + 0.866025*I of a^2-a+1
            >>> CC.exp_pi_i(QQ(1) / 3)
            ([0.500000000000000 +/- 3.94e-16] + [0.866025403784439 +/- 6.79e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_exp_pi_i, "exp_pi_i($x)")

    def log(ctx, x):
        """
        Natural logarithm:

            >>> QQ.log(1)
            0
            >>> QQ.log(2)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute log(x) in {Rational field (fmpq)} for {x = 2}
            >>> RR.log(2)
            [0.693147180559945 +/- 4.12e-16]
            >>> CC.log(1j)
            [1.570796326794897 +/- 5.54e-16]*I

        Matrix logarithms:

            >>> Mat(CC_ca)([[4,2],[2,4]]).log().det() == CC_ca.log(2)*CC_ca.log(6)
            True
            >>> Mat(QQ)([[1,1],[0,1]]).log()
            [[0, 1],
            [0, 0]]
            >>> Mat(QQ)([[0,1],[0,0]]).log()
            Traceback (most recent call last):
              ...
            FlintDomainError: log(x) is not an element of {Matrices (any shape) over Rational field (fmpq)} for {x = [[0, 1],
            [0, 0]]}
            >>> Mat(CC_ca)([[0,1],[0,0]]).log()
            Traceback (most recent call last):
              ...
            FlintDomainError: log(x) is not an element of {Matrices (any shape) over Complex numbers (ca)} for {x = [[0, 1],
            [0, 0]]}
            >>> Mat(CC_ca)([[0,0,1],[0,1,0],[1,0,0]]).log() / (CC_ca.pi() * CC_ca.i())
            [[0.500000 {1/2}, 0, -0.500000 {-1/2}],
            [0, 0, 0],
            [-0.500000 {-1/2}, 0, 0.500000 {1/2}]]
            >>> Mat(CC_ca)([[0,0,1],[0,1,0],[1,0,0]]).log().exp()
            [[0, 0, 1],
            [0, 1, 0],
            [1, 0, 0]]

        """
        return ctx._unary_op(x, libgr.gr_log, "log($x)")

    def log1p(ctx, x):
        """
            >>> RR.log1p(1)
            [0.693147180559945 +/- 4.12e-16]
            >>> CC.log1p(1j)
            ([0.346573590279973 +/- 4.20e-16] + [0.7853981633974483 +/- 7.66e-17]*I)
        """
        return ctx._unary_op(x, libgr.gr_log1p, "log1p($x)")

    def log_pi_i(ctx, x):
        """
            >>> QQbar.log_pi_i(-1j)
            -1/2
            >>> CC.log_pi_i(1j)
            [0.5000000000000000 +/- 7.07e-17]

        """
        return ctx._unary_op(x, libgr.gr_log_pi_i, "log_pi_i($x)")

    def log2(ctx, x):
        """
            >>> RR.log2(16)
            [4.00000000000000 +/- 1.45e-15]
        """
        return ctx._unary_op(x, libgr.gr_log2, "log2($x)")

    def log10(ctx, x):
        """
            >>> RR.log10(100)
            [2.000000000000000 +/- 7.72e-16]
        """
        return ctx._unary_op(x, libgr.gr_log10, "log10($x)")

    def sin(ctx, x):
        """
            >>> QQ.sin(0)
            0
            >>> QQ.sin(1)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute sin(x) in {Rational field (fmpq)} for {x = 1}
            >>> RR.sin(1)
            [0.841470984807897 +/- 6.08e-16]
        """
        return ctx._unary_op(x, libgr.gr_sin, "sin($x)")

    def cos(ctx, x):
        """
            >>> QQ.cos(0)
            1
            >>> QQ.cos(1)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute cos(x) in {Rational field (fmpq)} for {x = 1}
            >>> RR.cos(1)
            [0.540302305868140 +/- 4.59e-16]
        """
        return ctx._unary_op(x, libgr.gr_cos, "cos($x)")

    def sin_cos(ctx, x):
        """
            >>> RR.sin_cos(1)
            ([0.841470984807897 +/- 6.08e-16], [0.540302305868140 +/- 4.59e-16])
        """
        return ctx._unary_unary_op(x, libgr.gr_sin_cos, "sin_cos($x)")

    def tan(ctx, x):
        """
            >>> RR.tan(1)
            [1.557407724654902 +/- 3.26e-16]
            >>> CC.tan(1+1j)
            ([0.2717525853195117 +/- 9.11e-17] + [1.083923327338694 +/- 5.77e-16]*I)
            >>> x = RRser.gen(); RRser.tan(x)
            x + [0.333333333333333 +/- 3.99e-16]*x^3 + [0.1333333333333333 +/- 8.32e-17]*x^5 + O(x^6)
        """
        return ctx._unary_op(x, libgr.gr_tan, "tan($x)")

    def cot(ctx, x):
        """
            >>> RR.cot(1)
            [0.642092615934331 +/- 4.79e-16]
        """
        return ctx._unary_op(x, libgr.gr_cot, "cot($x)")

    def sec(ctx, x):
        """
            >>> RR.sec(1)
            [1.850815717680925 +/- 7.00e-16]
        """
        return ctx._unary_op(x, libgr.gr_sec, "sec($x)")

    def csc(ctx, x):
        """
            >>> RR.csc(1)
            [1.188395105778121 +/- 2.52e-16]
        """
        return ctx._unary_op(x, libgr.gr_csc, "csc($x)")

    def sin_pi(ctx, x):
        """
            >>> QQbar.sin_pi(QQ(1) / 3)
            Root a = 0.866025 of 4*a^2-3
        """
        return ctx._unary_op(x, libgr.gr_sin_pi, "sin_pi($x)")

    def cos_pi(ctx, x):
        """
            >>> QQbar.cos_pi(QQ(1) / 3)
            1/2
        """
        return ctx._unary_op(x, libgr.gr_cos_pi, "cos_pi($x)")

    def sin_cos_pi(ctx, x):
        return ctx._unary_unary_op(x, libgr.gr_sin_cos_pi, "sin_cos_pi($x)")

    def tan_pi(ctx, x):
        """
            >>> QQbar.tan_pi(QQ(1) / 3)
            Root a = 1.73205 of a^2-3
        """
        return ctx._unary_op(x, libgr.gr_tan_pi, "tan_pi($x)")

    def cot_pi(ctx, x):
        """
            >>> QQbar.cot_pi(QQ(1) / 3)
            Root a = 0.577350 of 3*a^2-1
        """
        return ctx._unary_op(x, libgr.gr_cot_pi, "cot_pi($x)")

    def sec_pi(ctx, x):
        """
            >>> QQbar.sec_pi(QQ(1) / 3)
            2
        """
        return ctx._unary_op(x, libgr.gr_sec_pi, "sec_pi($x)")

    def csc_pi(ctx, x):
        """
            >>> QQbar.csc_pi(QQ(1) / 3)
            Root a = 1.15470 of 3*a^2-4
        """
        return ctx._unary_op(x, libgr.gr_csc_pi, "csc_pi($x)")

    def sinc(ctx, x):
        """
            >>> RR.sinc(2)
            [0.4546487134128408 +/- 7.07e-17]
            >>> CC.sinc(1j)
            [1.175201193643801 +/- 6.61e-16]
        """
        return ctx._unary_op(x, libgr.gr_sinc, "sinc($x)")

    def sinc_pi(ctx, x):
        """
            >>> RR.sinc_pi(0.5)
            [0.636619772367581 +/- 4.04e-16]
            >>> CC.sinc_pi(1j)
            [3.67607791037498 +/- 3.11e-15]
        """
        return ctx._unary_op(x, libgr.gr_sinc_pi, "sinc_pi($x)")

    def sinh(ctx, x):
        """
            >>> RR.sinh(1)
            [1.175201193643801 +/- 6.18e-16]
            >>> CC.sinh(1+1j)
            ([0.634963914784736 +/- 4.68e-16] + [1.298457581415977 +/- 7.11e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_sinh, "sinh($x)")

    def cosh(ctx, x):
        """
            >>> RR.cosh(1)
            [1.543080634815244 +/- 5.28e-16]
            >>> CC.cosh(1+1j)
            ([0.833730025131149 +/- 5.04e-16] + [0.988897705762865 +/- 4.92e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_cosh, "cosh($x)")

    def sinh_cosh(ctx, x):
        """
            >>> RR.sinh_cosh(1)
            ([1.175201193643801 +/- 6.18e-16], [1.543080634815244 +/- 5.28e-16])
            >>> CC.sinh_cosh(1j)
            ([0.841470984807897 +/- 6.08e-16]*I, [0.540302305868140 +/- 4.59e-16])
        """
        return ctx._unary_unary_op(x, libgr.gr_sinh_cosh, "sinh_cos($x)")

    def tanh(ctx, x):
        """
            >>> RR.tanh(1)
            [0.761594155955765 +/- 2.81e-16]
            >>> CC.tanh(1+1j)
            ([1.083923327338694 +/- 5.77e-16] + [0.2717525853195117 +/- 9.11e-17]*I)
        """
        return ctx._unary_op(x, libgr.gr_tanh, "tanh($x)")

    def coth(ctx, x):
        """
            >>> RR.coth(1)
            [1.313035285499331 +/- 4.97e-16]
            >>> CC.coth(1+1j)
            ([0.868014142895925 +/- 1.88e-16] + [-0.2176215618544027 +/- 9.31e-17]*I)
        """
        return ctx._unary_op(x, libgr.gr_coth, "coth($x)")

    def sech(ctx, x):
        """
            >>> RR.sech(1)
            [0.648054273663885 +/- 4.67e-16]
        """
        return ctx._unary_op(x, libgr.gr_sech, "sech($x)")

    def csch(ctx, x):
        """
            >>> RR.csch(1)
            [0.850918128239321 +/- 5.70e-16]
        """
        return ctx._unary_op(x, libgr.gr_csch, "csch($x)")

    def asin(ctx, x):
        """
            >>> RR.asin(0.5)
            [0.523598775598299 +/- 3.79e-16]
            >>> CC.asin(1+1j)
            ([0.66623943249252 +/- 5.40e-15] + [1.06127506190504 +/- 5.04e-15]*I)
            >>> x = RRser.gen(); RRser.asin(x)
            x + [0.1666666666666667 +/- 7.04e-17]*x^3 + [0.0750000000000000 +/- 1.67e-17]*x^5 + O(x^6)
        """
        return ctx._unary_op(x, libgr.gr_asin, "asin($x)")

    def acos(ctx, x):
        """
            >>> RR.acos(0.5)
            [1.047197551196598 +/- 8.97e-16]
            >>> CC.acos(1+1j)
            ([0.90455689430238 +/- 2.07e-15] + [-1.06127506190504 +/- 5.04e-15]*I)
            >>> x = RRser.gen(); RRser.acos(x)
            [1.570796326794897 +/- 5.54e-16] - x + [-0.1666666666666667 +/- 7.04e-17]*x^3 + [-0.0750000000000000 +/- 1.67e-17]*x^5 + O(x^6)
        """
        return ctx._unary_op(x, libgr.gr_acos, "acos($x)")

    def atan(ctx, x):
        """
        Inverse tangent.

            >>> QQ.atan(0)
            0
            >>> RR.atan(1)
            [0.7853981633974483 +/- 7.66e-17]
            >>> RR_ca.atan(1)
            0.785398 {(a)/4 where a = 3.14159 [Pi]}
            >>> CC.atan(0.5j)
            [0.549306144334055 +/- 3.32e-16]*I
            >>> CC.atan(1j)
            Traceback (most recent call last):
              ...
            FlintDomainError: atan(x) is not an element of {Complex numbers (acb, prec = 53)} for {x = 1.000000000000000*I}
            >>> x = RRser.gen(); RRser.atan(x)
            x + [-0.3333333333333333 +/- 7.04e-17]*x^3 + [0.2000000000000000 +/- 4.45e-17]*x^5 + O(x^6)

        """
        return ctx._unary_op(x, libgr.gr_atan, "atan($x)")

    def atan2(ctx, y, x):
        """
            >>> RR.atan2(1,2)
            [0.4636476090008061 +/- 6.22e-17]
        """
        return ctx._binary_op(y, x, libgr.gr_atan2, "atan2($y, $x)")

    def acot(ctx, x):
        """
            >>> RR.acot(2)
            [0.4636476090008061 +/- 6.22e-17]
            >>> CC.acot(1+1j)
            ([0.553574358897045 +/- 5.74e-16] + [-0.402359478108525 +/- 4.46e-16]*I)
            >>> x = RRser.gen(); RRser.acot(1+x)
            [0.7853981633974483 +/- 7.66e-17] - 0.5000000000000000*x + 0.2500000000000000*x^2 + [-0.0833333333333333 +/- 4.26e-17]*x^3 + [0.02500000000000000 +/- 5.56e-18]*x^5 + O(x^6)
        """
        return ctx._unary_op(x, libgr.gr_acot, "acot($x)")

    def asec(ctx, x):
        """
            >>> RR.asec(2)
            [1.047197551196598 +/- 8.97e-16]
            >>> CC.asec(1+1j)
            ([1.118517879643706 +/- 7.12e-16] + [0.530637530952518 +/- 6.32e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_asec, "asec($x)")

    def acsc(ctx, x):
        """
            >>> RR.acsc(2)
            [0.523598775598299 +/- 3.79e-16]
            >>> CC.acsc(1+1j)
            ([0.452278447151191 +/- 5.95e-16] + [-0.530637530952518 +/- 6.32e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_acsc, "acsc($x)")

    def asin_pi(ctx, x):
        """
            >>> QQbar.asin_pi(QQ(1) / 2)
            1/6
        """
        return ctx._unary_op(x, libgr.gr_asin_pi, "asin_pi($x)")

    def acos_pi(ctx, x):
        """
            >>> QQbar.acos_pi(QQ(1) / 2)
            1/3
        """
        return ctx._unary_op(x, libgr.gr_acos_pi, "acos_pi($x)")

    def atan_pi(ctx, x):
        """
            >>> QQbar.atan_pi(QQbar(2).sqrt() - 1)
            1/8
        """
        return ctx._unary_op(x, libgr.gr_atan_pi, "atan_pi($x)")

    def acot_pi(ctx, x):
        """
            >>> QQbar.acot_pi(QQbar(2).sqrt() - 1)
            3/8
        """
        return ctx._unary_op(x, libgr.gr_acot_pi, "acot_pi($x)")

    def asec_pi(ctx, x):
        """
            >>> QQbar.asec_pi(2)
            1/3
        """
        return ctx._unary_op(x, libgr.gr_asec_pi, "asec_pi($x)")

    def acsc_pi(ctx, x):
        """
            >>> QQbar.acsc_pi(2)
            1/6
        """
        return ctx._unary_op(x, libgr.gr_acsc_pi, "acsc_pi($x)")

    def asinh(ctx, x):
        """
            >>> RR.asinh(1)
            [0.881373587019543 +/- 1.87e-16]
            >>> CC.asinh(1+1j)
            ([1.06127506190504 +/- 5.04e-15] + [0.66623943249252 +/- 5.40e-15]*I)
            >>> x = RRser.gen(); RRser.asinh(x)
            x + [-0.1666666666666667 +/- 7.04e-17]*x^3 + [0.0750000000000000 +/- 1.67e-17]*x^5 + O(x^6)
        """
        return ctx._unary_op(x, libgr.gr_asinh, "asinh($x)")

    def acosh(ctx, x):
        """
            >>> RR.acosh(2)
            [1.316957896924817 +/- 6.61e-16]
            >>> CC.acosh(1+1j)
            ([1.061275061905035 +/- 8.44e-16] + [0.904556894302381 +/- 8.22e-16]*I)
            >>> x = RRser.gen(); RRser.acosh(2+x)
            [1.316957896924817 +/- 6.61e-16] + [0.577350269189626 +/- 4.54e-16]*x + [-0.192450089729875 +/- 4.24e-16]*x^2 + [0.096225044864938 +/- 6.92e-16]*x^3 + [-0.058804194084128 +/- 7.16e-16]*x^4 + [0.040450157748779 +/- 4.91e-16]*x^5 + O(x^6)
        """
        return ctx._unary_op(x, libgr.gr_acosh, "acosh($x)")

    def atanh(ctx, x):
        """
            >>> RR.atanh(0.5)
            [0.549306144334055 +/- 3.32e-16]
            >>> CC.atanh(1+1j)
            ([0.4023594781085251 +/- 8.52e-17] + [1.017221967897851 +/- 4.37e-16]*I)
            >>> x = RRser.gen(); RRser.atanh(x)
            x + [0.3333333333333333 +/- 7.04e-17]*x^3 + [0.2000000000000000 +/- 4.45e-17]*x^5 + O(x^6)

        """
        return ctx._unary_op(x, libgr.gr_atanh, "atanh($x)")

    def acoth(ctx, x):
        """
            >>> RR.acoth(2)
            [0.549306144334055 +/- 3.32e-16]
            >>> CC.acoth(1+1j)
            ([0.4023594781085251 +/- 8.52e-17] + [-0.553574358897045 +/- 3.16e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_acoth, "acoth($x)")

    def asech(ctx, x):
        """
            >>> RR.asech(0.5)
            [1.316957896924817 +/- 6.61e-16]
            >>> CC.asech(1+1j)
            ([0.530637530952518 +/- 9.50e-16] + [-1.118517879643706 +/- 9.45e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_asech, "asech($x)")

    def acsch(ctx, x):
        """
            >>> RR.acsch(0.5)
            [1.443635475178810 +/- 5.32e-16]
            >>> CC.acsch(1+1j)
            ([0.530637530952518 +/- 5.66e-16] + [-0.452278447151191 +/- 7.14e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_acsch, "acsch($x)")

    def erf(ctx, x):
        """
            >>> RR.erf(1)
            [0.842700792949715 +/- 3.28e-16]
        """
        return ctx._unary_op(x, libgr.gr_erf, "erf($x)")

    def erfc(ctx, x):
        """
            >>> RR.erfc(1)
            [0.1572992070502851 +/- 3.71e-17]
        """
        return ctx._unary_op(x, libgr.gr_erfc, "erfc($x)")

    def erfi(ctx, x):
        """
            >>> RR.erfi(1)
            [1.650425758797543 +/- 4.58e-16]
        """
        return ctx._unary_op(x, libgr.gr_erfi, "erfi($x)")

    def erfcx(ctx, x):
        """
            >>> RR.erfcx(100)
            [0.00564161378298943 +/- 3.69e-18]
        """
        return ctx._unary_op(x, libgr.gr_erfcx, "erfcx($x)")

    def erfinv(ctx, x):
        """
            >>> RR.erfinv(0.5)
            [0.4769362762044698 +/- 7.79e-17]
        """
        return ctx._unary_op(x, libgr.gr_erfinv, "erfinv($x)")

    def erfcinv(ctx, x):
        """
            >>> RR.erfc(RR.erfcinv(0.25))
            [0.250000000000000 +/- 1.24e-16]
        """
        return ctx._unary_op(x, libgr.gr_erfcinv, "erfcinv($x)")

    def fresnel_s(ctx, x, normalized=False):
        """
            >>> RR.fresnel_s(1)
            [0.3102683017233811 +/- 2.67e-18]
            >>> RR.fresnel_s(1, normalized=True)
            [0.4382591473903548 +/- 9.24e-17]
        """
        return ctx._unary_op_with_flag(x, normalized, libgr.gr_fresnel_s, "fresnel_s($x)")

    def fresnel_c(ctx, x, normalized=False):
        """
            >>> RR.fresnel_c(1)
            [0.904524237900272 +/- 1.46e-16]
            >>> RR.fresnel_c(1, normalized=True)
            [0.779893400376823 +/- 3.59e-16]
        """
        return ctx._unary_op_with_flag(x, normalized, libgr.gr_fresnel_c, "fresnel_c($x)")

    def fresnel(ctx, x, normalized=False):
        """
            >>> RR.fresnel(1)
            ([0.3102683017233811 +/- 2.67e-18], [0.904524237900272 +/- 1.46e-16])
            >>> RR.fresnel(1, normalized=True)
            ([0.4382591473903548 +/- 9.24e-17], [0.779893400376823 +/- 3.59e-16])
        """
        return ctx._unary_unary_op_with_flag(x, normalized, libgr.gr_fresnel, "fresnel($x)")

    def gamma_upper(ctx, x, y, regularized=0):
        """
            >>> RR.gamma_upper(3, 4)
            [0.476206611107089 +/- 5.30e-16]
            >>> RR.gamma_upper(3, 4, regularized=True)
            [0.2381033055535443 +/- 8.24e-17]
        """
        return ctx._binary_op_with_flag(x, y, regularized, libgr.gr_gamma_upper, "gamma_upper($x, $y)")

    def gamma_lower(ctx, x, y, regularized=0):
        """
            >>> RR.gamma_lower(3, 4)
            [1.52379338889291 +/- 2.89e-15]
            >>> RR.gamma_lower(3, 4, regularized=True)
            [0.76189669444646 +/- 6.52e-15]
        """
        return ctx._binary_op_with_flag(x, y, regularized, libgr.gr_gamma_lower, "gamma_lower($x, $y)")

    def beta_lower(ctx, a, b, x, regularized=0):
        """
            >>> RR.beta_lower(2, 3, 0.5)
            [0.0572916666666667 +/- 6.08e-17]
            >>> RR.beta_lower(2, 3, 0.5, regularized=True)
            [0.687500000000000 +/- 5.24e-16]
        """
        return ctx._ternary_op_with_flag(a, b, x, regularized, libgr.gr_beta_lower, "beta_lower($a, $b, $x)")

    def exp_integral(ctx, x, y):
        """
            >>> RR.exp_integral(1, 2)
            [0.04890051070806 +/- 2.63e-15]
        """
        return ctx._binary_op(x, y, libgr.gr_exp_integral, "exp_integral($x, $y)")

    def exp_integral_ei(ctx, x):
        """
            >>> RR.exp_integral_ei(1)
            [1.89511781635594 +/- 5.11e-15]
        """
        return ctx._unary_op(x, libgr.gr_exp_integral_ei, "exp_integral_ei($x)")

    def sin_integral(ctx, x):
        """
            >>> RR.sin_integral(1)
            [0.946083070367183 +/- 1.35e-16]
        """
        return ctx._unary_op(x, libgr.gr_sin_integral, "sin_integral($x)")

    def cos_integral(ctx, x):
        """
            >>> RR.cos_integral(1)
            [0.3374039229009681 +/- 5.63e-17]
        """
        return ctx._unary_op(x, libgr.gr_cos_integral, "cos_integral($x)")

    def sinh_integral(ctx, x):
        """
            >>> RR.sinh_integral(1)
            [1.05725087537573 +/- 2.77e-15]
        """
        return ctx._unary_op(x, libgr.gr_sinh_integral, "sinh_integral($x)")

    def cosh_integral(ctx, x):
        """
            >>> RR.cosh_integral(1)
            [0.837866940980208 +/- 4.78e-16]
        """
        return ctx._unary_op(x, libgr.gr_cosh_integral, "cosh_integral($x)")

    def log_integral(ctx, x, offset=False):
        """
            >>> RR.log_integral(2)
            [1.04516378011749 +/- 4.01e-15]
            >>> RR.log_integral(2, offset=True)
            0
        """
        return ctx._unary_op_with_flag(x, offset, libgr.gr_log_integral, "log_integral($x)")

    def dilog(ctx, x):
        """
            >>> RR.dilog(1)
            [1.644934066848226 +/- 6.45e-16]
        """
        return ctx._unary_op(x, libgr.gr_dilog, "dilog($x)")

    def bessel_j(ctx, x, y):
        """
            >>> RR.bessel_j(2, 3)
            [0.486091260585891 +/- 4.75e-16]
        """
        return ctx._binary_op(x, y, libgr.gr_bessel_j, "bessel_j($n, $x)")

    def bessel_y(ctx, x, y):
        """
            >>> RR.bessel_y(2, 3)
            [-0.16040039348492 +/- 5.80e-15]
        """
        return ctx._binary_op(x, y, libgr.gr_bessel_y, "bessel_y($n, $x)")

    def bessel_i(ctx, x, y, scaled=False):
        """
            >>> RR.bessel_i(2, 3)
            [2.24521244092995 +/- 1.88e-15]
            >>> RR.bessel_i(2, 3, scaled=True)
            [0.111782545296958 +/- 2.09e-16]
        """
        if scaled:
            return ctx._binary_op(x, y, libgr.gr_bessel_i_scaled, "bessel_i($n, $x, scaled=True)")
        else:
            return ctx._binary_op(x, y, libgr.gr_bessel_i, "bessel_i($n, $x)")

    def bessel_k(ctx, x, y, scaled=False):
        """
            >>> RR.bessel_k(2, 3)
            [0.06151045847174 +/- 8.87e-15]
            >>> RR.bessel_k(2, 3, scaled=True)
            [1.235470584796 +/- 5.14e-13]
        """
        if scaled:
            return ctx._binary_op(x, y, libgr.gr_bessel_k_scaled, "bessel_k($n, $x, scaled=True)")
        else:
            return ctx._binary_op(x, y, libgr.gr_bessel_k, "bessel_k($n, $x)")

    def bessel_j_y(ctx, x, y):
        """
            >>> RR.bessel_j_y(1, 1)
            ([0.4400505857449335 +/- 5.91e-17], [-0.78121282130029 +/- 4.55e-15])
            >>> CC.bessel_j_y(1, 1j)
            ([0.565159103992485 +/- 1.89e-16]*I, ([-0.56515910399248 +/- 8.81e-15] + [0.38318604387456 +/- 8.19e-15]*I))
        """
        return ctx._binary_binary_op(x, y, libgr.gr_bessel_j_y, "bessel_j_y($n, $x)")

    def airy(ctx, x):
        """
            >>> RR.airy(1)
            ([0.1352924163128814 +/- 4.17e-17], [-0.1591474412967932 +/- 2.95e-17], [1.207423594952871 +/- 3.27e-16], [0.932435933392776 +/- 5.83e-16])
            >>> CC.airy(1)
            ([0.1352924163128814 +/- 4.17e-17], [-0.1591474412967932 +/- 2.95e-17], [1.207423594952871 +/- 3.27e-16], [0.932435933392776 +/- 5.83e-16])
        """
        return ctx._quaternary_unary_op(x, libgr.gr_airy, "airy($x)")

    def airy_ai(ctx, x):
        """
            >>> RR.airy_ai(1)
            [0.1352924163128814 +/- 4.17e-17]
        """
        return ctx._unary_op(x, libgr.gr_airy_ai, "airy_ai($x)")

    def airy_bi(ctx, x):
        """
            >>> RR.airy_bi(1)
            [1.207423594952871 +/- 3.27e-16]
        """
        return ctx._unary_op(x, libgr.gr_airy_bi, "airy_bi($x)")


    def airy_ai_prime(ctx, x):
        """
            >>> RR.airy_ai_prime(1)
            [-0.1591474412967932 +/- 2.95e-17]
        """
        return ctx._unary_op(x, libgr.gr_airy_ai_prime, "airy_ai_prime($x)")

    def airy_bi_prime(ctx, x):
        """
            >>> RR.airy_bi_prime(1)
            [0.932435933392776 +/- 5.83e-16]
        """
        return ctx._unary_op(x, libgr.gr_airy_bi_prime, "airy_bi_prime($x)")

    def airy_ai_zero(ctx, n):
        """
            >>> RR.airy_ai(RR.airy_ai_zero(1))
            [+/- 3.51e-16]
        """
        return ctx._unary_op_fmpz(n, libgr.gr_airy_ai_zero, "airy_ai_zero($n)")

    def airy_bi_zero(ctx, n):
        """
            >>> RR.airy_bi(RR.airy_bi_zero(1))
            [+/- 2.08e-16]
        """
        return ctx._unary_op_fmpz(n, libgr.gr_airy_bi_zero, "airy_bi_zero($n)")

    def airy_ai_prime_zero(ctx, n):
        """
            >>> RR.airy_ai_prime(RR.airy_ai_prime_zero(1))
            [+/- 1.44e-16]
        """
        return ctx._unary_op_fmpz(n, libgr.gr_airy_ai_prime_zero, "airy_ai_prime_zero($n)")

    def airy_bi_prime_zero(ctx, n):
        """
            >>> RR.airy_bi_prime(RR.airy_bi_prime_zero(1))
            [+/- 6.18e-16]
        """
        return ctx._unary_op_fmpz(n, libgr.gr_airy_bi_prime_zero, "airy_bi_prime_zero($n)")

    # todo: coulomb()

    def coulomb_f(ctx, x, y, z):
        """
            >>> CC.coulomb_f(2, 3, 4)
            [0.101631502833431 +/- 8.03e-16]
        """
        return ctx._ternary_op(x, y, z, libgr.gr_coulomb_f, "coulomb_f($x)")

    def coulomb_g(ctx, x, y, z):
        """
            >>> CC.coulomb_g(2, 3, 4)
            [5.371722466 +/- 6.15e-10]
        """
        return ctx._ternary_op(x, y, z, libgr.gr_coulomb_g, "coulomb_g($x)")

    def coulomb_hpos(ctx, x, y, z):
        """
            >>> CC.coulomb_hpos(2, 3, 4)
            ([5.371722466 +/- 6.15e-10] + [0.101631502833431 +/- 8.03e-16]*I)
        """
        return ctx._ternary_op(x, y, z, libgr.gr_coulomb_hpos, "coulomb_hpos($x)")

    def coulomb_hneg(ctx, x, y, z):
        """
            >>> CC.coulomb_hneg(2, 3, 4)
            ([5.371722466 +/- 6.15e-10] + [-0.101631502833431 +/- 8.03e-16]*I)
        """
        return ctx._ternary_op(x, y, z, libgr.gr_coulomb_hneg, "coulomb_hneg($x)")

    def chebyshev_t(ctx, n, x):
        """
        Chebyshev polynomial of the first kind.

            >>> [ZZ.chebyshev_t(n, 2) for n in range(5)]
            [1, 2, 7, 26, 97]
            >>> RR.chebyshev_t(0.5, 0.75)
            [0.935414346693485 +/- 5.18e-16]
            >>> ZZx.chebyshev_t(4, [0, 1])
            8*x^4-8*x^2+1
        """
        return ctx._binary_op_with_overloads(n, x, libgr.gr_chebyshev_t, fmpz_op=libgr.gr_chebyshev_t_fmpz, rstr="chebyshev_t($n, $x)")

    def chebyshev_u(ctx, n, x):
        """
        Chebyshev polynomial of the second kind.

            >>> [ZZ.chebyshev_u(n, 2) for n in range(5)]
            [1, 4, 15, 56, 209]
            >>> RR.chebyshev_u(0.5, 0.75)
            [1.33630620956212 +/- 2.68e-15]
            >>> ZZx.chebyshev_u(4, [0, 1])
            16*x^4-12*x^2+1
        """
        return ctx._binary_op_with_overloads(n, x, libgr.gr_chebyshev_u, fmpz_op=libgr.gr_chebyshev_u_fmpz, rstr="chebyshev_u($n, $x)")

    def jacobi_p(ctx, n, a, b, x):
        """
        Jacobi polynomial.

            >>> RR.jacobi_p(3, 1, 2, 4)
            [602.500000000000 +/- 3.28e-13]
        """
        return ctx._quaternary_op(n, a, b, x, libgr.gr_jacobi_p, rstr="jacobi_p($n, $a, $b, $x)")

    def gegenbauer_c(ctx, n, m, x):
        """
        Gegenbauer polynomial.

            >>> RR.gegenbauer_c(3, 2, 4)
            [2000.00000000000 +/- 3.60e-12]
        """
        return ctx._ternary_op(n, m, x, libgr.gr_gegenbauer_c, rstr="gegenbauer_c($n, $m, $x)")

    def laguerre_l(ctx, n, m, x):
        """
        Associated Laguerre polynomial (or Laguerre function).

            >>> RR.laguerre_l(3, 2, 4)
            [-0.66666666666667 +/- 5.71e-15]
        """
        return ctx._ternary_op(n, m, x, libgr.gr_laguerre_l, rstr="laguerre_l($n, $m, $x)")

    def hermite_h(ctx, n, x):
        """
        Hermite polynomial (Hermite function).

            >>> RR.hermite_h(3, 4)
            464.0000000000000
        """
        return ctx._binary_op(n, x, libgr.gr_hermite_h, rstr="hermite_h($n, $x)")

    def legendre_p(ctx, n, m, x, typ=0):
        """
        Associated Legendre function of the first kind.

            >>> RR.legendre_p(3, 1, 0.5)
            [-0.324759526419164 +/- 5.23e-16]
            >>> CC.legendre_p(3, 1, 0.5)
            [-0.324759526419164 +/- 5.23e-16]
            >>> CC.legendre_p(3, 1, 0.5, 1)
            [0.324759526419164 +/- 5.52e-16]*I
        """
        assert typ in (0, 1)
        return ctx._ternary_op_with_flag(n, m, x, typ, libgr.gr_legendre_p, rstr="legendre_p($n, $m, $x, $typ)")

    def legendre_q(ctx, n, m, x, typ=0):
        """
        Associated Legendre function of the second kind.

            >>> RR.legendre_q(3, 1, 0.5)
            [2.4918525917090 +/- 5.45e-14]
            >>> CC.legendre_q(3, 1, 0.5)
            [2.4918525917090 +/- 5.45e-14]
            >>> CC.legendre_q(3, 1, 0.5, 1)
            ([0.51013107119087 +/- 4.02e-15] + [-2.4918525917090 +/- 5.73e-14]*I)
        """
        assert typ in (0, 1)
        return ctx._ternary_op_with_flag(n, m, x, typ, libgr.gr_legendre_q, rstr="legendre_q($n, $m, $x, $typ)")

    def spherical_y(ctx, n, m, theta, phi):
        """
        Spherical harmonic.

            >>> CC.spherical_y(4, 3, 0.5, 0.75)
            ([0.076036396941350 +/- 2.18e-16] + [-0.094180781089734 +/- 4.96e-16]*I)
        """
        n = ctx._as_si(n)
        m = ctx._as_si(m)
        theta = ctx._as_elem(theta)
        phi = ctx._as_elem(phi)
        res = ctx._elem_type(context=ctx)
        status = libgr.gr_spherical_y_si(res._ref, n, m, theta._ref, phi._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "spherical_y($n, $m, $theta, $phi)", n, m, theta, phi)
        return res

    def legendre_p_root(ctx, n, k, weight=False):
        """
        Root of Legendre polynomial.
        With weight=True, also returns the corresponding weight for
        Gauss-Legendre quadrature.

            >>> RR.legendre_p_root(5, 1)
            [0.538469310105683 +/- 1.15e-16]
            >>> RR.legendre_p(5, 0, RR.legendre_p_root(5, 1))
            [+/- 8.15e-16]
            >>> RR.legendre_p_root(5, 1, weight=True)
            ([0.538469310105683 +/- 1.15e-16], [0.4786286704993664 +/- 7.10e-17])
        """
        n = ctx._as_si(n)
        k = ctx._as_si(k)
        if weight:
            res1 = ctx._elem_type(context=ctx)
            res2 = ctx._elem_type(context=ctx)
            status = libgr.gr_legendre_p_root_ui(res1._ref, res2._ref, n, k, ctx._ref)
        else:
            res1 = ctx._elem_type(context=ctx)
            res2 = None
            status = libgr.gr_legendre_p_root_ui(res1._ref, res2, n, k, ctx._ref)
        if status:
            _handle_error(ctx, status, "legendre_p_root($n, $k)", n, k)
        if weight:
            return (res1, res2)
        else:
            return res1

    def hypgeom_0f1(ctx, a, z, regularized=False):
        """
        Hypergeometric function 0F1, optionally regularized.

            >>> RR.hypgeom_0f1(3, 4)
            [3.21109468764205 +/- 5.00e-15]
            >>> RR.hypgeom_0f1(3, 4, regularized=True)
            [1.60554734382103 +/- 5.20e-15]
            >>> CC.hypgeom_0f1(1, 2+2j)
            ([2.435598449671389 +/- 7.27e-16] + [4.43452765355337 +/- 4.91e-15]*I)
        """
        flags = int(regularized)
        return ctx._binary_op_with_flag(a, z, flags, libgr.gr_hypgeom_0f1, rstr="hypgeom_0f1($a, $x)")

    def hypgeom_1f1(ctx, a, b, z, regularized=False):
        """
        Hypergeometric function 1F1, optionally regularized.

            >>> RR.hypgeom_1f1(3, 4, 5)
            [60.504568913851 +/- 3.82e-13]
            >>> RR.hypgeom_1f1(3, 4, 5, regularized=True)
            [10.0840948189752 +/- 3.31e-14]
        """
        flags = int(regularized)
        return ctx._ternary_op_with_flag(a, b, z, flags, libgr.gr_hypgeom_1f1, rstr="hypgeom_1f1($a, $b, $x)")

    def hypgeom_u(ctx, a, b, z):
        """
        Hypergeometric function U.

            >>> RR.hypgeom_u(1, 2, 3)
            [0.3333333333333333 +/- 7.04e-17]
        """
        flags = 0
        return ctx._ternary_op_with_flag(a, b, z, flags, libgr.gr_hypgeom_u, rstr="hypgeom_u($a, $b, $x)")

    def hypgeom_2f1(ctx, a, b, c, z, regularized=False):
        """
        Hypergeometric function 2F1, optionally regularized.

            >>> RR.hypgeom_2f1(1, 2, 3, -4)
            [0.29882026094574 +/- 8.48e-15]
            >>> RR.hypgeom_2f1(1, 2, 3, -4, regularized=True)
            [0.14941013047287 +/- 4.24e-15]
        """
        flags = int(regularized)
        return ctx._quaternary_op_with_flag(a, b, c, z, flags, libgr.gr_hypgeom_2f1, rstr="hypgeom_2f1($a, $b, $c, $x)")

    def hypgeom_pfq(ctx, a, b, z, regularized=False):
        """
        Generalized hypergeometric function, optionally regularized.

            >>> RR.hypgeom_pfq([1,2], [3,4], 0.5)
            [1.09002619782383 +/- 4.32e-15]
            >>> RR.hypgeom_pfq([1, 2], [3, 4], 0.5, regularized=True)
            [0.090835516485319 +/- 2.36e-16]
            >>> CC.hypgeom_pfq([1,2], [3,4], 0.5+0.5j)
            ([1.08239550393928 +/- 2.16e-15] + [0.096660812453003 +/- 5.55e-16]*I)
            >>> x = CCser.gen()
            >>> CCser.hypgeom_pfq([1, 2+x], [3+x], CC(-0.5)+x)
            [0.75627913513468 +/- 8.06e-15] + [0.32390309682858 +/- 8.68e-15]*x + [0.2403726319128 +/- 9.25e-14]*x^2 + [0.11701345954 +/- 2.11e-12]*x^3 + [0.0806045925 +/- 6.82e-11]*x^4 + [0.043967479 +/- 9.09e-10]*x^5 + O(x^6)
            >>> CCser.hypgeom_pfq([2+x], [3+x], CC(-0.5)+x)
            [0.72163208344840 +/- 3.62e-15] + [0.41823982602421 +/- 4.49e-15]*x + [0.24476273640228 +/- 2.45e-15]*x^2 + [0.0532046593049 +/- 2.33e-14]*x^3 + [0.021808666094 +/- 4.11e-13]*x^4 + [0.000742930119 +/- 5.65e-13]*x^5 + O(x^6)

        """
        a = ctx._as_vec(a)
        b = ctx._as_vec(b)
        z = ctx._as_elem(z)
        res = ctx._elem_type(context=ctx)
        flags = int(regularized)
        status = libgr.gr_hypgeom_pfq(res._ref, a._ref, b._ref, z._ref, flags, ctx._ref)
        if status:
            _handle_error(ctx, status, "hypgeom_pfq($a, $b, $x)", a, b, z)
        return res

    def fac(ctx, x):
        """
        Factorial.

            >>> ZZ.fac(10)
            3628800
            >>> ZZ.fac(-1)
            Traceback (most recent call last):
              ...
            FlintDomainError: fac(x) is not an element of {Integer ring (fmpz)} for {x = -1}

        Real and complex factorials extend using the gamma function:

            >>> RR.fac(10**20)
            [1.93284951431010e+1956570551809674817245 +/- 3.03e+1956570551809674817230]
            >>> RR.fac(0.5)
            [0.886226925452758 +/- 1.78e-16]
            >>> CC.fac(1+1j)
            ([0.652965496420167 +/- 6.21e-16] + [0.343065839816545 +/- 5.38e-16]*I)

        Factorials mod N:

            >>> ZZmod(10**7 + 19).fac(10**7)
            2343096

        More tests:

            >>> RF.fac(10**6)
            8.263931688331239e+5565708
            >>> RF.fac(10**20)
            1.932849514310098e+1956570551809674817245

        """
        return ctx._unary_op_with_fmpz_fmpq_overloads(x, libgr.gr_fac, op_fmpz=libgr.gr_fac_fmpz, rstr="fac($x)")

    def fac_vec(ctx, length):
        """
        Vector of factorials.

            >>> ZZ.fac_vec(10)
            [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880]
            >>> QQ.fac_vec(10) / 3
            [1/3, 1/3, 2/3, 2, 8, 40, 240, 1680, 13440, 120960]
            >>> ZZmod(7).fac_vec(10)
            [1, 1, 2, 6, 3, 1, 6, 0, 0, 0]
            >>> sum(RR.fac_vec(100))
            [9.427862397658e+155 +/- 3.19e+142]

        """
        return ctx._op_vec_len(length, libgr.gr_fac_vec, "fac_vec($length)")

    def rfac(ctx, x):
        """
        Reciprocal factorial.

            >>> QQ.rfac(5)
            1/120
            >>> ZZ.rfac(-2)
            0
            >>> ZZ.rfac(2)
            Traceback (most recent call last):
              ...
            FlintDomainError: rfac(x) is not an element of {Integer ring (fmpz)} for {x = 2}
            >>> RR.rfac(0.5)
            [1.128379167095513 +/- 7.02e-16]

        """
        return ctx._unary_op_with_fmpz_fmpq_overloads(x, libgr.gr_rfac, op_fmpz=libgr.gr_rfac_fmpz, rstr="rfac($x)")

    def rfac_vec(ctx, length):
        """
        Vector of reciprocal factorials.

            >>> QQ.rfac_vec(8)
            [1, 1, 1/2, 1/6, 1/24, 1/120, 1/720, 1/5040]
            >>> ZZmod(7).rfac_vec(7)
            [1, 1, 4, 6, 5, 1, 6]
            >>> ZZmod(7).rfac_vec(8)
            Traceback (most recent call last):
              ...
            FlintDomainError: rfac_vec(length) is not an element of {Integers mod 7 (_gr_nmod)} for {length = 8}
            >>> sum(RR.rfac_vec(20))
            [2.71828182845904 +/- 8.66e-15]
        """
        return ctx._op_vec_len(length, libgr.gr_rfac_vec, "rfac_vec($length)")

    def rising(ctx, x, n):
        """
        Rising factorial.

            >>> [ZZ.rising(3, k) for k in range(5)]
            [1, 3, 12, 60, 360]
            >>> ZZx.rising(ZZx([0,1]), 5)
            x^5+10*x^4+35*x^3+50*x^2+24*x
            >>> RR.rising(1, 10**7)
            [1.202423400515903e+65657059 +/- 5.57e+65657043]
        """
        return ctx._binary_op_with_overloads(x, n, libgr.gr_rising, op_ui=libgr.gr_rising_ui, rstr="rising($x, $n)")

    def falling(ctx, x, n):
        """
        Falling factorial.

            >>> [ZZ.falling(3, k) for k in range(5)]
            [1, 3, 6, 6, 0]
            >>> ZZx.falling(ZZx([0,1]), 5)
            x^5-10*x^4+35*x^3-50*x^2+24*x
            >>> RR.log(RR.falling(RR.pi(), 10**7))
            [151180898.7174084 +/- 9.72e-8]
            >>> RR.falling(10.5, 3.5)
            [2360.99664364330 +/- 4.00e-12]

        """
        return ctx._binary_op_with_overloads(x, n, libgr.gr_falling, op_ui=libgr.gr_falling_ui, rstr="falling($x, $n)")

    def bin(ctx, x, y):
        """
        Binomial coefficient.

            >>> [ZZ.bin(5, k) for k in range(7)]
            [1, 5, 10, 10, 5, 1, 0]
            >>> RR.bin(100000, 50000)
            [2.52060836892200e+30100 +/- 5.36e+30085]
            >>> ZZmod(1000).bin(10000, 3000)
            200
            >>> ZZp64.bin(100000, 50000)
            5763493550349629692
            >>> ZZp64.bin(10**30, 2)
            998763921924463582
            >>> RR.bin(1.5, 0.75)
            [1.57378746535479 +/- 5.62e-15]

        """
        try:
            x = ctx._as_ui(x)
            y = ctx._as_ui(y)
            return ctx._op_uiui(x, y, libgr.gr_bin_uiui, "bin($x, $y)")
        except:
            return ctx._binary_op_with_overloads(x, y, libgr.gr_bin, op_ui=libgr.gr_bin_ui, rstr="bin($x, $y)")

    def bin_vec(ctx, n, length=None):
        """
        Vector of binomial coefficients, optionally truncated to specified length.

            >>> ZZ.bin_vec(8)
            [1, 8, 28, 56, 70, 56, 28, 8, 1]
            >>> ZZmod(5).bin_vec(8)
            [1, 3, 3, 1, 0, 1, 3, 3, 1]
            >>> ZZ.bin_vec(0)
            [1]
            >>> ZZ.bin_vec(1000, 3)
            [1, 1000, 499500]
            >>> ZZ.bin_vec(4, 8)
            [1, 4, 6, 4, 1, 0, 0, 0]
            >>> QQ.bin_vec(QQ(1)/2, 5)
            [1, 1/2, -1/8, 1/16, -5/128]
            >>> ZZmod(7).bin_vec(10)
            [1, 3, 3, 1, 0, 0, 0, 1, 3, 3, 1]
            >>> ZZmod(7).bin_vec(3)
            [1, 3, 3, 1]
            >>> ZZmod(7).bin_vec(10, 1)
            [1]
        """
        try:
            n = ctx._as_ui(n)
        except:
            return ctx._op_vec_arg_len(n, length, libgr.gr_bin_vec, "bin_vec($n, $length)")
        if length is None:
            length = n + 1
        return ctx._op_vec_ui_len(n, length, libgr.gr_bin_ui_vec, "bin_vec($n, $length)")

    def gamma(ctx, x):
        """
            >>> RR.gamma(10)
            362880.0000000000
            >>> RR.gamma(0.5)
            [1.772453850905516 +/- 3.41e-16]
            >>> RR.gamma(QQ(1) / 3)
            [2.678938534707747 +/- 8.99e-16]
            >>> CC.gamma(1+1j) / CC.gamma(1j)
            ([+/- 6.32e-16] + [1.0000000000000 +/- 1.03e-15]*I)
        """
        return ctx._unary_op_with_fmpz_fmpq_overloads(x, libgr.gr_gamma, op_fmpz=libgr.gr_gamma_fmpz, op_fmpq=libgr.gr_gamma_fmpq, rstr="gamma($x)")

    def lgamma(ctx, x):
        """
            >>> RR.lgamma(10)
            [12.80182748008147 +/- 2.69e-15]
            >>> CC.lgamma(10j)
            ([-15.94031728124131 +/- 6.90e-15] + [12.23211664743500 +/- 4.89e-15]*I)
        """
        return ctx._unary_op(x, libgr.gr_lgamma, "lgamma($x)")

    def rgamma(ctx, x):
        """
            >>> RR.rgamma(10)
            [2.755731922398589e-6 +/- 5.96e-22]
            >>> CC.rgamma(10+1j)
            ([-1.83246026966323e-6 +/- 5.08e-21] + [-2.25314671311995e-6 +/- 5.78e-21]*I)
        """
        return ctx._unary_op(x, libgr.gr_rgamma, "lgamma($x)")

    def digamma(ctx, x):
        """
            >>> RR.digamma(2)
            [0.4227843350984671 +/- 4.84e-17]
            >>> CC.digamma(2j)
            ([0.714591515373977 +/- 6.06e-16] + [1.820807282642230 +/- 3.65e-16]*I)
        """
        return ctx._unary_op(x, libgr.gr_digamma, "digamma($x)")

    def doublefac(ctx, x):
        """
        Double factorial (semifactorial).

            >>> [ZZ.doublefac(n) for n in range(10)]
            [1, 1, 2, 3, 8, 15, 48, 105, 384, 945]
            >>> RR.doublefac(2.5)
            [2.40706945611604 +/- 5.54e-15]
            >>> CC.doublefac(1+1j)
            ([0.250650779545753 +/- 7.56e-16] + [0.100474421235437 +/- 4.14e-16]*I)
        """
        return ctx._unary_op_with_fmpz_fmpq_overloads(x, libgr.gr_doublefac, op_ui=libgr.gr_doublefac_ui, rstr="doublefac($x)")

    def harmonic(ctx, x):
        """
        Harmonic numbers.

            >>> [QQ.harmonic(n) for n in range(6)]
            [0, 1, 3/2, 11/6, 25/12, 137/60]
            >>> RR.harmonic(10**9)
            [21.30048150234794 +/- 8.48e-15]
            >>> ZZp64.harmonic(1000)
            6514760847963681162
            >>> RR.harmonic(10.5)
            [2.97545479443731 +/- 5.16e-15]
            >>> RR.harmonic(15092688622113788323693563264538101449859497)
            [100.000000000000 +/- 4.35e-14]
        """
        return ctx._unary_op_with_fmpz_fmpq_overloads(x, libgr.gr_harmonic, op_ui=libgr.gr_harmonic_ui, rstr="harmonic($x)")

    def beta(ctx, x, y):
        """
        Beta function.

            >>> RR.beta(3, 4.5)
            [0.01243201243201243 +/- 6.93e-18]
            >>> CC.beta(1j, 1+1j)
            ([-1.18807306241087 +/- 5.32e-15] + [-1.31978426013907 +/- 4.09e-15]*I)
        """
        return ctx._binary_op(y, x, libgr.gr_beta, "beta($x, $y)")

    def barnes_g(ctx, x):
        """
        Barnes G-function.

            >>> RR.barnes_g(7)
            34560.00000000000
            >>> CC.barnes_g(1+2j)
            ([0.54596949228965 +/- 7.69e-15] + [-3.98421873125106 +/- 8.76e-15]*I)
        """
        return ctx._unary_op(x, libgr.gr_barnes_g, "barnes_g($x)")

    def log_barnes_g(ctx, x):
        """
        Logarithmic Barnes G-function.

            >>> RR.log_barnes_g(100)
            [15258.0613921488 +/- 3.87e-11]
            >>> CC.log_barnes_g(10+20j)
            ([-452.057343313397 +/- 6.85e-13] + [121.014356688943 +/- 2.52e-13]*I)
        """
        return ctx._unary_op(x, libgr.gr_log_barnes_g, "log_barnes_g($x)")

    def zeta(ctx, s):
        """
        Riemann zeta function.

            >>> RR.zeta(2)
            [1.644934066848226 +/- 4.57e-16]
            >>> CC.zeta(1+1j)
            ([0.5821580597520036 +/- 5.17e-17] + [-0.9268485643308071 +/- 2.75e-17]*I)
        """
        return ctx._unary_op(s, libgr.gr_zeta, "zeta($s)")

    def hurwitz_zeta(ctx, s, a):
        """
        Hurwitz zeta function.

            >>> RR.hurwitz_zeta(2, 2)
            [0.6449340668482264 +/- 3.72e-17]
            >>> CC.hurwitz_zeta(1j, 1)
            ([0.0033002236853241 +/- 2.42e-17] + [-0.4181554491413217 +/- 4.51e-17]*I)
        """
        return ctx._binary_op(s, a, libgr.gr_hurwitz_zeta, "hurwitz_zeta($s, $a)")

    def stieltjes(ctx, n, a=1):
        """
        Stieltjes constant.

            >>> CC.stieltjes(1)
            [-0.0728158454836767 +/- 2.78e-17]
            >>> CC.stieltjes(1, a=0.5)
            [-1.353459680804942 +/- 7.22e-16]
        """
        n = ctx._as_fmpz(n)
        a = ctx._as_elem(a)
        res = ctx._elem_type(context=ctx)
        libgr.gr_stieltjes.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p)
        status = libgr.gr_stieltjes(res._ref, n._ref, a._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "stieltjes($n, $a)", n, a)
        return res

    def polylog(ctx, s, z):
        """
        Polylogarithm.

            >>> CC.polylog(2, -1)
            [-0.822467033424113 +/- 3.22e-16]
        """
        return ctx._binary_op(s, z, libgr.gr_polylog, "polylog($s, $z)")

    def polygamma(ctx, s, z):
        """
        Polygamma function.

            >>> CC.polygamma(2, 3)
            [-0.1541138063191886 +/- 7.16e-17]
        """
        return ctx._binary_op(s, z, libgr.gr_polygamma, "polygamma($s, $z)")

    def lerch_phi(ctx, z, s, a):
        """
            >>> CC.lerch_phi(2, 3, 4)
            ([-0.00213902437921 +/- 1.70e-15] + [-0.04716836434127 +/- 5.28e-15]*I)
        """
        return ctx._ternary_op(z, s, a, libgr.gr_lerch_phi, "lerch_phi($z, $s, $a)")

    def dirichlet_eta(ctx, x):
        """
        Dirichlet eta function.

            >>> CC.dirichlet_eta(1)
            [0.6931471805599453 +/- 6.93e-17]
            >>> CC.dirichlet_eta(2)
            [0.822467033424113 +/- 2.36e-16]
        """
        return ctx._unary_op(x, libgr.gr_dirichlet_eta, "dirichlet_eta($x)")

    def riemann_xi(ctx, x):
        """
        Riemann xi function.

            >>> s = 2+3j; CC.riemann_xi(s); CC.riemann_xi(1-s)
            ([0.41627125989962 +/- 4.65e-15] + [0.08882330496564 +/- 1.43e-15]*I)
            ([0.41627125989962 +/- 4.65e-15] + [0.08882330496564 +/- 1.43e-15]*I)
        """
        return ctx._unary_op(x, libgr.gr_riemann_xi, "riemann_xi($x)")

    def lambertw(ctx, x, k=None):
        """
            >>> RR.lambertw(1)
            [0.567143290409784 +/- 2.72e-16]
            >>> RR.lambertw(-0.25)
            [-0.3574029561813889 +/- 5.91e-17]
            >>> RR.lambertw(-0.25, -1)
            [-2.153292364110349 +/- 8.59e-16]
            >>> CC.lambertw(-1)
            ([-0.318131505204764 +/- 1.92e-16] + [1.337235701430689 +/- 5.99e-16]*I)
            >>> CC.lambertw(1, 5)
            ([-3.398692196764719 +/- 6.76e-16] + [29.73131070782852 +/- 7.03e-15]*I)
        """
        if k is None:
            return ctx._unary_op(x, libgr.gr_lambertw, "lambertw($x)")
        else:
            return ctx._binary_op_fmpz(x, k, libgr.gr_lambertw_fmpz, "lambertw($x, $k)")

    def bernoulli(ctx, n):
        """
        Bernoulli number `B_n` as an element of this domain.

            >>> QQ.bernoulli(10)
            5/66
            >>> RR.bernoulli(10)
            [0.0757575757575757 +/- 5.97e-17]

            >>> ZZ.bernoulli(0)
            1
            >>> ZZ.bernoulli(1)
            Traceback (most recent call last):
              ...
            FlintDomainError: bernoulli(n) is not an element of {Integer ring (fmpz)} for {n = 1}

        Huge Bernoulli numbers can be computed numerically:

            >>> RR.bernoulli(10**20)
            [-1.220421181609039e+1876752564973863312289 +/- 4.69e+1876752564973863312273]
            >>> RF.bernoulli(10**20)
            -1.220421181609039e+1876752564973863312289
            >>> QQ.bernoulli(10**20)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute bernoulli(n) in {Rational field (fmpq)} for {n = 100000000000000000000}

        """
        return ctx._op_fmpz(n, libgr.gr_bernoulli_fmpz, "bernoulli($n)")

    def bernoulli_vec(ctx, length):
        """
        Vector of Bernoulli numbers.

            >>> QQ.bernoulli_vec(12)
            [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66, 0]
            >>> CC_ca.bernoulli_vec(5)
            [1, -0.500000 {-1/2}, 0.166667 {1/6}, 0, -0.0333333 {-1/30}]
            >>> sum(RR.bernoulli_vec(100))
            [1.127124216595034e+76 +/- 6.74e+60]
            >>> sum(RF.bernoulli_vec(100))
            1.127124216595034e+76
            >>> sum(CC.bernoulli_vec(100))
            [1.127124216595034e+76 +/- 6.74e+60]

        """
        return ctx._op_vec_len(length, libgr.gr_bernoulli_vec, "bernoulli_vec($length)")

    def eulernum(ctx, n):
        """
        Euler number `E_n` as an element of this domain.

            >>> ZZ.eulernum(10)
            -50521
            >>> RR.eulernum(10)
            -50521.00000000000

        Huge Euler numbers can be computed numerically:

            >>> RR.eulernum(10**20)
            [4.346791453661149e+1936958564106659551331 +/- 8.35e+1936958564106659551315]
            >>> RF.eulernum(10**20)
            4.346791453661149e+1936958564106659551331
            >>> ZZ.eulernum(10**20)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute eulernum(n) in {Integer ring (fmpz)} for {n = 100000000000000000000}

        """
        return ctx._op_fmpz(n, libgr.gr_eulernum_fmpz, "eulernum($n)")

    def eulernum_vec(ctx, length):
        """
        Vector of Euler numbers.

            >>> ZZ.eulernum_vec(12)
            [1, 0, -1, 0, 5, 0, -61, 0, 1385, 0, -50521, 0]
            >>> QQ.eulernum_vec(12) / 3
            [1/3, 0, -1/3, 0, 5/3, 0, -61/3, 0, 1385/3, 0, -50521/3, 0]
            >>> sum(RR.eulernum_vec(100))
            [-7.23465655613392e+134 +/- 3.20e+119]
            >>> sum(RF.eulernum_vec(100))
            -7.234656556133921e+134
        """
        return ctx._op_vec_len(length, libgr.gr_eulernum_vec, "eulernum_vec($length)")

    def fib(ctx, n):
        """
        Fibonacci number `F_n` as an element of this domain.

            >>> ZZ.fib(10)
            55
            >>> RR.fib(10)
            55.00000000000000
            >>> ZZ.fib(-10)
            -55

        Huge Fibonacci numbers can be computed numerically and in modular arithmetic:

            >>> RR.fib(10**20)
            [3.78202087472056e+20898764024997873376 +/- 4.02e+20898764024997873361]
            >>> RF.fib(10**20)
            3.782020874720557e+20898764024997873376
            >>> F = FiniteField_fq(17, 1)
            >>> n = 10**20; F.fib(n); F.fib(n-1) + F.fib(n-2)
            13
            13

        """
        return ctx._op_fmpz(n, libgr.gr_fib_fmpz, "fib($n)")

    def fib_vec(ctx, length):
        """
        Vector of Fibonacci numbers.

            >>> ZZ.fib_vec(10)
            [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
            >>> QQ.fib_vec(10) / 3
            [0, 1/3, 1/3, 2/3, 1, 5/3, 8/3, 13/3, 7, 34/3]
            >>> sum(RR.fib_vec(100))            # doctest: +ELLIPSIS
            [5.7314784401...e+20 +/- ...]
            >>> sum(RF.fib_vec(100))
            5.731478440138172e+20
        """
        return ctx._op_vec_len(length, libgr.gr_fib_vec, "fib($length)")

    def stirling_s1u(ctx, n, k):
        """
        Unsigned Stirling number of the first kind.

            >>> ZZ.stirling_s1u(5, 2)
            50
            >>> QQ.stirling_s1u(5, 2)
            50
            >>> ZZ.stirling_s1u(50, 21)
            33187391298039120738041153829116024033357291261862000
            >>> RR.stirling_s1u(50, 21)
            [3.318739129803912e+52 +/- 8.66e+36]
        """
        return ctx._op_uiui(n, k, libgr.gr_stirling_s1u_uiui, "stirling_s1u($n, $k)")

    def stirling_s1(ctx, n, k):
        """
        Signed Stirling number of the first kind.

            >>> ZZ.stirling_s1(5, 2)
            -50
            >>> QQ.stirling_s1(5, 2)
            -50
            >>> RR.stirling_s1(5, 2)
            -50.00000000000000
        """
        return ctx._op_uiui(n, k, libgr.gr_stirling_s1_uiui, "stirling_s1($n, $k)")

    def stirling_s2(ctx, n, k):
        """
        Stirling number of the second kind.

            >>> ZZ.stirling_s2(5, 2)
            15
            >>> QQ.stirling_s2(5, 2)
            15
            >>> RR.stirling_s2(5, 2)
            15.00000000000000
            >>> RR.stirling_s2(50, 20)
            [7.59792160686099e+45 +/- 5.27e+30]
        """
        return ctx._op_uiui(n, k, libgr.gr_stirling_s2_uiui, "stirling_s2($n, $k)")

    def stirling_s1u_vec(ctx, n, length=None):
        """
        Vector of unsigned Stirling numbers of the first kind,
        optionally truncated to specified length.

            >>> ZZ.stirling_s1u_vec(5)
            [0, 24, 50, 35, 10, 1]
            >>> QQ.stirling_s1u_vec(5) / 3
            [0, 8, 50/3, 35/3, 10/3, 1/3]
            >>> RR.stirling_s1u_vec(5, 3)
            [0, 24.00000000000000, 50.00000000000000]
        """
        if length is None:
            length = n + 1
        return ctx._op_vec_ui_len(n, length, libgr.gr_stirling_s1u_ui_vec, "stirling_s1u_vec($n, $length)")

    def stirling_s1_vec(ctx, n, length=None):
        """
        Vector of signed Stirling numbers of the first kind,
        optionally truncated to specified length.

            >>> ZZ.stirling_s1_vec(5)
            [0, 24, -50, 35, -10, 1]
            >>> QQ.stirling_s1_vec(5) / 3
            [0, 8, -50/3, 35/3, -10/3, 1/3]
            >>> RR.stirling_s1_vec(5, 3)
            [0, 24.00000000000000, -50.00000000000000]
        """
        if length is None:
            length = n + 1
        return ctx._op_vec_ui_len(n, length, libgr.gr_stirling_s1_ui_vec, "stirling_s1_vec($n, $length)")

    def stirling_s2_vec(ctx, n, length=None):
        """
        Vector of Stirling numbers of the second kind,
        optionally truncated to specified length.

            >>> ZZ.stirling_s2_vec(5)
            [0, 1, 15, 25, 10, 1]
            >>> QQ.stirling_s2_vec(5) / 3
            [0, 1/3, 5, 25/3, 10/3, 1/3]
            >>> RR.stirling_s2_vec(5, 3)
            [0, 1, 15.00000000000000]
        """
        if length is None:
            length = n + 1
        return ctx._op_vec_ui_len(n, length, libgr.gr_stirling_s2_ui_vec, "stirling_s2_vec($n, $length)")

    def bellnum(ctx, n):
        """
        Bell number `E_n` as an element of this domain.

            >>> ZZ.bellnum(10)
            115975
            >>> RR.bellnum(10)
            115975.0000000000
            >>> ZZp64.bellnum(10000)
            355901145009109239
            >>> ZZmod(1000).bellnum(10000)
            635

        Huge Bell numbers can be computed numerically:

            >>> RR.bellnum(10**20)
            [5.38270113176282e+1794956117137290721328 +/- 5.44e+1794956117137290721313]
            >>> ZZ.bellnum(10**20)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute bellnum(n) in {Integer ring (fmpz)} for {n = 100000000000000000000}
        """
        return ctx._op_fmpz(n, libgr.gr_bellnum_fmpz, "bellnum($n)")

    def bellnum_vec(ctx, length):
        """
        Vector of Bell numbers.

            >>> ZZ.bellnum_vec(10)
            [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]
            >>> QQ.bellnum_vec(10) / 3
            [1/3, 1/3, 2/3, 5/3, 5, 52/3, 203/3, 877/3, 1380, 7049]
            >>> RR.bellnum_vec(100).sum()
            [1.67618752079292e+114 +/- 4.30e+99]
            >>> RF.bellnum_vec(100).sum()
            1.676187520792924e+114
            >>> ZZmod(10000).bellnum_vec(10000).sum()
            7337

        """
        return ctx._op_vec_len(length, libgr.gr_bellnum_vec, "bellnum_vec($length)")

    def partitions(ctx, n):
        """
        Partition function `p(n)` as an element of this domain.

            >>> ZZ.partitions(10)
            42
            >>> QQ.partitions(10) / 5
            42/5
            >>> RR.partitions(10)
            42.00000000000000
            >>> RR.partitions(10**20)
            [1.838176508344883e+11140086259 +/- 8.18e+11140086243]
        """
        return ctx._op_fmpz(n, libgr.gr_partitions_fmpz, "partitions($n)")

    def partitions_vec(ctx, length):
        """
        Vector of partition numbers.

            >>> ZZ.partitions_vec(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
            >>> QQ.partitions_vec(10) / 3
            [1/3, 1/3, 2/3, 1, 5/3, 7/3, 11/3, 5, 22/3, 10]
            >>> ZZmod(10).partitions_vec(10)
            [1, 1, 2, 3, 5, 7, 1, 5, 2, 0]
            >>> sum(ZZmod(10).partitions_vec(100))
            6
            >>> sum(RR.partitions_vec(100))
            1452423276.000000
        """
        return ctx._op_vec_len(length, libgr.gr_partitions_vec, "partitions($length)")

    def zeta_zero(ctx, n):
        """
        Zero of the Riemann zeta function.

            >>> CC.zeta_zero(1)
            (0.5000000000000000 + [14.13472514173469 +/- 4.71e-15]*I)
            >>> CC.zeta_zero(2)
            (0.5000000000000000 + [21.02203963877155 +/- 6.02e-15]*I)
        """
        return ctx._unary_op_fmpz(n, libgr.gr_zeta_zero, "zeta_zero($n)")

    def zeta_zeros(ctx, num, start=1):
        """
        Zeros of the Riemann zeta function.

            >>> [x.im() for x in CC.zeta_zeros(4)]
            [[14.13472514173469 +/- 4.71e-15], [21.02203963877155 +/- 6.02e-15], [25.01085758014569 +/- 7.84e-15], [30.42487612585951 +/- 5.96e-15]]
            >>> [x.im() for x in CC.zeta_zeros(2, start=100)]
            [[236.5242296658162 +/- 3.51e-14], [237.7698204809252 +/- 5.29e-14]]
        """
        return ctx._op_vec_fmpz_len(start, num, libgr.gr_zeta_zero_vec, "zeta_zeros($n)")

    def zeta_nzeros(ctx, t):
        """
        Number of zeros of Riemann zeta function up to given height.

            >>> CC.zeta_nzeros(100)
            29.00000000000000
        """
        return ctx._unary_op(t, libgr.gr_zeta_nzeros, "zeta_nzeros($t)")

    def dirichlet_l(ctx, s, chi):
        """
        Dirichlet L-function with character chi.

            >>> CC.dirichlet_l(2, DirichletGroup(1)(1))
            [1.644934066848226 +/- 4.57e-16]
            >>> RR.dirichlet_l(2, DirichletGroup(4)(3))
            [0.915965594177219 +/- 2.68e-16]
            >>> CC.dirichlet_l(2+3j, DirichletGroup(7)(3))
            ([1.273313649440491 +/- 9.69e-16] + [-0.074323294425594 +/- 6.96e-16]*I)
        """
        s = ctx._as_elem(s)
        assert isinstance(chi, dirichlet_char)
        res = ctx._elem_type(context=ctx)
        libgr.gr_dirichlet_l.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p)
        G = libgr.gr_ctx_data_as_ptr(chi.parent()._ref)
        status = libgr.gr_dirichlet_l(res._ref, G, chi._ref, s._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "dirichlet_l($s, $chi)", s, chi)
        return res

    def hardy_theta(ctx, s, chi=None):
        """
        Hardy theta function.

            >>> CC.hardy_theta(10)
            [-3.06707439628989 +/- 6.66e-15]
            >>> RR.hardy_theta(2)
            [-2.525910918816132 +/- 9.34e-16]
            >>> CC.hardy_theta(10, DirichletGroup(4)(3))
            [4.64979557270698 +/- 4.41e-15]
        """
        s = ctx._as_elem(s)
        res = ctx._elem_type(context=ctx)
        libgr.gr_dirichlet_hardy_theta.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p)
        if chi is None:
            chi_ref = G = None
        else:
            assert isinstance(chi, dirichlet_char)
            G = libgr.gr_ctx_data_as_ptr(chi.parent()._ref)
            chi_ref = chi._ref
        status = libgr.gr_dirichlet_hardy_theta(res._ref, G, chi_ref, s._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "hardy_theta($s, $chi)", s, chi)
        return res

    def hardy_z(ctx, s, chi=None):
        """
        Hardy Z-function.

            >>> CC.hardy_z(2)
            [-0.539633125646145 +/- 8.59e-16]
            >>> RR.hardy_z(2)
            [-0.539633125646145 +/- 8.59e-16]
            >>> CC.hardy_z(2, DirichletGroup(4)(3))
            [1.15107760668266 +/- 5.01e-15]
        """
        s = ctx._as_elem(s)
        res = ctx._elem_type(context=ctx)
        libgr.gr_dirichlet_hardy_z.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p)
        if chi is None:
            chi_ref = G = None
        else:
            assert isinstance(chi, dirichlet_char)
            G = libgr.gr_ctx_data_as_ptr(chi.parent()._ref)
            chi_ref = chi._ref
        status = libgr.gr_dirichlet_hardy_z(res._ref, G, chi_ref, s._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "hardy_z($s, $chi)", s, chi)
        return res

    def dirichlet_chi(ctx, n, chi):
        """
        Value of the Dirichlet character chi(n).

            >>> chi = DirichletGroup(5)(3)
            >>> [CC.dirichlet_chi(n, chi) for n in range(5)]
            [0, 1, -1.000000000000000*I, 1.000000000000000*I, -1]
        """
        n = ctx._as_fmpz(n)
        assert isinstance(chi, dirichlet_char)
        res = ctx._elem_type(context=ctx)
        libgr.gr_dirichlet_chi_fmpz.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p)
        G = libgr.gr_ctx_data_as_ptr(chi.parent()._ref)
        status = libgr.gr_dirichlet_chi_fmpz(res._ref, G, chi._ref, n._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "dirichlet_chi($n, $chi)", n, chi)
        return res

    def dirichlet_chi_vec(ctx, chi, n):
        """
        Vector of values of the given Dirichlet character.

            >>> CC.dirichlet_chi_vec(DirichletGroup(4)(3), 5)
            [0, 1, 0, -1, 0]
        """
        n = ctx._as_si(n)
        assert n >= 0
        assert n <= HUGE_LENGTH
        assert isinstance(chi, dirichlet_char)
        libgr.gr_dirichlet_chi_vec.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.c_void_p)
        G = libgr.gr_ctx_data_as_ptr(chi.parent()._ref)
        res = Vec(ctx)()
        assert not libgr.gr_vec_set_length(res._ref, n, ctx._ref)
        status = libgr.gr_dirichlet_chi_vec(libgr.gr_vec_entry_ptr(res._ref, 0, ctx._ref), G, chi._ref, n, ctx._ref)
        if status:
            _handle_error(ctx, status, "dirichlet_chi_vec($chi, $n)", chi, n)
        return res

    def modular_j(ctx, tau):
        """
        j-invariant j(tau).

            >>> CC.modular_j(1j)
            [1728.0000000000 +/- 5.10e-11]
        """
        return ctx._unary_op(tau, libgr.gr_modular_j, "modular_j($tau)")

    def modular_lambda(ctx, tau):
        """
        Modular lambda function lambda(tau).

            >>> CC.modular_lambda(1j)
            [0.50000000000000 +/- 2.16e-15]
        """
        return ctx._unary_op(tau, libgr.gr_modular_lambda, "modular_lambda($tau)")

    def modular_delta(ctx, tau):
        """
        Modular discriminant delta(tau).

            >>> CC.modular_delta(1j)
            [0.0017853698506421 +/- 6.01e-17]
        """
        return ctx._unary_op(tau, libgr.gr_modular_delta, "modular_delta($tau)")

    def dedekind_eta(ctx, tau):
        """
        Dedekind eta function eta(tau).

            >>> CC.dedekind_eta(1j)
            [0.768225422326057 +/- 9.03e-16]
        """
        return ctx._unary_op(tau, libgr.gr_dedekind_eta, "dedekind_eta($tau)")

    def hilbert_class_poly(ctx, D, x):
        """
        Hilbert class polynomial H_D(x) evaluated at x.

            >>> ZZx.hilbert_class_poly(-20, ZZx.gen())
            x^2-1264000*x-681472000
            >>> CC.hilbert_class_poly(-20, 1+1j)
            (-682736000.0000000 - 1263998.000000000*I)
            >>> ZZx.hilbert_class_poly(-21, ZZx.gen())
            Traceback (most recent call last):
              ...
            FlintDomainError: hilbert_class_poly(D, x) is not an element of {Polynomials over integers (fmpz_poly)} for {D = -21}, {x = x}
        """
        D = ctx._as_si(D)
        x = ctx._as_elem(x)
        res = ctx._elem_type(context=ctx)
        libgr.gr_hilbert_class_poly.argtypes = (ctypes.c_void_p, c_slong, ctypes.c_void_p, ctypes.c_void_p)
        status = libgr.gr_hilbert_class_poly(res._ref, D, x._ref, ctx._ref)
        if status:
            _handle_error(ctx, status, "hilbert_class_poly($D, $x)", D, x)
        return res

    def eisenstein_g(ctx, n, tau):
        """
        Eisenstein series G_n(tau).

            >>> CC.eisenstein_g(2, 1j)
            [3.14159265358979 +/- 8.71e-15]
            >>> CC.eisenstein_g(4, 1j); RR.gamma(0.25)**8 / (960 * RR.pi()**2)
            [3.1512120021539 +/- 3.41e-14]
            [3.15121200215390 +/- 7.72e-15]

        """
        return ctx._ui_binary_op(n, tau, libgr.gr_eisenstein_g, "eisenstein_g($n, $tau)")

    def eisenstein_e(ctx, n, tau):
        """
        Eisenstein series E_n(tau).

            >>> CC.eisenstein_e(2, 1j)
            [0.95492965855137 +/- 3.85e-15]
            >>> CC.eisenstein_e(4, 1j); 3*RR.gamma(0.25)**8/(64*RR.pi()**6)
            [1.4557628922687 +/- 1.32e-14]
            [1.45576289226871 +/- 3.76e-15]

        """
        return ctx._ui_binary_op(n, tau, libgr.gr_eisenstein_e, "eisenstein_e($n, $tau)")

    def eisenstein_g_vec(ctx, tau, n):
        """
        Vector of Eisenstein series [G_4(tau), G_6(tau), ...].
        Note that G_2(tau) is omitted.

            >>> CC.eisenstein_g_vec(1j, 3)
            [[3.1512120021539 +/- 3.41e-14], [+/- 4.40e-14], [4.255773035365 +/- 2.12e-13]]
        """
        return ctx._op_vec_arg_len(tau, n, libgr.gr_eisenstein_g_vec, "eisenstein_g_vec($tau, $n)")

    def agm(ctx, x, y=None):
        """
        Arithmetic-geometric mean.

            >>> RR.agm(2)
            [1.45679103104691 +/- 3.98e-15]
            >>> RR.agm(2, 3)
            [2.47468043623630 +/- 4.68e-15]
        """
        if y is None:
            return ctx._unary_op(x, libgr.gr_agm1, "agm1($x)")
        else:
            return ctx._binary_op(x, y, libgr.gr_agm, "agm($x, $y)")

    def elliptic_k(ctx, m):
        return ctx._unary_op(m, libgr.gr_elliptic_k, "elliptic_k($m)")

    def elliptic_e(ctx, m):
        return ctx._unary_op(m, libgr.gr_elliptic_e, "elliptic_e($m)")

    def elliptic_pi(ctx, n, m):
        return ctx._binary_op(n, m, libgr.gr_elliptic_pi, "elliptic_pi($n, $m)")

    def elliptic_f(ctx, phi, m, pi=0):
        return ctx._binary_op_with_flag(phi, m, pi, libgr.gr_elliptic_f, "elliptic_f($phi, $m, $pi)")

    def elliptic_e_inc(ctx, phi, m, pi=0):
        return ctx._binary_op_with_flag(phi, m, pi, libgr.gr_elliptic_e_inc, "elliptic_e_inc($phi, $m, $pi)")

    def elliptic_pi_inc(ctx, n, phi, m, pi=0):
        return ctx._ternary_op_with_flag(n, phi, m, pi, libgr.gr_elliptic_pi_inc, "elliptic_pi_inc($n, $phi, $m, $pi)")

    def carlson_rc(ctx, x, y, flags=0):
        return ctx._binary_op_with_flag(x, y, flags, libgr.gr_carlson_rc, "carlson_rc($x, $y)")

    def carlson_rf(ctx, x, y, z, flags=0):
        return ctx._ternary_op_with_flag(x, y, z, flags, libgr.gr_carlson_rf, "carlson_rf($x, $y, $z)")

    def carlson_rg(ctx, x, y, z, flags=0):
        return ctx._ternary_op_with_flag(x, y, z, flags, libgr.gr_carlson_rg, "carlson_rg($x, $y, $z)")

    def carlson_rd(ctx, x, y, z, flags=0):
        return ctx._ternary_op_with_flag(x, y, z, flags, libgr.gr_carlson_rd, "carlson_rd($x, $y, $z)")

    def carlson_rj(ctx, x, y, z, w, flags=0):
        return ctx._quaternary_op_with_flag(x, y, z, w, flags, libgr.gr_carlson_rd, "carlson_rj($x, $y, $z, $w)")

    def jacobi_theta(ctx, z, tau):
        """
        Simultaneous computation of the four Jacobi theta functions.

            >>> CC.jacobi_theta(0.125, 1j)
            ([0.347386687929454 +/- 3.21e-16], [0.843115469091413 +/- 8.18e-16], [1.061113709291166 +/- 5.74e-16], [0.938886290708834 +/- 3.52e-16])
        """
        return ctx._quaternary_binary_op(z, tau, libgr.gr_jacobi_theta, "jacobi_theta($z, $tau)")

    def jacobi_theta_1(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_1(0.125, 1j)
            [0.347386687929454 +/- 3.21e-16]
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_1, "jacobi_theta_1($z, $tau)")

    def jacobi_theta_2(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_2(0.125, 1j)
            [0.843115469091413 +/- 8.18e-16]
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_2, "jacobi_theta_2($z, $tau)")

    def jacobi_theta_3(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_3(0.125, 1j)
            [1.061113709291166 +/- 5.74e-16]
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_3, "jacobi_theta_3($z, $tau)")

    def jacobi_theta_4(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_4(0.125, 1j)
            [0.938886290708834 +/- 3.52e-16]
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_4, "jacobi_theta_4($z, $tau)")

    def elliptic_invariants(ctx, tau):
        """
            >>> g2, g3 = CC.elliptic_invariants(1j)
            >>> CC.weierstrass_p_prime(0.25, 1j)**2; 4*CC.weierstrass_p(0.25, 1j)**3 - g2*CC.weierstrass_p(0.25, 1j) - g3
            [15152.862386715 +/- 7.03e-10]
            [15152.86238672 +/- 5.09e-9]
        """
        return ctx._unary_unary_op(tau, libgr.gr_elliptic_invariants, "elliptic_invariants($tau)")

    def elliptic_roots(ctx, tau):
        """
            >>> e1, e2, e3 = CC.elliptic_roots(1j)
            >>> g2, g3 = CC.elliptic_invariants(1j)
            >>> 4*e1**3 - g2*e1 - g3
            [+/- 3.12e-11]
            >>> 4*e2**3 - g2*e2 - g3
            [+/- 8.29e-12]
            >>> 4*e3**3 - g2*e3 - g3
            [+/- 3.14e-11]
        """
        return ctx._ternary_unary_op(tau, libgr.gr_elliptic_roots, "elliptic_roots($tau)")

    def weierstrass_p(ctx, z, tau):
        return ctx._binary_op(z, tau, libgr.gr_weierstrass_p, "weierstrass_p($z, $tau)")

    def weierstrass_p_prime(ctx, z, tau):
        return ctx._binary_op(z, tau, libgr.gr_weierstrass_p_prime, "weierstrass_p_prime($z, $tau)")

    def weierstrass_p_inv(ctx, z, tau):
        """
        Inverse Weierstrass elliptic function.

            >>> CC.weierstrass_p(CC.weierstrass_p_inv(0.5, 1j), 1j)
            ([0.50000000000 +/- 4.61e-12] + [+/- 6.98e-12]*I)
        """
        return ctx._binary_op(z, tau, libgr.gr_weierstrass_p_inv, "weierstrass_p_inv($z, $tau)")

    def weierstrass_zeta(ctx, z, tau):
        return ctx._binary_op(z, tau, libgr.gr_weierstrass_zeta, "weierstrass_zeta($z, $tau)")

    def weierstrass_sigma(ctx, z, tau):
        return ctx._binary_op(z, tau, libgr.gr_weierstrass_sigma, "weierstrass_sigma($z, $tau)")


def _gr_set_int(self, val):
    if WORD_MIN <= val <= WORD_MAX:
        status = libgr.gr_set_si(self._ref, val, self._ctx)
    else:
        n = fmpz_struct()
        nref = ctypes.byref(n)
        libflint.fmpz_init(nref)
        libflint.fmpz_set_str(nref, ctypes.c_char_p(str(val).encode('ascii')), 10)
        status = libgr.gr_set_fmpz(self._ref, nref, self._ctx)
        libflint.fmpz_clear(nref)
    return status

class gr_elem:
    """
    Base class for elements.
    """

    @staticmethod
    def _default_context():
        return None

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

    def __init__(self, val=None, context=None, random=False):
        """
            >>> ZZ(QQ(1))
            1
            >>> ZZ(QQ(1) / 3)
            Traceback (most recent call last):
              ...
            FlintDomainError: 1/3 is not defined in Integer ring (fmpz)
        """
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
                status = _gr_set_int(self, val)
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
            elif typ is fexpr:
                # XXX
                xvec = yvec = Vec(context)()
                status = libgr.gr_set_fexpr(self._ref, xvec._ref, yvec._ref, val._ref, self._ctx)
            elif hasattr(val, "_gr_elem_"):
                val = val._gr_elem_(context)
                assert val.parent() is context
                status = libgr.gr_set_other(self._ref, val._ref, val._ctx, self._ctx)
            elif typ.__name__ == "mpz":
                status = _gr_set_int(self, int(val))
            else:
                status = GR_UNABLE
            if status:
                if status & GR_UNABLE: raise FlintUnableError(f"unable to create element of {self.parent()} from {val} of type {type(val)}")
                if status & GR_DOMAIN: raise FlintDomainError(f"{val} is not defined in {self.parent()}")
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

    def nstr(self, n):
        """
        Return a string representation of this element, where
        real and complex numbers may be rounded to n digits.

            >>> RR.pi().nstr(10)
            '3.141592654'
            >>> CC(1+1j).exp().nstr(10)
            '(1.468693940 + 2.287355287*I)'
        """
        arr = ctypes.c_char_p()
        n = self._ctx_python._as_si(n)
        if libgr.gr_get_str_n(ctypes.byref(arr), self._ref, n, self._ctx) != GR_SUCCESS:
            raise NotImplementedError
        try:
            return ctypes.cast(arr, ctypes.c_char_p).value.decode("ascii")
        finally:
            libflint.flint_free(arr)

    def nprint(self, n):
        """
        Print a string representation of this element, where
        real and complex numbers may be rounded to *n* digits.

            >>> RR.pi().nprint(10)
            3.141592654
            >>> CC(1+1j).exp().nprint(10)
            (1.468693940 + 2.287355287*I)
        """
        print(self.nstr(n))

    def fexpr(self, serialize=True):
        res = fexpr()
        if serialize:
            status = libflint.gr_get_fexpr_serialize(res._ref, self._ref, self._ctx)
        else:
            status = libflint.gr_get_fexpr(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "fexpr($x)")
        return res

    def latex(self):
        return self.fexpr(serialize=False).latex()

    def _repr_latex_(self):
        return "$$" + self.latex() + "$$"

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
            _handle_error(self.parent(), status, rstr, self, other)
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
            _handle_error(self.parent(), status, rstr, self, other)
            # if status & GR_UNABLE: raise NotImplementedError(f"unable to compute {rstr} for x = {self}, y = {other} over {self.parent()}")
            # if status & GR_DOMAIN: raise ValueError(f"{rstr} is not defined for x = {self}, y = {other} over {self.parent()}")
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
            _handle_error(self.parent(), status, rstr, self)
        return res

    @staticmethod
    def _unary_op_get_fmpz(self, op, rstr):
        res = ZZ()
        status = op(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, rstr, self)
        return res

    @staticmethod
    def _binary_op_fmpz(self, other, op, rstr):
        other = ZZ(other)
        elem_type = type(self)
        res = elem_type(context=self._ctx_python)
        status = op(res._ref, self._ref, other._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, rstr, self, other)
        return res

    @staticmethod
    def _constant(self, op, rstr):
        elem_type = type(self)
        res = elem_type(context=self._ctx_python)
        status = op(res._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, rstr)
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
        return self._binary_op2(self, other, _add_methods, "$x + $y")

    def __radd__(self, other):
        return self._binary_op2(other, self, _add_methods, "$x + $y")

    def __sub__(self, other):
        return self._binary_op2(self, other, _sub_methods, "$x - $y")

    def __rsub__(self, other):
        return self._binary_op2(other, self, _sub_methods, "$x - $y")

    def __mul__(self, other):
        return self._binary_op2(self, other, _mul_methods, "$x * $y")

    def __rmul__(self, other):
        return self._binary_op2(other, self, _mul_methods, "$x * $y")

    def __truediv__(self, other):
        return self._binary_op2(self, other, _div_methods, "$x / $y")

    def __rtruediv__(self, other):
        return self._binary_op2(other, self, _div_methods, "$x / $y")

    def __pow__(self, other):
        return self._binary_op2(self, other, _pow_methods, "$x ** $y")

    def __rpow__(self, other):
        return self._binary_op2(other, self, _pow_methods, "$x ** $y")

    def __floordiv__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_div, "$x // $y")

    def __rfloordiv__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_div, "$x // $y")

    def __mod__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_rem, "$x % $y")

    def __rmod__(self, other):
        return self._binary_op(self, other, libgr.gr_euclidean_rem, "$x % $y")

    def is_zero(self):
        return self._unary_predicate(self, libgr.gr_is_zero, "is_zero")

    def is_one(self):
        return self._unary_predicate(self, libgr.gr_is_one, "is_one")

    def is_neg_one(self):
        return self._unary_predicate(self, libgr.gr_is_neg_one, "is_neg_one")

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
            >>> pi = CC_ca.pi(); i = CC_ca.i(); x = PolynomialRing(CC_ca).gen(); (x**2 + pi**2).gcd(x+i*pi)
            (3.14159*I {a*b where a = 3.14159 [Pi], b = I [b^2+1=0]}) + x
            >>> QQx([1,1,2,-1,3]).gcd(QQx([1,-1,1]))
            1 - x + x^2
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
        exponents = VecZZ()
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
            >>> QQbar(3).is_square()
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
            if status & GR_UNABLE: raise NotImplementedError(f"x = {self} is not a float in {self.parent()}")
            if status & GR_DOMAIN: raise ValueError(f"x = {self} is not a float in {self.parent()}")
        return c[0]

    # todo
    def __complex__(self):
        return float(self.re()) + float(self.im()) * 1j

    def inv(self):
        """
        Multiplicative inverse of this element.

            >>> QQ(3).inv()
            1/3
            >>> QQ(0).inv()
            Traceback (most recent call last):
              ...
            FlintDomainError: inv(x) is not an element of {Rational field (fmpq)} for {x = 0}
        """
        return self._unary_op(self, libgr.gr_inv, "inv($x)")

    def sqrt(self):
        """
        Square root of this element.

            >>> ZZ(4).sqrt()
            2
            >>> ZZ(2).sqrt()
            Traceback (most recent call last):
              ...
            FlintDomainError: sqrt(x) is not an element of {Integer ring (fmpz)} for {x = 2}
            >>> QQbar(2).sqrt()
            Root a = 1.41421 of a^2-2
            >>> (QQ(25)/16).sqrt()
            5/4
            >>> QQbar(-1).sqrt()
            Root a = 1.00000*I of a^2+1
            >>> RR(-1).sqrt()
            Traceback (most recent call last):
              ...
            FlintDomainError: sqrt(x) is not an element of {Real numbers (arb, prec = 53)} for {x = -1}
            >>> RF(-1).sqrt()
            nan

        """
        return self._unary_op(self, libgr.gr_sqrt, "sqrt($x)")

    def rsqrt(self):
        """
        Reciprocal square root of this element.

            >>> QQ(25).rsqrt()
            1/5
        """
        return self._unary_op(self, libgr.gr_rsqrt, "rsqrt($x)")

    def numerator(self):
        r"""
        Numerator of this element.

            >>> (QQ(-2) / 3).numerator()
            -2
            >>> ZZ(5).numerator()
            5
        """
        return self._unary_op(self, libgr.gr_numerator, "numerator($x)")

    def denominator(self):
        r"""
        Denominator of this element.

            >>> (QQ(-2) / 3).denominator()
            3
            >>> ZZ(5).denominator()
            1

        Depending on the ring, the denominator need not be minimal.
        This is currently not the case for algebraic numbers:

            >>> a = ((2 + QQbar.i()) / 4)
            >>> a.numerator()
            Root a = 8.00000 + 4.00000*I of a^2-16*a+80
            >>> a.denominator()
            16

        """
        return self._unary_op(self, libgr.gr_denominator, "denominator($x)")

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
        return self._unary_op(self, libgr.gr_floor, "floor($x)")

    def ceil(self):
        r"""
        Ceiling function: closest integer in the direction of `+\infty`.

            >>> (QQ(3) / 2).ceil()
            2
        """
        return self._unary_op(self, libgr.gr_ceil, "ceil($x)")

    def trunc(self):
        r"""
        Truncate to integer: closest integer in the direction of zero.

            >>> (QQ(3) / 2).trunc()
            1
        """
        return self._unary_op(self, libgr.gr_trunc, "trunc($x)")

    def nint(self):
        r"""
        Nearest integer function: nearest integer, rounding to
        even on a tie.

            >>> (QQ(3) / 2).nint()
            2
        """
        return self._unary_op(self, libgr.gr_nint, "nint($x)")

    def abs(self):
        return self._unary_op(self, libgr.gr_abs, "abs($x)")

    def conj(self):
        """
        Complex conjugate.

            >>> QQbar.i().conj()
            Root a = -1.00000*I of a^2+1
            >>> CC(-2).log().conj()
            ([0.693147180559945 +/- 4.12e-16] + [-3.141592653589793 +/- 3.39e-16]*I)
            >>> QQ(3).conj()
            3
        """
        return self._unary_op(self, libgr.gr_conj, "conj($x)")

    def re(self):
        """
        Real part.

            >>> QQ(1).re()
            1
            >>> (QQbar(-1) ** (QQ(1) / 3)).re()
            1/2
        """
        return self._unary_op(self, libgr.gr_re, "re($x)")

    def im(self):
        """
        Imaginary part.

            >>> QQ(1).im()
            0
            >>> (QQbar(-1) ** (QQ(1) / 3)).im()
            Root a = 0.866025 of 4*a^2-3
        """
        return self._unary_op(self, libgr.gr_im, "im($x)")

    def sgn(self):
        """
        Sign function.

            >>> QQ(-5).sgn()
            -1
            >>> CC(-10).sqrt().sgn()
            1.000000000000000*I
        """
        return self._unary_op(self, libgr.gr_sgn, "sgn($x)")

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
        return self._unary_op(self, libgr.gr_csgn, "csgn($x)")

    def arg(self):
        """
        >>> RR(1).arg()
        0
        >>> RR(-1).arg()
        [3.141592653589793 +/- 3.39e-16]
        >>> CC(1j).arg()
        [1.570796326794897 +/- 5.54e-16]
        >>> CC_ca.i().arg()
        1.57080 {(a)/2 where a = 3.14159 [Pi]}
        >>> RR("+/- 0.1").arg()
        [+/- 3.15]
        """
        return self._unary_op(self, libgr.gr_arg, "arg($x)")

    def mul_2exp(self, other):
        """
        Exact multiplication by a dyadic number `2^y`.

            >>> QQ(3).mul_2exp(5)
            96
            >>> QQ(3).mul_2exp(-5)
            3/32
            >>> ZZ(100).mul_2exp(-2)
            25
            >>> ZZ(100).mul_2exp(-3)
            Traceback (most recent call last):
              ...
            FlintDomainError: mul_2exp(x, y) is not an element of {Integer ring (fmpz)} for {x = 100}, {y = -3}
        """
        return self._binary_op_fmpz(self, other, libgr.gr_mul_2exp_fmpz, "mul_2exp($x, $y)")

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
            FlintUnableError: failed to compute exp(x) in {Rational field (fmpq)} for {x = 1}
            >>> QQser.gen().exp()
            1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 + (1/120)*x^5 + O(x^6)

        """
        return self._unary_op(self, libgr.gr_exp, "exp($x)")

    def expm1(self):
        """
        Exponential function minus 1.

            >>> RR("1e-10").expm1()
            [1.000000000050000e-10 +/- 3.86e-26]
            >>> CC(RR("1e-10")).expm1()
            [1.000000000050000e-10 +/- 3.86e-26]
            >>> RF("1e-10").expm1()
            1.000000000050000e-10
            >>> CF(RF("1e-10")).expm1()
            1.000000000050000e-10
            >>> QQ(0).expm1()
            0
            >>> QQ(1).expm1()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute expm1(x) in {Rational field (fmpq)} for {x = 1}
            >>> PowerSeriesModRing(RR, 4).gen().expm1()
            x + 0.5000000000000000*x^2 + [0.1666666666666667 +/- 7.04e-17]*x^3 (mod x^4)
            >>> (PowerSeriesModRing(RR, 2).gen() + 1).expm1()
            [1.718281828459045 +/- 5.41e-16] + [2.718281828459045 +/- 5.41e-16]*x (mod x^2)

        """
        return self._unary_op(self, libgr.gr_expm1, "expm1($x)")

    def exp2(self):
        """
        Exponential function with base 2.

            >>> QQ(5).exp2()
            32
            >>> RF(0.5).exp2()
            1.414213562373095
        """
        return self._unary_op(self, libgr.gr_exp2, "exp2($x)")

    def exp10(self):
        """
        Exponential function with base 10.

            >>> QQ(5).exp2()
            32
            >>> RF(0.5).exp10()
            3.162277660168380
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.exp()
            1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 (mod x^5)

        """
        return self._unary_op(self, libgr.gr_exp10, "exp10($x)")

    def log(self):
        """
        Natural logarithm.

            >>> QQ(1).log()
            0
            >>> QQ(2).log()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute log(x) in {Rational field (fmpq)} for {x = 2}
            >>> RF(2).log()
            0.6931471805599453
            >>> QQser(QQx([1, 1])).log()
            x + (-1/2)*x^2 + (1/3)*x^3 + (-1/4)*x^4 + (1/5)*x^5 + O(x^6)
            >>> QQser(1).log()
            0
            >>> x = PowerSeriesModRing(QQ, 5).gen(); (1+x).log()
            x + (-1/2)*x^2 + (1/3)*x^3 + (-1/4)*x^4 (mod x^5)

        """
        return self._unary_op(self, libgr.gr_log, "log($x)")

    def log1p(self):
        """
        Natural logarithm with one added to the argument.

            >>> QQ(0).log1p()
            0
            >>> RF(-0.5).log1p()
            -0.6931471805599453
            >>> RR(1).log1p()
            [0.693147180559945 +/- 4.12e-16]
        """
        return self._unary_op(self, libgr.gr_log1p, "log1p($x)")

    def sin(self):
        return self._unary_op(self, libgr.gr_sin, "sin($x)")

    def cos(self):
        return self._unary_op(self, libgr.gr_cos, "cos($x)")

    def tan(self):
        """
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.tan()
            x + (1/3)*x^3 (mod x^5)
        """
        return self._unary_op(self, libgr.gr_tan, "tan($x)")

    def sinh(self):
        return self._unary_op(self, libgr.gr_sinh, "sinh($x)")

    def cosh(self):
        return self._unary_op(self, libgr.gr_cosh, "cosh($x)")

    def tanh(self):
        return self._unary_op(self, libgr.gr_tanh, "tanh($x)")

    def asin(self):
        """
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.asin()
            x + (1/6)*x^3 (mod x^5)
        """
        return self._unary_op(self, libgr.gr_asin, "asin($x)")

    def acos(self):
        """
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.acos()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute acos(x) in {Power series over Rational field (fmpq) mod x^5} for {x = x (mod x^5)}
            >>> x = PowerSeriesModRing(RR_ca, 5).gen(); x.acos()
            (1.57080 {(a)/2 where a = 3.14159 [Pi]}) - x + (-0.166667 {-1/6})*x^3 (mod x^5)
        """
        return self._unary_op(self, libgr.gr_acos, "acos($x)")

    def atan(self):
        """
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.atan()
            x + (-1/3)*x^3 (mod x^5)
        """
        return self._unary_op(self, libgr.gr_atan, "atan($x)")

    def asinh(self):
        """
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.asinh()
            x + (-1/6)*x^3 (mod x^5)
        """
        return self._unary_op(self, libgr.gr_asinh, "asinh($x)")

    def acosh(self):
        """
            >>> x = PowerSeriesModRing(CC, 3).gen(); x.acosh()
            ([1.570796326794897 +/- 5.54e-16]*I) + (-1.000000000000000*I)*x (mod x^3)
        """
        return self._unary_op(self, libgr.gr_acosh, "acosh($x)")

    def atanh(self):
        """
            >>> x = PowerSeriesModRing(QQ, 5).gen(); x.atanh()
            x + (1/3)*x^3 (mod x^5)
        """
        return self._unary_op(self, libgr.gr_atanh, "atanh($x)")


    def exp_pi_i(self):
        r"""
        `\exp(\pi i x)` evaluated at self.

            >>> (QQbar(1) / 3).exp_pi_i()
            Root a = 0.500000 + 0.866025*I of a^2-a+1
            >>> (QQbar(2).sqrt()).exp_pi_i()
            Traceback (most recent call last):
              ...
            FlintDomainError: exp_pi_i(x) is not an element of {Complex algebraic numbers (qqbar)} for {x = Root a = 1.41421 of a^2-2}
        """
        return self._unary_op(self, libgr.gr_exp_pi_i, "exp_pi_i($x)")

    def log_pi_i(self):
        r"""
        `\log(x) / (\pi i)` evaluated at self.

            >>> (QQbar(-1) ** (QQbar(7) / 5)).log_pi_i()
            -3/5
            >>> (QQbar(1) / 2).log_pi_i()
            Traceback (most recent call last):
              ...
            FlintDomainError: log_pi_i(x) is not an element of {Complex algebraic numbers (qqbar)} for {x = 1/2}
        """
        return self._unary_op(self, libgr.gr_log_pi_i, "log_pi_i($x)")

    def sin_pi(self):
        r"""
        `\sin(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).sin_pi()
            Root a = 0.866025 of 4*a^2-3
        """
        return self._unary_op(self, libgr.gr_sin_pi, "sin_pi($x)")

    def cos_pi(self):
        r"""
        `\cos(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).cos_pi()
            1/2
        """
        return self._unary_op(self, libgr.gr_cos_pi, "cos_pi($x)")

    def tan_pi(self):
        r"""
        `\tan(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).tan_pi()
            Root a = 1.73205 of a^2-3
        """
        return self._unary_op(self, libgr.gr_tan_pi, "tan_pi($x)")

    def cot_pi(self):
        r"""
        `\cot(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).cot_pi()
            Root a = 0.577350 of 3*a^2-1
        """
        return self._unary_op(self, libgr.gr_cot_pi, "cot_pi($x)")

    def sec_pi(self):
        r"""
        `\sec(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).sec_pi()
            2
        """
        return self._unary_op(self, libgr.gr_sec_pi, "sec_pi($x)")

    def csc_pi(self):
        r"""
        `\csc(\pi x)` evaluated at self.

            >>> (QQbar(1) / 3).csc_pi()
            Root a = 1.15470 of 3*a^2-4
        """
        return self._unary_op(self, libgr.gr_csc_pi, "csc_pi($x)")

    def asin_pi(self):
        return self._unary_op(self, libgr.gr_asin_pi, "asin_pi($x)")

    def acos_pi(self):
        return self._unary_op(self, libgr.gr_acos_pi, "acos_pi($x)")

    def atan_pi(self):
        return self._unary_op(self, libgr.gr_atan_pi, "atan_pi($x)")

    def acot_pi(self):
        return self._unary_op(self, libgr.gr_acot_pi, "acot_pi($x)")

    def asec_pi(self):
        return self._unary_op(self, libgr.gr_asec_pi, "asec_pi($x)")

    def acsc_pi(self):
        return self._unary_op(self, libgr.gr_acsc_pi, "acsc_pi($x)")

    def erf(self):
        return self._unary_op(self, libgr.gr_erf, "erf($x)")

    def erfi(self):
        return self._unary_op(self, libgr.gr_erfi, "erfi($x)")

    def erfc(self):
        return self._unary_op(self, libgr.gr_erfc, "erfc($x)")

    def gamma(self):
        return self._unary_op(self, libgr.gr_gamma, "gamma($x)")

    def lgamma(self):
        return self._unary_op(self, libgr.gr_lgamma, "lgamma($x)")

    def rgamma(self):
        return self._unary_op(self, libgr.gr_rgamma, "lgamma($x)")

    def digamma(self):
        return self._unary_op(self, libgr.gr_digamma, "digamma($x)")

    def zeta(self):
        return self._unary_op(self, libgr.gr_zeta, "zeta($x)")


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

class QQbarField(gr_ctx):

    def set_limits(self, deg_limit=-1, bits_limit=-1):
        """
        Set evaluation limits preventing the creation of excessively
        large (degree or bits) algebraic numbers.
        Warning: currently not all methods respect these limits.

            >>> QQbar
            Complex algebraic numbers (qqbar)
            >>> QQbar.set_limits(deg_limit=6, bits_limit=100)
            >>> QQbar
            Complex algebraic numbers (qqbar), deg_limit = 6, bits_limit = 100
            >>> QQbar(2).sqrt() + QQbar(3).sqrt() + QQbar(5).sqrt()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute x + y in {Complex algebraic numbers (qqbar), deg_limit = 6, bits_limit = 100} for {x = Root a = 3.14626 of a^4-10*a^2+1}, {y = Root a = 2.23607 of a^2-5}
            >>> (QQbar(2).sqrt() + 1) ** 100
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute x ** y in {Complex algebraic numbers (qqbar), deg_limit = 6, bits_limit = 100} for {x = Root a = 2.41421 of a^2-2*a-1}, {y = 100}
            >>> QQbar.set_limits(deg_limit=-1, bits_limit=-1)
            >>> (QQbar(2).sqrt() + 1) ** 100
            Root a = 1.89482e+38 of a^2-189482250299273866835746159841800035874*a+1

        """
        libgr._gr_ctx_qqbar_set_limits(self._ref, deg_limit, bits_limit)
        self._str = None

class ComplexAlgebraicField_qqbar(QQbarField):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_qqbar(self._ref)
        self._elem_type = qqbar

class RealAlgebraicField_qqbar(QQbarField):
    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_real_qqbar(self._ref)
        self._elem_type = qqbar

class gr_arb_ctx(gr_ctx):
    pass


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

class ComplexExtended_ca(gr_ctx_ca):
    def __init__(self, **kwargs):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_complex_extended_ca(self._ref)
        self._elem_type = ca
        self._set_options(kwargs)


class PolynomialRing_gr_poly(gr_ctx):
    def __init__(self, coefficient_ring, var=None):
        assert isinstance(coefficient_ring, gr_ctx)
        gr_ctx.__init__(self)
        #if libgr.gr_ctx_is_ring(coefficient_ring._ref) != T_TRUE:
        #    raise ValueError("coefficient structure must be a ring")
        libgr.gr_ctx_init_gr_poly(self._ref, coefficient_ring._ref)
        coefficient_ring._refcount += 1
        self._coefficient_ring = coefficient_ring
        self._elem_type = gr_poly

        if var is not None:
            self._set_gen_name(var)

    def __del__(self):
        self._coefficient_ring._decrement_refcount()


class PowerSeriesRing_gr_series(gr_ctx):

    def __init__(self, coefficient_ring, prec=6, var=None):
        """
            >>> PowerSeriesRing(QQ, 5, var="y")
            Power series over Rational field (fmpq) with precision O(y^5)
        """
        assert isinstance(coefficient_ring, gr_ctx)
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_gr_series(self._ref, coefficient_ring._ref, prec)
        coefficient_ring._refcount += 1
        self._coefficient_ring = coefficient_ring
        self._elem_type = gr_series
        if var is not None:
            self._set_gen_name(var)

    def __del__(self):
        self._coefficient_ring._decrement_refcount()

class PowerSeriesModRing_gr_poly(gr_ctx):
    """
        >>> x = PowerSeriesModRing(ZZ, 3).gen()
        >>> (1+x)**1000
        1 + 1000*x + 499500*x^2 (mod x^3)
        >>> ((1+x)**2).sqrt() == (1+x)
        True
        >>> Rxy = PowerSeriesModRing(PowerSeriesModRing(RR, 2, "x"), 2, "y")
        >>> x, y = Rxy.gens(recursive=True)
        >>> (1+x+y).exp().log() - (1+x+y)
        ([+/- 3.89e-16] + [+/- 3.32e-16]*x (mod x^2)) + ([+/- 3.32e-16] + [+/- 6.63e-16]*x (mod x^2))*y (mod y^2)
    """

    def __init__(self, coefficient_ring, mod, var=None):
        assert isinstance(coefficient_ring, gr_ctx)
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_series_mod_gr_poly(self._ref, coefficient_ring._ref, mod)
        coefficient_ring._refcount += 1
        self._coefficient_ring = coefficient_ring
        self._elem_type = gr_poly
        if var is not None:
            self._set_gen_name(var)

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
    """
    Wrapper around the qqbar type, representing an algebraic number.

        >>> (qqbar(2).sqrt() / qqbar(-2).sqrt()) ** 2
        -1
        >>> qqbar(0.5) == qqbar(1) / 2
        True
        >>> qqbar(0.1) == qqbar(1) / 10
        False
        >>> qqbar(3+4j)
        Root a = 3.00000 + 4.00000*I of a^2-6*a+25
        >>> qqbar(3+4j).root(5)
        Root a = 1.35607 + 0.254419*I of a^10-6*a^5+25
        >>> qqbar(3+4j).root(5) ** 5
        Root a = 3.00000 + 4.00000*I of a^2-6*a+25

    The constructor can evaluate fexpr symbolic expressions
    provided that the expressions are constant and composed strictly
    of algebraic-valued basic operations applied to algebraic numbers.

        >>> fexpr.inject()
        >>> qqbar(Pow(0, 0))
        1
        >>> qqbar(Sqrt(2) * Abs(1+1j) + (+Re(3-4j)) + (-Im(5+6j)))
        -1
        >>> qqbar((Floor(Sqrt(1000)) + Ceil(Sqrt(1000)) + Sign(1+1j) / Sign(1-1j) + Csgn(1j) + Conjugate(1j)) ** Div(-1, 3))
        1/4
        >>> [qqbar(RootOfUnity(3)), qqbar(RootOfUnity(3,2))]
        [Root a = -0.500000 + 0.866025*I of a^2+a+1, Root a = -0.500000 - 0.866025*I of a^2+a+1]
        >>> qqbar(Decimal("0.125")) == qqbar(125)/1000
        True
        >>> qqbar(Decimal("-2.7e5")) == -270000
        True

    """

    _struct_type = qqbar_struct

    @staticmethod
    def _default_context():
        return QQbar

    # todo: generic
    def root(self, n):
        return self ** (qqbar(1) / n)

    def fexpr(self, formula=True, root_index=False, serialize=False,
            gaussians=True, quadratics=True, cyclotomics=True, cubics=True,
            quartics=True, quintics=True, depression=True, deflation=True,
            separation=True):
        """
        """
        res = fexpr()
        if serialize:
            libcalcium.qqbar_get_fexpr_repr(res, self)
            return res
        if formula:
            flags = 0
            if gaussians: flags |= 1
            if quadratics: flags |= 2
            if cyclotomics: flags |= 4
            if cubics: flags |= 8
            if quartics: flags |= 16
            if quintics: flags |= 32
            if depression: flags |= 64
            if deflation: flags |= 128
            if separation: flags |= 256
            if libcalcium.qqbar_get_fexpr_formula(res, self, flags):
                return res
        if root_index:
            libcalcium.qqbar_get_fexpr_root_indexed(res, self)
            return res
        libcalcium.qqbar_get_fexpr_root_nearest(res, self)
        return res

    def fexpr_repr(self):
        """
        """
        res = fexpr()
        libcalcium.qqbar_get_fexpr_repr(res, self)
        return res

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
    pass

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


@functools.cache
def get_nfloat_class(prec):
    n = (prec + FLINT_BITS - 1) // FLINT_BITS
    prec = n * FLINT_BITS

    class _nfloat_struct(ctypes.Structure):
        _fields_ = [('val', c_ulong * (n + 2))]

    _nfloat_struct.__qualname__ = _nfloat_struct.__name__ = ("nfloat" + str(prec) + "_struct")

    class _nfloat_class(gr_elem):
        _struct_type = _nfloat_struct

        @staticmethod
        def _default_context():
            raise NotImplementedError

    _nfloat_class.__qualname__ = _nfloat_class.__name__ = ("nfloat" + str(prec))

    return _nfloat_class

class RealFloat_nfloat(gr_ctx):
    def __init__(self, prec=128):
        gr_ctx.__init__(self)
        libflint.nfloat_ctx_init(self._ref, prec, 0)
        self._elem_type = get_nfloat_class(prec)

@functools.cache
def get_nfloat_complex_class(prec):
    n = (prec + FLINT_BITS - 1) // FLINT_BITS
    prec = n * FLINT_BITS

    class _nfloat_complex_struct(ctypes.Structure):
        _fields_ = [('val', c_ulong * (2 * (n + 2)))]

    _nfloat_complex_struct.__qualname__ = _nfloat_complex_struct.__name__ = ("nfloat" + str(prec) + "_complex_struct")

    class _nfloat_complex_class(gr_elem):
        _struct_type = _nfloat_complex_struct

        @staticmethod
        def _default_context():
            raise NotImplementedError

    _nfloat_complex_class.__qualname__ = _nfloat_complex_class.__name__ = ("nfloat" + str(prec) + "_complex")

    return _nfloat_complex_class

class ComplexFloat_nfloat_complex(gr_ctx):
    def __init__(self, prec=128):
        gr_ctx.__init__(self)
        libflint.nfloat_complex_ctx_init(self._ref, prec, 0)
        self._elem_type = get_nfloat_complex_class(prec)




class IntegersMod_nmod(gr_ctx):
    def __init__(self, n, n_is_prime=None):
        n = self._as_ui(n)
        assert n >= 1
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_nmod(self._ref, n)
        self._elem_type = nmod
        if n_is_prime is not None:
            libgr.gr_ctx_set_is_field(self, T_TRUE if n_is_prime else T_FALSE)

class nmod(gr_elem):
    _struct_type = nmod_struct


class IntegersMod_mpn_mod(gr_ctx):
    def __init__(self, n, n_is_prime=None):
        n = self._as_fmpz(n)
        # todo: error handling (must handle cleanup when ctx has not been initialized
        assert n >= (1 << FLINT_BITS) and n < (1 << (8 * FLINT_BITS))
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_mpn_mod(self._ref, n._ref)
        self._elem_type = mpn_mod
        if n_is_prime is not None:
            libgr.gr_ctx_set_is_field(self, T_TRUE if n_is_prime else T_FALSE)

class mpn_mod(gr_elem):
    _struct_type = mpn_mod_struct


class IntegersMod_fmpz_mod(gr_ctx):
    def __init__(self, n, n_is_prime=None):
        n = self._as_fmpz(n)
        assert n >= 1
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fmpz_mod(self._ref, n._ref)
        self._elem_type = fmpz_mod
        if n_is_prime is not None:
            libgr.gr_ctx_set_is_field(self, T_TRUE if n_is_prime else T_FALSE)

class fmpz_mod(gr_elem):
    _struct_type = fmpz_struct


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
    def __init__(self, p, n, var=None):
        gr_ctx.__init__(self)
        p = ZZ(p)
        n = int(n)
        assert p.is_prime()
        assert n >= 1
        if var is not None:
            var = ctypes.c_char_p(str(var).encode('ascii'))
        libgr.gr_ctx_init_fq(self._ref, p._ref, n, var)
        self._elem_type = fq

class FiniteField_fq_nmod(FiniteField_base):
    def __init__(self, p, n, var=None):
        gr_ctx.__init__(self)
        p = self._as_ui(p)
        n = int(n)
        assert ZZ(p).is_prime()
        assert n >= 1
        if var is not None:
            var = ctypes.c_char_p(str(var).encode('ascii'))
        libgr.gr_ctx_init_fq_nmod(self._ref, p, n, var)
        self._elem_type = fq_nmod

class FiniteField_fq_zech(FiniteField_base):
    def __init__(self, p, n, var=None):
        gr_ctx.__init__(self)
        p = self._as_ui(p)
        n = int(n)
        assert ZZ(p).is_prime()
        assert n >= 1
        if var is not None:
            var = ctypes.c_char_p(str(var).encode('ascii'))
        libgr.gr_ctx_init_fq_zech(self._ref, p, n, var)
        self._elem_type = fq_zech


class fq_elem(gr_elem):

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

    def pth_root(self):        return self._unary_op(self, libgr.gr_fq_pth_root, "pth_root")


class fq(fq_elem):
    _struct_type = fq_struct

class fq_nmod(fq_elem):
    _struct_type = fq_nmod_struct

class fq_zech(fq_elem):
    _struct_type = fq_zech_struct


class NumberField_nf(gr_ctx):
    def __init__(self, pol, var=None):
        pol = ZZx(pol)
        # assert pol.is_irreducible()
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_nf_fmpz_poly(self._ref, pol._ref)
        self._elem_type = nf_elem
        if var is not None:
            self._set_gen_name(var)

class nf_elem(gr_elem):
    _struct_type = nf_elem_struct



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

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

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

    def is_monic(self):
        """
            >>> RRx([2,3,4]).is_monic()
            False
            >>> RRx([2,3,1]).is_monic()
            True
            >>> RRx([]).is_monic()
            False

        """
        R = self.parent()._coefficient_ring
        truth = libgr.gr_poly_is_monic(self._ref, R._ref)
        def op(*args):
            return truth
        return gr_elem._unary_predicate(self, op, "is_monic")

    def monic(self):
        """
        Return self rescaled to a monic polynomial.

            >>> f = RRx([1,RR.pi()])
            >>> f.monic()
            [0.318309886183791 +/- 4.43e-16] + x
            >>> RRx([]).monic()   # the zero polynomial cannot be made monic
            Traceback (most recent call last):
              ...
            ValueError
            >>> (f - f).monic()   # unknown whether it is the zero polynomial
            Traceback (most recent call last):
              ...
            NotImplementedError

        """
        Rx = self.parent()
        R = Rx._coefficient_ring
        res = Rx()
        status = libgr.gr_poly_make_monic(res._ref, self._ref, R._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def derivative(self):
        Rx = self.parent()
        R = Rx._coefficient_ring
        res = Rx()
        status = libgr.gr_poly_derivative(res._ref, self._ref, R._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def integral(self):
        Rx = self.parent()
        R = Rx._coefficient_ring
        res = Rx()
        status = libgr.gr_poly_integral(res._ref, self._ref, R._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def roots(self, domain=None):
        """
        Computes the roots in the coefficient ring of this polynomial,
        returning a tuple (``roots``, ``multiplicities``).
        If the ring is not algebraically closed, the sum of multiplicities
        can be smaller than the degree of the polynomial.
        If ``domain`` is given, returns roots in that ring instead.

            >>> (ZZx([3,2]) * ZZx([15,1])**2 * ZZx([-10,1])).roots()
            ([10, -15], [1, 2])
            >>> ZZx([1]).roots()
            ([], [])

        We consider roots of the zero polynomial to be ill-defined:

            >>> ZZx([]).roots()
            Traceback (most recent call last):
              ...
            ValueError

        We construct an integer polynomial with rational, real algebraic
        and complex algebraic roots and extract its roots over
        different domains:

            >>> f = ZZx([-2,0,1]) * ZZx([1, 0, 1]) * ZZx([3, 2])**2
            >>> f.roots()   # integer roots (there are none)
            ([], [])
            >>> f.roots(domain=QQ)    # rational roots
            ([-3/2], [2])
            >>> f.roots(domain=AA)     # real algebraic roots
            ([Root a = 1.41421 of a^2-2, Root a = -1.41421 of a^2-2, -3/2], [1, 1, 2])
            >>> f.roots(domain=QQbar)     # complex algebraic roots
            ([Root a = 1.00000*I of a^2+1, Root a = -1.00000*I of a^2+1, Root a = 1.41421 of a^2-2, Root a = -1.41421 of a^2-2, -3/2], [1, 1, 1, 1, 2])
            >>> f.roots(domain=RR)      # real ball roots
            ([[-1.414213562373095 +/- 4.91e-17], [1.414213562373095 +/- 4.91e-17], -1.500000000000000], [1, 1, 2])
            >>> f.roots(domain=CC)      # complex ball roots
            ([[-1.414213562373095 +/- 4.91e-17], [1.414213562373095 +/- 4.91e-17], ([+/- 1.45e-19] + [1.000000000000000 +/- 2.7e-19]*I), ([+/- 1.45e-19] + [-1.000000000000000 +/- 2.7e-19]*I), -1.500000000000000], [1, 1, 1, 1, 2])
            >>> f.roots(RF)     # real floating-point roots
            ([-1.414213562373095, 1.414213562373095, -1.500000000000000], [1, 1, 2])
            >>> f.roots(CF)     # complex floating-point roots
            ([-1.414213562373095, 1.414213562373095, 1.000000000000000*I, -1.000000000000000*I, -1.500000000000000], [1, 1, 1, 1, 2])

        Calcium examples/tests:

            >>> PolynomialRing(CC_ca)([2,11,20,12]).roots()
            ([-0.666667 {-2/3}, -0.500000 {-1/2}], [1, 2])
            >>> PolynomialRing(RR_ca)([1,-1,0,1]).roots()
            ([-1.32472 {a where a = -1.32472 [a^3-a+1=0]}], [1])
            >>> PolynomialRing(CC_ca)([1,-1,0,1]).roots()
            ([-1.32472 {a where a = -1.32472 [a^3-a+1=0]}, 0.662359 + 0.562280*I {a where a = 0.662359 + 0.562280*I [a^3-a+1=0]}, 0.662359 - 0.562280*I {a where a = 0.662359 - 0.562280*I [a^3-a+1=0]}], [1, 1, 1])

        """
        Rx = self.parent()
        R = Rx._coefficient_ring
        mult = VecZZ()
        if domain is None:
            roots = Vec(R)()
            status = libgr.gr_poly_roots(roots._ref, mult._ref, self._ref, 0, R._ref)
        else:
            C = domain
            roots = Vec(C)()
            status = libgr.gr_poly_roots_other(roots._ref, mult._ref, self._ref, R._ref, 0, C._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (roots, mult)

    def _series_op(self, n, op, rstr):
        Rx = self.parent()
        R = Rx._coefficient_ring
        res = Rx()
        n = int(n)
        status = op(res._ref, self._ref, n, R._ref)
        if status:
            return _handle_error(Rx, status, rstr, self, n)
        return res

    def _series_op_fmpz_fmpq_overloads(self, other, n, op, op_fmpz, op_fmpq, rstr):
        Rx = self.parent()
        R = Rx._coefficient_ring
        res = Rx()
        n = int(n)
        # todo: variants
        other = QQ(other)
        status = op_fmpq(res._ref, self._ref, other._ref, n, R._ref)
        if status:
            return _handle_error(Rx, status, rstr, self, other, n)
        return res

    def _series_binary_op(self, other, n, op, rstr):
        Rx = self.parent()
        R = Rx._coefficient_ring
        # fixme
        other = Rx(other)
        res = Rx()
        n = int(n)
        status = op(res._ref, self._ref, other._ref, n, R._ref)
        if status:
            return _handle_error(Rx, status, rstr, self, n)
        return res

    def inv_series(self, n):
        """
        Reciprocal of this polynomial viewed as a power series,
        truncated to length n.

            >>> ZZx([1,2,3]).inv_series(10)
            22*x^9+73*x^8-56*x^7+13*x^6+10*x^5-11*x^4+4*x^3+x^2-2*x+1
            >>> ZZx([2,3,4]).inv_series(5)
            Traceback (most recent call last):
              ...
            FlintDomainError: f.inv_series(n) is not an element of {Polynomials over integers (fmpz_poly)} for {f = 4*x^2+3*x+2}, {n = 5}
            >>> QQx([2,3,4]).inv_series(5)
            (1/2) + (-3/4)*x + (1/8)*x^2 + (21/16)*x^3 + (-71/32)*x^4
        """
        return self._series_op(n, libgr.gr_poly_inv_series, "$f.inv_series($n)")

    def div_series(self, other, n):
        return self._series_binary_op(other, n, libgr.gr_poly_div_series, "$f.div_series(%g, $n)")

    def log_series(self, n):
        """
        Logarithm of this polynomial viewed as a power series,
        truncated to length n.

            >>> QQx([1,1]).log_series(8)
            x + (-1/2)*x^2 + (1/3)*x^3 + (-1/4)*x^4 + (1/5)*x^5 + (-1/6)*x^6 + (1/7)*x^7
            >>> RRx([2,1]).log_series(3)
            [0.693147180559945 +/- 4.12e-16] + 0.5000000000000000*x - 0.1250000000000000*x^2
            >>> RRx([0,0]).log_series(3)
            Traceback (most recent call last):
              ...
            FlintDomainError: f.log_series(n) is not an element of {Ring of polynomials over Real numbers (arb, prec = 53)} for {f = 0}, {n = 3}
        """
        return self._series_op(n, libgr.gr_poly_log_series, "$f.log_series($n)")

    def exp_series(self, n):
        """
        Exponential of this polynomial viewed as a power series,
        truncated to length n.

            >>> QQx([0,1]).exp_series(8)
            1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 + (1/120)*x^5 + (1/720)*x^6 + (1/5040)*x^7
            >>> QQx([1,1]).exp_series(2)
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute f.exp_series(n) in {Ring of polynomials over Rational field (fmpq)} for {f = 1 + x}, {n = 2}
            >>> RRx([1,1]).exp_series(2)
            [2.718281828459045 +/- 5.41e-16] + [2.718281828459045 +/- 5.41e-16]*x
            >>> RRx([2,3]).log_series(3).exp_series(3)
            [2.000000000000000 +/- 6.97e-16] + [3.00000000000000 +/- 1.61e-15]*x + [+/- 1.49e-15]*x^2
        """
        return self._series_op(n, libgr.gr_poly_exp_series, "$f.exp_series($n)")

    def pow_series(self, other, n):
        """
        Power of this polynomial viewed as a power series,
        truncated to length n.

            >>> QQx([4,3,2]).pow_series(QQ(1) / 2, 6)
            2 + (3/4)*x + (23/64)*x^2 + (-69/512)*x^3 + (299/16384)*x^4 + (2277/131072)*x^5
            >>> (QQx([4,3,2]) ** 2).pow_series(QQ(1) / 2, 6)
            4 + 3*x + 2*x^2
        """
        # todo
        return self._series_op_fmpz_fmpq_overloads(other, n, None, None, libgr.gr_poly_pow_series_fmpq_recurrence, "$f.pow_series($g, $n)")

    def atan_series(self, n):
        """
        Inverse tangent of this polynomial viewed as a power series,
        truncated to length n.

            >>> f = PolynomialRing(CC_ca)([2,3,4])
            >>> 2*f.atan_series(5) - ((2*f).div_series(1-f**2, 5)).atan_series(5) == CC_ca.pi()
            True
        """
        return self._series_op(n, libgr.gr_poly_atan_series, "$f.atan_series($n)")

    def atanh_series(self, n):
        return self._series_op(n, libgr.gr_poly_atanh_series, "$f.atanh_series($n)")

    def sqrt_series(self, n):
        return self._series_op(n, libgr.gr_poly_sqrt_series, "$f.sqrt_series($n)")

    def rsqrt_series(self, n):
        return self._series_op(n, libgr.gr_poly_rsqrt_series, "$f.rsqrt_series($n)")



class gr_series(gr_elem):
    _struct_type = gr_series_struct

    def __init__(self, val=None, error=None, context=None, random=False):
        gr_elem.__init__(self, val, context, random)
        if error is not None:
            raise NotImplementedError


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
    """
    Group of Dirichlet characters of given modulus.

        >>> G = DirichletGroup(10)
        >>> G.q
        10
        >>> len(G)
        4
        >>> [G(i) for i in [1,3,7,9]]
        [chi_10(1, .), chi_10(3, .), chi_10(7, .), chi_10(9, .)]
    """

    def __init__(self, q, **kwargs):
        # todo: automatic range checking with ctypes int -> c_ulong cast?
        if q <= 0:
            raise ValueError(f"modulus must not be zero")
        if q > UWORD_MAX:
            raise NotImplementedError(f"only word-size moduli are supported")
        gr_ctx.__init__(self)
        status = libgr.gr_ctx_init_dirichlet_group(self._ref, q)
        if status & GR_UNABLE: raise NotImplementedError(f"modulus with prime factor p > 10^16 is not currently supported")
        if status & GR_DOMAIN: raise ValueError(f"modulus must not be zero")
        self._elem_type = dirichlet_char
        self.q = int(q)    # for easy access

    def __len__(self):
        libarb.dirichlet_group_size.restype = c_slong
        libarb.dirichlet_group_size.argtypes = (ctypes.c_void_p,)
        return libarb.dirichlet_group_size(libgr.gr_ctx_data_as_ptr(self._ref))

    def __call__(self, n):
        n = int(n)
        assert 1 <= n <= max(self.q, 2) - 1
        assert ZZ(n).gcd(self.q) == 1
        x = dirichlet_char(context=self)
        libarb.dirichlet_char_log.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_ulong)
        libarb.dirichlet_char_log(x._ref, libgr.gr_ctx_data_as_ptr(self._ref), n)
        return x


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

    def norm_max(self):
        """
            >>> Mat(RR)([[1,2,3],[4,5,6],[7,8,9]]).norm_max()
            9.000000000000000
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        status = libgr.gr_mat_norm_max(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res


    def norm_1(self):
        """
            >>> Mat(RR)([[1,2,3],[4,5,6],[7,8,9]]).norm_1()
            18.00000000000000
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        status = libgr.gr_mat_norm_1(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def norm_inf(self):
        """
            >>> Mat(RR)([[1,2,3],[4,5,6],[7,8,9]]).norm_inf()
            24.00000000000000
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        status = libgr.gr_mat_norm_inf(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def norm_frobenius(self):
        """
            >>> Mat(RR)([[1,2,3],[4,5,6],[7,8,9]]).norm_frobenius()
            [16.88194301613413 +/- 3.73e-15]
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        status = libgr.gr_mat_norm_frobenius(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res


    def nullspace(self):
        """
        Right kernel (nullspace) of this matrix.

            >>> M = Mat(QQ)([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
            >>> X = M.nullspace()
            >>> X
            [[1],
            [-2],
            [1]]
            >>> M * X
            [[0],
            [0],
            [0]]
        """
        element_ring = self.parent()._element_ring
        X = self.parent()(0, 0)
        status = libgr.gr_mat_nullspace(X._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return X

    def det(self, algorithm=None):
        """
        Determinant of this matrix.

            >>> MatZZ(3, 3, ZZ.fac_vec(9)).det()
            233280
            >>> MatRR(3, 3, ZZ.fac_vec(9)).det()
            233280.0000000000
            >>> MatRR(3, 3, ZZ.fac_vec(9)).det(algorithm="lu")
            [233280.000000000 +/- 2.67e-10]
        """
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

    def trace(self):
        """
            >>> MatZZ([[3,4],[5,6]]).trace()
            9
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        status = libgr.gr_mat_trace(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def rank(self):
        """
            >>> MatZZ([[1,2,3],[4,5,6],[7,8,9]]).rank()
            2
            >>> Mat(CC_ca)([[1, 0, 0], [0, 1-(CC_ca(2)**-10).exp(), 0]]).rank()
            2
            >>> Mat(CC_ca)([[1, 0, 0], [0, 1-(CC_ca(2)**-10000).exp(), 0]]).rank()
            Traceback (most recent call last):
              ...
            NotImplementedError
        """
        element_ring = self.parent()._element_ring
        r = (ctypes.c_long * 1)()
        status = libgr.gr_mat_rank(r, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return ZZ(r[0])

    def solve(self, B):
        """
        Solves `AX = B` where `A` is given by self.
        Allows the system to be singular, undetermined, or
        overdetermined. If there are multiple solutions, an arbitrary
        solution is returned.

        This function currently only makes sense over fields.

            >>> A = MatQQ([[1,2,0], [0,1,0], [2,-2,0]])
            >>> B = MatQQ([[9], [2], [6]])
            >>> A.nonsingular_solve(B)
            Traceback (most recent call last):
              ...
            ValueError
            >>> A.solve(B)
            [[5],
            [2],
            [0]]
            >>> X = A.solve(B)
            >>> X
            [[5],
            [2],
            [0]]
            >>> A * X == B
            True

        """
        r = self.nrows()
        c = self.ncols()
        if r != c or r != B.nrows():
            raise ValueError
        element_ring = self.parent()._element_ring
        X = self.parent()(r, B.ncols())
        status = libgr.gr_mat_solve_field(X._ref, self._ref, B._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return X

    def nonsingular_solve(self, B, algorithm=None):
        """
        Proves invertibility of A (self) over the corresponding fraction field
        and solves `AX = B`.

            >>> A = MatQQ([[1,2],[3,4]])
            >>> B = MatQQ([[4],[5]])
            >>> X = A.nonsingular_solve(B)
            >>> A * X == B
            True

        The optional algorithm can be "lu" or "fflu".

            >>> MatZZ([[3,5],[1,2]]).nonsingular_solve(MatZZ([[1],[2]]), algorithm="fflu")
            [[-8],
            [5]]
        """
        r = self.nrows()
        c = self.ncols()
        if r != c or r != B.nrows():
            raise ValueError
        element_ring = self.parent()._element_ring
        X = self.parent()(r, B.ncols())
        if algorithm is None:
            status = libgr.gr_mat_nonsingular_solve(X._ref, self._ref, B._ref, element_ring._ref)
        elif algorithm == "lu":
            status = libgr.gr_mat_nonsingular_solve_lu(X._ref, self._ref, B._ref, element_ring._ref)
        elif algorithm == "fflu":
            status = libgr.gr_mat_nonsingular_solve_fflu(X._ref, self._ref, B._ref, element_ring._ref)
        else:
            raise ValueError("unknown algorithm")
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return X

    def nonsingular_solve_den(self, B):
        """
        Proves invertibility of A (self) over the corresponding fraction field
        and solves `A(X/d) = B`.

            >>> A = MatZZ([[3,4],[5,8]]); B = MatZZ([[1],[1]])
            >>> X, d = A.nonsingular_solve_den(B)
            >>> X
            [[4],
            [-2]]
            >>> d
            4
            >>> A*X == B*d
            True
        """
        r = self.nrows()
        c = self.ncols()
        if r != c or r != B.nrows():
            raise ValueError
        element_ring = self.parent()._element_ring
        X = self.parent()(r, B.ncols())
        den = element_ring()
        status = libgr.gr_mat_nonsingular_solve_den(X._ref, den._ref, self._ref, B._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return X, den

    def pascal(self, triangular=0):
        """
        Returns a Pascal matrix of the same shape.

            >>> MatZZ(4,5).pascal()
            [[1, 1, 1, 1, 1],
            [1, 2, 3, 4, 5],
            [1, 3, 6, 10, 15],
            [1, 4, 10, 20, 35]]
            >>> MatZZ(4,5).pascal(1)
            [[1, 1, 1, 1, 1],
            [0, 1, 2, 3, 4],
            [0, 0, 1, 3, 6],
            [0, 0, 0, 1, 4]]
            >>> MatZZ(4,5).pascal(-1)
            [[1, 0, 0, 0, 0],
            [1, 1, 0, 0, 0],
            [1, 2, 1, 0, 0],
            [1, 3, 3, 1, 0]]
        """
        element_ring = self.parent()._element_ring
        res = self.parent()(self.nrows(), self.ncols())
        status = libgr.gr_mat_pascal(res._ref, triangular, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def stirling(self, kind=0):
        """
        Returns a Stirling matrix of the same shape.

            >>> MatZZ(4,5).stirling()
            [[1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 1, 1, 0, 0],
            [0, 2, 3, 1, 0]]
            >>> MatZZ(4,5).stirling(1)
            [[1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, -1, 1, 0, 0],
            [0, 2, -3, 1, 0]]
            >>> MatZZ(4,5).stirling(2)
            [[1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 1, 1, 0, 0],
            [0, 1, 3, 1, 0]]
        """
        element_ring = self.parent()._element_ring
        res = self.parent()(self.nrows(), self.ncols())
        status = libgr.gr_mat_stirling(res._ref, kind, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def hilbert(self):
        """
        Returns a Hilbert matrix of the same shape.

            >>> MatQQ(2,3).hilbert()
            [[1, 1/2, 1/3],
            [1/2, 1/3, 1/4]]
        """
        element_ring = self.parent()._element_ring
        res = self.parent()(self.nrows(), self.ncols())
        status = libgr.gr_mat_hilbert(res._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def hadamard(self):
        """
        Returns a Hadamard matrix of the same shape.

            >>> MatZZ(4,4).hadamard()
            [[1, 1, 1, 1],
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, -1, -1, 1]]
            >>> MatZZ(3,3).hadamard()
            Traceback (most recent call last):
              ...
            ValueError

        """
        element_ring = self.parent()._element_ring
        res = self.parent()(self.nrows(), self.ncols())
        status = libgr.gr_mat_hadamard(res._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def charpoly(self, R=None, algorithm=None):
        """
        Characteristic polynomial of this matrix.

            >>> MatZZ([[1,0,1],[0,0,0],[1,0,1]]).charpoly()
            -2*x^2 + x^3
            >>> MatRR([[1,0,1],[0,0,0],[1,0,1]]).charpoly()
            -2.000000000000000*x^2 + x^3
            >>> Mat(CC_ca)([[5,CC_ca.pi()],[1,-1]]).charpoly()
            (-8.14159 {-a-5 where a = 3.14159 [Pi]}) - 4*x + x^2
        """
        mat_ring = self.parent()
        element_ring = mat_ring._element_ring
        poly_ring = R
        if poly_ring is None:
            poly_ring = PolynomialRing_gr_poly(element_ring)
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

    def minpoly(self, R=None):
        """
        Minimal polynomial of this matrix.
        This currently only makes sense over fields.

            >>> A = MatrixRing(QQ,3)([[1,0,1],[0,0,0],[1,0,1]])
            >>> A.minpoly()
            -2*x + x^2
            >>> A.minpoly()(A)
            [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]
        """
        mat_ring = self.parent()
        element_ring = mat_ring._element_ring
        poly_ring = R
        if poly_ring is None:
            poly_ring = PolynomialRing_gr_poly(element_ring)
        poly_element_ring = poly_ring._coefficient_ring
        assert element_ring is poly_element_ring
        res = poly_ring()
        status = libgr.gr_mat_minpoly_field(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def transpose(self):
        """
            >>> MatZZ(3,4,range(12)).transpose()
            [[0, 4, 8],
            [1, 5, 9],
            [2, 6, 10],
            [3, 7, 11]]
        """
        r = self.nrows()
        c = self.ncols()
        element_ring = self.parent()._element_ring
        res = gr_mat(c, r, context=self.parent())
        status = libgr.gr_mat_transpose(res._ref, self._ref, element_ring._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def is_scalar(self):
        """
        Return whether this matrix is a scalar matrix.
        """
        R = self.parent()._element_ring
        truth = libgr.gr_mat_is_scalar(self._ref, R._ref)
        def op(*args):
            return truth
        return gr_elem._unary_predicate(self, op, "is_scalar")

    def is_diagonal(self):
        """
        Return whether this matrix is a diagonal matrix.
        """
        R = self.parent()._element_ring
        truth = libgr.gr_mat_is_diagonal(self._ref, R._ref)
        def op(*args):
            return truth
        return gr_elem._unary_predicate(self, op, "is_diagonal")

    def is_upper_triangular(self):
        """
        Return whether this matrix is upper triangular.
        """
        R = self.parent()._element_ring
        truth = libgr.gr_mat_is_upper_triangular(self._ref, R._ref)
        def op(*args):
            return truth
        return gr_elem._unary_predicate(self, op, "is_upper_triangular")

    def is_lower_triangular(self):
        """
        Return whether this matrix is lower triangular.
        """
        R = self.parent()._element_ring
        truth = libgr.gr_mat_is_lower_triangular(self._ref, R._ref)
        def op(*args):
            return truth
        return gr_elem._unary_predicate(self, op, "is_lower_triangular")

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

    def is_hessenberg(self):
        """
        Return whether this matrix is in upper Hessenberg form.
        """
        R = self.parent()._element_ring
        truth = libgr.gr_mat_is_hessenberg(self._ref, R._ref)
        def op(*args):
            return truth
        return gr_elem._unary_predicate(self, op, "is_hessenberg")

    def eigenvalues(self, domain=None):
        """
        Computes the eigenvalues in the coefficient ring of this matrix,
        returning a tuple (``eigenvalues``, ``multiplicities``).
        If the ring is not algebraically closed, the sum of multiplicities
        can be smaller than the dimension of the matrix.
        If ``domain`` is given, returns eigenvalues in that ring instead.

            >>> Mat(ZZ)([[1,2],[3,4]]).eigenvalues()
            ([], [])
            >>> Mat(ZZ)([[1,2],[3,-4]]).eigenvalues()
            ([2, -5], [1, 1])
            >>> Mat(ZZ)([[1,2],[3,4]]).eigenvalues(domain=QQbar)
            ([Root a = 5.37228 of a^2-5*a-2, Root a = -0.372281 of a^2-5*a-2], [1, 1])
            >>> Mat(ZZ)([[1,2],[3,4]]).eigenvalues(domain=RR)
            ([[-0.3722813232690143 +/- 3.01e-17], [5.372281323269014 +/- 3.32e-16]], [1, 1])
            >>> Mat(QQbar)([[1, 0, QQbar.i()], [0, 0, 1], [1, 1, 1]]).eigenvalues()
            ([Root a = 1.94721 + 0.604643*I of a^6-4*a^5+4*a^4+2*a^3-3*a^2+1, Root a = 0.654260 - 0.430857*I of a^6-4*a^5+4*a^4+2*a^3-3*a^2+1, Root a = -0.601467 - 0.173786*I of a^6-4*a^5+4*a^4+2*a^3-3*a^2+1], [1, 1, 1])
            >>> Mat(ZZi)([[1, 0, ZZi.i()], [0, 0, 1], [1, 1, 1]]).eigenvalues(domain=QQbar)
            ([Root a = 1.94721 + 0.604643*I of a^6-4*a^5+4*a^4+2*a^3-3*a^2+1, Root a = 0.654260 - 0.430857*I of a^6-4*a^5+4*a^4+2*a^3-3*a^2+1, Root a = -0.601467 - 0.173786*I of a^6-4*a^5+4*a^4+2*a^3-3*a^2+1], [1, 1, 1])

        The matrix must be square:

            >>> Mat(ZZ)([[1,2,3],[4,5,6]]).eigenvalues()
            Traceback (most recent call last):
              ...
            ValueError

        """
        Rmat = self.parent()
        R = Rmat._element_ring
        mult = VecZZ()
        if domain is None:
            roots = Vec(R)()
            status = libgr.gr_mat_eigenvalues(roots._ref, mult._ref, self._ref, 0, R._ref)
        else:
            C = domain
            roots = Vec(C)()
            status = libgr.gr_mat_eigenvalues_other(roots._ref, mult._ref, self._ref, R._ref, 0, C._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (roots, mult)

    def diagonalization(self):
        """
        Matrix diagonalization: returns (D, L, R) where D is a vector
        of eigenvalues, LAR = diag(D) and LR = 1.

            >>> A = Mat(QQ)([[1,2],[-1,4]])
            >>> D, L, R = A.diagonalization()
            >>> L*A*R
            [[3, 0],
            [0, 2]]
            >>> D
            [3, 2]
            >>> L*R
            [[1, 0],
            [0, 1]]

            >>> A = Mat(CC)([[1,2],[-1,4]])
            >>> D, L, R = A.diagonalization()
            >>> D
            [([2.00000000000000 +/- 1.86e-15] + [+/- 1.86e-15]*I), ([3.00000000000000 +/- 2.90e-15] + [+/- 1.86e-15]*I)]
            >>> L*A*R
            [[([2.00000000000 +/- 1.10e-12] + [+/- 1.08e-12]*I), ([+/- 1.44e-12] + [+/- 1.42e-12]*I)],
            [([+/- 9.77e-13] + [+/- 9.63e-13]*I), ([3.00000000000 +/- 1.27e-12] + [+/- 1.25e-12]*I)]]
            >>> L*R
            [[([1.00000000000 +/- 3.26e-13] + [+/- 3.20e-13]*I), ([+/- 3.73e-13] + [+/- 3.67e-13]*I)],
            [([+/- 2.77e-13] + [+/- 2.73e-13]*I), ([1.00000000000 +/- 3.17e-13] + [+/- 3.13e-13]*I)]]

            >>> A = Mat(CF)([[1,2],[-1,4]])
            >>> D, L, R = A.diagonalization()
            >>> D
            [2.000000000000000, 3.000000000000000]
            >>> L*A*R
            [[2.000000000000000, -1.655022760610928e-16],
            [0, 3.000000000000000]]
            >>> L*R
            [[0.9999999999999998, -8.275113803054644e-17],
            [0, 1.000000000000000]]

            >>> M = Mat(CC_ca)
            >>> A = M([[1,2],[3,4]])
            >>> D, L, R = A.diagonalization()
            >>> D
            [5.37228 {(a+5)/2 where a = 5.74456 [a^2-33=0]}, -0.372281 {(-a+5)/2 where a = 5.74456 [a^2-33=0]}]
            >>> R * M([[D[0], 0], [0, D[1]]]) * L
            [[1, 2],
            [3, 4]]

        A diagonalizable matrix without distinct eigenvalues:

            >>> A = M([[-1,3,-1],[-3,5,-1],[-3,3,1]])
            >>> D, L, R = A.diagonalization()
            >>> D
            [1, 2, 2]
            >>> L
            [[3, -3, 1],
            [-3, 4, -1],
            [-3, 3, 0]]
            >>> R
            [[1, 1, -0.333333 {-1/3}],
            [1, 1, 0],
            [1, 0, 1]]
            >>> R * M([[D[0],0,0],[0,D[1],0],[0,0,D[2]]]) * L == A
            True

        """
        Rmat = self.parent()
        C = Rmat._element_ring
        D = Vec(C)()
        n = self.nrows()
        L = gr_mat(n, n, context=self.parent())
        R = gr_mat(n, n, context=self.parent())
        status = libgr.gr_mat_diagonalization(D._ref, L._ref, R._ref, self._ref, 0, C._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (D, L, R)


    #def __getitem__(self, i):
    #    pass



libgr.gr_mat_entry_ptr.argtypes = (ctypes.c_void_p, c_slong, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_mat_entry_ptr.restype = ctypes.POINTER(ctypes.c_char)

libgr.gr_vec_entry_ptr.restype = ctypes.POINTER(ctypes.c_char)


# todo singleton/cached domains (also for matrices, etc...)
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
        """
            >>> VecZZ(range(3, 20, 3))
            [3, 6, 9, 12, 15, 18]
        """
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
                elif isinstance(val, range):
                    start = val.start
                    step = val.step
                    n = len(val)
                    # todo: watch for slong -> int
                    status = libgr._gr_vec_check_resize(self._ref, n, self._ctx)
                    if not status:
                        start = element_ring(start)
                        step = element_ring(step)
                        iptr = libgr.gr_vec_entry_ptr(self._ref, 0, element_ring._ref)
                        status = libgr._gr_vec_step(iptr, start._ref, step._ref, n, element_ring._ref)
                if status:
                    if status & GR_UNABLE: raise NotImplementedError
                    if status & GR_DOMAIN: raise ValueError

    def __len__(self):
        return self._data.length

    def __getitem__(self, i):
        i = int(i)
        if not 0 <= i < len(self):
            raise IndexError
        element_ring = self.parent()._element_ring
        res = element_ring()
        iptr = libgr.gr_vec_entry_ptr(self._ref, i, res._ctx)
        status = libgr.gr_set(res._ref, iptr, res._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def __setitem__(self, i, v):
        i = int(i)
        if not 0 <= i < len(self):
            raise IndexError
        element_ring = self.parent()._element_ring
        # todo: avoid copy
        x = element_ring(v)
        iptr = libgr.gr_vec_entry_ptr(self._ref, i, x._ctx)
        status = libgr.gr_set(iptr, x._ref, x._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return x

    def sum(self):
        """
        Sum of the elements in this vector.

            >>> VecZZ(list(range(1,101))).sum()
            5050
            >>> VecZZ([]).sum()
            0
            >>> Vec(ZZmod(100))(list(range(1,101))).sum()
            50
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        ptr = libgr.gr_vec_entry_ptr(self._ref, 0, res._ctx)
        status = libgr._gr_vec_sum(res._ref, ptr, len(self), res._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    def product(self):
        """
        Product of the elements in this vector.

            >>> VecZZ(list(range(1,11))).product()
            3628800
            >>> VecZZ([]).product()
            1
            >>> Vec(ZZmod(103))(list(range(1,101))).product()
            51

        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        ptr = libgr.gr_vec_entry_ptr(self._ref, 0, res._ctx)
        status = libgr._gr_vec_product(res._ref, ptr, len(self), res._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res


class fmpz_poly(gr_poly):
    _struct_type = fmpz_poly_struct

    @staticmethod
    def _default_context():
        return ZZx_fmpz_poly

    def __init__(self, val=None, context=None):
        if isinstance(val, (list, tuple)):
            gr_elem.__init__(self, ZZx_gr_poly(val), context)
        else:
            gr_elem.__init__(self, val, context)

class fmpq_poly(gr_elem):
    _struct_type = fmpq_poly_struct

    @staticmethod
    def _default_context():
        return QQx_fmpq_poly

    def __init__(self, val=None, context=None):
        if isinstance(val, (list, tuple)):
            gr_elem.__init__(self, QQx_gr_poly(val), context)
        else:
            gr_elem.__init__(self, val, context)


class PolynomialRing_fmpz_poly(gr_ctx):

    def __init__(self, var=None):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fmpz_poly(self._ref)
        self._elem_type = fmpz_poly
        if var is not None:
            self._set_gen_name(var)

    @property
    def _coefficient_ring(self):
        return ZZ


class PolynomialRing_fmpq_poly(gr_ctx):

    def __init__(self, var=None):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fmpq_poly(self._ref)
        self._elem_type = fmpq_poly
        if var is not None:
            self._set_gen_name(var)

ZZx_fmpz_poly = PolynomialRing_fmpz_poly()
QQx_fmpq_poly = PolynomialRing_fmpq_poly()


class fmpz_poly(gr_poly):
    _struct_type = fmpz_poly_struct

    @staticmethod
    def _default_context():
        return ZZx_fmpz_poly

    def __init__(self, val=None, context=None):
        if isinstance(val, (list, tuple)):
            gr_elem.__init__(self, ZZx_gr_poly(val), context)
        else:
            gr_elem.__init__(self, val, context)


class fmpz_mpoly(gr_elem):
    _struct_type = fmpz_mpoly_struct

class PolynomialRing_fmpz_mpoly(gr_ctx):

    def __init__(self, nvars, vars=None):
        gr_ctx.__init__(self)
        nvars = gr_ctx._as_si(nvars)
        assert nvars >= 0
        libgr.gr_ctx_init_fmpz_mpoly(self._ref, nvars, 0)
        self._elem_type = fmpz_mpoly
        if vars is not None:
            assert len(vars) == nvars
            self._set_gen_names(vars)

    @property
    def _coefficient_ring(self):
        return ZZ



class gr_mpoly(gr_elem):
    _struct_type = gr_mpoly_struct

class PolynomialRing_gr_mpoly(gr_ctx):
    def __init__(self, coefficient_ring, nvars, vars=None):
        assert isinstance(coefficient_ring, gr_ctx)
        gr_ctx.__init__(self)

        gr_ctx.__init__(self)
        nvars = gr_ctx._as_si(nvars)
        assert nvars >= 0
        libgr.gr_ctx_init_gr_mpoly(self._ref, coefficient_ring._ref, nvars, 0)
        self._elem_type = fmpz_mpoly

        coefficient_ring._refcount += 1
        self._coefficient_ring = coefficient_ring
        self._elem_type = gr_mpoly

        if vars is not None:
            assert len(vars) == nvars
            self._set_gen_names(vars)

    def __del__(self):
        self._coefficient_ring._decrement_refcount()




class fmpz_mpoly_q(gr_elem):
    _struct_type = fmpz_mpoly_q_struct

class FractionField_fmpz_mpoly_q(gr_ctx):

    def __init__(self, nvars, vars=None):
        gr_ctx.__init__(self)
        nvars = gr_ctx._as_si(nvars)
        assert nvars >= 0
        libgr.gr_ctx_init_fmpz_mpoly_q(self._ref, nvars, 0)
        self._elem_type = fmpz_mpoly_q

        if vars is not None:
            assert len(vars) == nvars
            self._set_gen_names(vars)

    @property
    def _coefficient_ring(self):
        return QQ








class fexpr(gr_elem):

    _struct_type = fexpr_struct

    @staticmethod
    def inject(vars=False):
        """
        Inject all builtin symbol names into the calling namespace.
        For interactive use only!

            >>> fexpr.inject()
            >>> n = fexpr("n")
            >>> Sum(Sin(Pi*n/3)/Factorial(n), For(n,0,Infinity))
            Sum(Div(Sin(Div(Mul(Pi, n), 3)), Factorial(n)), For(n, 0, Infinity))

        """
        from inspect import currentframe
        frame = currentframe().f_back
        num = libflint.fexpr_builtin_length()
        for i in range(num):
            # memory leak
            symbol_name = libflint.fexpr_builtin_name(i)
            symbol_name = symbol_name.decode('ascii')
            if not symbol_name[0].islower():
                frame.f_globals[symbol_name] = fexpr(symbol_name)
        if vars:
            def inject_vars(string):
                for s in string.split():
                    for symbol_name in [s, s + "_"]:
                        frame.f_globals[symbol_name] = fexpr(symbol_name)
            inject_vars("""a b c d e f g h i j k l m n o p q r s t u v w x y z""")
            inject_vars("""A B C D E F G H I J K L M N O P Q R S T U V W X Y Z""")
            inject_vars("""alpha beta gamma delta epsilon zeta eta theta iota kappa lamda mu nu xi pi rho sigma tau phi chi psi omega ell varphi vartheta""")
            inject_vars("""Alpha Beta GreekGamma Delta Epsilon Zeta Eta Theta Iota Kappa Lamda Mu Nu Xi GreekPi Rho Sigma Tau Phi Chi Psi Omega""")
        del frame

    def inject_dict():
        vs = {}
        num = libflint.fexpr_builtin_length()
        for i in range(num):
            # memory leak
            symbol_name = libflint.fexpr_builtin_name(i)
            symbol_name = symbol_name.decode('ascii')
            if not symbol_name[0].islower():
                vs[symbol_name] = fexpr(symbol_name)
        def inject_vars(string):
            for s in string.split():
                for symbol_name in [s, s + "_"]:
                    vs[symbol_name] = fexpr(symbol_name)
        inject_vars("""a b c d e f g h i j k l m n o p q r s t u v w x y z""")
        inject_vars("""A B C D E F G H I J K L M N O P Q R S T U V W X Y Z""")
        inject_vars("""alpha beta gamma delta epsilon zeta eta theta iota kappa lamda mu nu xi pi rho sigma tau phi chi psi omega ell varphi vartheta""")
        inject_vars("""Alpha Beta GreekGamma Delta Epsilon Zeta Eta Theta Iota Kappa Lamda Mu Nu Xi GreekPi Rho Sigma Tau Phi Chi Psi Omega""")
        return vs

    def builtins():
        num = libflint.fexpr_builtin_length()
        names = []
        for i in range(num):
            # memory leak
            symbol_name = libflint.fexpr_builtin_name(i)
            symbol_name = symbol_name.decode('ascii')
            names.append(symbol_name)
        return names

    # todo: deduplicate conversion code
    def __init__(self, val=None, context=None):
        if context is None:
            context = FExpressions
        gr_elem.__init__(self, None, context)
        #self._data = fexpr_struct()
        #self._ref = ctypes.byref(self._data)
        #libflint.fexpr_init(self)
        if val is not None:
            typ = type(val)
            if typ is int:
                b = sys.maxsize
                if -b <= val <= b:
                    libflint.fexpr_set_si(self, val)
                else:
                    n = fmpz_struct()
                    nref = ctypes.byref(n)
                    libflint.fmpz_init(nref)
                    libflint.fmpz_set_str(nref, ctypes.c_char_p(str(val).encode('ascii')), 10)
                    libflint.fexpr_set_fmpz(self, nref)
                    libflint.fmpz_clear(nref)
            elif typ is str:
                if val[0] == "'" or val[0] == '"':
                    libflint.fexpr_set_string(self, val[1:-1].encode('ascii'))
                else:
                    libflint.fexpr_set_symbol_str(self, val.encode('ascii'))
            elif typ is float:
                libflint.fexpr_set_d(self, val)
            elif typ is complex:
                libflint.fexpr_set_re_im_d(self, val.real, val.imag)
            elif typ is bool:
                if val:
                    libflint.fexpr_set_symbol_str(self, ("True").encode('ascii'))
                else:
                    libflint.fexpr_set_symbol_str(self, ("False").encode('ascii'))
            elif issubclass(typ, gr_elem):
                status = libflint.gr_get_fexpr(self, val._ref, val._ctx)
                if status:
                    _handle_error(context, status, "fexpr($x)", val)
            #elif typ is qqbar:
            #    #libflint.qqbar_get_fexpr_repr(self, val, val._ctx)
            #    tmp = val.fexpr()
            #    libflint.fexpr_set(self, tmp)
            #elif typ is ca:
            #    libflint.ca_get_fexpr(self, val, 0, val._ctx)
            #elif typ is ca_mat:
            #    libflint.ca_mat_get_fexpr(self, val, 0, val._ctx)
            #elif typ is ca_poly:
            #    libflint.ca_poly_get_fexpr(self, val, 0, val._ctx)
            elif typ is tuple:
                tmp = fexpr("Tuple")(*val)         # todo: create without copying
                libflint.fexpr_set(self, tmp)
            elif typ is list:
                tmp = fexpr("List")(*val)
                libflint.fexpr_set(self, tmp)
            elif typ is set:
                tmp = fexpr("Set")(*val)
                libflint.fexpr_set(self, tmp)
            else:
                raise TypeError

    def __del__(self):
        libflint.fexpr_clear(self)

    @property
    def _as_parameter_(self):
        return self._ref

    @staticmethod
    def from_param(arg):
        return arg

    def __repr__(self):
        ptr = libflint.fexpr_get_str(self)
        try:
            return ctypes.cast(ptr, ctypes.c_char_p).value.decode("ascii")
        finally:
            libflint.flint_free(ptr)

    def latex(self):
        ptr = libflint.fexpr_get_str_latex(self, 0)
        try:
            return ctypes.cast(ptr, ctypes.c_char_p).value.decode()
        finally:
            libflint.flint_free(ptr)

    def _repr_latex_(self):
        return "$$" + self.latex() + "$$"

    def nwords(self):
        return libflint.fexpr_size(self)

    def size_bytes(self):
        return libflint.fexpr_size_bytes(self)

    def allocated_bytes(self):
        return libflint.fexpr_allocated_bytes(self)

    def num_leaves(self):
        return libflint.fexpr_num_leaves(self)

    def depth(self):
        return libflint.fexpr_depth(self)

    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        if libflint.fexpr_equal(self, other):
            return True
        return False

    def is_atom(self):
        return bool(libflint.fexpr_is_atom(self))

    def is_atom_integer(self):
        return bool(libflint.fexpr_is_integer(self))

    def is_symbol(self):
        return bool(libflint.fexpr_is_symbol(self))

    def head(self):
        if libflint.fexpr_is_atom(self):
            return None
        res = fexpr()
        libflint.fexpr_func(res, self)
        return res

    def nargs(self):
        # todo: long
        if self.is_atom():
            return None
        return libflint.fexpr_nargs(self)

    def args(self):
        if libflint.fexpr_is_atom(self):
            return None
        n = self.nargs()
        args = [fexpr() for i in range(n)]
        for i in range(n):
            libflint.fexpr_arg(args[i], self, i)
        return tuple(args)

    def __hash__(self):
        return libflint.fexpr_hash(self)

    def __call__(self, *args):
        args2 = []
        for arg in args:
            tp = type(arg)
            if tp is not fexpr:
                if tp is str:
                    arg = "'" + arg + "'"
                arg = fexpr(arg)
            args2.append(arg)
        n = len(args2)
        res = fexpr()
        if n == 0:
            libflint.fexpr_call0(res, self)
        elif n == 1:
            libflint.fexpr_call1(res, self, args2[0])
        elif n == 2:
            libflint.fexpr_call2(res, self, args2[0], args2[1])
        elif n == 3:
            libflint.fexpr_call3(res, self, args2[0], args2[1], args2[2])
        elif n == 4:
            libflint.fexpr_call4(res, self, args2[0], args2[1], args2[2], args2[3])
        else:
            vec = libflint.flint_malloc(n * ctypes.sizeof(fexpr_struct))
            vec = ctypes.cast(vec, ctypes.POINTER(fexpr_struct))
            for i in range(n):
                vec[i] = args2[i]._data
            libflint.fexpr_call_vec(res, self, vec, n)
            libflint.flint_free(vec)
        return res

    def contains(self, x):
        """
        Check if *x* appears exactly as a subexpression in *self*.

            >>> f = fexpr("f"); x = fexpr("x"); y = fexpr("y")
            >>> (f(x+1).contains(f), f(x+1).contains(x), f(x+1).contains(y))
            (True, True, False)
            >>> (f(x+1).contains(1), f(x+1).contains(2))
            (True, False)
            >>> (f(x+1).contains(x+1), f(x+1).contains(f(x+1)))
            (True, True)
        """
        if type(x) is not fexpr:
            x = fexpr(x)
        if libflint.fexpr_contains(self, x):
            return True
        return False

    def replace(self, old, new=None):
        """
        Replace subexpression.

            >>> f = fexpr("f"); x = fexpr("x"); y = fexpr("y")
            >>> f(x+1, x-1).replace(x, y)
            f(Add(y, 1), Sub(y, 1))
            >>> f(x+1, x-1).replace(x+1, y-1)
            f(Sub(y, 1), Sub(x, 1))
            >>> f(x+1, x-1).replace(f, f+1)
            Add(f, 1)(Add(x, 1), Sub(x, 1))
            >>> f(x+1, x-1).replace(x+2, y)
            f(Add(x, 1), Sub(x, 1))
        """
        # todo: dict replacement
        if type(old) is not fexpr:
            old = fexpr(old)
        if type(new) is not fexpr:
            new = fexpr(new)
        res = fexpr()
        libflint.fexpr_replace(res, self, old, new)
        return res

    def __add__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_add(res, self, other)
        return res

    def __radd__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_add(res, other, self)
        return res

    def __sub__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_sub(res, self, other)
        return res

    def __rsub__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_sub(res, other, self)
        return res

    def __mul__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_mul(res, self, other)
        return res

    def __rmul__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_mul(res, other, self)
        return res

    def __truediv__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_div(res, self, other)
        return res

    def __rtruediv__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_div(res, other, self)
        return res

    def __pow__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_pow(res, self, other)
        return res

    def __rpow__(self, other):
        if type(self) is not type(other):
            try:
                other = fexpr(other)
            except TypeError:
                return NotImplemented
        res = fexpr()
        libflint.fexpr_pow(res, other, self)
        return res

    # def __floordiv__(self, other):
    #     return (self / other).floor()
    # def __rfloordiv__(self, other):
    #     return (other / self).floor()

    def __bool__(self):
        return True

    def __abs__(self):
        return fexpr("Abs")(self)

    def __neg__(self):
        res = fexpr()
        libflint.fexpr_neg(res, self)
        return res

    def __pos__(self):
        return fexpr("Pos")(self)

    def expanded_normal_form(self):
        """
        Converts this expression to expanded normal form as
        a formal rational function of its non-arithmetic subexpressions.

            >>> x = fexpr("x"); y = fexpr("y")
            >>> (x / x**2).expanded_normal_form()
            Div(1, x)
            >>> (((x ** 0) + 3) ** 5).expanded_normal_form()
            1024
            >>> ((x+y+1)**3 - (y+1)**3 - (x+y)**3 - (x+1)**3).expanded_normal_form()
            Add(Mul(-1, Pow(x, 3)), Mul(6, x, y), Mul(-1, Pow(y, 3)), -1)
            >>> (1/((1/y + 1/x))).expanded_normal_form()
            Div(Mul(x, y), Add(x, y))
            >>> (((x+y)**5 * (x-y)) / (x**2 - y**2)).expanded_normal_form()
            Add(Pow(x, 4), Mul(4, Pow(x, 3), y), Mul(6, Pow(x, 2), Pow(y, 2)), Mul(4, x, Pow(y, 3)), Pow(y, 4))
            >>> (1 / (x - x)).expanded_normal_form()
            Traceback (most recent call last):
              ...
            ValueError: expanded_normal_form: overflow, formal division by zero or unsupported expression
        """
        res = fexpr()
        if not libflint.fexpr_expanded_normal_form(res, self, 0):
            raise ValueError("expanded_normal_form: overflow, formal division by zero or unsupported expression")
        return res

    def nstr(self, n=16):
        """
        Evaluates this expression numerically using Arb, returning
        a decimal string correct within 1 ulp in the last output digit.
        Attempts to obtain *n* digits (but the actual output accuracy
        may be lower).

            >>> Exp = fexpr("Exp"); Exp(1).nstr()
            '2.718281828459045'
            >>> Pi = fexpr("Pi"); Pi.nstr(30)
            '3.14159265358979323846264338328'
            >>> Log = fexpr("Log"); Log(-2).nstr()
            '0.6931471805599453 + 3.141592653589793*I'
            >>> Im = fexpr("Im")
            >>> Im(Log(2)).nstr()   # exact zero
            '0'

        Here the imaginary part is zero, but Arb is not able to
        compute so exactly. The output ``0e-N``
        indicates only that the absolute value is bounded by ``1e-N``:

            >>> Exp(Log(-2)).nstr()
            '-2.000000000000000 + 0e-22*I'
            >>> Im(Exp(Log(-2))).nstr()
            '0e-731'

        The algorithm fails if the expression or any subexpression
        is not a finite complex number:

            >>> Log(0).nstr()
            Traceback (most recent call last):
              ...
            ValueError: nstr: unable to evaluate to a number

        Expressions must be constant:

            >>> fexpr("x").nstr()
            Traceback (most recent call last):
              ...
            ValueError: nstr: unable to evaluate to a number

        """
        ptr = libflint.fexpr_get_decimal_str(self, n, 0)
        try:
            s = ctypes.cast(ptr, ctypes.c_char_p).value.decode("ascii")
            if s == "?":
                raise ValueError("nstr: unable to evaluate to a number")
            return s
        finally:
            libflint.flint_free(ptr)

class Expressions_fexpr(gr_ctx):

    def __init__(self):
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_fexpr(self._ref)
        self._elem_type = fexpr

FExpressions = Expressions_fexpr()


libflint.fexpr_builtin_name.restype = ctypes.c_char_p
libflint.fexpr_set_symbol_str.argtypes = ctypes.c_void_p, ctypes.c_char_p
libflint.fexpr_get_str.restype = ctypes.c_void_p
libflint.fexpr_get_str_latex.restype = ctypes.c_void_p
libflint.fexpr_set_si.argtypes = fexpr, ctypes.c_long
libflint.fexpr_set_d.argtypes = fexpr, ctypes.c_double
libflint.fexpr_set_re_im_d.argtypes = fexpr, ctypes.c_double, ctypes.c_double
libflint.fexpr_get_decimal_str.restype = ctypes.c_void_p





# todo: def .one()


PolynomialRing = PolynomialRing_gr_poly
PowerSeriesRing = PowerSeriesRing_gr_series
PowerSeriesModRing = PowerSeriesModRing_gr_poly

NumberField = NumberField_nf

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

def ZZmod(n):
    # todo: selection
    return IntegersMod_nmod(n)

ZZp16 = ZZmod((1 << 15) + 3)
ZZp32 = ZZmod((1 << 31) + 11)
ZZp63 = ZZmod((1 << 62) + 135)
ZZp64 = ZZmod((1 << 63) + 29)

VecZZ = Vec(ZZ)
VecQQ = Vec(QQ)
VecRR = Vec(RR)
VecCC = Vec(CC)
VecRF = Vec(RF)
VecCF = Vec(CF)

MatZZ = Mat(ZZ)
MatQQ = Mat(QQ)
MatRR = Mat(RR)
MatCC = Mat(CC)
MatRF = Mat(RF)
MatCF = Mat(CF)

ZZx_gr_poly = PolynomialRing_gr_poly(ZZ)
ZZx = ZZx_fmpz_poly
QQx = QQx_gr_poly = PolynomialRing_gr_poly(QQ)
RRx = RRx_arb = PolynomialRing_gr_poly(RR_arb)
CCx = CCx_acb = PolynomialRing_gr_poly(CC_acb)
RRx_ca = PolynomialRing_gr_poly(RR_ca)
CCx_ca = PolynomialRing_gr_poly(CC_ca)

ZZser = PowerSeriesRing(ZZ)
QQser = PowerSeriesRing(QQ)
RRser = RRser_arb = PowerSeriesRing(RR_arb)
CCser = CCser_acb = PowerSeriesRing(CC_acb)
RRser_ca = PowerSeriesRing(RR_ca)
CCser_ca = PowerSeriesRing(CC_ca)

# QQx = QQx_fmpq_poly

ModularGroup = ModularGroup_psl2z
DirichletGroup = DirichletGroup_dirichlet_char

PSL2Z = ModularGroup()
SymmetricGroup = SymmetricGroup_perm

def timing(f, *args, **kwargs):
    once = kwargs.get('once')
    if 'once' in kwargs:
        del kwargs['once']
    if args or kwargs:
        if len(args) == 1 and not kwargs:
            arg = args[0]
            g = lambda: f(arg)
        else:
            g = lambda: f(*args, **kwargs)
    else:
        g = f
    from timeit import default_timer as clock
    t1=clock(); v=g(); t2=clock(); t=t2-t1
    if t > 0.05 or once:
        return t
    for i in range(3):
        t1=clock();
        # Evaluate multiple times because the timer function
        # has a significant overhead
        g();g();g();g();g();g();g();g();g();g()
        t2=clock()
        t=min(t,(t2-t1)/10)
    return t

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

def test_polynomial():
    poly_types = [ZZx_fmpz_poly, ZZx_gr_poly, QQx_fmpq_poly, QQx_gr_poly]
    for A in poly_types:
        for B in poly_types + [VecZZ, VecQQ]:
            for C in poly_types:
                assert A(B([1,2,3])) == C([1,2,3])

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

    MatZZ = Mat(ZZ)
    A = MatZZ([[1, 2, 3], [0, 4, 5], [0, 0, 6]])
    assert A.is_upper_triangular()
    assert not A.is_lower_triangular()
    assert A.transpose().is_lower_triangular()
    assert not A.transpose().is_upper_triangular()
    A = MatZZ([[1, 2, 3], [1, 4, 5], [0, 5, 6]])
    assert A.is_hessenberg()
    assert not A.transpose().is_hessenberg()
    assert not A.is_diagonal()
    assert MatZZ([[1, 0, 0], [0, 2, 0], [0, 0, 3]]).is_diagonal()

    assert not A.is_scalar()
    assert not MatZZ([[1,0],[0,2]]).is_scalar()
    assert MatZZ([[1,0],[0,1]]).is_scalar()

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

    assert RR(3).nint() == 3
    assert 2.9 < RR("3.5 +/- 0.1").nint() < 4.1

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
        x = R(3) / 2 + R.i()
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
    x = QQ.bernoulli(50)
    sign, primes, exponents = x.factor()
    assert (sign * (primes ** exponents)).product() == x

def test_qqbar():
    a = (-23 + 5*ZZi.i())
    assert ZZi(QQbar(a**2).sqrt()) == -a
    x = (1 + 100*qqbar(2)**(-1000)).root(100)
    y = (1 + 101*qqbar(2)**(-1000)).root(101)
    assert x > y

    assert QQ(ComplexAlgebraicField_qqbar()(ca(-3).sqrt())**2) == -3

def test_nf():
    a = NumberField_nf(ZZx([1,2,3])).gen()
    assert (a+5)**(-1) * (2*a+10) == 2

def test_ca():
    R = ComplexField_ca()
    X = ComplexExtended_ca()
    assert X.inf() == -X.neg_inf()
    assert 1 / X(0) == X.uinf()
    assert 0 / X(0) == X.undefined()
    assert X(0).log() == X.neg_inf()
    assert X.inf().exp() == X.inf()
    assert raises(lambda: R(X.inf()), FlintDomainError)
    assert raises(lambda: R(X.undefined()), FlintDomainError)
    assert raises(lambda: R(X.unknown()), FlintUnableError)
    assert X(R(1)) == 1
    assert R(X(1)) == 1

def test_ca_more():
    sqrt = CC_ca.sqrt
    exp = CC_ca.exp
    log = CC_ca.log
    tan = CC_ca.tan
    tanh = CC_ca.tanh
    sin = CC_ca.sin
    cos = CC_ca.cos
    acos = CC_ca.acos
    sqrt = CC_ca.sqrt
    arg = CC_ca.arg
    atan = CC_ca.atan
    re = CC_ca.re
    im = CC_ca.im
    floor = CC_ca.floor
    ceil = CC_ca.ceil
    pi = CC_ca.pi()
    i = CC_ca.i()
    e = CC_ca.exp(1)
    erf = CC_ca.erf
    erfc = CC_ca.erfc
    erfi = CC_ca.erfi

    Sqrt = fexpr("Sqrt")

    def gd(x):
        return 2*atan(exp(x))-pi/2

    # todo: proper gr versions
    def cosh(x):
        y = exp(x)
        return (y + 1/y)/2

    def sinh(x):
        y = exp(x)
        return (y - 1/y)/2

    def tanh(x):
        return sinh(x)/cosh(x)

    def gamma(x):
        return CC_ca(x).gamma()

    assert floor(sqrt(2)) == 1
    assert ceil(sqrt(2)) == 2

    assert (sqrt(2)**sqrt(2))**sqrt(2) == 2
    assert (sqrt(-2)**sqrt(2))**sqrt(2) == -2
    assert (sqrt(3)**sqrt(3))**sqrt(3) == 3*sqrt(3)
    assert sqrt(-pi)**2 == -pi

    assert log(1+pi) - log(pi) - log(1+1/pi) == 0
    assert log(log(-log(log(exp(exp(-exp(exp(3)))))))) == 3

    assert exp(pi*i) + 1 == 0
    assert exp(pi*i) == -1
    assert exp(log(2)*log(3)) > 2
    assert e**2 == exp(2)

    assert erf(2*log(sqrt(ca(1)/2-sqrt(2)/4))+log(4)) - erf(log(2-sqrt(2))) == 0
    assert 1-erf(pi)-erfc(pi) == 0
    assert erf(sqrt(2))**2 + erfi(sqrt(-2))**2 == 0

    assert sin(gd(1)) == tanh(1)
    assert tan(gd(1)) == sinh(1)
    assert sin(gd(sqrt(2))) == tanh(sqrt(2))
    assert tan(gd(1)/2) - tanh(ca(1)/2) == 0

    assert gamma(1) == 1
    assert gamma(ca(1)/2) == sqrt(pi)
    assert gamma(sqrt(2)*sqrt(3)) == gamma(sqrt(6))
    assert gamma(pi+1)/gamma(pi) == pi
    assert gamma(pi)/gamma(pi-1) == pi-1
    assert log(gamma(pi+1)) - log(gamma(pi)) - log(pi) == 0
    assert log(gamma(-pi+1)) - log(gamma(-pi)) - log(pi) == pi * i

    assert raises(lambda: gamma(0), FlintDomainError)
    assert 1/ComplexExtended_ca()(0).gamma() == 0

    assert ca(qqbar(sqrt(2))) == sqrt(2)
    assert ca(fexpr(ca(Sqrt(2)))) == sqrt(2)

    assert sqrt(2)*(1+i)/2 * pi - exp(pi*i/4) * pi == 0
    assert (sqrt(3) + i)/2  * pi - exp(pi*i/6) * pi == 0
    assert arg(sqrt(-pi*i)) == -pi/4
    assert (pi + sqrt(2) + sqrt(3)) / (pi + sqrt(5 + 2*sqrt(6))) == 1
    assert log(1/exp(sqrt(2)+1)) == -sqrt(2)-1
    assert abs(exp(sqrt(1+i))) == exp(re(sqrt(1+i)))
    assert tan(pi*sqrt(2))*tan(pi*sqrt(3)) == (cos(pi*sqrt(5-2*sqrt(6))) - cos(pi*sqrt(5+2*sqrt(6))))/(cos(pi*sqrt(5-2*sqrt(6))) + cos(pi*sqrt(5+2*sqrt(6))))
    assert log(exp(i) / exp(-i)) == 2*i

    v = cos(acos(sqrt(2) - sqrt(3))/3)
    assert 1 - 90*v**2 + 321*v**4 - 592*v**6 + 864*v**8 - 768*v**10 + 256*v**12 == 0

    def expect_not_implemented(f):
        try:
            v = f()
            assert v
        except NotImplementedError:
            return
        raise AssertionError

    expect_not_implemented(lambda: acos(cos(1)) == 1)
    expect_not_implemented(lambda: acos(cos(sqrt(2) - 1)) == sqrt(2) - 1)
    expect_not_implemented(lambda: tan(sqrt(pi*2))*tan(sqrt(pi*3)) == \
        (cos(sqrt(pi*(5-2*sqrt(6)))) - cos(sqrt(pi*(5+2*sqrt(6)))))/(cos(sqrt(pi*(5-2*sqrt(6)))) + cos(sqrt(pi*(5+2*sqrt(6))))))
    expect_not_implemented(lambda: sqrt(exp(2*sqrt(2)) + exp(-2*sqrt(2)) - 2) == (exp(2*sqrt(2))-1)/sqrt(exp(2*sqrt(2))))

    # Some examples from Stoutemyer
    z = pi
    assert str(i*((3-5*i)*z+1)/(((5+3*i)*z+i)*z)) == "0.318310 {(1)/(a) where a = 3.14159 [Pi], b = I [b^2+1=0]}"
    assert str(ca(-1)**(ca(1)/8) * sqrt(i+1) / ca(2)**(ca(3)/4) + i*exp(i*pi/2)) == "-0.500000 + 0.500000*I {(a-1)/2 where a = I [a^2+1=0]}"
    assert sin(2*atan(z)) + (15*sqrt(3)+26)**(ca(1)/3) == 2*z/(z**2+1) + sqrt(3) + 2

def test_ca_comparisons():

    R = ComplexExtended_ca()

    inf = R.inf()
    uinf = R.uinf()
    undefined = R.undefined()
    i = R.i()
    pi = R.pi()
    sqrt = R.sqrt
    exp = R.exp
    sin = R.sin
    cos = R.cos
    asin = R.asin

    assert inf >= inf
    assert inf >= -inf
    assert inf >= 3
    assert not -inf >= inf
    assert -inf >= -inf
    assert not -inf >= 3

    assert not inf > inf
    assert inf > -inf
    assert inf > 3
    assert not -inf > inf
    assert not -inf > -inf
    assert not -inf >= 3

    assert inf <= inf
    assert not inf <= -inf
    assert not inf <= 3
    assert -inf <= inf
    assert -inf <= -inf
    assert -inf <= 3

    assert not inf < inf
    assert not inf < -inf
    assert not inf < 3
    assert -inf < inf
    assert not -inf < -inf
    assert -inf < 3

    # as currently implemented, these comparisons raise
    #for a in [inf, -inf, ca(3)]:
    #    for b in [uinf, undefined, inf*i, 2+i]:
    #        assert not a < b
    #        assert not a <= b
    #        assert not a > b
    #        assert not a >= b
    #        assert not b < a
    #        assert not b <= a
    #        assert not b > a
    #        assert not b >= a

    assert pi + sqrt(2) + sqrt(3) <= pi + sqrt(5 + 2 * sqrt(6))
    assert not pi + sqrt(2) + sqrt(3) < pi + sqrt(5 + 2 * sqrt(6))
    assert pi + sqrt(2) + sqrt(3) >= pi + sqrt(5 + 2 * sqrt(6))
    assert not pi + sqrt(2) + sqrt(3) > pi + sqrt(5 + 2 * sqrt(6))
    assert pi + exp(-1000) > pi
    assert pi - exp(-1000) < pi
    assert sin(1) < 1
    assert sin(1) <= sqrt(1 - cos(1)**2)
    assert asin(ca(1)/10) > ca(1)/10

def test_latex():

    latex_test_cases = [
        ("""f(0)""", "f(0)"),
        ("""f("Hello, world!")""", "f(\\text{``Hello, world!''})"),
        ("f", r"f"),
        ("f_", r"f"),
        ("f_(0)", r"f_{0}"),
        ("f()", r"f()"),
        ("Add(Add(Add(f(a, b), c_(n)), f_(x, y)), f_())", r"f(a, b) + c_{n} + f_{x, y} + f_{}"),
        ("f(alpha, beta, chi, delta, ell, epsilon, eta)", r"f(\alpha, \beta, \chi, \delta, \ell, \varepsilon, \eta)"),
        ("f(gamma, iota, kappa, lamda, mu, nu, omega, phi)", r"f(\gamma, \iota, \kappa, \lambda, \mu, \nu, \omega, \phi)"),
        ("f(pi, rho, sigma, tau, theta, varphi, vartheta, xi, zeta)", r"f(\pi, \rho, \sigma, \tau, \theta, \varphi, \vartheta, \xi, \zeta)"),
        ("f(Delta, GreekGamma, GreekPi, Lamda, Omega, Phi, Psi, Sigma, Theta, Xi)", r"f(\Delta, \Gamma, \Pi, \Lambda, \Omega, \Phi, \Psi, \Sigma, \Theta, \Xi)"),
        ("f(alpha_, beta_, chi_, delta_, ell_, epsilon_, eta_)", r"f(\alpha, \beta, \chi, \delta, \ell, \varepsilon, \eta)"),
        ("f(gamma_, iota_, kappa_, lamda_, mu_, nu_, omega_, phi_)", r"f(\gamma, \iota, \kappa, \lambda, \mu, \nu, \omega, \phi)"),
        ("f(pi_, rho_, sigma_, tau_, theta_, varphi_, vartheta_, xi_, zeta)", r"f(\pi, \rho, \sigma, \tau, \theta, \varphi, \vartheta, \xi, \zeta)"),
        ("f(Delta_, GreekGamma_, GreekPi_, Lamda_, Omega_, Phi_, Psi_, Sigma_, Theta_, Xi_)", r"f(\Delta, \Gamma, \Pi, \Lambda, \Omega, \Phi, \Psi, \Sigma, \Theta, \Xi)"),
        ("f(alpha(x), beta(x), chi(x), delta(x), ell(x), epsilon(x), eta(x))", r"f\!\left(\alpha(x), \beta(x), \chi(x), \delta(x), \ell(x), \varepsilon(x), \eta(x)\right)"),
        ("f(gamma(x), iota(x), kappa(x), lamda(x), mu(x), nu(x), omega(x), phi(x))", r"f\!\left(\gamma(x), \iota(x), \kappa(x), \lambda(x), \mu(x), \nu(x), \omega(x), \phi(x)\right)"),
        ("f(pi(x), rho(x), sigma(x), tau(x), theta(x), varphi(x), vartheta(x), xi(x), zeta(x))", r"f\!\left(\pi(x), \rho(x), \sigma(x), \tau(x), \theta(x), \varphi(x), \vartheta(x), \xi(x), \zeta(x)\right)"),
        ("f(Delta(x), GreekGamma(x), GreekPi(x), Lamda(x), Omega(x), Phi(x), Psi(x), Sigma(x), Theta(x), Xi(x))", r"f\!\left(\Delta(x), \Gamma(x), \Pi(x), \Lambda(x), \Omega(x), \Phi(x), \Psi(x), \Sigma(x), \Theta(x), \Xi(x)\right)"),
        ("f(alpha_(n), beta_(n), chi_(n), delta_(n), ell_(n), epsilon_(n), eta_(n))", r"f\!\left(\alpha_{n}, \beta_{n}, \chi_{n}, \delta_{n}, \ell_{n}, \varepsilon_{n}, \eta_{n}\right)"),
        ("f(gamma_(n), iota_(n), kappa_(n), lamda_(n), mu_(n), nu_(n), omega_(n), phi_(n))", r"f\!\left(\gamma_{n}, \iota_{n}, \kappa_{n}, \lambda_{n}, \mu_{n}, \nu_{n}, \omega_{n}, \phi_{n}\right)"),
        ("f(pi_(n), rho_(n), sigma_(n), tau_(n), theta_(n), varphi_(n), vartheta_(n), xi_(n), zeta_(n))", r"f\!\left(\pi_{n}, \rho_{n}, \sigma_{n}, \tau_{n}, \theta_{n}, \varphi_{n}, \vartheta_{n}, \xi_{n}, \zeta_{n}\right)"),
        ("f(Delta_(n), GreekGamma_(n), GreekPi_(n), Lamda_(n), Omega_(n), Phi_(n), Psi_(n), Sigma_(n), Theta_(n), Xi_(n))", r"f\!\left(\Delta_{n}, \Gamma_{n}, \Pi_{n}, \Lambda_{n}, \Omega_{n}, \Phi_{n}, \Psi_{n}, \Sigma_{n}, \Theta_{n}, \Xi_{n}\right)"),
        ("f(a)(b)", r"f(a)(b)"),
        ("Mul(f, g)(x)", r"\left(f g\right)(x)"),
        ("Add(f, g)(x)", r"\left(f + g\right)(x)"),
        ("c_(m, n, p)(x, y, z)", r"c_{m, n, p}(x, y, z)"),
        ("f_(Div(-3, 2))(Div(-3, 2))", r"f_{-3 / 2}\!\left(-\frac{3}{2}\right)"),
        ("Mul(2, Pi, NumberI)", r"2 \pi i"),
        ("Mul(-2, Pi, NumberI)", r"-2 \pi i"),
        ("Mul(1)", r"1"),
        ("Mul(-1)", r"-1"),
        ("Mul(1, x, y)", r" x y"),
        ("Mul(-1, x, y)", r"- x y"),
        ("Mul(1, 2, 3)", r"1 \cdot 2 \cdot 3"),
        ("Mul(-1, 2, 3)", r"-1 \cdot 2 \cdot 3"),
        ("Mul(-1, -1)", r"-1 \cdot \left(-1\right)"),
        ("Add(2, x)", r"2 + x"),
        ("Sub(2, x)", r"2 - x"),
        ("Add(2, Neg(x))", r"2 + \left(-x\right)"),
        ("Sub(2, Neg(x))", r"2 - \left(-x\right)"),
        ("Add(2, Sub(x, y))", r"2 + \left(x - y\right)"),
        ("Sub(2, Add(x, y))", r"2 - \left(x + y\right)"),
        ("Sub(2, Sub(x, y))", r"2 - \left(x - y\right)"),
        ("Add(-3, -4, -5)", r"-3-4-5"),
        ("Add(-3, Mul(-4, x), -5)", r"-3-4 x-5"),
        ("Add(-3, Mul(-4, x), Pos(5))", r"-3-4 x+5"),
        ("Add(-3, Div(Mul(-4, x), 7), -5)", r"-3-\frac{4 x}{7}-5"),
        ("Add(-3, Div(Mul(4, x), 7), -5)", r"-3 + \frac{4 x}{7}-5"),
        ("Mul(0, 1)", r"0 \cdot 1"),
        ("Mul(3, Pow(2, n))", r"3 \cdot {2}^{n}"),
        ("Mul(3, Pow(-1, n))", r"3 \cdot {\left(-1\right)}^{n}"),
        ("Mul(-3, Pow(-1, n))", r"-3 \cdot {\left(-1\right)}^{n}"),
        ("Mul(-1, -2, -3)", r"-1 \cdot \left(-2\right) \cdot \left(-3\right)"),
        ("Div(-1, 3)", r"-\frac{1}{3}"),
        ("Div(Mul(-5, Pi), 3)", r"-\frac{5 \pi}{3}"),
        ("Div(Neg(Mul(5, Pi)), 3)", r"-\frac{5 \pi}{3}"),
        ("Div(Add(Add(Mul(-5, Pow(x, 2)), Mul(4, x)), 1), Add(Mul(3, x), y))", r"\frac{-5 {x}^{2} + 4 x + 1}{3 x + y}"),
        ("Pow(2, n)", r"{2}^{n}"),
        ("Pow(-1, n)", r"{\left(-1\right)}^{n}"),
        ("Pow(10, Pow(10, -10))", r"{10}^{{10}^{-10}}"),
        ("Pow(2, Div(-1, 3))", r"{2}^{-1 / 3}"),
        ("Pow(2, Div(-1, Mul(3, n)))", r"{2}^{-1 / \left(3 n\right)}"),
        ("Equal(Add(Pow(Sin(x), 2), Pow(Cos(x), 2)), 1)", r"\sin^{2}\!\left(x\right) + \cos^{2}\!\left(x\right) = 1"),
        ("Set(Set(), Set(1), Set(1, 2, 3))", r"\left\{\left\{\right\}, \left\{1\right\}, \left\{1, 2, 3\right\}\right\}"),
        ("Tuple(Tuple(), Tuple(1), Tuple(1, 2, 3))", r"\left(\left(\right), \left(1\right), \left(1, 2, 3\right)\right)"),
        ("List(List(), List(1), List(1, 2, 3))", r"\left[\left[\right], \left[1\right], \left[1, 2, 3\right]\right]"),
        ("Set(f(x), For(x, CC))", r"\left\{ f(x) : x \in \mathbb{C} \right\}"),
        ("Set(f(x), For(x, CC), NotEqual(x, 0))", r"\left\{ f(x) : x \in \mathbb{C}\,\mathbin{\operatorname{and}}\, x \ne 0 \right\}"),
        ("List(Floor(x), Ceil(x), Abs(x), RealAbs(x), Conjugate(z), Sqrt(x))", r"\left[\left\lfloor x \right\rfloor, \left\lceil x\right\rceil, \left|x\right|, \left|x\right|, \overline{z}, \sqrt{x}\right]"),
        ("List(Floor(Div(1, 2)), Ceil(Div(1, 2)), Abs(Div(1, 2)), RealAbs(Div(1, 2)), Conjugate(Div(1, 2)), Sqrt(Div(1, 2)))", r"\left[\left\lfloor \frac{1}{2} \right\rfloor, \left\lceil \frac{1}{2}\right\rceil, \left|\frac{1}{2}\right|, \left|\frac{1}{2}\right|, \overline{\frac{1}{2}}, \sqrt{\frac{1}{2}}\right]"),
        ("And(Equal(Length(Tuple(1, 2, 3)), 3), Equal(Cardinality(Set()), 0))", r"\# \left(1, 2, 3\right) = 3 \;\mathbin{\operatorname{and}}\; \# \left\{\right\} = 0"),
        ("Tuple(Parentheses(x), Brackets(x), Braces(x), AngleBrackets(x))", r"\left(\left(x\right), \left[x\right], \left\{x\right\}, \left\langle x\right\rangle\right)"),
        ("Tuple(Parentheses(Div(1, 2)), Brackets(Div(1, 2)), Braces(Div(1, 2)), AngleBrackets(Div(1, 2)))", r"\left(\left(\frac{1}{2}\right), \left[\frac{1}{2}\right], \left\{\frac{1}{2}\right\}, \left\langle \frac{1}{2}\right\rangle\right)"),
        ("Concatenation(A, B)", r"A  \,^\frown  B"),
        ("Equal(Concatenation(Tuple(a, b), Tuple(c, d, e), Tuple()), Tuple(a, b, c, d, e))", r"\left(a, b\right)  \,^\frown  \left(c, d, e\right)  \,^\frown  \left(\right) = \left(a, b, c, d, e\right)"),
        ("Equal(f(x), Cases(Case(y, P(x)), Case(Neg(y), Q(x))))", r"f(x) = \begin{cases} y, & P(x)\\-y, & Q(x)\\ \end{cases}"),
        ("Equal(f(x), Cases(Case(y, P(x)), Case(Neg(y), Q(x)), Case(0, Otherwise)))", r"f(x) = \begin{cases} y, & P(x)\\-y, & Q(x)\\0, & \text{otherwise}\\ \end{cases}"),
        ("And(Equal(True, Not(False)), Equal(False, Not(True)))", r"\operatorname{True} = \operatorname{not} \operatorname{False} \;\mathbin{\operatorname{and}}\; \operatorname{False} = \operatorname{not} \operatorname{True}"),
        ("All(Greater(x, 0), For(x, S))", r"x > 0 \;\text{ for all } x \in S"),
        ("All(Greater(x, 0), For(x, S), P(x))", r"x > 0 \;\text{ for all } x \in S \text{ with } P(x)"),
        ("Exists(Greater(x, 0), For(x, S))", r"x > 0 \;\text{ for some } x \in S"),
        ("Exists(Greater(x, 0), For(x, S), P(x))", r"x > 0 \;\text{ for some } x \in S \text{ with } P(x)"),
        ("Logic(All(Greater(x, 0), For(x, S)))", r"\forall x \in S : \, x > 0"),
        ("Logic(All(Greater(x, 0), For(x, S), P(x)))", r"\forall x \in S, \,P(x) : \, x > 0"),
        ("Logic(Exists(Greater(x, 0), For(x, S)))", r"\exists x \in S : \, x > 0"),
        ("Logic(Exists(Greater(x, 0), For(x, S), P(x)))", r"\exists x \in S, \,P(x) : \, x > 0"),
        ("Or(Q, And(P, Q, Not(P), Or(Q, P), Not(Or(Q, P))))", r"Q \;\mathbin{\operatorname{or}}\; \left(P \;\mathbin{\operatorname{and}}\; Q \;\mathbin{\operatorname{and}}\; \operatorname{not} P \;\mathbin{\operatorname{and}}\; \left(Q \;\mathbin{\operatorname{or}}\; P\right) \;\mathbin{\operatorname{and}}\; \operatorname{not} \,\left(Q \;\mathbin{\operatorname{or}}\; P\right)\right)"),
        ("Logic(Or(Q, And(P, Q, Not(P), Or(Q, P), Not(Or(Q, P)))))", r"Q \,\lor\, \left(P \,\land\, Q \,\land\, \neg P \,\land\, \left(Q \,\lor\, P\right) \,\land\, \neg \left(Q \,\lor\, P\right)\right)"),
        ("Equivalent(A, B)", r"A \iff B"),
        ("Equivalent(Not(Equal(x, y)), NotEqual(x, y))", r"\left(\operatorname{not} \,\left(x = y\right)\right) \iff \left(x \ne y\right)"),
        ("Implies(P, And(R, S))", r"P \;\implies\; \left(R \;\mathbin{\operatorname{and}}\; S\right)"),
        ("Implies(Element(x, QQ), Element(x, RR))", r"x \in \mathbb{Q} \;\implies\; x \in \mathbb{R}"),
        ("And(Less(x, y), Less(x, y, z), LessEqual(x, y), LessEqual(x, y, z))", r"x < y \;\mathbin{\operatorname{and}}\; x < y < z \;\mathbin{\operatorname{and}}\; x \le y \;\mathbin{\operatorname{and}}\; x \le y \le z"),
        ("And(Greater(x, y), Greater(x, y, z), GreaterEqual(x, y), GreaterEqual(x, y, z))", r"x > y \;\mathbin{\operatorname{and}}\; x > y > z \;\mathbin{\operatorname{and}}\; x \ge y \;\mathbin{\operatorname{and}}\; x \ge y \ge z"),
        ("Subset(Primes, NN, ZZ, QQ, RR, CC)", r"\mathbb{P} \subset \mathbb{N} \subset \mathbb{Z} \subset \mathbb{Q} \subset \mathbb{R} \subset \mathbb{C}"),
        ("Subset(QQ, AlgebraicNumbers, CC)", r"\mathbb{Q} \subset \overline{\mathbb{Q}} \subset \mathbb{C}"),
        ("SubsetEqual(S, QQ)", r"S \subseteq \mathbb{Q}"),
        ("NotElement(123456789012345678901234567890, SetMinus(QQ, ZZ))", r"123456789012345678901234567890 \notin \mathbb{Q} \setminus \mathbb{Z}"),
        ("KroneckerDelta(x, Div(1, 2))", r"\delta_{(x,1 / 2)}"),
        ("Set(Interval(a, b), OpenInterval(a, b), ClosedOpenInterval(a, b), OpenClosedInterval(a, b))", r"\left\{\left[a, b\right], \left(a, b\right), \left[a, b\right), \left(a, b\right]\right\}"),
        ("Set(Interval(a, Div(1, 2)), OpenInterval(a, Div(1, 2)), ClosedOpenInterval(a, Div(1, 2)), OpenClosedInterval(a, Div(1, 2)))", r"\left\{\left[a, 1 / 2\right], \left(a, 1 / 2\right), \left[a, 1 / 2\right), \left(a, 1 / 2\right]\right\}"),
        ("Set(RealBall(m, r), OpenRealBall(m, r))", r"\left\{\left[m \pm r\right], \left(m \pm r\right)\right\}"),
        ("Set(ClosedComplexDisk(m, r), OpenComplexDisk(m, r))", r"\left\{\overline{D}(m, r), D(m, r)\right\}"),
        ("Set(Undefined, UnsignedInfinity, Pos(Infinity), Neg(Infinity))", r"\left\{\mathfrak{u}, \hat{\infty}, +\infty, -\infty\right\}"),
        ("Equal(RealSignedInfinities, Set(Pos(Infinity), Neg(Infinity)))", r"\{\pm \infty\} = \left\{+\infty, -\infty\right\}"),
        ("Equal(ComplexSignedInfinities, Set(Mul(Exp(Mul(NumberI, theta)), Infinity), For(theta, OpenClosedInterval(Neg(Pi), Pi))))", r"\{[e^{i \theta}] \infty\} = \left\{ e^{i \theta} \cdot \infty : \theta \in \left(-\pi, \pi\right] \right\}"),
        ("Equal(RealInfinities, Union(RealSignedInfinities, Set(UnsignedInfinity)))", r"\{\hat{\infty}, \pm \infty\} = \{\pm \infty\} \cup \left\{\hat{\infty}\right\}"),
        ("Equal(ComplexInfinities, Union(ComplexSignedInfinities, Set(UnsignedInfinity)))", r"\{\hat{\infty}, [e^{i \theta}] \infty\} = \{[e^{i \theta}] \infty\} \cup \left\{\hat{\infty}\right\}"),
        ("Equal(ExtendedRealNumbers, Union(RR, RealSignedInfinities))", r"\overline{\mathbb{R}} = \mathbb{R} \cup \{\pm \infty\}"),
        ("Equal(SignExtendedComplexNumbers, Union(CC, ComplexSignedInfinities))", r"\overline{\mathbb{C}}_{[e^{i \theta}] \infty} = \mathbb{C} \cup \{[e^{i \theta}] \infty\}"),
        ("Equal(ProjectiveRealNumbers, Union(RR, Set(UnsignedInfinity)))", r"\hat{\mathbb{R}}_{\infty} = \mathbb{R} \cup \left\{\hat{\infty}\right\}"),
        ("Equal(ProjectiveComplexNumbers, Union(CC, Set(UnsignedInfinity)))", r"\hat{\mathbb{C}}_{\infty} = \mathbb{C} \cup \left\{\hat{\infty}\right\}"),
        ("Equal(RealSingularityClosure, Union(RR, Set(UnsignedInfinity), RealSignedInfinities, Set(Undefined)))", r"\overline{\mathbb{R}}_{\text{Sing}} = \mathbb{R} \cup \left\{\hat{\infty}\right\} \cup \{\pm \infty\} \cup \left\{\mathfrak{u}\right\}"),
        ("Equal(ComplexSingularityClosure, Union(CC, Set(UnsignedInfinity), ComplexSignedInfinities, Set(Undefined)))", r"\overline{\mathbb{C}}_{\text{Sing}} = \mathbb{C} \cup \left\{\hat{\infty}\right\} \cup \{[e^{i \theta}] \infty\} \cup \left\{\mathfrak{u}\right\}"),
        ("Set(Mul(Infinity, Infinity), Mul(Mul(a, b), Infinity), Mul(NumberI, Infinity), Mul(Infinity, NumberI))", r"\left\{\infty \cdot \infty, a b \cdot \infty, i \cdot \infty, \infty i\right\}"),
        ("ArgMin(Add(f(x), g(x)), For(x, RR), NotEqual(x, 0))", r"\mathop{\operatorname{arg\,min}\,}\limits_{x \in \mathbb{R},\,x \ne 0} \left[f(x) + g(x)\right]"),
        ("List(ArgMin(f(x), For(x, S)), ArgMax(f(x), For(x, S)), ArgMin(f(x), For(x, S), P(x)), ArgMax(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{arg\,min}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,max}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,min}\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\operatorname{arg\,max}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
        ("List(Minimum(f(x), For(x, S)), Maximum(f(x), For(x, S)), Minimum(f(x), For(x, S), P(x)), Maximum(f(x), For(x, S), P(x)))", r"\left[\mathop{\min\,}\limits_{x \in S} f(x), \mathop{\max\,}\limits_{x \in S} f(x), \mathop{\min\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\max\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
        ("List(ArgMinUnique(f(x), For(x, S)), ArgMaxUnique(f(x), For(x, S)), ArgMinUnique(f(x), For(x, S), P(x)), ArgMaxUnique(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{arg\,min*}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,max*}\,}\limits_{x \in S} f(x), \mathop{\operatorname{arg\,min*}\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\operatorname{arg\,max*}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
        ("List(Infimum(f(x), For(x, S)), Supremum(f(x), For(x, S)), Infimum(f(x), For(x, S), P(x)), Supremum(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{inf}\,}\limits_{x \in S} f(x), \mathop{\operatorname{sup}\,}\limits_{x \in S} f(x), \mathop{\operatorname{inf}\,}\limits_{x \in S,\,P(x)} f(x), \mathop{\operatorname{sup}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
        ("List(Solutions(Q(x), For(x, S)), Zeros(f(x), For(x, S)), Solutions(Q(x), For(x, S), P(x)), Zeros(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{solutions}\,}\limits_{x \in S} Q(x), \mathop{\operatorname{zeros}\,}\limits_{x \in S} f(x), \mathop{\operatorname{solutions}\,}\limits_{x \in S,\,P(x)} Q(x), \mathop{\operatorname{zeros}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
        ("List(UniqueSolution(Q(x), For(x, S)), UniqueZero(f(x), For(x, S)), UniqueSolution(Q(x), For(x, S), P(x)), UniqueZero(f(x), For(x, S), P(x)))", r"\left[\mathop{\operatorname{solution*}\,}\limits_{x \in S} Q(x), \mathop{\operatorname{zero*}\,}\limits_{x \in S} f(x), \mathop{\operatorname{solution*}\,}\limits_{x \in S,\,P(x)} Q(x), \mathop{\operatorname{zero*}\,}\limits_{x \in S,\,P(x)} f(x)\right]"),
        ("Sum(f(n) + g(n), For(n, a, b))", r"\sum_{n=a}^{b} \left(f(n) + g(n)\right)"),
        ("Sum(f(n), For(n, ZZ))", r"\sum_{n  \in \mathbb{Z}} f(n)"),
        ("Sum(f(n), For(n, ZZ), NotEqual(n, 0))", r"\sum_{\textstyle{n  \in \mathbb{Z} \atop n \ne 0}} f(n)"),
        ("Sum(f(n), For(n, a, b), NotEqual(n, 0))", r"\sum_{\textstyle{n=a \atop n \ne 0}}^{b} f(n)"),
        ("Sum(f(n), For(n, a, b))", r"\sum_{n=a}^{b} f(n)"),
        ("Product(f(n) + g(n), For(n, a, b))", r"\prod_{n=a}^{b} \left(f(n) + g(n)\right)"),
        ("Product(f(n), For(n, NN))", r"\prod_{n  \in \mathbb{N}} f(n)"),
        ("Product(f(n), For(n, NN), NotEqual(g(n), 0))", r"\prod_{\textstyle{n  \in \mathbb{N} \atop g(n) \ne 0}} f(n)"),
        ("Product(f(n), For(n, a, b), NotEqual(n, 0))", r"\prod_{\textstyle{n=a \atop n \ne 0}}^{b} f(n)"),
        ("Product(f(n), For(n, a, b))", r"\prod_{n=a}^{b} f(n)"),
        ("Equal(Set(f(n), For(n, ZZ)), Union(Set(f(n), For(n, ZZ), IsEven(n)), Set(f(n), For(n, ZZ), IsOdd(n))))", r"\left\{ f(n) : n \in \mathbb{Z} \right\} = \left\{ f(n) : n \in \mathbb{Z}\,\mathbin{\operatorname{and}}\, n \text{ even} \right\} \cup \left\{ f(n) : n \in \mathbb{Z}\,\mathbin{\operatorname{and}}\, n \text{ odd} \right\}"),
        ("Equal(Primes, Set(p, For(p, NN), IsPrime(p)))", r"\mathbb{P} = \left\{ p : p \in \mathbb{N}\,\mathbin{\operatorname{and}}\, p \text{ prime} \right\}"),
        ("Equal(Sum(f(n), Element(n, ZZ)), Add(Sum(f(n), Element(n, ZZ), IsOdd(n)), Sum(f(n), Element(n, ZZ), IsEven(n))))", r"\sum_{n  \in \mathbb{Z}} f(n) = \sum_{\textstyle{n  \in \mathbb{Z} \atop n \text{ odd}}} f(n) + \sum_{\textstyle{n  \in \mathbb{Z} \atop n \text{ even}}} f(n)"),
        ("Set(DivisorSum(f(d), For(d, n)), DivisorSum(f(d), For(d, n), IsOdd(d)), DivisorSum(Add(f(d), g(d)), For(d, n)))", r"\left\{\sum_{d \mid n} f(d), \sum_{d \mid n,\, d \text{ odd}} f(d), \sum_{d \mid n} \left(f(d) + g(d)\right)\right\}"),
        ("Set(DivisorProduct(f(d), For(d, n)), DivisorProduct(f(d), For(d, n), IsOdd(d)), DivisorProduct(Add(f(d), g(d)), For(d, n)))", r"\left\{\prod_{d \mid n} f(d), \prod_{d \mid n,\, d \text{ odd}} f(d), \prod_{d \mid n} \left(f(d) + g(d)\right)\right\}"),
        ("Set(PrimeSum(f(p), For(p)), PrimeSum(f(p), For(p), NotElement(p, S)), PrimeProduct(f(p), For(p)), PrimeProduct(f(p), For(p), NotElement(p, S)))", r"\left\{\sum_{p} f(p), \sum_{p \notin S} f(p), \prod_{p} f(p), \prod_{p \notin S} f(p)\right\}"),
        ("Integral(f(x), For(x, -Infinity, Infinity))", r"\int_{-\infty}^{\infty} f(x) \, dx"),
        ("Integral(f(x), For(x, RR))", r"\int_{x \in \mathbb{R}} f(x) \, dx"),
        ("Integral(f(x) + g(x) / h(x), For(x, a, b))", r"\int_{a}^{b} \left(f(x) + \frac{g(x)}{h(x)}\right) \, dx"),
        ("Set(Derivative(f(x_), For(x_, x)), Derivative(f(x_), For(x_, Div(x, y))), Derivative(Gamma(x_), For(x_, 1)))", r"\left\{f'\!\left(x\right), f'\!\left(\frac{x}{y}\right), \Gamma'\!\left(1\right)\right\}"),
        ("Set(Derivative(f(x_), For(x_, x, 0)), Derivative(f(x_), For(x_, x, 1)), Derivative(f(x_), For(x_, x, 2)), Derivative(f(x_), For(x_, x, 3)), Derivative(f(x_), For(x_, x, 4)), Derivative(f(x_), For(x_, x, n)), Derivative(f(x_), For(x_, x, Add(Mul(2, n), 3))))", r"\left\{{f}^{(0)}\!\left(x\right), f'\!\left(x\right), f''\!\left(x\right), f'''\!\left(x\right), {f}^{(4)}\!\left(x\right), {f}^{(n)}\!\left(x\right), {f}^{(2 n + 3)}\!\left(x\right)\right\}"),
        ("Set(Derivative(f(Add(x, 1)), For(x, x)), Derivative(f(Add(x, 1)), For(x, x, 0)), Derivative(f(Add(x, 1)), For(x, x, 1)), Derivative(f(Add(x, 1)), For(x, x, n)))", r"\left\{\frac{d}{d x}\, f\!\left(x + 1\right), \frac{d^{0}}{{d x}^{0}}\, f\!\left(x + 1\right), \frac{d}{d x}\, f\!\left(x + 1\right), \frac{d^{n}}{{d x}^{n}}\, f\!\left(x + 1\right)\right\}"),
        ("Set(Derivative(Add(f(x), g(x)), For(x, Add(y, 3))), Derivative(Add(f(x), g(x)), For(x, Add(y, 3), 5)))", r"\left\{\left[\frac{d}{d x}\, \left[f(x) + g(x)\right] \right]_{x = y + 3}, \left[\frac{d^{5}}{{d x}^{5}}\, \left[f(x) + g(x)\right] \right]_{x = y + 3}\right\}"),
        ("Set(RealDerivative(f(x), For(x, 1)), ComplexDerivative(f(x), For(x, 1)), ComplexBranchDerivative(f(x), For(x, 1)), MeromorphicDerivative(f(x), For(x, 1)))", r"\left\{f'\!\left(1\right), f'\!\left(1\right), f'\!\left(1\right), f'\!\left(1\right)\right\}"),
        ("Set(Limit(f(x), For(x, a)), Limit(f(x), For(x, a), P(x)))", r"\left\{\lim_{x \to a} f(x), \lim_{x \to a,\,P(x)} f(x)\right\}"),
        ("Set(Limit(f(x), For(x, a)), RealLimit(f(x), For(x, a)), ComplexLimit(f(x), For(x, a)), MeromorphicLimit(f(x), For(x, a)))", r"\left\{\lim_{x \to a} f(x), \lim_{x \to a} f(x), \lim_{x \to a} f(x), \lim_{x \to a} f(x)\right\}"),
        ("Set(LeftLimit(f(x), For(x, 0)), RightLimit(f(x), For(x, 0)))", r"\left\{\lim_{x \to {0}^{-}} f(x), \lim_{x \to {0}^{+}} f(x)\right\}"),
        ("Set(SequenceLimit(f(n), For(n, Infinity)), SequenceLimitInferior(f(n), For(n, Infinity)), SequenceLimitSuperior(f(n), For(n, Infinity)))", r"\left\{\lim_{n \to \infty} f(n), \liminf_{n \to \infty} f(n), \limsup_{n \to \infty} f(n)\right\}"),
        ("Sub(Limit(Add(f(x), g(x)), For(x, a)), Limit(Sub(f(x), g(x)), For(x, a)))", r"\lim_{x \to a} \left[f(x) + g(x)\right] - \lim_{x \to a} \left[f(x) - g(x)\right]"),
        ("Divides(GCD(a, b), LCM(a, b))", r"\gcd(a, b) \mid \operatorname{lcm}(a, b)"),
        ("Set(Exp(x), Exp(Div(3, 2)), Exp(Add(Neg(Pow(x, 2)), x)), Exp(Abs(Im(z))), Exp(Div(3, Add(2, x))), Exp(Sin(x)))", r"\left\{e^{x}, e^{3 / 2}, e^{-{x}^{2} + x}, e^{\left|\operatorname{Im}(z)\right|}, \exp\!\left(\frac{3}{2 + x}\right), \exp\!\left(\sin(x)\right)\right\}"),
        ("Add(Sin(x), Cos(x), Tan(x), Cot(x), Sec(x), Csc(x))", r"\sin(x) + \cos(x) + \tan(x) + \cot(x) + \sec(x) + \csc(x)"),
        ("Add(Sinh(x), Cosh(x), Tanh(x), Coth(x), Sech(x), Csch(x))", r"\sinh(x) + \cosh(x) + \tanh(x) + \coth(x) + \operatorname{sech}(x) + \operatorname{csch}(x)"),
        ("Add(Asin(x), Acos(x), Atan(x), Acot(x), Asec(x), Acsc(x))", r"\operatorname{asin}(x) + \operatorname{acos}(x) + \operatorname{atan}(x) + \operatorname{acot}(x) + \operatorname{asec}(x) + \operatorname{acsc}(x)"),
        ("Add(Asinh(x), Acosh(x), Atanh(x), Acoth(x), Asech(x), Acsch(x))", r"\operatorname{asinh}(x) + \operatorname{acosh}(x) + \operatorname{atanh}(x) + \operatorname{acoth}(x) + \operatorname{asech}(x) + \operatorname{acsch}(x)"),
        ("Exp(Neg(Euler))", r"e^{-\gamma}"),
        ("Set(Re(z), Im(z), Atan2(y, x))", r"\left\{\operatorname{Re}(z), \operatorname{Im}(z), \operatorname{atan2}(y, x)\right\}"),
        ("Add(NumberE, GoldenRatio, CatalanConstant)", r"e + \varphi + G"),
        ("Add(Sinc(x), Pow(Sinc(x), 2))", r"\operatorname{sinc}(x) + \operatorname{sinc}^{2}\!\left(x\right)"),
        ("AGM(a, b)", r"\operatorname{agm}(a, b)"),
        ("And(Equal(LogBarnesG(z), Log(BarnesG(z))), Equal(LogGamma(z), Log(Gamma(z))))", r"\log G(z) = \log\!\left(G(z)\right) \;\mathbin{\operatorname{and}}\; \log \Gamma(z) = \log\!\left(\Gamma(z)\right)"),
        ("DirichletL(s, chi)", r"L(s, \chi)"),
        ("DirichletLambda(s, chi)", r"\Lambda(s, \chi)"),
        ("Implies(GeneralizedRiemannHypothesis, RiemannHypothesis)", r"\operatorname{GRH} \;\implies\; \operatorname{RH}"),
        ("Set(ModularJ(tau), ModularLambda(tau), JacobiTheta(n, z, tau))", r"\left\{j(\tau), \lambda(\tau), \theta_{n}\!\left(z, \tau\right)\right\}"),
        ("Set(WeierstrassP(z, tau), WeierstrassSigma(z, tau), WeierstrassZeta(z, tau))", r"\left\{\wp(z, \tau), \sigma(z, \tau), \zeta(z, \tau)\right\}"),
        ("Mul(ChebyshevT(n, x), ChebyshevU(n, x))", r"T_{n}\!\left(x\right) U_{n}\!\left(x\right)"),
        ("Add(FresnelC(z), FresnelS(z))", r"C(z) + S(z)"),
        ("Div(EisensteinE(Mul(2, n), tau), EisensteinG(Mul(2, n), tau))", r"\frac{E_{2 n}\!\left(\tau\right)}{G_{2 n}\!\left(\tau\right)}"),
        ("Equal(Div(IncompleteBeta(z, a, b), IncompleteBetaRegularized(z, a, b)), BetaFunction(a, b))", r"\frac{\mathrm{B}_{z}\!\left(a, b\right)}{I_{z}\!\left(a, b\right)} = \mathrm{B}(a, b)"),
        ("Set(PolyLog(s, z), HurwitzZeta(s, z), LerchPhi(z, s, a))", r"\left\{\operatorname{Li}_{s}\!\left(z\right), \zeta(s, z), \Phi(z, s, a)\right\}"),
        ("Equal(PartitionsP(n), Mul(Div(1, n), Sum(Mul(DivisorSigma(1, Sub(n, k)), PartitionsP(k)), For(k, 0, Sub(n, 1)))))", r"p(n) = \frac{1}{n} \sum_{k=0}^{n - 1} \sigma_{1}\!\left(n - k\right) p(k)"),
        ("MultiZetaValue(a, b, c)", r"\zeta(a, b, c)"),
        ("RiemannXi(s)", r"\xi(s)"),
        ("Mul(LiouvilleLambda(n), EulerPhi(n), MoebiusMu(n))", r"\lambda(n) \varphi(n) \mu(n)"),
        ("BetaFunction(a, b)", r"\mathrm{B}(a, b)"),
        ("PrimePi(x)", r"\pi(x)"),
        ("Equal(Min(a, b), Neg(Max(Neg(a), Neg(b))))", r"\min(a, b) = -\max\!\left(-a, -b\right)"),
        ("Equal(Arg(z), Div(Pi, 2))", r"\arg(z) = \frac{\pi}{2}"),
        ("NotEqual(Csgn(z), Sign(z))", r"\operatorname{csgn}(z) \ne \operatorname{sgn}(z)"),
        ("Add(Factorial(0), Factorial(1), Div(1, Factorial(-3)), Factorial(Div(1, 2)), Factorial(Factorial(n)), DoubleFactorial(n))", r"0! + 1! + \frac{1}{\left(-3\right)!} + \left(\frac{1}{2}\right)! + \left(n!\right)! + n!!"),
        ("List(Binomial(x, n), RisingFactorial(x, n), FallingFactorial(x, n), StirlingCycle(x, n), StirlingS1(x, n), StirlingS2(x, n))", r"\left[{x \choose n}, \left(x\right)_{n}, \left(x\right)^{\underline{n}}, \left[{x \atop n}\right], s\!\left(x, n\right), \left\{{x \atop n}\right\}\right]"),
        ("Add(BellNumber(5), BernoulliB(5), EulerE(5), Fibonacci(5), HarmonicNumber(5), Prime(5), RiemannZetaZero(5))", r"\operatorname{B}_{5} + B_{5} + E_{5} + F_{5} + H_{5} + p_{5} + \rho_{5}"),
        ("List(LegendreSymbol(p, q), JacobiSymbol(p, q), KroneckerSymbol(p, q))", r"\left[\left(\frac{p}{q}\right), \left(\frac{p}{q}\right), \left(\frac{p}{q}\right)\right]"),
        ("Add(ExpIntegralEi(x), ExpIntegralE(n, x), SinIntegral(x), SinhIntegral(x), CosIntegral(x), CoshIntegral(x), LogIntegral(x))", r"\operatorname{Ei}(x) + E_{n}\!\left(x\right) + \operatorname{Si}(x) + \operatorname{Shi}(x) + \operatorname{Ci}(x) + \operatorname{Chi}(x) + \operatorname{li}(x)"),
        ("Mul(BesselJ(nu, z), BesselI(nu, z), BesselY(nu, z), BesselK(nu, z))", r"J_{\nu}\!\left(z\right) I_{\nu}\!\left(z\right) Y_{\nu}\!\left(z\right) K_{\nu}\!\left(z\right)"),
        ("Equal(AiryAi(AiryAiZero(n)), AiryBi(AiryBiZero(n)), 0)", r"\operatorname{Ai}\!\left(a_{n}\right) = \operatorname{Bi}\!\left(b_{n}\right) = 0"),
        ("Equal(BesselJ(nu, BesselJZero(nu, n)), BesselY(nu, BesselYZero(nu, n)), 0)", r"J_{\nu}\!\left(j_{\nu, n}\right) = Y_{\nu}\!\left(y_{\nu, n}\right) = 0"),
        ("Equal(RiemannZeta(s), Mul(Mul(Mul(Mul(2, Pow(Mul(2, Pi), Sub(s, 1))), Sin(Div(Mul(Pi, s), 2))), Gamma(Sub(1, s))), RiemannZeta(Sub(1, s))))", r"\zeta(s) = 2 {\left(2 \pi\right)}^{s - 1} \sin\!\left(\frac{\pi s}{2}\right) \Gamma\!\left(1 - s\right) \zeta\!\left(1 - s\right)"),
        ("Pow(Div(Pow(DedekindEta(Mul(2, tau)), 2), Mul(DedekindEta(tau), DedekindEta(Mul(4, tau)))), 24)", r"{\left(\frac{\eta^{2}\!\left(2 \tau\right)}{\eta(\tau) \eta\!\left(4 \tau\right)}\right)}^{24}"),
        ("Mul(Mul(Erf(z), Erfc(z)), Erfi(z))", r"\operatorname{erf}(z) \operatorname{erfc}(z) \operatorname{erfi}(z)"),
        ("Mul(EllipticK(m), EllipticE(m), EllipticPi(n, m))", r"K(m) E(m) \Pi(n, m)"),
        ("Mul(IncompleteEllipticE(z, m), IncompleteEllipticF(z, m), IncompleteEllipticPi(n, z, m))", r"E(z, m) F(z, m) \Pi(n, z, m)"),
        ("Add(CarlsonRF(x, y, z), CarlsonRG(x, y, z), CarlsonRJ(x, y, z, w), CarlsonRD(x, y, z), CarlsonRC(x, y))", r"R_F(x, y, z) + R_G(x, y, z) + R_J(x, y, z, w) + R_D(x, y, z) + R_C(x, y)"),
        ("Mul(Hypergeometric0F1(b, z), Hypergeometric0F1Regularized(b, z))", r"\,{}_0F_1(b, z) \,{}_0{\textbf F}_1(b, z)"),
        ("Mul(Hypergeometric1F1(a, b, z), Hypergeometric1F1Regularized(a, b, z))", r"\,{}_1F_1(a, b, z) \,{}_1{\textbf F}_1(a, b, z)"),
        ("Hypergeometric2F0(a, b, z)", r"\,{}_2F_0(a, b, z)"),
        ("Mul(HypergeometricU(a, b, z), HypergeometricUStar(a, b, z))", r"U(a, b, z) U^{*}(a, b, z)"),
        ("Mul(Hypergeometric2F1(a, b, c, z), Hypergeometric2F1Regularized(a, b, c, z))", r"\,{}_2F_1(a, b, c, z) \,{}_2{\textbf F}_1(a, b, c, z)"),
        ("Mul(Hypergeometric1F2(a, b, c, z), Hypergeometric1F2Regularized(a, b, c, z))", r"\,{}_1F_2(a, b, c, z) \,{}_1{\textbf F}_2(a, b, c, z)"),
        ("Mul(Hypergeometric2F2(a, b, c, d, z), Hypergeometric2F2Regularized(a, b, c, d, z))", r"\,{}_2F_2(a, b, c, d, z) \,{}_2{\textbf F}_2(a, b, c, d, z)"),
        ("Mul(Hypergeometric3F2(a, b, c, d, e, z), Hypergeometric3F2Regularized(a, b, c, d, e, z))", r"\,{}_3F_2(a, b, c, d, e, z) \,{}_3{\textbf F}_2(a, b, c, d, e, z)"),
        ("(Hypergeometric2F1Regularized(Div(-1,4),Div(1,4),1/2, (x-1)/2)**2)", r"{\left(\,{}_2{\textbf F}_1\!\left(-\frac{1}{4}, \frac{1}{4}, \frac{1}{2}, \frac{x - 1}{2}\right)\right)}^{2}"),
        ("Matrix(List(List(a, b, c), List(d, e, f), List(g, h, 0)))", r"\displaystyle{\begin{pmatrix}a & b & c \\d & e & f \\g & h & 0\end{pmatrix}}"),
        ("Matrix2x2(a, b, c, d)", r"\displaystyle{\begin{pmatrix}a & b \\ c & d\end{pmatrix}}"),
        ("Set(RowMatrix(), RowMatrix(a), RowMatrix(a, b), RowMatrix(a, b, c), RowMatrix(a, b, c, d))", r"\left\{\displaystyle{\begin{pmatrix}\end{pmatrix}}, \displaystyle{\begin{pmatrix}a\end{pmatrix}}, \displaystyle{\begin{pmatrix}a & b\end{pmatrix}}, \displaystyle{\begin{pmatrix}a & b & c\end{pmatrix}}, \displaystyle{\begin{pmatrix}a & b & c & d\end{pmatrix}}\right\}"),
        ("Set(ColumnMatrix(), ColumnMatrix(a), ColumnMatrix(a, b), ColumnMatrix(a, b, c), ColumnMatrix(a, b, c, d))", r"\left\{\displaystyle{\begin{pmatrix}\end{pmatrix}}, \displaystyle{\begin{pmatrix}a\end{pmatrix}}, \displaystyle{\begin{pmatrix}a \\ b\end{pmatrix}}, \displaystyle{\begin{pmatrix}a \\ b \\ c\end{pmatrix}}, \displaystyle{\begin{pmatrix}a \\ b \\ c \\ d\end{pmatrix}}\right\}"),
        ("Set(DiagonalMatrix(), DiagonalMatrix(a), DiagonalMatrix(a, b), DiagonalMatrix(a, b, c), DiagonalMatrix(a, b, c, d))", r"\left\{\displaystyle{\begin{pmatrix}\end{pmatrix}}, \displaystyle{\begin{pmatrix}a\end{pmatrix}}, \displaystyle{\begin{pmatrix}a &  \\  & b\end{pmatrix}}, \displaystyle{\begin{pmatrix}a &  &  \\  & b &  \\  &  & c\end{pmatrix}}, \displaystyle{\begin{pmatrix}a &  &  &  \\  & b &  &  \\  &  & c &  \\  &  &  & d\end{pmatrix}}\right\}"),
        ("Matrix(c_(m, n), For(m, 1, N), For(n, 1, 10))", r"\displaystyle{\begin{pmatrix} c_{1, 1} & c_{1, 2} & \cdots & c_{1, 10} \\ c_{2, 1} & c_{2, 2} & \cdots & c_{2, 10} \\ \vdots & \vdots & \ddots & \vdots \\ c_{N, 1} & c_{N, 2} & \cdots & c_{N, 10} \end{pmatrix}}"),
        ("Matrix(ShowExpandedNormalForm(Div(1, Sub(Add(m, n), 1))), For(m, 1, 10), For(n, 1, 10))", r"\displaystyle{\begin{pmatrix} 1 & \frac{1}{2} & \cdots & \frac{1}{10} \\ \frac{1}{2} & \frac{1}{3} & \cdots & \frac{1}{11} \\ \vdots & \vdots & \ddots & \vdots \\ \frac{1}{10} & \frac{1}{11} & \cdots & \frac{1}{19} \end{pmatrix}}"),
        ("Add(ZeroMatrix(2), IdentityMatrix(2), HilbertMatrix(2))", r"0_{2} + I_{2} + H_{2}"),
        ("Set(SpecialLinearGroup(n, ZZ), GeneralLinearGroup(n, ZZ))", r"\left\{\operatorname{SL}_{n}\!\left(\mathbb{Z}\right), \operatorname{GL}_{n}\!\left(\mathbb{Z}\right)\right\}"),
        ("Equal(One(QQ), 1)", r"1_{\mathbb{Q}} = 1"),
        ("Equal(Zero(QQ), 0)", r"0_{\mathbb{Q}} = 0"),
        ("List(Polynomials(QQ, x), Polynomials(QQ, x, y), Polynomials(QQ, Tuple()), Polynomials(QQ, Tuple(x)), Polynomials(QQ, Tuple(x, y)))", r"\left[\mathbb{Q}[x], \mathbb{Q}[x, y], \mathbb{Q}[], \mathbb{Q}[x], \mathbb{Q}[x, y]\right]"),
        ("List(Polynomials(QQ, x), PolynomialFractions(QQ, x), FormalPowerSeries(QQ, x), FormalLaurentSeries(QQ, x), FormalPuiseuxSeries(QQ, x))", r"\left[\mathbb{Q}[x], \mathbb{Q}(x), \mathbb{Q}[[x]], \mathbb{Q}(\!(x)\!), \mathbb{Q}\!\left\langle\!\left\langle x \right\rangle\!\right\rangle\right]"),
        ("Set(IntegersGreaterEqual(0), IntegersGreaterEqual(n), IntegersLessEqual(0), IntegersLessEqual(n))", r"\left\{\mathbb{Z}_{\ge 0}, \mathbb{Z}_{\ge n}, \{0, -1, \ldots\}, \mathbb{Z}_{\le n}\right\}"),
        ("List(Range(a, b), Range(1, b), Range(-3, 5))", r"\left[\{a, a + 1, \ldots, b\}, \{1, 2, \ldots, b\}, \{-3, -2, \ldots, 5\}\right]"),
        ("CongruentMod(f(n), 0, p)", r"f(n) \equiv 0 \pmod {p }"),
        ("PrimitiveReducedPositiveIntegralBinaryQuadraticForms(D)", r"\mathcal{Q}^{*}_{D}"),
        ("Set(EllipticRootE(1, tau), EllipticRootE(2, tau), EllipticRootE(3, tau))", r"\left\{e_{1}\!\left(\tau\right), e_{2}\!\left(\tau\right), e_{3}\!\left(\tau\right)\right\}"),
        ("GaussSum(n, chi)", r"G_{n}\!\left(\chi\right)"),
        ("Set(GlaisherConstant, KhinchinConstant)", r"\left\{A, K\right\}"),
        ('Decimal("0.3141")', r"0.3141"),
        ('Decimal("0.3141e-27")', r"0.3141 \cdot 10^{-27}"),
        ("Set(DigammaFunction(z), DigammaFunction(z, 1), DigammaFunction(z, n))", r"\left\{\psi(z), \psi'\!\left(z\right), {\psi}^{(n)}\!\left(z\right)\right\}"),
        ("Set(AiryAi(z, 1), AiryAi(z, 2), AiryBi(z, n), BesselJ(n, z, 1), BesselY(n, z, 2), BesselK(n, z, Add(Mul(3, r), 1)))", r"\left\{\operatorname{Ai}'\!\left(z\right), \operatorname{Ai}''\!\left(z\right), {\operatorname{Bi}}^{(n)}\!\left(z\right), J'_{n}\!\left(z\right), Y''_{n}\!\left(z\right), {K}^{(3 r + 1)}_{n}\!\left(z\right)\right\}"),
        ("Set(HankelH1(n, z), HankelH2(n, z))", r"\left\{H^{(1)}_{n}\!\left(z\right), H^{(2)}_{n}\!\left(z\right)\right\}"),
        ("Element(DirichletCharacter(q, k), PrimitiveDirichletCharacters(q))", r"\chi_{q \, . \, k} \in G^{\text{Primitive}}_{q}"),
        ("DirichletCharacter(q, k, n)", r"\chi_{q \, . \, k}(n)"),
        ("JacobiTheta(3, z, tau, 2)", r"\theta''_{3}\!\left(z, \tau\right)"),
        ("Set(RiemannZeta(s, 1), RiemannZeta(s, r))", r"\left\{\zeta'\!\left(s\right), {\zeta}^{(r)}\!\left(s\right)\right\}"),
        ("Subscript(x, y)", r"{x}_{y}"),
        ("Set(BernsteinEllipse(r), UnitCircle, PSL2Z, AGMSequence(n, a, b), CarlsonHypergeometricR, CarlsonHypergeometricT)", r"\left\{\mathcal{E}_{r}, \mathbb{T}, \operatorname{PSL}_2(\mathbb{Z}), \operatorname{agm}_{n}\!\left(a, b\right), R, T\right\}"),
        ("Mul(Mul(Pow(Fibonacci(n), 2), Pow(x_(a), 2)), Pow(alpha_(n), 2))", r"F_{n}^{2} x_{a}^{2} \alpha_{n}^{2}"),
        ("Set(Derivative(ChebyshevT(n, x_), For(x_, x)), Derivative(ChebyshevT(n, x_), For(x_, x, 2)), Derivative(ChebyshevT(n, x_), For(x_, x, 4)))", r"\left\{T'_{n}\!\left(x\right), T''_{n}\!\left(x\right), {T}^{(4)}_{n}\!\left(x\right)\right\}"),
        ("Poles(Gamma(z), For(z, CC))", r"\mathop{\operatorname{poles}\,}\limits_{z \in \mathbb{C}} \Gamma(z)"),
        ("Equal(Item(Tuple(a, b, c), 2), b)", r"{\left(a, b, c\right)}_{2} = b"),
        ("Set(Tuple(n, For(n, a, b)), List(Pow(Neg(n), 2), For(n, 1, 100)), Set(f_(n), For(n, 0, N)))", r"\left\{\left(a, a + 1, \ldots, b\right), \left[{\left(-1\right)}^{2}, {\left(-2\right)}^{2}, \ldots, {\left(-100\right)}^{2}\right], \left\{f_{0}, f_{1}, \ldots, f_{N}\right\}\right\}"),
        ("Set(Tuple(1, 0, Repeat(3, N)), Tuple(1, 0, Repeat(1, 2, 3, N)))", r"\left\{\left(1, 0, \underbrace{3, \ldots, 3}_{N \text{ times}}\right), \left(1, 0, \underbrace{1, 2, 3, \ldots, 1, 2, 3}_{\left(1, 2, 3\right) \; N \text{ times}}\right)\right\}"),
        ("Tuple(Sub(A, 2), Sub(A, 1), Step(n, For(n, A, B)), Add(B, 1), Add(B, 2))", r"\left(A - 2, A - 1, A, A + 1, \ldots, B, B + 1, B + 2\right)"),
        ("Lattice(1, tau)", r"\Lambda_{(1, \tau)}"),
        ("DiscreteLog(n, 2, q)", r"(\epsilon : {2}^{\epsilon} \equiv n \text{ mod }q)"),
        ("AsymptoticTo(f(n), g(n), n, Infinity)", r"f(n) \sim g(n), \; n \to \infty"),
        ("Set(CoulombF(l, eta, z), CoulombG(l, eta, z), CoulombH(1, l, eta, z), CoulombH(-1, l, eta, z), CoulombH(omega, l, eta, z))", r"\left\{F_{l,\eta}(z), G_{l,\eta}(z), H^{+}_{l,\eta}(z), H^{-}_{l,\eta}(z), H^{\omega}_{l,\eta}(z)\right\}"),
        ("Set(Matrices(CC, n), Matrices(CC, n, m))", r"\left\{\operatorname{M}_{n}(\mathbb{C}), \operatorname{M}_{n \times m}(\mathbb{C})\right\}"),
        ('Set(SloaneA(40, n), SloaneA(12345, n), SloaneA("A553322", n))', r"\left\{\text{A000040}\!\left(n\right), \text{A012345}\!\left(n\right), \text{A553322}\!\left(n\right)\right\}"),
        ('And(EqualAndElement(x, Pi, RR), EqualNearestDecimal(Pi, Decimal("3.14"), 3))', r"x = \pi \in \mathbb{R} \;\mathbin{\operatorname{and}}\; \pi = 3.14 \;\, {\scriptstyle (\text{nearest } 3 \text{ digits})}"),
        ("Less(0, Same(Div(1, Sqrt(2)), Div(Sqrt(2), 2)), 1)", r"0 < \frac{1}{\sqrt{2}} = \frac{\sqrt{2}}{2} < 1"),
        ("Fun(x, Pow(x, 2))", r"x \mapsto {x}^{2}"),
        ("GeneralizedBernoulliB(n, chi)", r"B_{n, \chi}"),
        ("HurwitzZeta(s, a, 2)", r"\zeta''\!\left(s, a\right)"),
        ("Set(StieltjesGamma(n), StieltjesGamma(n, a))", r"\left\{\gamma_{n}, \gamma_{n}\!\left(a\right)\right\}"),
        ("And(IsHolomorphicOn(f(z), For(z, CC)), IsMeromorphicOn(g(z), For(z, CC)))", r"f(z) \text{ is holomorphic on } z \in \mathbb{C} \;\mathbin{\operatorname{and}}\; g(z) \text{ is meromorphic on } z \in \mathbb{C}"),
        ("Set(StirlingSeriesRemainder(N, z), LogBarnesGRemainder(N, z))", r"\left\{R_{N}\!\left(z\right), R_{N}\!\left(z\right)\right\}"),
        ("AnalyticContinuation(f(z), For(z, a, b))", r"\mathop{\text{Continuation}}\limits_{\displaystyle{z: a \rightsquigarrow b}} \, f(z)"),
        ("AnalyticContinuation(f(z), For(z, CurvePath(g(t), For(t, a, b))))", r"\mathop{\text{Continuation}}\limits_{\displaystyle{z: \left(g(t),\, t : a \rightsquigarrow b\right)}} \, f(z)"),
        ("BernoulliPolynomial(n, x)", r"B_{n}\!\left(x\right)"),
        ("Call(f, x)", r"f(x)"),
        ("CallIndeterminate(f, x, v)", r"f(v)"),
        ("CartesianProduct(ZZ, QQ)", r"\mathbb{Z} \times \mathbb{Q}"),
        ("CartesianPower(RR, 3)", r"{\mathbb{R}}^{3}"),
        ("Characteristic(R)", r"\operatorname{char}(R)"),
        ("Coefficient(f, x, 2)", r"[{x}^{2}] f"),
        ("ComplexZeroMultiplicity(f(z), For(z, a))", r"\mathop{\operatorname{ord}}\limits_{z=a} f(z)"),
        ("Residue(f(z), For(z, a))", r"\mathop{\operatorname{res}}\limits_{z=a} f(z)"),
        ("ConreyGenerator(q)", r"g_{q}"),
        ("CoulombC(l, eta)", r"C_{l}\!\left(\eta\right)"),
        ("CoulombSigma(l, eta)", r"\sigma_{l}\!\left(\eta\right)"),
        ("Cyclotomic(n, x)", r"\Phi_{n}\!\left(x\right)"),
        ("DedekindEtaEpsilon(a, b, c, d)", r"\varepsilon(a, b, c, d)"),
        ("DedekindSum(a, b)", r"s(a, b)"),
        ("Where(f(Add(x, 1)), Def(f(t), Div(1, t)))", r"f\!\left(x + 1\right)\; \text{ where } f(t) = \frac{1}{t}"),
        ("Det(Matrix2x2(a, b, c, d))", r"\operatorname{det} \displaystyle{\begin{pmatrix}a & b \\ c & d\end{pmatrix}}"),
        ("Det(A)", r"\operatorname{det}(A)"),
        ("f(a, b, Ellipsis, z)", r"f(a, b, \ldots, z)"),
        ("Equal(DirichletL(DirichletLZero(n, chi), chi), 0)", r"L\!\left(\rho_{n, \chi}, \chi\right) = 0"),
        ("DirichletGroup(q)", r"G_{q}"),
        ("Set(ModularGroupFundamentalDomain, ModularLambdaFundamentalDomain)", r"\left\{\mathcal{F}, \mathcal{F}_{\lambda}\right\}"),
        ("Path(a, b, c, d)", r"a \rightsquigarrow b \rightsquigarrow c \rightsquigarrow d"),
        ("PolynomialDegree(f)", r"\deg(f)"),
        ("RootOfUnity(5)", r"\zeta_{5}"),
        ("SL2Z", r"\operatorname{SL}_2(\mathbb{Z})"),
        ("UpperHalfPlane", r"\mathbb{H}"),
        ("XGCD(m, n)", r"\operatorname{xgcd}(m, n)"),
        ("JacobiThetaQ(3, z, q)", r"\theta_{3}\!\left(z, q\right)"),
        ("KeiperLiLambda(n)", r"\lambda_{n}"),
        ("Set(LambertW(z), LambertW(z, n), LambertW(z, n, 1), LambertW(z, n, r))", r"\left\{W(z), W_{n}(z), W'_{n}(z), {W}^{(r)}_{n}(z)\right\}"),
        ("LandauG(n)", r"g(n)"),
        ("SquaresR(k, n)", r"r_{k}\!\left(n\right)"),
        ("ModularGroupAction(gamma, tau)", r"\gamma \circ \tau"),
        ("Mod(n, p)", r"n \bmod p"),
        ("LessEqual(0, Step(f(n), For(n, a, b)), 1)", r"0 \le f(a) \le f\!\left(a + 1\right) \le \ldots \le f(b) \le 1"),
        ("LessEqual(0, Step(f(n), For(n, 1, b)), 1)", r"0 \le f(1) \le f(2) \le \ldots \le f(b) \le 1"),
        ("Add(LowerGamma(s, z), UpperGamma(s, z))", r"\gamma(s, z) + \Gamma(s, z)"),
        ("DigammaFunctionZero(n)", r"x_{n}"),
        ("Set(EulerPolynomial(n, x), HermiteH(n, x), HilbertClassPolynomial(n, x))", r"\left\{E_{n}\!\left(x\right), H_{n}\!\left(x\right), H_{n}\!\left(x\right)\right\}"),
        ("JacobiThetaEpsilon(j, a, b, c, d)", r"\varepsilon_{j}\!\left(a, b, c, d\right)"),
        ("JacobiThetaPermutation(j, a, b, c, d)", r"S_{j}\!\left(a, b, c, d\right)"),
        ("SymmetricPolynomial(k, List(X_(t), For(t, 1, n)))", r"e_{k}\!\left(\left[X_{1}, X_{2}, \ldots, X_{n}\right]\right)"),
        ("Intersection(A, B)", r"A \cap B"),
        ("RealAlgebraicNumbers", r"\overline{\mathbb{Q}}_{\mathbb{R}}"),
    ]

    def test_latex(fexpr):
        namespace = fexpr.inject_dict()
        for formula, expected in latex_test_cases:
            expr = eval(formula, globals(), namespace)
            latex = expr.latex()
            if latex != expected:
                raise AssertionError("%s:  got '%s', expected '%s'" % (formula, latex, expected))

    test_latex(fexpr)

    def latex_report(fexpr):
        namespace = fexpr.inject_dict()
        formulas = [eval(formula, globals(), namespace) for formula, expected in latex_test_cases]
        from os.path import expanduser
        from time import time
        fp = open(expanduser("~/Desktop/latex_report.html"), "w")
        fp.write("""
    <!DOCTYPE html>
    <html>
    <head>
    <title>fexpr to LaTeX test sheet</title>
    <meta http-equiv="Content-Type" content="text/html;charset=utf-8" >
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
    tt { padding: 0.1em; background-color: #f8f8f8; border:1px solid #eee; }
    table { border-collapse:collapse; margin: 1em; }
    table, th, td { border: 1px solid #aaa; }
    th, td { padding:0.1em 0.3em 0.1em 0.3em; }
    table { width: 95%; }
    .katex { font-size: 1.1em !important; } 
    .katex-display { margin:0.1em; padding:0.1em; }
    </style>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css" integrity="sha384-AfEj0r4/OFrOo5t7NnNe46zW/tFgW6x/bCJG8FqQCEo3+Aro6EYUG4+cU+KJWu/X" crossorigin="anonymous">
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.js" integrity="sha384-g7c+Jr9ZivxKLnZTDUhnkOnsh30B4H0rpLUpJ4jAIKs4fnJI+sEnkvrMWph2EDg4" crossorigin="anonymous"></script>
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/contrib/auto-render.min.js" integrity="sha384-mll67QQFJfxn0IYznZYonOWZ644AWYC+Pt2cHqMaRhXVrursRwvLnLaebdGIlYNa" crossorigin="anonymous"
        onload="renderMathInElement(document.body);"></script>
    <script>
      document.addEventListener("DOMContentLoaded", function() {
          renderMathInElement(document.body, {
              delimiters: [
                {left: "$$", right: "$$", display: true},
                {left: "$", right: "$", display: false}
              ]
          });
      });
    </script>
    </head>
    <body>
    """)
        output = [formula.latex() for formula in formulas]
        one_big = fexpr("BigLatex")(*formulas)
        t1 = time()
        one_big_latex = one_big.latex()
        t2 = time()
        fp.write("""<h1>fexpr to LaTeX test sheet</h1>""")
        fp.write("""<p>Converted %i formulas (%i leaves, %i bytes) to LaTeX in %f seconds.</p>""" % (len(formulas), one_big.num_leaves(), one_big.size_bytes(), (t2-t1)))
        fp.write("""<table>""")
        fp.write("""<tr><th>fexpr</th> <th>Generated LaTeX</th> <th>KaTeX display</th>""")
        for formula, latex in zip(formulas, output):
            fp.write("""<tr>""")
            fp.write("""<td><tt>%s</tt></td>""" % formula)
            fp.write("""<td><tt>%s</tt></td>""" % latex)
            fp.write("""<td>$$%s$$</td>""" % latex)
            fp.write("""</tr>""")
        fp.write("""</table>""")

        fp.write("""<br/><p>Untested builtins:</p> <p><tt>""")
        s = str(one_big)
        for c in '-+()_,"':
            s = s.replace(c, " ")
        used = set(s.split())
        builtins = [name.strip("_") for name in fexpr.builtins()]
        unused = [name for name in builtins if name not in used]
        for name in unused:
            fp.write(name)
            fp.write(" ")
        fp.write("""</tt></p>""")
        fp.write("""</body></html>""")
        fp.close()

    # latex_report(fexpr)



def test_ca_trigonometric():

    def expect_not_implemented(f):
        try:
            v = f()
            assert v
        except NotImplementedError:
            return
        raise AssertionError

    sqrt = CC_ca.sqrt
    pi = CC_ca.pi()
    xsin = CC_ca.sin
    xcos = CC_ca.cos
    xtan = CC_ca.tan

    a = 1+sqrt(2)
    b = 2+sqrt(2)

    assert xsin(a)**2 + xcos(a)**2 == 1
    assert xsin(-a)**2 + xcos(a)**2 == 1
    assert xsin(a) == -xsin(-a)
    assert xcos(a) == xcos(-a)
    assert xtan(a) == -xtan(-a)
    assert xsin(a+2*pi) == xsin(a)
    assert xcos(a+2*pi) == xcos(a)
    assert xtan(a+pi) == xtan(a)
    assert xsin(a+pi) == -xsin(a)
    assert xcos(a+pi) == -xcos(a)
    assert xtan(a+pi/2) == -1/xtan(a)
    assert xsin(a+pi/2) == xcos(a)
    assert xcos(a+pi/2) == -xsin(a)
    assert xsin(a-pi/2) == -xcos(a)
    assert xcos(a-pi/2) == xsin(a)
    assert xtan(a+pi/4) == (xtan(a)+1)/(1-xtan(a))
    assert xtan(a-pi/4) == (xtan(a)-1)/(1+xtan(a))
    assert xsin(pi/2 - a) == xcos(a)
    assert xcos(pi/2 - a) == xsin(a)
    assert xtan(pi/2 - a) == 1 / xtan(a)
    assert xsin(pi - a) == xsin(a)
    assert xcos(pi - a) == -xcos(a)
    assert xtan(pi - a) == -xtan(a)
    assert xsin(2*pi - a) == -xsin(a)
    assert xcos(2*pi - a) == xcos(a)
    assert xsin(a+b) == xsin(a)*xcos(b) + xcos(a)*xsin(b)
    assert xsin(a-b) == xsin(a)*xcos(b) - xcos(a)*xsin(b)
    assert xcos(a+b) == xcos(a)*xcos(b) - xsin(a)*xsin(b)
    assert xcos(a-b) == xcos(a)*xcos(b) + xsin(a)*xsin(b)
    assert xtan(a+b) == (xtan(a)+xtan(b)) / (1 - xtan(a)*xtan(b))
    assert xtan(a-b) == (xtan(a)-xtan(b)) / (1 + xtan(a)*xtan(b))
    assert xsin(2*a) == 2*xsin(a)*xcos(a)
    assert xsin(2*a) == 2*xtan(a)/(1+xtan(a)**2)
    assert xcos(2*a) == xcos(a)**2 - xsin(a)**2
    assert xcos(2*a) == 2*xcos(a)**2 - 1
    assert xcos(2*a) == 1 - 2*xsin(a)**2
    assert xcos(2*a) == (1 - xtan(a)**2) / (1 + xtan(a)**2)
    assert xtan(2*a) == (2*xtan(a)) / (1 - xtan(a)**2)

    assert xsin(3*a) == 3*xsin(a) - 4*xsin(a)**3
    assert xcos(3*a) == 4*xcos(a)**3 - 3*xcos(a)
    assert xtan(3*a) == (3*xtan(a) - xtan(a)**3) / (1 - 3*xtan(a)**2)
    assert xsin(a/2)**2 == (1-xcos(a))/2
    assert xcos(a/2)**2 == (1+xcos(a))/2
    assert xtan((a-b)/2) == (xsin(a) - xsin(b)) / (xcos(a) + xcos(b))

    assert 2*xcos(a)*xcos(b) == xcos(a-b) + xcos(a+b)
    assert 2*xsin(a)*xsin(b) == xcos(a-b) - xcos(a+b)
    assert 2*xsin(a)*xcos(b) == xsin(a+b) + xsin(a-b)
    assert 2*xcos(a)*xsin(b) == xsin(a+b) - xsin(a-b)

    assert xsin(a) + xsin(b) == 2*xsin((a+b)/2)*xcos((a-b)/2)
    assert xsin(a) - xsin(b) == 2*xsin((a-b)/2)*xcos((a+b)/2)

    assert xcos(a) + xcos(b) == 2*xcos((a+b)/2)*xcos((a-b)/2)
    assert xcos(a) - xcos(b) == -2*xsin((a+b)/2)*xsin((a-b)/2)
    assert xsin(a) == sqrt(1 - xcos(a)**2)

    for N in range(1,17):
        assert sum(xcos(n*a) for n in range(1,N+1)) == xsin((N+ca(1)/2)*a)/(2*xsin(a/2)) - ca(1)/2

    assert xcos(a) == -sqrt(1 - xsin(a)**2)
    assert xsin(a/2) == sqrt((1-xcos(a))/2)

    expect_not_implemented(lambda: xsin(3*a) == 4*xsin(a)*xsin(pi/3-a)*xsin(pi/3+a))
    assert xtan((a+b)/2) == (xsin(a) + xsin(b)) / (xcos(a) + xcos(b))
    assert xtan((a+b)/2) == (xsin(a) + xsin(b)) / (xcos(a) + xcos(b))
    assert xtan(a)*xtan(b) == ((xcos(a-b)-xcos(a+b))/(xcos(a-b)+xcos(a+b)))


def test_arb():
    a = arb(2.5)
    assert a  == arb("2.5")
    b = acb(2.5)
    assert a == b
    c = acb(2.5+1j)
    assert c == b + 1j
    assert raises(lambda: arb(2.5+1j), ValueError)
    assert acb(3+1j) == acb(ZZi(3+1j))
    assert arb(ZZi(3)) == 3
    assert raises(lambda: arb(ZZi(2.5+1j)), ValueError)

def test_vec():
    a = VecZZ([1,2,3])
    b = VecQQ([2,3,4])
    assert a[0] == 1
    assert a[2] == 3
    assert raises(lambda: a[-1], IndexError)
    assert raises(lambda: a[3], IndexError)
    assert a + a == VecZZ([2,4,6])
    assert a + b == VecQQ([3,5,7])
    assert b + a == VecQQ([3,5,7])
    assert a + ZZ(1) == VecZZ([2,3,4])
    assert ZZ(1) + a == VecZZ([2,3,4])
    assert b + ZZ(1) == VecQQ([3,4,5])
    assert ZZ(1) + b == VecQQ([3,4,5])
    assert b ** -5 == 1 / b ** 5
    assert raises(lambda: Vec(ZZi).i(), ValueError)
    i = ZZi.i()
    V = Vec(ZZi,3)
    assert V.i()  == V([i,i,1j])

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

    assert ZZp64(QQ(1) / 3) * 3 == ZZp64(1)
    assert ZZp64(QQ(1)) ** (QQ(1) / 2) == 1
    assert ZZp64(QQ(5)) ** (QQ(5)) == 3125
    assert ZZp32(10001).sqrt() ** 2 == 10001

    assert abs(VecZZ([-3,2,5])) == [3, 2, 5]

    b, t = PolynomialRing(PowerSeriesModRing(ZZ, 6, var="b"), "t").gens(recursive=True)
    assert (5+2*b+3*t)**5 / (5+2*b+3*t)**5 == 1

    assert CCser(1+ZZser.gen()) == 1 + RRser.gen()

def test_float():
    assert RF(5).mul_2exp(-1) == RF(2.5)
    assert CF(2+3j).mul_2exp(-1) == CF(1+1.5j)

def test_special():
    a = ZZ.fib_vec(100)
    for i in range(100):
        assert ZZ.fib(i) == a[i]
    F = FiniteField_fq(17, 1)
    for i in range(-10,10):
        assert QQ.fib(i) == QQ.fib(i-1) + QQ.fib(i-2)
        assert F.fib(i) == F.fib(i-1) + F.fib(i-2)

def test_mpoly():
    ZZxyz = PolynomialRing_fmpz_mpoly(3)
    x, y, z = ZZxyz.gens()
    f = (-72 * (1+x)**2 * (y+z+1))
    c, fac, exp = f.factor()
    assert c == -72
    assert ((fac, exp) == ([1+x, y+z+1], [2, 1])) or ((fac, exp) == ([y+z+1, 1+x], [1, 2]))
    assert f.gcd(-100-100*x) == 4+4*x

    assert str(PolynomialRing_fmpz_mpoly(2).gens()) == '[x1, x2]'
    assert str(PolynomialRing_fmpz_mpoly(2, ["a", "b"]).gens()) == '[a, b]'
    assert str(PolynomialRing_gr_mpoly(ZZi, 2).gens()) == '[x1, x2]'

    I, x, y, z = PolynomialRing_gr_mpoly(ZZi, 3, ["x", "y", "z"]).gens(recursive=True)
    assert str(x) == "x"
    assert str(y) == "y"
    assert str(z) == "z"
    assert str(x-y) == "x - y"
    assert str(x+2*y) == "x + 2*y"
    assert str(x-2*y) == "x - 2*y"
    assert str(-x) == "-x"
    assert str(-3*x) == "-3*x"
    assert str(x+1) == "x + 1"
    assert str(x-1) == "x - 1"
    assert str(x+2) == "x + 2"
    assert str(x-2) == "x - 2"
    assert str(x*0) == "0"
    assert str(x**0) == "1"
    assert str(-x**0) == "-1"
    assert str(-2*x**0) == "-2"
    assert str(x*y*z) == "x*y*z"
    assert str(x*y**2*z) == "x*y^2*z"
    assert str(3*x*y**2*z) == "3*x*y^2*z"
    assert str((1+I)*x + I*y) == "(1+I)*x + I*y"
    assert str((1+I)*x - I*y) == "(1+I)*x - I*y"
    assert str(x+1+I) == "x + (1+I)"

    x, y = PolynomialRing_gr_mpoly(ZZx, 1, ["y"]).gens(recursive=True)
    assert str(x+1) == "(x+1)"
    assert str((x+1)*y) == "(x+1)*y"
    assert str((x+1)*y + (x+2)) == "(x+1)*y + (x+2)"
    assert str((x+1)*y**2 + (x+2)*y)
    assert str((x+1)*y**2 - (x+2)*y) == "(x+1)*y^2 + (-x-2)*y"

    assert str(FiniteField_fq(3, 2, "c").gen()) == "c"
    assert str(FiniteField_fq_nmod(3, 2, "d").gen()) == "d"
    assert str(FiniteField_fq_zech(3, 2, "e").gen()) == "e^1"

def test_mpoly_q():
    assert str(FractionField_fmpz_mpoly_q(2).gens()) == '[x1, x2]'
    assert str(FractionField_fmpz_mpoly_q(2, ["a", "b"]).gens()) == '[a, b]'

def test_set_str():
    assert RR("1/4") == RR(1)/4
    v = RR("pi ^ (1 / 2)")
    assert 1.77 < v < 1.78
    assert CC("1/2 + i/4") == 0.5+0.25j
    assert CC("(-1)^(1/2)") == 1j
    assert CF("(-1)^(1/2)") == 1j
    assert abs(RF("2^(1/2)") ** 2) - 2 < 1e-14
    assert CC_ca("i*pi/2").exp() == CC_ca.i()

    assert ZZ("1 + 2^10") == 1025
    x = ZZx.gen(); R = PolynomialRing_gr_mpoly(NumberField(x**3+x+1), 3, ["x", "y", "z"])
    a, x, y, z = R.gens(recursive=True)
    assert R("((a-x+1)*y + (a^2+1)*y^2 + z)^2") == ((a-x+1)*y + (a**2+1)*y**2 + z)**2
    x = ZZx.gen()
    R = NumberField(x**2+1, "b")
    assert R("b-1") == R.gen()-1

    R = FractionField_fmpz_mpoly_q(2, ["x", "y"])
    x, y = R.gens()
    assert R("(4+4*x-y*(-4))^2 / (1+x+y) / 16") == 1+x+y

    assert RRx("1 +/- 0") == RR(1)

def test_qqbar_roots():
    for R in [ZZ, QQ, ZZi, QQbar, AA, QQbar_ca, AA_ca, RR_ca, CC_ca]:
        Rx = PolynomialRing(R)
        assert Rx([-2,0,1]).roots(domain=AA) == ([AA(2).sqrt(), -AA(2).sqrt()], [1, 1])
        assert Rx([2,0,1]).roots(domain=AA) == ([], [])
        assert Rx([2,0,1]).roots(domain=QQbar) == ([QQbar(-2).sqrt(), -QQbar(-2).sqrt()], [1, 1])
        assert (Rx([-2,0,1]) ** 2).roots(domain=AA) == ([AA(2).sqrt(), -AA(2).sqrt()], [2, 2])
    Rx = PolynomialRing(QQbar, "x")
    x = Rx.gen()
    g = -QQbar(3).sqrt() + x
    f = 2 + QQbar(2).sqrt()*x + x**2
    h = g**2 * f
    ((r1, r2, r3), (e1, e2, e3)) = h.roots(domain=QQbar)
    assert (x-r1)**e1 * (x-r2)**e2 * (x-r3)**e3 == h

def test_qqbar_sage_bug_37927():
    # check that the example in https://github.com/sagemath/sage/issues/37927
    # works with our implementation of qqbar
    for R in [QQbar, QQbar_ca]:
        I = R.i()
        v1 = -R.i()
        v2 = -R(2).sqrt()
        M = Mat(R)([[0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                           [-4, 2*v1, 1, 64, -32*v1, -16, 8*v1, 4, -2*v1, -1],
                           [4*v1, 1, 0, -192*v1, -80, 32*v1, 12, -4*v1, -1, 0],
                           [2, 0, 0, -480, 160*v1, 48, -12*v1, -2, 0, 0],
                           [-4, 2*I, 1, 64, -32*I, -16, 8*I, 4, -2*I, -1],
                           [4*I, 1, 0, -192*I, -80, 32*I, 12, -4*I, -1, 0],
                           [2, 0, 0, -480, 160*I, 48, -12*I, -2, 0, 0],
                           [0, 0, 0, 8, 4*v2, 4, 2*v2, 2, v2, 1],
                           [0, 0, 0, 24*v2, 20, 8*v2, 6, 2*v2, 1, 0],
                           [0, 0, 0, 8, 4*v2, 4, 2*v2, 2, -v2, 1],
                           [0, 0, 0, 24*v2, 20, 8*v2, 6, 2*v2, 1, 0],
                           [0, 0, 0, -4096, -1024*I, 256, 64*I, -16, -4*I, 1],
                           [0, 0, 0, -4096, 1024*I, 256, -64*I, -16, 4*I, 1]])
        X = M.nullspace()
        v = Mat(R)(10, 1, [-108, 0, 0, 1, 0, 12, 0, -60, 0, 64])
        assert (M * v).is_zero()

def test_ca_notebook_examples():
    # algebraic number identity
    NumberI = fexpr("NumberI")
    Sqrt = fexpr("Sqrt")
    Div = fexpr("Div")
    I = NumberI
    lhs = Sqrt(36 + 3*(-54+35*I*Sqrt(3))**Div(1,3)*3**Div(1,3) + \
                117/(-162+105*I*Sqrt(3))**Div(1,3))/3 + \
                Sqrt(5)*(1296*I+840*Sqrt(3)-35*3**Div(5,6)*(-54+35*I*Sqrt(3))**Div(1,3)-\
                54*I*(-162+105*I*Sqrt(3))**Div(1,3)+13*I*(-162+105*I*Sqrt(3))**Div(2,3))/(5*(162*I+105*Sqrt(3)))
    rhs = Sqrt(5) + Sqrt(7)
    i = CC_ca.i()
    pi = CC_ca.pi()
    exp = CC_ca.exp
    log = CC_ca.log
    sqrt = CC_ca.sqrt
    assert qqbar(lhs) == qqbar(rhs)
    assert ca(lhs) == ca(rhs)
    assert fexpr(ca(lhs) - ca(rhs)) == fexpr(0)
    # misc
    assert fexpr(exp(pi) * exp(-pi + log(2))) == fexpr(2)
    assert i**i - exp(pi / ((sqrt(-2)**sqrt(2)) ** sqrt(2))) == 0
    assert log(sqrt(2)+sqrt(3)) / log(5 + 2*sqrt(6)) == ca(1)/2
    assert ca(10)**-30 < (640320**3 + 744)/exp(pi*sqrt(163)) - 1 < ca(10)**-29
    M = Mat(CC_ca, 2, 2)
    A = M([[5, pi], [1, -1]])**4
    assert A.charpoly()(A) == M([[0,0],[0,0]])
    # comparison with higher precision
    ctx = ComplexField_ca()
    ctx._set_options({"prec_limit":65536})
    eps = ca(10, context=ctx) ** (-10000)
    assert (eps.exp() == 1) == False

    assert ZZ("0.000") == 0
    assert ZZ("3.0") == 3
    assert ZZ("-0.03e+2") == -3
    assert ZZ("0e-100000000000000") == 0
    assert raises(lambda: ZZ("0.1"), FlintUnableError)

    assert str(RR("(+/- 1e-3) + 0.1")) == "[0.1 +/- 1.01e-3]"
    assert str(RR("1/2 +/- 1/100000")) == "[0.5000 +/- 1.01e-5]"
    assert str(CC("(1+i) +/- (1e-5 + 1e-7*i)")) == "([1.0000 +/- 1.01e-5] + [1.000000 +/- 1.01e-7]*I)"
    assert str(CC("-1e23 +/- -5e12")) == "[-1.000000000e+23 +/- 5.01e+12]"

    assert RR("0.5 + 1") == 1.5
    assert QQ("0.01") == QQ(1) / 100
    assert QQ("1.01") == QQ(101) / 100
    assert QQ("1.01e+1 + 1") == QQ(111) / 10
    assert QQ("1.01e1 + 1") == QQ(111) / 10
    assert QQ("1.01e-1 + 1") == QQ(1101) / 1000
    assert RRx("0.75") == RRx(3)/4
    assert CCx("0.75") == RRx(3)/4

    assert str(RRx("(0.5 +/- 1e-10) + (0.6 +/- 1e-11)*x")) == "[0.500000000 +/- 1.01e-10] + [0.6000000000 +/- 1.01e-11]*x"
    assert str(RRx("(1 + x + x^2) +/- 0.003")) == "[1.00 +/- 3.01e-3] + x + x^2"
    assert str(RRx("(1 + x + x^2) +/- (0.003 + 0.004*x + 0.005*x^2 + 0.006*x^3)")) == "[1.00 +/- 3.01e-3] + [1.00 +/- 4.01e-3]*x + [1.00 +/- 5.01e-3]*x^2 + [+/- 6.01e-3]*x^3"
    assert str(RRx("1 + (+/- 1e-10)*x")) == "1 + [+/- 1.01e-10]*x"

    assert str(CCx("x +/- 1e-6*I*x")) == "(1.000000000000000 + [+/- 1.01e-6]*I)*x"

    with optimistic_logic:
        RRx(str(RRx("1+2*x")/3)) == RRx([1,2])/3

def test_gr_series():

    x = QQser.gen()
    # default prec is 6
    O6 = x**6
    On = lambda n: QQser(PowerSeriesRing(QQ, n, "x").gen() ** n)
    O0 = On(0)
    O1 = On(1)
    O2 = On(2)
    O3 = On(3)
    O4 = On(4)
    O5 = On(5)

    assert x / x == 1
    assert (2 * x) / (3 * x) == QQ(2) / 3
    assert (2 + 2*x) / (1 + x) == 2
    assert str(x / (x.exp() - 1)) == "1 + (-1/2)*x + (1/12)*x^2 + (-1/720)*x^4 + O(x^5)"

    assert str(x**6 / 1) == "0 + O(x^6)"
    assert str(x**6 / x) == "0 + O(x^5)"
    assert str(x**6 / x**5) == "0 + O(x^1)"
    assert str(x**6 / x**5) == "0 + O(x^1)"

    assert str(1 / (1 + x)) == "1 - x + x^2 - x^3 + x^4 - x^5 + O(x^6)"
    assert str(1 / (1 + x + O5)) == "1 - x + x^2 - x^3 + x^4 + O(x^5)"
    assert str((1 + O5) / (1 + x)) == "1 - x + x^2 - x^3 + x^4 + O(x^5)"
    assert str((1 + O5) / (1 + x + O4)) == "1 - x + x^2 - x^3 + O(x^4)"
    assert str((1 + O4) / (1 + x + O5)) == "1 - x + x^2 - x^3 + O(x^4)"
    assert str(x / (x + x**2)) == "1 - x + x^2 - x^3 + x^4 - x^5 + O(x^6)"
    assert str(x / (x + x**2 + O5)) == "1 - x + x^2 - x^3 + O(x^4)"
    assert str((x + O5) / (x + x**2)) == "1 - x + x^2 - x^3 + O(x^4)"
    assert str((x + O5) / (x + x**2 + O4)) == "1 - x + x^2 + O(x^3)"
    assert str((x + O4) / (x + x**2 + O5)) == "1 - x + x^2 + O(x^3)"

    assert str(O5 / 1) == "0 + O(x^5)"
    assert str(O5 / x) == "0 + O(x^4)"
    assert str(O5 / x**4) == "0 + O(x^1)"
    assert str(O5 / x**5) == "0 + O(x^0)"
    assert raises(lambda: O5 / x**6, FlintUnableError)

    assert raises(lambda: (0 * x) / 0, FlintDomainError)
    assert raises(lambda: x / 0, FlintDomainError)

    assert raises(lambda: O0 / 0, FlintDomainError)
    assert raises(lambda: O1 / 0, FlintDomainError)

    assert raises(lambda: 0 / O0, FlintUnableError)
    assert raises(lambda: 0 / O0, FlintUnableError)
    assert raises(lambda: 0 / O1, FlintUnableError)
    assert raises(lambda: 0 / O2, FlintUnableError)

    assert raises(lambda: O0 / O0, FlintUnableError)
    assert raises(lambda: O1 / O0, FlintUnableError)
    assert raises(lambda: O0 / O1, FlintUnableError)

    assert raises(lambda: 1 / O0, FlintUnableError)
    assert raises(lambda: 1 / O1, FlintDomainError)
    assert raises(lambda: 1 / O2, FlintDomainError)

    assert raises(lambda: x / O0, FlintUnableError)
    assert raises(lambda: x / O1, FlintUnableError)
    assert raises(lambda: x / O2, FlintDomainError)

    assert raises(lambda: x**2 / O0, FlintUnableError)
    assert raises(lambda: x**2 / O1, FlintUnableError)
    assert raises(lambda: x**2 / O2, FlintUnableError)
    assert raises(lambda: (x**0) / 0, FlintDomainError)

    assert raises(lambda: x**3 / O2, FlintUnableError)
    assert raises(lambda: x**3 / O3, FlintUnableError)

    R3 = PowerSeriesModRing(QQ, 3)
    assert R3(3 + O4) == R3(3)
    assert R3(3 + O3) == R3(3)
    assert raises(lambda: R3(3 + O2), FlintUnableError)

    R2 = PowerSeriesModRing(QQ, 2)
    R2b = PowerSeriesModRing(QQ, 2)
    assert R2(R3(5)) == 5
    assert R2(2) + R2b(3) == 5
    assert raises(lambda: R3(R2(5)), FlintDomainError)

    R = RRser
    x = R.gen()
    a = R(RR("0 +/- 1e-10"))
    On = lambda n: QQser(PowerSeriesRing(QQ, n, "x").gen() ** n)
    O6 = x**6
    O0 = On(0)
    O1 = On(1)
    O2 = On(2)
    O3 = On(3)
    O4 = On(4)
    O5 = On(5)

    assert raises(lambda: a == 0, Undecidable)
    assert raises(lambda: a * x == 0, Undecidable)
    assert not (a + x == 0)
    assert (a + x != 0)

    assert raises(lambda: 1 / a, FlintUnableError)
    assert raises(lambda: 1 / (a * x), FlintDomainError)
    assert raises(lambda: (a * x) / (a * x**2), FlintUnableError)
    assert raises(lambda: (a * x + x**3) / (a * x**2), FlintUnableError)
    assert raises(lambda: (x**3) / (a * x**4), FlintDomainError)
    assert raises(lambda: (a * x + x**3) / (a * x**2), FlintUnableError)

    x = PowerSeriesRing(ZZmod(1)).gen()
    assert x + x == 0
    assert x - x == 0
    assert x * x == 0
    assert x / x == 0


    R = PowerSeriesModRing(QQ, 6)
    x = R.gen()

    assert x**6 == 0
    assert x**5 != 0
    assert str(1 / (1 + x)) == "1 - x + x^2 - x^3 + x^4 - x^5 (mod x^6)"

    # Deflating quotients are nonunique and not supported by / by default
    assert raises(lambda: x / x, FlintDomainError)
    # assert x / x == 1
    # assert str(x / (x + x**2)) == "1 - x + x^2 - x^3 + x^4 - x^5 (mod x^6)"
    # assert str(x / (x.exp() - 1)) == "1 + (-1/2)*x + (1/12)*x^2 + (-1/720)*x^4 + (1/720)*x^5 (mod x^6)"

    assert raises(lambda: x / 0, FlintDomainError)

def test_integers_mod():
    R = IntegersMod_mpn_mod(10**20 + 1)
    c = ZZ(2) ** 4321
    assert R(3) * c == R(c) * 3
    assert R(3) + c == R(c) + 3
    assert R(3) - c == -(R(c) - 3)
    assert IntegersMod_mpn_mod(10**20)(IntegersMod_mpn_mod(10**20)(17)) == 17
    assert IntegersMod_mpn_mod(10**20)(IntegersMod_fmpz_mod(10**20)(17)) == 17
    assert raises(lambda: IntegersMod_mpn_mod(10**20)(IntegersMod_fmpz_mod(10**20 + 1)(17)), NotImplementedError)
    assert raises(lambda: IntegersMod_mpn_mod(10**20)(IntegersMod_fmpz_mod(10**50)(17)), NotImplementedError)
    assert raises(lambda: IntegersMod_mpn_mod(10**20)(IntegersMod_mpn_mod(10**20 + 1)(17)), NotImplementedError)
    assert raises(lambda: IntegersMod_mpn_mod(10**20)(IntegersMod_mpn_mod(10**50)(17)), NotImplementedError)

def test_nfloat():
    R = RealFloat_nfloat(128)
    R2 = RealField_arb(192)
    R3 = RealFloat_nfloat(256)
    tol = R(2.0**(-120))
    tol2 = R2(2.0**(-120))
    assert R2(R(3)) == 3
    assert R(R2(3)) == 3
    assert R2(R(0)) == 0
    assert R(R2(0)) == 0
    assert abs(R(R.pi() - R2.pi())) < tol
    assert abs(R2(R.pi() - R2.pi())) < tol
    assert abs(R(R.pi() - R2.pi())) < tol2
    assert abs(R2(R.pi() - R2.pi())) < tol2
    assert abs(R(QQ(1)/3) - QQ(1)/3) < tol
    assert R(ZZ(5)) == 5
    assert R(-ZZ(5)) == -5
    c = ZZ(5)**100
    assert abs(R2(R(c)) - c) < tol * c
    assert abs(R2(R(-c)) - (-c)) < tol * c
    assert R(R2(5)) == R(5)
    assert str(R(0)) == '0'
    assert str(R(1) / 4) == '0.250000000000000000000000000000000000000'
    assert R("0.25") == 0.25
    assert R(-0.25) == -0.25
    assert abs(R(R.pi() - R3.pi())) < tol
    assert abs(R3(R.pi() - R3.pi())) < tol
    assert R3(R(0)) == 0
    assert R3(R(-1)) == -1
    assert R(R3(0)) == 0
    assert R(-R3(1)) == -1
    assert R(3) <= R(3)
    assert R(-3) <= R(3)
    assert not (R(3) < R(3))
    assert not (R(3) <= R(-3))
    assert R(0) <= R(3)
    assert R(0) <= R(0)
    assert not (R(0) < R(0))
    assert not (R(0) <= R(-3))

def test_gen_name():
    for R in [NumberField(ZZx.gen() ** 2 + 1, "b"),
              PolynomialRing_fmpz_poly("b"),
              PolynomialRing_fmpq_poly("b"),
              PolynomialRing_gr_poly(QQbar, "b"),
              PowerSeriesRing_gr_series(ZZ, var="b"),
              PowerSeriesModRing_gr_poly(ZZ, 3, "b"),
              FiniteField_fq(3, 2, "b"),
              FiniteField_fq_nmod(3, 2, "b"),
              FiniteField_fq_zech(3, 2, "b")]:
        assert str(R.gen()) in ["b", "b^1", "b (mod b^3)"]
        R._set_gen_name("c")
        assert str(R.gen()) in ["c", "c^1", "c (mod c^3)"]
        R._set_gen_names(["d"])
        assert str(R.gen()) in ["d", "d^1", "d (mod d^3)"]

if __name__ == "__main__":
    from time import time
    print("Testing flint_ctypes")
    print("----------------------------------------------------------")
    for fname in dir():
        if fname.startswith("test_"):
            print(fname + "...", end="")
            import sys
            sys.stdout.flush()
            t1 = time()
            globals()[fname]()
            t2 = time()
            print("PASS", end="     ")
            print("%.2f" % (t2-t1))
    print("----------------------------------------------------------")
    import doctest
    __r = doctest.testmod(optionflags=(doctest.FAIL_FAST | doctest.ELLIPSIS), verbose=False)[0]
    if __r:
        sys.exit(__r)
    print("----------------------------------------------------------")
