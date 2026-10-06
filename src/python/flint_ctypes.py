import ctypes
import math
import random
import weakref
import ctypes.util
import sys
import functools

if sys.platform == "darwin":
    libflint = ctypes.CDLL("libflint.dylib")
else:
    libflint = ctypes.CDLL("libflint.so")
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

class radix_integer_struct(ctypes.Structure):
    _fields_ = [('d', ctypes.c_void_p),
                ('alloc', c_slong),
                ('size', c_slong)]

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

class fmpq_mpoly_struct(ctypes.Structure):
    _fields_ = [('content', fmpq_struct),
                ('zpoly', fmpz_mpoly_struct)]

class fmpz_mpoly_q_struct(ctypes.Structure):
    _fields_ = [('num', fmpz_mpoly_struct),
                ('den', fmpz_mpoly_struct)]


class padic_radix_struct(ctypes.Structure):
    _fields_ = [('u', radix_integer_struct),
                ('v', c_slong),
                ('N', c_slong)]

class decfloat_struct(ctypes.Structure):
    _fields_ = [('m', radix_integer_struct),
                ('exp', fmpz_struct)]

class decmag_struct(ctypes.Structure):
    _fields_ = [('m', c_ulong),
                ('exp', fmpz_struct)]

class decball_struct(ctypes.Structure):
    _fields_ = [('mid', decfloat_struct),
                ('rad', decmag_struct)]

class deccfloat_struct(ctypes.Structure):
    _fields_ = [('re', decfloat_struct),
                ('im', decfloat_struct)]

class deccball_struct(ctypes.Structure):
    _fields_ = [('re', decball_struct),
                ('im', decball_struct)]

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

class gr_tower_lazy_elem_struct(ctypes.Structure):
    # a lazy tower element: pointer to the tower, level, versions, the
    # representation tag and a union of an fmpq, an fmpz_mpoly_q with its
    # context pointer, or the dense form (17 words in total); allocated
    # with slack (the size is checked when a context is created)
    _fields_ = [('data', c_slong * 24)]

class gr_tower_field_struct(ctypes.Structure):
    # an element of the top field of a fixed tower: a gr_poly (nested
    # representation), or an element of the base field for the trivial
    # tower (an fmpq, or a rational function: an fmpz_mpoly_q, 12 words);
    # allocated with slack
    _fields_ = [('data', c_slong * 14)]

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
    _fields_ = [('n', c_ulong)]

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
    __module__ = Exception.__module__

class gr_ctx_struct(ctypes.Structure):
    _fields_ = [('content', ctypes.c_char * libgr.gr_ctx_sizeof_ctx())]


libflint.flint_malloc.restype = ctypes.c_void_p
libflint.flint_free.argtypes = (ctypes.c_void_p,)
libflint.fmpz_set_str.argtypes = ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int
libflint.fmpz_get_str.argtypes = ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(fmpz_struct)
libflint.fmpz_get_str.restype = ctypes.c_void_p
libflint.n_is_prime.argtypes = (c_ulong,)

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
libgr.gr_derivative_gen.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))

libgr.gr_set_d.argtypes = (ctypes.c_void_p, ctypes.c_double, ctypes.POINTER(gr_ctx_struct))

libgr.gr_set_str.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str_n.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmp.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmpabs.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_heap_clear.argtypes = (ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_tower_heap_init.argtypes = (ctypes.POINTER(gr_ctx_struct),)
libgr.gr_tower_heap_init.restype = ctypes.c_void_p
libgr.gr_tower_heap_clear.argtypes = (ctypes.c_void_p,)
libgr.gr_tower_get_str.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p)
libgr.gr_tower_degree.argtypes = (ctypes.c_void_p,)
libgr.gr_tower_degree.restype = c_slong
libgr.gr_tower_length_si.argtypes = (ctypes.c_void_p,)
libgr.gr_tower_length_si.restype = c_slong
libgr.gr_tower_num_gens_si.argtypes = (ctypes.c_void_p,)
libgr.gr_tower_num_gens_si.restype = c_slong
libgr.gr_tower_lazy_ctx_num_towers.restype = c_slong
libgr.gr_tower_lazy_ctx_num_towers.argtypes = (ctypes.c_void_p,)
libgr.gr_tower_gen_name.argtypes = (ctypes.c_void_p, c_slong)
libgr.gr_tower_gen_name.restype = ctypes.c_void_p
libgr.gr_tower_gen_get.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong)
libgr.gr_tower_adjoin_root_ui.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_ulong, ctypes.c_char_p)
libgr.gr_tower_adjoin_qqbar.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p)
libgr.gr_tower_adjoin_root_of_unity.argtypes = (ctypes.c_void_p, c_ulong, ctypes.c_char_p)
libgr.gr_tower_adjoin_algebraic.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_char_p)
libgr.gr_tower_adjoin_pi.argtypes = (ctypes.c_void_p, ctypes.c_char_p)
libgr.gr_tower_adjoin_exp.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p)
libgr.gr_tower_adjoin_log.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p)
libgr.gr_ctx_init_tower_field.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_void_p)
libgr.gr_ctx_init_tower_lazy.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.POINTER(gr_ctx_struct), ctypes.c_int)
libgr.gr_ctx_init_tower_lazy_view.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.POINTER(gr_ctx_struct), ctypes.c_int)
libgr.gr_tower_lazy_ctx_set_print.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_int, c_slong)
libgr.gr_tower_lazy_ctx_print_flags.argtypes = (ctypes.POINTER(gr_ctx_struct),)
libgr.gr_tower_lazy_ctx_set_gen_flags.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_int)
libgr.gr_tower_lazy_ctx_gen_flags.argtypes = (ctypes.POINTER(gr_ctx_struct),)
libgr.gr_tower_lazy_ctx_set_option.argtypes = (ctypes.POINTER(gr_ctx_struct), c_slong, c_slong)
libgr.gr_tower_lazy_ctx_set_option.restype = ctypes.c_int
libgr.gr_tower_option_name.argtypes = (c_slong,)
libgr.gr_tower_option_name.restype = ctypes.c_char_p
libgr.gr_tower_option_find.argtypes = (ctypes.c_char_p,)
libgr.gr_tower_option_find.restype = c_slong
libgr.gr_tower_option_default.argtypes = (c_slong,)
libgr.gr_tower_option_default.restype = c_slong
libgr.gr_tower_lazy_ctx_get_option.argtypes = (ctypes.POINTER(gr_ctx_struct), c_slong)
libgr.gr_tower_lazy_ctx_get_option.restype = c_slong
libgr.gr_tower_lazy_get_tower.argtypes = (ctypes.POINTER(c_slong), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_tower_lazy_get_tower.restype = ctypes.c_void_p

libgr.gr_ctx_init_nmod.argtypes = (ctypes.POINTER(gr_ctx_struct), c_ulong)
libgr.gr_ctx_init_dirichlet_group.argtypes = (ctypes.POINTER(gr_ctx_struct), c_ulong)
libgr.gr_ctx_init_padic_radix.argtypes = (ctypes.POINTER(gr_ctx_struct), c_ulong, c_slong, c_slong, ctypes.c_int)

libgr._gr_ctx_init_decimal.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_int, ctypes.c_uint, c_slong, ctypes.c_int, ctypes.c_int)
libgr.decimal_ctx_set_prec.argtypes = (ctypes.POINTER(gr_ctx_struct), c_slong)
libgr.decimal_ctx_set_rnd.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_int)
libgr.decimal_ctx_set_rad_prec.argtypes = (ctypes.POINTER(gr_ctx_struct), c_slong)
libgr.decimal_ctx_set_exp_limits.argtypes = (ctypes.POINTER(gr_ctx_struct), c_slong, c_slong)
libgr.decimal_ctx_set_flags.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_int)
libgr.decimal_ctx_get_prec.restype = c_slong
libgr.decimal_ctx_get_rad_prec.restype = c_slong
libgr.decimal_ctx_get_limb_digits.restype = c_slong
libgr.decimal_ctx_get_exp_limits.argtypes = (ctypes.POINTER(c_slong), ctypes.POINTER(c_slong), ctypes.POINTER(gr_ctx_struct))
libgr.decfloat_set_round.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.c_int, ctypes.POINTER(gr_ctx_struct))
libgr.decfloat_digits.restype = c_slong
libgr.decfloat_limbs.restype = c_slong
libgr.decfloat_get_digit_si.argtypes = (ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.decfloat_get_digit_si.restype = c_ulong
libgr.decfloat_set_digit_si.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, c_ulong, ctypes.POINTER(gr_ctx_struct))
libgr.decfloat_get_sci_exp_si.argtypes = (ctypes.POINTER(c_slong), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.decfloat_get_str_sci.restype = ctypes.c_void_p
libgr.decball_set_round.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.decball_add_error_10exp_si.argtypes = (ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.decball_rel_accuracy_digits.restype = c_slong
libgr.decimal_ctx_set_rnd_im.argtypes = (ctypes.POINTER(gr_ctx_struct), ctypes.c_int)
libgr.deccfloat_set_round.argtypes = (ctypes.c_void_p, ctypes.c_void_p, c_slong, ctypes.c_int, ctypes.c_int, ctypes.POINTER(gr_ctx_struct))
libgr.deccball_add_error_decmag.argtypes = (ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.deccball_rel_accuracy_digits.restype = c_slong

_add_methods = [libgr.gr_add, libgr.gr_add_si, libgr.gr_add_fmpz, libgr.gr_add_other, libgr.gr_other_add]
_sub_methods = [libgr.gr_sub, libgr.gr_sub_si, libgr.gr_sub_fmpz, libgr.gr_sub_other, libgr.gr_other_sub]
_mul_methods = [libgr.gr_mul, libgr.gr_mul_si, libgr.gr_mul_fmpz, libgr.gr_mul_other, libgr.gr_other_mul]
_div_methods = [libgr.gr_div, libgr.gr_div_si, libgr.gr_div_fmpz, libgr.gr_div_other, libgr.gr_other_div]
_divexact_methods = [libgr.gr_divexact, libgr.gr_divexact_si, libgr.gr_divexact_fmpz, libgr.gr_divexact_other, libgr.gr_other_divexact]
_pow_methods = [libgr.gr_pow, libgr.gr_pow_si, libgr.gr_pow_fmpz, libgr.gr_pow_other, libgr.gr_other_pow]

libgr._gr_series_set_error.argtypes = (ctypes.c_void_p, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr._gr_series_get_error.restype = c_slong



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
        libgr.gr_ctx_uninitialized(self._ref)
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
            >>> ZZser.is_integral_domain()
            True
            >>> ZZmod(4).is_integral_domain()
            False
            >>> PowerSeriesRing(ZZmod(4)).is_integral_domain()
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

    def is_pretend_field(self):
        """
        Return whether this structure is a field or pretends to be one
        (see :meth:`set_is_pretend_field`).

            >>> QQ.is_pretend_field()
            True
            >>> ZZ.is_pretend_field()
            False

        """
        return self._ctx_predicate(libflint.gr_ctx_is_pretend_field, "is_pretend_field")

    def set_is_pretend_field(self, flag=True):
        """
        Make this ring compute as if it were a field: an operation which
        meets a nonzero non-invertible element raises
        ``FlintUnableError`` after recording a zero divisor, which
        :meth:`recover_zero_divisor` returns.

            >>> R = IntegersMod_fmpz_mod(91)
            >>> R.set_is_pretend_field()
            >>> R.is_pretend_field()
            True
            >>> R(14).inv()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute inv(x) in {Integers mod 91 (fmpz)} for {x = 14}
            >>> R.recover_zero_divisor()
            7
            >>> R(5).inv()
            73

        """
        status = libflint.gr_ctx_set_is_pretend_field(self._ref, T_TRUE if flag else T_FALSE)
        if status:
            _handle_error(self, status, "set_is_pretend_field")

    def recover_zero_divisor(self):
        """
        Return the zero divisor recorded by a ring pretending to be a
        field (see :meth:`set_is_pretend_field`).
        """
        return self._constant(self, libflint.gr_ctx_recover_zero_divisor, "recover_zero_divisor")

    def is_rational_vector_space(self):
        """
        Return whether this structure is a vector space over the rational numbers.

            >>> QQx.is_rational_vector_space()
            True
            >>> ZZi.is_rational_vector_space()
            False

        """
        return self._ctx_predicate(libflint.gr_ctx_is_rational_vector_space, "is_rational_vector_space")

    def is_real_vector_space(self):
        """
        Return whether this structure is a vector space over the real numbers.

            >>> QQx.is_real_vector_space()
            False
            >>> RRx.is_real_vector_space()
            True
            >>> CC.is_real_vector_space()
            True

        """
        return self._ctx_predicate(libflint.gr_ctx_is_real_vector_space, "is_real_vector_space")

    def is_complex_vector_space(self):
        """
        Return whether this structure is a vector space over the complex numbers.

            >>> QQx.is_complex_vector_space()
            False
            >>> CC.is_complex_vector_space()
            True
            >>> CCx.is_complex_vector_space()
            True

        """
        return self._ctx_predicate(libflint.gr_ctx_is_complex_vector_space, "is_complex_vector_space")

    def is_exact(self):
        """
        Return whether elements of this structure are represented exactly.

            >>> QQ.is_exact()
            True
            >>> RR.is_exact()
            False
            >>> RealFloat_decfloat(None).is_exact(), RealFloat_decfloat(10).is_exact()
            (True, False)

        """
        return self._ctx_predicate(libflint.gr_ctx_is_exact, "is_exact")

    def is_canonical(self):
        """
        Return whether equal elements of this structure are guaranteed to
        have a unique representation (so that structural equality is
        mathematical equality).

            >>> QQ.is_canonical()
            True
            >>> RR.is_canonical()
            False
            >>> RealFloat_decfloat(10).is_canonical(), ComplexFloat_deccfloat(10).is_canonical()
            (True, True)

        """
        return self._ctx_predicate(libflint.gr_ctx_is_canonical, "is_canonical")


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

    def O(ctx, x, exp):
        """
        Create a big-O error term. The base must be a generator in
        a power series ring or the base of a p-adic ring.

            >>> ZZser
            Power series over Integer ring (fmpz) with precision O(x^6)
            >>> x = ZZser.gen()
            >>> ZZser.O(x, 4)
            0 + O(x^4)
            >>> 3*x + 5*x**4 + ZZser.O(x, 3)
            3*x + O(x^3)
            >>> QQser("1/3 + x/5 - O(x^2)")
            (1/3) + (1/5)*x + O(x^2)
            >>> (2**32 - 1) * (1 + Qp_padic_radix(2).O(2, 8))
            255 + O(2^8)

        Constant input is ambiguous. Use the two-argument version instead:

        """
        return ctx._binary_op_fmpz(x, exp, libgr.gr_big_o_base_fmpz, "big_o($x, $y)")

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
            >>> PowerSeriesRing(CC).i()
            (1.000000000000000*I)

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
            >>> PowerSeriesRing(CC).pi()
            [3.141592653589793 +/- 3.39e-16]
            >>> PowerSeriesRing(CC, prec=0).pi()
            0 + O(x^0)

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
            [3.141592653589793 +/- 3.39e-16]
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
            FlintDomainError: log(x) is not an element of {Rational field (fmpq)} for {x = 2}
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
            >>> Mat(QQ)([[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,4,-6,4]]).log()
            [[-11/6, 3, -3/2, 1/3],
            [-1/3, -1/2, 1, -1/6],
            [1/6, -1, 1/2, 1/3],
            [-1/3, 3/2, -3, 11/6]]
            >>> _.exp()
            [[0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [-1, 4, -6, 4]]

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
            x + [0.3333333333333333 +/- 7.04e-17]*x^3 + [0.1333333333333333 +/- 5.37e-17]*x^5 + O(x^6)
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
            [3.676077910374977 +/- 9.92e-16]
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
            >>> CC.erf(1j)
            [1.650425758797543 +/- 4.58e-16]*I
            >>> RRser.erf(RRser("1+x", error=2))
            [0.842700792949715 +/- 3.28e-16] + [0.415107497420595 +/- 6.47e-16]*x + O(x^2)
            >>> RRser.erf(RRser("1+x", error=1))
            [0.842700792949715 +/- 3.28e-16] + O(x^1)
            >>> RRser.erf(RRser("1+x", error=0))
            0 + O(x^0)
            >>> QQser.erf(QQser("1+x"))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute erf(x) in {Power series over Rational field (fmpq) with precision O(x^6)} for {x = 1 + x}
        """
        return ctx._unary_op(x, libgr.gr_erf, "erf($x)")

    def erfc(ctx, x):
        """
            >>> RR.erfc(1)
            [0.1572992070502851 +/- 3.71e-17]
            >>> CC.erfc(1j)
            (1.000000000000000 + [-1.650425758797543 +/- 2.45e-16]*I)
            >>> RRser.erfc(RRser("1+x", error=2))
            [0.1572992070502851 +/- 3.71e-17] + [-0.415107497420595 +/- 6.47e-16]*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_erfc, "erfc($x)")

    def erfi(ctx, x):
        """
            >>> RR.erfi(1)
            [1.650425758797543 +/- 4.58e-16]
            >>> CC.erfi(1j)
            [0.842700792949715 +/- 3.28e-16]*I
            >>> RRser.erfi(RRser("1+x", error=2))
            [1.650425758797543 +/- 4.58e-16] + [3.06725258552748 +/- 6.09e-15]*x + O(x^2)
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
            >>> CC.fresnel_s(1j)
            [-0.3102683017233811 +/- 2.67e-18]*I
            >>> CC.fresnel_s(1j, normalized=True)
            [-0.4382591473903548 +/- 9.24e-17]*I
            >>> RRser.fresnel_s(RRser("1+x", error=2))
            [0.3102683017233811 +/- 2.67e-18] + [0.841470984807897 +/- 6.08e-16]*x + O(x^2)
            >>> RRser.fresnel_s(RRser("1+x", error=2), normalized=True)
            [0.4382591473903548 +/- 9.24e-17] + x + O(x^2)
            >>> CCser.fresnel_s(CCser("I+x", error=2))
            ([-0.3102683017233811 +/- 2.67e-18]*I) + [-0.841470984807897 +/- 6.08e-16]*x + O(x^2)
            >>> CCser.fresnel_s(CCser("I+x", error=2), normalized=True)
            ([-0.4382591473903548 +/- 9.24e-17]*I) - x + O(x^2)
            >>> RRser.fresnel_c(1)
            [0.904524237900272 +/- 1.46e-16]
        """
        return ctx._unary_op_with_flag(x, normalized, libgr.gr_fresnel_s, "fresnel_s($x)")

    def fresnel_c(ctx, x, normalized=False):
        """
            >>> RR.fresnel_c(1)
            [0.904524237900272 +/- 1.46e-16]
            >>> RR.fresnel_c(1, normalized=True)
            [0.779893400376823 +/- 3.59e-16]
            >>> CC.fresnel_c(1j)
            [0.904524237900272 +/- 1.46e-16]*I
            >>> CC.fresnel_c(1j, normalized=True)
            [0.779893400376823 +/- 3.59e-16]*I
            >>> RRser.fresnel_c(RRser("1+x", error=2))
            [0.904524237900272 +/- 1.46e-16] + [0.540302305868140 +/- 4.59e-16]*x + O(x^2)
            >>> RRser.fresnel_c(RRser("1+x", error=2), normalized=True)
            [0.779893400376823 +/- 3.59e-16] + O(x^2)
            >>> CCser.fresnel_c(CCser("I+x", error=2))
            ([0.904524237900272 +/- 1.46e-16]*I) + [0.540302305868140 +/- 4.59e-16]*x + O(x^2)
            >>> CCser.fresnel_c(CCser("I+x", error=2), normalized=True)
            ([0.779893400376823 +/- 3.59e-16]*I) + O(x^2)
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
            >>> CC.gamma_upper(3, 1)
            [1.839397205857211 +/- 7.30e-16]
            >>> CCser.gamma_upper(3, CCser("1+I*x", error=2))
            [1.839397205857211 +/- 7.30e-16] + ([-0.3678794411714423 +/- 7.79e-17]*I)*x + O(x^2)

        """
        return ctx._binary_op_with_flag(x, y, regularized, libgr.gr_gamma_upper, "gamma_upper($x, $y)")

    def gamma_lower(ctx, x, y, regularized=0):
        """
            >>> RR.gamma_lower(3, 4)
            [1.52379338889291 +/- 2.89e-15]
            >>> RR.gamma_lower(3, 4, regularized=True)
            [0.76189669444646 +/- 6.52e-15]
            >>> RR.gamma_lower(3, 1)
            [0.160602794142788 +/- 5.34e-16]
            >>> RRser.gamma_lower(RRser("3+x"), RRser("1+x", error=1))
            [0.160602794142788 +/- 5.34e-16] + O(x^1)
            >>> RRser.gamma_lower(RRser("3+x"), RRser("1+x", error=2))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute gamma_lower(x, y) in {Power series over Real numbers (arb, prec = 53) with precision O(x^6)} for {x = 3.000000000000000 + x}, {y = 1 + x + O(x^2)}
            >>> RRser.gamma_lower(RRser("3"), RRser("1+x", error=2))
            [0.160602794142788 +/- 5.34e-16] + [0.3678794411714423 +/- 7.79e-17]*x + O(x^2)
            >>> CCser.gamma_lower(3, CCser("1+I*x", error=2))
            [0.160602794142788 +/- 5.34e-16] + ([0.3678794411714423 +/- 7.79e-17]*I)*x + O(x^2)

        """
        return ctx._binary_op_with_flag(x, y, regularized, libgr.gr_gamma_lower, "gamma_lower($x, $y)")

    def beta_lower(ctx, a, b, x, regularized=0):
        """
            >>> RR.beta_lower(2, 3, 0.5)
            [0.0572916666666667 +/- 6.08e-17]
            >>> RR.beta_lower(2, 3, 0.5, regularized=True)
            [0.687500000000000 +/- 5.24e-16]
            >>> RRser.beta_lower(2, 3, RRser("1+x"))
            [0.0833333333333333 +/- 4.26e-17] + [0.3333333333333333 +/- 7.04e-17]*x^3 + 0.2500000000000000*x^4 + O(x^6)
            >>> CCser.beta_lower(2, 3, CCser("1+I*x"))
            [0.0833333333333333 +/- 4.26e-17] + ([-0.3333333333333333 +/- 7.04e-17]*I)*x^3 + 0.2500000000000000*x^4 + O(x^6)
            >>> CCser.beta_lower(2, RRser("1+x"), RRser("1+x"))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute beta_lower(a, b, x) in {Power series over Complex numbers (acb, prec = 53) with precision O(x^6)} for {a = 2.000000000000000}, {b = 1 + x}, {x = 1 + x}
            >>> CCser.beta_lower(RRser("1+x"), 2, RRser("1+x"))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute beta_lower(a, b, x) in {Power series over Complex numbers (acb, prec = 53) with precision O(x^6)} for {a = 1 + x}, {b = 2.000000000000000}, {x = 1 + x}
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
            [1.895117816355937 +/- 6.91e-16]
            >>> RRser.exp_integral_ei(RRser("1+x", error=2))
            [1.895117816355937 +/- 6.91e-16] + [2.718281828459045 +/- 5.41e-16]*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_exp_integral_ei, "exp_integral_ei($x)")

    def sin_integral(ctx, x):
        """
            >>> RR.sin_integral(1)
            [0.946083070367183 +/- 1.35e-16]
            >>> CC.sin_integral(CC("(5 +/- 1e-10)*I"))
            [20.09321183 +/- 5.79e-9]*I
            >>> CC.sin_integral(CC("10 + (1 +/- 1e-10)*I"))
            ([1.7002629761 +/- 6.33e-11] + [-0.0667638998 +/- 2.17e-11]*I)
            >>> RRser.sin_integral(RRser("1+x", error=2))
            [0.946083070367183 +/- 1.35e-16] + [0.8414709848078965 +/- 3.37e-17]*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_sin_integral, "sin_integral($x)")

    def cos_integral(ctx, x):
        """
            >>> RR.cos_integral(1)
            [0.3374039229009681 +/- 5.63e-17]
            >>> RRser.cos_integral(RRser("1+x", error=2))
            [0.3374039229009681 +/- 5.63e-17] + [0.540302305868140 +/- 4.59e-16]*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_cos_integral, "cos_integral($x)")

    def sinh_integral(ctx, x):
        """
            >>> RR.sinh_integral(1)
            [1.057250875375728 +/- 6.29e-16]
            >>> RRser.sinh_integral(RRser("1+x", error=2))
            [1.057250875375728 +/- 6.29e-16] + [1.175201193643801 +/- 6.61e-16]*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_sinh_integral, "sinh_integral($x)")

    def cosh_integral(ctx, x):
        """
            >>> RR.cosh_integral(1)
            [0.837866940980208 +/- 3.15e-16]
            >>> RRser.cosh_integral(RRser("1+x", error=2))
            [0.837866940980208 +/- 3.15e-16] + [1.543080634815244 +/- 5.28e-16]*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_cosh_integral, "cosh_integral($x)")

    def log_integral(ctx, x, offset=False):
        """
            >>> RR.log_integral(2)
            [1.045163780117493 +/- 9.71e-16]
            >>> RR.log_integral(2, offset=True)
            0
            >>> RRser.log_integral(RRser("2+x", error=2))
            [1.045163780117493 +/- 9.71e-16] + [1.442695040888963 +/- 8.70e-16]*x + O(x^2)
            >>> CCser.log_integral(CCser("2+x", error=2))
            [1.045163780117493 +/- 9.71e-16] + [1.442695040888963 +/- 8.70e-16]*x + O(x^2)
            >>> CCser.log_integral(CCser("2+x", error=2), offset=True)
            [1.442695040888963 +/- 8.70e-16]*x + O(x^2)
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
            >>> sum(CC.bessel_j(0, k) for k in range(100))
            [1.419207859380 +/- 2.35e-13]
            >>> w = CC.exp_pi_i(QQ(1)/100); sum(CC.bessel_j(0, w*k) for k in range(101))
            ([0.78030446659 +/- 2.44e-12] + [0.62756344686 +/- 5.49e-12]*I)
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
            >>> CC.airy_ai(1j)
            ([0.3314933054321412 +/- 8.35e-17] + [-0.3174498589684437 +/- 9.89e-17]*I)
            >>> RRser.airy_ai(RRser(1))
            [0.1352924163128814 +/- 4.17e-17]
            >>> RRser.airy_ai(RRser("1+x", error=2))
            [0.1352924163128814 +/- 4.17e-17] + [-0.1591474412967932 +/- 2.95e-17]*x + O(x^2)
            >>> CCser.airy_ai(CCser("1+I*x", error=2))
            [0.1352924163128814 +/- 4.17e-17] + ([-0.1591474412967932 +/- 2.95e-17]*I)*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_airy_ai, "airy_ai($x)")

    def airy_bi(ctx, x):
        """
            >>> RR.airy_bi(1)
            [1.207423594952871 +/- 3.27e-16]
            >>> CC.airy_bi(1j)
            ([0.648858208330395 +/- 2.42e-16] + [0.3449586347680483 +/- 8.89e-17]*I)
            >>> RRser.airy_bi(RRser("1+x", error=2))
            [1.207423594952871 +/- 3.27e-16] + [0.932435933392776 +/- 5.83e-16]*x + O(x^2)
            >>> CCser.airy_bi(CCser("1+I*x", error=2))
            [1.207423594952871 +/- 3.27e-16] + ([0.932435933392776 +/- 5.83e-16]*I)*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_airy_bi, "airy_bi($x)")


    def airy_ai_prime(ctx, x):
        """
            >>> RR.airy_ai_prime(1)
            [-0.1591474412967932 +/- 2.95e-17]
            >>> CC.airy_ai_prime(1j)
            ([-0.4324926598418071 +/- 8.56e-17] + [0.0980478562292432 +/- 3.82e-17]*I)
            >>> RRser.airy_ai_prime(RRser("1+x", error=2))
            [-0.1591474412967932 +/- 2.95e-17] + [0.1352924163128814 +/- 4.17e-17]*x + O(x^2)
            >>> CCser.airy_ai_prime(CCser("1+I*x", error=2))
            [-0.1591474412967932 +/- 2.95e-17] + ([0.1352924163128814 +/- 4.17e-17]*I)*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_airy_ai_prime, "airy_ai_prime($x)")

    def airy_bi_prime(ctx, x):
        """
            >>> RR.airy_bi_prime(1)
            [0.932435933392776 +/- 5.83e-16]
            >>> CC.airy_bi_prime(1j)
            ([0.1350266467108190 +/- 6.84e-17] + [-0.1288373867812549 +/- 6.33e-17]*I)
            >>> RRser.airy_bi_prime(RRser("1+x", error=2))
            [0.932435933392776 +/- 5.83e-16] + [1.207423594952871 +/- 3.27e-16]*x + O(x^2)
            >>> CCser.airy_bi_prime(CCser("1+I*x", error=2))
            [0.932435933392776 +/- 5.83e-16] + ([1.207423594952871 +/- 3.27e-16]*I)*x + O(x^2)
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
            [2.49185259170895 +/- 9.81e-15]
            >>> CC.legendre_q(3, 1, 0.5)
            [2.49185259170895 +/- 9.81e-15]
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
            ([+/- 6.32e-16] + [1.00000000000000 +/- 1.03e-15]*I)
            >>> RRser.gamma(RRser("1+x",error=2))
            1 + [-0.577215664901533 +/- 3.58e-16]*x + O(x^2)
            >>> CCser.gamma(CCser("1+I*x",error=2))
            ([1.00000000000000 +/- 3.36e-16] + [+/- 3.63e-21]*I) + ([+/- 2.61e-20] + [-0.57721566490153 +/- 4.09e-15]*I)*x + O(x^2)
            >>> RRser("x").gamma()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute gamma(x) in {Power series over Real numbers (arb, prec = 53) with precision O(x^6)} for {x = x}
            >>> CCser("x").gamma()
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute gamma(x) in {Power series over Complex numbers (acb, prec = 53) with precision O(x^6)} for {x = x}
        """
        return ctx._unary_op_with_fmpz_fmpq_overloads(x, libgr.gr_gamma, op_fmpz=libgr.gr_gamma_fmpz, op_fmpq=libgr.gr_gamma_fmpq, rstr="gamma($x)")

    def lgamma(ctx, x):
        """
            >>> RR.lgamma(10)
            [12.80182748008147 +/- 2.69e-15]
            >>> CC.lgamma(10j)
            ([-15.94031728124131 +/- 6.90e-15] + [12.23211664743500 +/- 4.89e-15]*I)
            >>> RRser.lgamma(RRser("1+x",error=2))
            [-0.577215664901533 +/- 3.58e-16]*x + O(x^2)
            >>> CCser.lgamma(CCser("1+I*x",error=2))
            ([-0.5772156649015329 +/- 4.84e-17]*I)*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_lgamma, "lgamma($x)")

    def rgamma(ctx, x):
        """
            >>> RR.rgamma(10)
            [2.755731922398589e-6 +/- 5.96e-22]
            >>> CC.rgamma(10+1j)
            ([-1.83246026966323e-6 +/- 5.08e-21] + [-2.25314671311995e-6 +/- 5.78e-21]*I)
            >>> RRser.rgamma(RRser("1+x",error=2))
            1 + [0.577215664901533 +/- 3.58e-16]*x + O(x^2)
            >>> CCser.rgamma(CCser("1+I*x",error=2))
            ([1.000000000000000 +/- 3.30e-16] + [+/- 3.63e-21]*I) + ([+/- 2.61e-20] + [0.57721566490153 +/- 4.16e-15]*I)*x + O(x^2)
        """
        return ctx._unary_op(x, libgr.gr_rgamma, "lgamma($x)")

    def digamma(ctx, x):
        """
            >>> RR.digamma(2)
            [0.4227843350984671 +/- 4.84e-17]
            >>> CC.digamma(2j)
            ([0.714591515373977 +/- 6.06e-16] + [1.820807282642230 +/- 3.65e-16]*I)
            >>> RRser.digamma(RRser("1+x",error=2))
            [-0.5772156649015329 +/- 9.00e-17] + [1.644934066848226 +/- 4.71e-16]*x + O(x^2)
            >>> CCser.digamma(CCser("1+I*x",error=2))
            ([-0.5772156649015329 +/- 5.63e-17] + [+/- 7.87e-20]*I) + ([+/- 1.44e-19] + [1.644934066848226 +/- 6.75e-16]*I)*x + O(x^2)
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
            >>> CCser.hurwitz_zeta(CCser("0.5+x", error=2), CCser("0.5"))
            [-0.6048986434216304 +/- 3.21e-17] + [-3.056337630862498 +/- 5.73e-16]*x + O(x^2)
            >>> CCser.hurwitz_zeta(CCser("0.5+x", error=2), CCser("0.5+x"))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute hurwitz_zeta(s, a) in {Power series over Complex numbers (acb, prec = 53) with precision O(x^6)} for {s = 0.5000000000000000 + x + O(x^2)}, {a = 0.5000000000000000 + x}

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
            >>> CCser.polylog(CCser("0.5+x", error=2), CCser("0.5"))
            [0.80612672304285 +/- 5.56e-15] + ([-0.2905212160378 +/- 4.15e-14] + [+/- 3.02e-14]*I)*x + O(x^2)
            >>> CCser.polylog(CCser("0.5+x", error=2), CCser("0.5+x"))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute polylog(s, z) in {Power series over Complex numbers (acb, prec = 53) with precision O(x^6)} for {s = 0.5000000000000000 + x + O(x^2)}, {z = 0.5000000000000000 + x}

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
            >>> sum(RR.bernoulli_vec(100)).nprint(16)
            1.127124216595034e+76
            >>> sum(RF.bernoulli_vec(100))
            1.127124216595034e+76
            >>> sum(CC.bernoulli_vec(100)).nprint(16)
            1.127124216595034e+76

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
            >>> CC.dirichlet_l(2, DirichletGroup(4)(3))
            [0.915965594177219 +/- 2.68e-16]
            >>> CCser.dirichlet_l(CCser("2+x", error=2), DirichletGroup(4)(3))
            [0.915965594177219 +/- 5.98e-16] + [0.08158073611659 +/- 4.36e-15]*x + O(x^2)
            >>> QQser.dirichlet_l(2)
            Traceback (most recent call last):
              File "<stdin>", line 1, in <module>
            TypeError: gr_ctx.dirichlet_l() missing 1 required positional argument: 'chi'
            >>> QQser.dirichlet_l(2, DirichletGroup(1)(1))
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute dirichlet_l(s, chi) in {Power series over Rational field (fmpq) with precision O(x^6)} for {s = 2}, {chi = chi_1(1, .)}

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
            >>> CCser.hardy_theta(CCser("10+x", error=2), DirichletGroup(4)(3))
            [4.64979557270698 +/- 4.41e-15] + [0.925292493992482 +/- 5.80e-16]*x + O(x^2)
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
            >>> CCser.hardy_z(CCser("2+x", error=2), DirichletGroup(4)(3))
            [1.15107760668266 +/- 5.01e-15] + [0.36975256259156 +/- 6.01e-15]*x + O(x^2)
            >>> CCser.hardy_z(CCser(2), DirichletGroup(4)(3))
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
            [[3.1512120021539 +/- 3.41e-14], [+/- 4.40e-14], [4.2557730353652 +/- 9.85e-14]]
        """
        return ctx._op_vec_arg_len(tau, n, libgr.gr_eisenstein_g_vec, "eisenstein_g_vec($tau, $n)")

    def agm(ctx, x, y=None):
        """
        Arithmetic-geometric mean.

            >>> RR.agm(2)
            [1.456791031046907 +/- 9.72e-16]
            >>> RR.agm(2, 3)
            [2.47468043623630 +/- 4.68e-15]
            >>> CCser.agm(CCser("1+x", error=3))
            1 + 0.5000000000000000*x - 0.06250000000000000*x^2 + O(x^3)
            >>> CCser.agm(CCser("2+x", error=2))
            [1.456791031046907 +/- 4.54e-16] + [0.4257908959543789 +/- 6.82e-17]*x + O(x^2)
        """
        if y is None:
            return ctx._unary_op(x, libgr.gr_agm1, "agm1($x)")
        else:
            return ctx._binary_op(x, y, libgr.gr_agm, "agm($x, $y)")

    def elliptic_k(ctx, m):
        """
            >>> CC.elliptic_k(0.5)
            [1.85407467730137 +/- 3.40e-15]
            >>> CCser.elliptic_k(CCser("0.5 + x", error=2))
            [1.85407467730137 +/- 3.43e-15] + [0.84721308479398 +/- 5.48e-15]*x + O(x^2)
        """
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
            >>> CCser.jacobi_theta_1(CCser("0.125+x", error=2), CCser("I"))
            ([0.347386687929454 +/- 3.21e-16] + [+/- 1.18e-16]*I) + ([2.64053630032731 +/- 2.35e-15] + [+/- 2.68e-16]*I)*x + O(x^2)

        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_1, "jacobi_theta_1($z, $tau)")

    def jacobi_theta_2(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_2(0.125, 1j)
            [0.843115469091413 +/- 8.18e-16]
            >>> CCser.jacobi_theta_2(CCser("0.125+x", error=2), CCser("I"))
            ([0.843115469091413 +/- 8.18e-16] + [+/- 7.86e-17]*I) + ([-1.11111761491452 +/- 5.95e-15] + [+/- 4.05e-16]*I)*x + O(x^2)
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_2, "jacobi_theta_2($z, $tau)")

    def jacobi_theta_3(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_3(0.125, 1j)
            [1.061113709291166 +/- 5.74e-16]
            >>> CCser.jacobi_theta_3(CCser("0.125+x", error=2), CCser("I"))
            ([1.061113709291166 +/- 5.74e-16] + [+/- 3.09e-17]*I) + ([-0.38407640677719 +/- 2.70e-15] + [+/- 2.18e-16]*I)*x + O(x^2)
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_3, "jacobi_theta_3($z, $tau)")

    def jacobi_theta_4(ctx, z, tau):
        """
        Jacobi theta function.

            >>> CC.jacobi_theta_4(0.125, 1j)
            [0.938886290708834 +/- 3.52e-16]
            >>> CCser.jacobi_theta_4(CCser("0.125+x", error=2), CCser("I"))
            ([0.938886290708834 +/- 3.52e-16] + [+/- 3.09e-17]*I) + ([0.38390111383116 +/- 3.67e-15] + [+/- 2.18e-16]*I)*x + O(x^2)
        """
        return ctx._binary_op(z, tau, libgr.gr_jacobi_theta_4, "jacobi_theta_4($z, $tau)")

    def elliptic_invariants(ctx, tau):
        """
            >>> g2, g3 = CC.elliptic_invariants(1j)
            >>> CC.weierstrass_p_prime(0.25, 1j)**2; 4*CC.weierstrass_p(0.25, 1j)**3 - g2*CC.weierstrass_p(0.25, 1j) - g3
            [15152.862386715 +/- 7.03e-10]
            [15152.862386715 +/- 9.62e-10]
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
        """
            >>> CC.weierstrass_p(CC("0.5"), CC("I"))
            [6.875185818020 +/- 4.48e-13]
            >>> CCser.weierstrass_p(CCser("0.5+x", error=3), CCser("I"))
            [6.875185818020 +/- 4.53e-13] + [+/- 2.91e-19]*x + [94.53636006462 +/- 4.09e-12]*x^2 + O(x^3)

        """
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
    def _binary_op_si(self, other, op, rstr):
        other = int(other)
        elem_type = type(self)
        res = elem_type(context=self._ctx_python)
        status = op(res._ref, self._ref, other, self._ctx)
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

    def overlaps(self, other):
        self, other = gr_elem._binary_coercion(self, other)
        truth = libgr.gr_equal(self._ref, other._ref, self._ctx)
        if truth == T_FALSE: return False
        return True

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

    def derivative_gen(self, i):
        """
        Derivative with respect to ith generator.

            >>> ZZx("3*x^2").derivative_gen(0)
            6*x
            >>> ZZx("3*x^2").derivative_gen(1)
            Traceback (most recent call last):
              ...
            FlintDomainError: ...

            >>> RRx("x^3/3").derivative_gen(0)
            [1.00000000000000 +/- 3.89e-16]*x^2
            >>> RRx("x").derivative_gen(1)
            Traceback (most recent call last):
              ...
            FlintDomainError: ...

            >>> x, y = PolynomialRing_fmpz_mpoly(2, ["x", "y"]).gens()
            >>> ((3*x + 5*y)**2).derivative_gen(0)
            18*x+30*y
            >>> ((3*x + 5*y)**2).derivative_gen(1)
            30*x+50*y
            >>> x.derivative_gen(2)
            Traceback (most recent call last):
              ...
            FlintDomainError: ...

            >>> x, y = PolynomialRing_fmpq_mpoly(2, ["x", "y"]).gens()
            >>> ((3*x + 5*y/2)**2).derivative_gen(0)
            18*x + 15*y
            >>> ((3*x + 5*y/2)**2).derivative_gen(1)
            15*x + 25/2*y
            >>> x.derivative_gen(2)
            Traceback (most recent call last):
              ...
            FlintDomainError: ...

        """
        return self._binary_op_si(self, i, libgr.gr_derivative_gen, "derivative_gen")


    def divexact(self, other):
        """
        Assert that other divides self exactly in the ring and compute
        the quotient, allowing the implementation to skip checks
        for correctness. This will often be faster than performing a checked
        division with the / operator.

        Exact division is useful, for example, when manipulating
        inexact power series over non-fields, where the default checked
        division cannot tell whether terms beyond the O(x^n) error term
        would divide if the denominator is not a unit:

            >>> x = ZZser.gen()
            >>> A = 7 + 3*x + 5*x**2
            >>> B = 12 + 2*x + 19*x**3 + 3*x**4
            >>> A * B
            84 + 50*x + 66*x^2 + 143*x^3 + 78*x^4 + 104*x^5 + O(x^6)
            >>> (A * B) / B
            Traceback (most recent call last):
              ...
            FlintUnableError: failed to compute x / y in {Power series over Integer ring (fmpz) with precision O(x^6)} for {x = 84 + 50*x + 66*x^2 + 143*x^3 + 78*x^4 + 104*x^5 + O(x^6)}, {y = 12 + 2*x + 19*x^3 + 3*x^4}
            >>> (A * B).divexact(B)
            7 + 3*x + 5*x^2 + O(x^6)
        """
        return self._binary_op2(self, other, _divexact_methods, "$x / $y (exact division)")

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
            >>> Mat(QQ)([[1,0,1],[1,1,0],[0,0,1]]).sqrt()
            [[1, 0, 1/2],
            [1/2, 1, -1/8],
            [0, 0, 1]]
            >>> Mat(QQbar)([[1,2,3],[4,5,6],[7,8,9]]).sqrt()
            [[Root a = 0.449756 + 0.762279*I of 132*a^4+100*a^2+81, Root a = 0.552622 + 0.206796*I of 99*a^4-52*a^2+12, Root a = 0.655487 - 0.348687*I of 1188*a^4-732*a^2+361],
            [Root a = 1.01852 + 0.0841514*I of 33*a^4-68*a^2+36, Root a = 1.25147 + 0.0228291*I of 99*a^4-310*a^2+243, Root a = 1.48442 - 0.0384931*I of 297*a^4-1308*a^2+1444],
            [Root a = 1.58729 - 0.593976*I of 12*a^4-52*a^2+99, Root a = 1.95032 - 0.161138*I of 9*a^4-68*a^2+132, Root a = 2.31335 + 0.271701*I of 108*a^4-1140*a^2+3179]]
            >>> _ ** 2
            [[1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]]
            >>> QQser(4).sqrt()
            2
            >>> QQser(3).sqrt()
            Traceback (most recent call last):
              ...
            FlintDomainError: sqrt(x) is not an element of {Power series over Rational field (fmpq) with precision O(x^6)} for {x = 3}
        """
        return self._unary_op(self, libgr.gr_sqrt, "sqrt($x)")

    def rsqrt(self):
        """
        Reciprocal square root of this element.

            >>> QQ(25).rsqrt()
            1/5
            >>> Mat(QQ)([[1,0,1],[1,1,0],[0,0,1]]).rsqrt()
            [[1, 0, -1/2],
            [-1/2, 1, 3/8],
            [0, 0, 1]]
            >>> PowerSeriesRing(ZZmod(1))().rsqrt()
            0
            >>> QQser(4).rsqrt()
            (1/2)
            >>> QQser(3).rsqrt()
            Traceback (most recent call last):
              ...
            FlintDomainError: rsqrt(x) is not an element of {Power series over Rational field (fmpq) with precision O(x^6)} for {x = 3}
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
            FlintDomainError: exp(x) is not an element of {Rational field (fmpq)} for {x = 1}
            >>> QQser.gen().exp()
            1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 + (1/120)*x^5 + O(x^6)

        """
        return self._unary_op(self, libgr.gr_exp, "exp($x)")

    def expm1(self):
        """
        Exponential function minus 1.

            >>> RR("1e-10").expm1()
            [1.000000000050000e-10 +/- 1.69e-26]
            >>> CC(RR("1e-10")).expm1()
            [1.000000000050000e-10 +/- 1.69e-26]
            >>> RF("1e-10").expm1()
            1.000000000050000e-10
            >>> CF(RF("1e-10")).expm1()
            1.000000000050000e-10
            >>> QQ(0).expm1()
            0
            >>> QQ(1).expm1()
            Traceback (most recent call last):
              ...
            FlintDomainError: expm1(x) is not an element of {Rational field (fmpq)} for {x = 1}
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
            FlintDomainError: log(x) is not an element of {Rational field (fmpq)} for {x = 2}
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

class IntegerRing_radix_integer(gr_ctx):
    """
        >>> ZZdec = IntegerRing_radix_integer(10)
        >>> ZZ(5)**100 + 123
        7888609052210118054117285652827862296732064351090230047702789306640748
        >>> ZZdec(5)**100 + 123
        7888609052210118054117285652827862296732064351090230047702789306640748
        >>> (ZZdec(10) // 3, ZZdec(10) % 3, ZZdec(-10) // 3, ZZdec(-10) % 3)
        (3, 1, -4, 2)

        >>> ZZ7 = IntegerRing_radix_integer(7, 10)
        >>> ZZ7
        Integers in radix 7^10 (radix_integer)
        >>> ZZ7(5**100) + 123
        194 * 7^80 + 171307382 * 7^70 + 127454085 * 7^60 + 241336803 * 7^50 + 220718973 * 7^40 + 16449109 * 7^30 + 123134395 * 7^20 + 38207260 * 7^10 + 51613155
        >>> ZZ7(str(_)) == ZZ(str(_))
        True

    """

    def __init__(self, base, exp=0):
        assert base >= 2 and base <= UWORD_MAX
        assert exp >= 0 and exp <= FLINT_BITS
        assert base**exp <= UWORD_MAX
        gr_ctx.__init__(self)
        libgr.gr_ctx_init_radix_integer(self._ref, base, exp)
        self._elem_type = radix_integer

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

GR_TOWER_MERGE_EXPRESS = 1
GR_TOWER_LAZY_REAL = 2
GR_TOWER_LAZY_ALGEBRAIC = 4

GR_TOWER_PRINT_NUMERIC = 1
GR_TOWER_PRINT_SYMBOLIC = 2
GR_TOWER_PRINT_DEFS = 4
GR_TOWER_GENS_SPLIT_IMAGINARY = 1
GR_TOWER_GENS_COMPOSITE_ROOTS = 2
GR_TOWER_GENS_COMPOSITE_RADICALS = 4

# tuning options of lazy tower fields (GR_TOWER_OPT_*, in order), from the library
def _gr_tower_option_names():
    names = []
    k = 0
    while True:
        name = libgr.gr_tower_option_name(k)
        if name is None:
            return names
        names.append(name.decode('ascii'))
        k += 1
GR_TOWER_OPTION_NAMES = _gr_tower_option_names()
GR_TOWER_TRIG_EXPONENTIAL = 0
GR_TOWER_TRIG_TANGENT = 1

class gr_tower_lazy_ctx(gr_ctx):
    r"""
    Base class of the lazy tower fields (``gr_ctx_init_tower_lazy``):
    elements live in towers of algebraic and transcendental extensions
    of QQ which grow as needed, with complete zero tests (Richardson's
    algorithm, assuming Schanuel's conjecture for termination).

    Every extension element created during the lifetime of the context
    gets a unique symbol (``a1``, ``a2``, ... for algebraic definitions,
    ``t1``, ``t2``, ... for exponentials and logarithms, and ``pi``,
    ``i``), listed by :meth:`gens`:

        >>> C = ComplexField_tower()
        >>> [C(2).sqrt(), C(3).sqrt(), C(1).exp().exp(), C(2).log().log()]
        [a1 {a1 = sqrt(2)}, a2 {a2 = sqrt(3)}, t2 {t1 = exp(1); t2 = exp(t1)}, t4 {t3 = log(2); t4 = log(t3)}]
        >>> C.gens()
        [a1 {a1 = sqrt(2)}, a2 {a2 = sqrt(3)}, t1 {t1 = exp(1)}, t2 {t1 = exp(1); t2 = exp(t1)}, t3 {t3 = log(2)}, t4 {t3 = log(2); t4 = log(t3)}]

    The printing of elements is controlled by :meth:`set_print`: a
    numerical value, the expression in the generators, and the
    definitions of the generators involved, in any combination:

        >>> x = 2 + C(2).exp() + 3*C(2).exp().exp()
        >>> x
        3*t6+t5+2 {t5 = exp(2); t6 = exp(t5)}
        >>> C.set_print(numeric=True, symbolic=True, defs=False); x
        4863.92 {3*t6+t5+2}
        >>> C.set_print(numeric=True, symbolic=False); x
        4863.92
        >>> C.set_print(numeric=True, symbolic=True, defs=True, digits=10); x
        4863.923032 {3*t6+t5+2 where t5 = exp(2); t6 = exp(t5)}
        >>> C.set_print(); x
        3*t6+t5+2 {t5 = exp(2); t6 = exp(t5)}

    The printed form (with the definitions) can be read back, in the
    same context or in another one -- the definitions are evaluated in
    order, and the names are local to the string -- which is also how
    elements convert between lazy contexts:

        >>> D = ComplexField_tower()
        >>> D("3*t2+t1+2 {t1 = exp(2); t2 = exp(t1)}") == 2 + D(2).exp() + 3*D(2).exp().exp()
        True
        >>> D(x)
        3*t2+t1+2 {t1 = exp(2); t2 = exp(t1)}
        >>> r = PolynomialRing(C)([-1, -1, 0, 0, 0, 1]).roots()[0][1]; r     # doctest: +ELLIPSIS
        a... {a... = root(-1 - a... + a...^5, 0.181232 + 1.08395*i)}
        >>> D(str(r)) == D(r), D(r) ** 5 - D(r) - 1
        (True, 0)
        >>> RealField_tower()(C(2).sqrt()), RealAlgebraicField_tower()(C(2).sqrt() + C(3).sqrt())
        (a1 {a1 = sqrt(2)}, a1+a2 {a1 = sqrt(2); a2 = sqrt(3)})

    Symbolic expressions (:meth:`gr_elem.fexpr`) and LaTeX are available
    for elements:

        >>> x.fexpr()
        Where(Add(Mul(3, t_6), t_5, 2), Def(t_5, Exp(2)), Def(t_6, Exp(t_5)))
        >>> C(2).sqrt().latex()
        'a_{1}\\; \\text{ where } a_{1} = \\sqrt{2}'
        >>> C.acos(C(1)/2), C.log((C(3).sqrt() + C.i())/2)
        (pi/3, pi*i/6)

    Special functions apply their functional equations and special
    values before creating generators (for canonical arguments only),
    and Lambert W values take part in the relation search:

        >>> C = ComplexField_tower()
        >>> (C(1)/3).gamma() * (C(2)/3).gamma() == 2*C.pi()/C(3).sqrt()
        True
        >>> C(2).sqrt().gamma()
        t2*a2-t2 {a2 = sqrt(2); t2 = gamma(a2-1)}
        >>> (C(1)/4).digamma()
        (-pi-6*t4-2*t3)/2 {t3 = euler; t4 = log(2)}
        >>> C.lambertw(3*C(3).exp()), C.polylog(2, C(1)/2)
        (3, (pi^2-6*t4^2)/12 {t4 = log(2)})
        >>> C.lambertw(C(1)).exp() * C.lambertw(C(1))
        1
        >>> C.polygamma(1, C(3)/4), C.polylog(2, C.i())
        (pi^2-8*t10 {t10 = catalan}, (-pi^2+48*t10*i)/48 {t10 = catalan})
    """
    _flags = GR_TOWER_MERGE_EXPRESS

    def __init__(self, flags=None, parent=None, **options):
        """
        With *parent* (a lazy field), the new field is a view of it: it
        shares its towers and elements (conversions between the two are
        free) and restricts it further to the real or algebraic
        subfield. Keyword arguments set tuning options (``set_option``).
        """
        gr_ctx.__init__(self)
        self._base = QQ
        if flags is None:
            flags = self._flags
        if parent is None:
            libgr.gr_ctx_init_tower_lazy(self._ref, QQ._ref, flags)
        else:
            libgr.gr_ctx_init_tower_lazy_view(self._ref, parent._ref, flags & (GR_TOWER_LAZY_REAL | GR_TOWER_LAZY_ALGEBRAIC))
            self._parent = parent
        self._elem_type = gr_tower_lazy
        assert libgr.gr_ctx_sizeof_elem(self._ref) <= ctypes.sizeof(gr_tower_lazy_elem_struct)
        for name, value in options.items():
            self.set_option(name, value)

    def _view_class(self, flags):
        return {0: ComplexField_tower, GR_TOWER_LAZY_REAL: RealField_tower,
                GR_TOWER_LAZY_ALGEBRAIC: ComplexAlgebraicField_tower,
                GR_TOWER_LAZY_REAL | GR_TOWER_LAZY_ALGEBRAIC: RealAlgebraicField_tower}[flags]

    def real_field(self):
        """
        The real subfield of this field, as a view sharing its elements.

            >>> C = ComplexField_tower(); R = C.real_field(); R
            Real field (lazy towers)
            >>> x = C(2).sqrt() + C.pi()
            >>> R(x)
            pi+a1 {a1 = sqrt(2)}
            >>> R(C.i())     # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
              ...
            FlintDomainError
        """
        flags = libgr.gr_tower_lazy_ctx_field_flags(self._ref) | GR_TOWER_LAZY_REAL
        return self._view_class(flags)(parent=self)

    def algebraic_field(self):
        """
        The algebraic subfield of this field, as a view sharing its
        elements.
        """
        flags = libgr.gr_tower_lazy_ctx_field_flags(self._ref) | GR_TOWER_LAZY_ALGEBRAIC
        return self._view_class(flags)(parent=self)

    def stats(self):
        """
        Prints statistics about the towers of the context (with the
        environment variable ``GR_TOWER_STATS_VERBOSE``, the sizes of
        their moduli).
        """
        libgr.gr_tower_lazy_ctx_stats(self._ref)

    def num_towers(self):
        """
        The number of towers of the context. Towers in which no element
        lives any more are collected (except those holding canonical
        definitions: pi, roots of unity, radicals of integers), so this
        stays bounded in long sessions:

            >>> C = ComplexField_tower()
            >>> x = PolynomialRing(C).gen()
            >>> for k in range(200):
            ...     r = (x**3 - C(2).sqrt()*x - k).roots()
            >>> C.num_towers() < 20
            True
        """
        return libgr.gr_tower_lazy_ctx_num_towers(self._ref)

    def set_print(self, numeric=False, symbolic=True, defs=True, digits=6):
        """
        Sets the printing of elements: a numerical value with the
        given number of digits, the symbolic expression in the
        generators, and the definitions of the generators involved.
        Called without arguments, restores the default (expression with
        definitions).
        """
        flags = 0
        if numeric: flags |= GR_TOWER_PRINT_NUMERIC
        if symbolic: flags |= GR_TOWER_PRINT_SYMBOLIC
        if defs: flags |= GR_TOWER_PRINT_DEFS
        libgr.gr_tower_lazy_ctx_set_print(self._ref, flags, digits)

    @property
    def print_flags(self):
        return libgr.gr_tower_lazy_ctx_print_flags(self._ref)

    def set_gens(self, split_imaginary=False, composite_roots=False, composite_radicals=False):
        """
        Sets how new generators are chosen. By default, the square
        root of a negative rational number is a single generator
        sqrt(-A) (except sqrt(-1) = i and sqrt(-3) = 2*exp(2*pi*i/3)+1).
        With ``split_imaginary=True``, it is instead i*sqrt(A), whose
        generators are shared with the real square roots:

            >>> C = ComplexField_tower()
            >>> C(-163).sqrt(), C(-12).sqrt()
            (a1 {a1 = sqrt(-163)}, 4*a2+2 {a2 = exp(2*pi*i/3)})
            >>> C.set_gens(split_imaginary=True)
            >>> C(-163).sqrt()
            a3*i {a3 = sqrt(163)}

        By default, roots of unity of composite orders are products of
        roots of unity of prime power orders, and square roots of
        positive rationals are products of square roots of primes. With
        ``composite_roots=True``, a root of unity is a power of one
        generator exp(2*pi*i/N) for N the least common multiple of the
        orders requested (Calcium's choice: faster arithmetic within one
        cyclotomic field, slower when many orders occur; the same as
        ``set_option("cyclotomic_degree_limit", 128)``), and with ``composite_radicals=True`` the square root of
        A B^2 (A squarefree) is B sqrt(A):

            >>> C = ComplexField_tower()
            >>> C.set_gens(composite_roots=True, composite_radicals=True)
            >>> C.i() * C("exp(2*pi*i/15)")
            a1^15+a1^13-a1^7-a1^5+a1 {a1 = exp(2*pi*i/60)}
            >>> C(24).sqrt()
            2*a2 {a2 = sqrt(6)}

        Elements created before the change are unaffected (and remain
        compatible with those created after it).
        """
        flags = 0
        if split_imaginary: flags |= GR_TOWER_GENS_SPLIT_IMAGINARY
        if composite_roots: flags |= GR_TOWER_GENS_COMPOSITE_ROOTS
        if composite_radicals: flags |= GR_TOWER_GENS_COMPOSITE_RADICALS
        libgr.gr_tower_lazy_ctx_set_gen_flags(self._ref, flags)

    @property
    def gen_flags(self):
        return libgr.gr_tower_lazy_ctx_gen_flags(self._ref)

    @staticmethod
    def _option_index(name):
        if isinstance(name, int):
            return name
        try:
            return GR_TOWER_OPTION_NAMES.index(name)
        except ValueError:
            raise ValueError("unknown option: %s" % name)

    def set_option(self, name, value):
        r"""
        Sets a tuning option, shared with the views of this field (see
        ``options`` for the names). It affects elements created
        afterwards. For example, ``cyclotomic_degree_limit`` (default 0)
        is the largest degree phi(N) of a root of unity exp(2*pi*i/N)
        of composite order N used as one generator (the roots of unity
        requested are its powers); beyond it, roots of unity of prime
        power orders are used. ``trig_form`` selects the form of sin,
        cos and tan of real arguments in the complex field:
        GR_TOWER_TRIG_EXPONENTIAL (through exp) or GR_TOWER_TRIG_TANGENT
        (through tangents, as in the real field). Elements of number
        fields of one generator of degree up to ``dense_form_degree_limit``
        (default 2^20) also have a dense form (a polynomial in the
        generator, multiplied with FLINT's polynomial arithmetic), on which
        arithmetic needs no lock; polynomials longer than 256 with fewer
        than one nonzero coefficient in ``dense_form_sparsity`` (default 8)
        stay sparse. With ``primitive_degree_limit`` (default 0), merged
        number fields of several steps up to that degree become
        QQ(theta) for a primitive element (worthwhile for small degrees
        only):

            >>> C = ComplexField_tower(cyclotomic_degree_limit=16)
            >>> C.i() * C("exp(2*pi*i/15)")
            a1^15+a1^13-a1^7-a1^5+a1 {a1 = exp(2*pi*i/60)}
            >>> C.get_option("cyclotomic_degree_limit")
            16
            >>> C.set_option("trig_form", GR_TOWER_TRIG_TANGENT)
            >>> C.cos(C.pi() / 5)
            (-a2^2+7)/8 {a2 = tan(pi/5)}
            >>> D = ComplexField_tower(primitive_degree_limit=4)
            >>> D.sqrt(2) + D.sqrt(3)
            a3 {a3 = root(1 - 10*a3^2 + a3^4, 3.14626)}
        """
        if libgr.gr_tower_lazy_ctx_set_option(self._ref, self._option_index(name), value) != GR_SUCCESS:
            raise ValueError("value %s out of range for option %s" % (value, name))

    def get_option(self, name):
        return libgr.gr_tower_lazy_ctx_get_option(self._ref, self._option_index(name))

    @property
    def options(self):
        r"""
        The tuning options and their values:

            >>> C = ComplexField_tower()
            >>> C.options["cyclotomic_degree_limit"], C.options["prec_limit"]
            (0, 256)
        """
        return {name: self.get_option(k) for k, name in enumerate(GR_TOWER_OPTION_NAMES)}


class ComplexField_tower(gr_tower_lazy_ctx):
    r"""
    Field of exactly represented complex numbers (algebraic numbers,
    exponentials and logarithms and what is built from them), with
    specific embeddings in the complex numbers.

        >>> C = ComplexField_tower()
        >>> C
        Complex field (lazy towers)
        >>> C.exp(C.pi() * C.i()) + 1
        0
        >>> C.sqrt(2) * C.sqrt(3) == C.sqrt(6)
        True
        >>> (C.log(3) - C.log(2)) - C.log(C(3)/2)
        0
        >>> C.sqrt(C(1)/2)
        a1/2 {a1 = sqrt(2)}
        >>> C.sqrt(C(1)/2) ** 2
        1/2
        >>> C(-1).sqrt(), C.i()
        (i, i)
    """
    _flags = GR_TOWER_MERGE_EXPRESS

class RealField_tower(gr_tower_lazy_ctx):
    r"""
    The real subfield of :class:`ComplexField_tower`: an operation whose
    result is not real fails with a domain error.

        >>> R = RealField_tower()
        >>> R
        Real field (lazy towers)
        >>> R(2).sqrt() + R(2).exp()
        t1+a1 {a1 = sqrt(2); t1 = exp(2)}
        >>> R(-2).sqrt()
        Traceback (most recent call last):
          ...
        FlintDomainError: sqrt(x) is not an element of {Real field (lazy towers)} for {x = -2}
        >>> R.i()
        Traceback (most recent call last):
          ...
        FlintDomainError: i is not an element of {Real field (lazy towers)}
        >>> R.log(-1)
        Traceback (most recent call last):
          ...
        FlintDomainError: log(x) is not an element of {Real field (lazy towers)} for {x = -1}
        >>> R.asin(2)
        Traceback (most recent call last):
          ...
        FlintDomainError: asin(x) is not an element of {Real field (lazy towers)} for {x = 2}
        >>> R.acos(R(1)/2)
        pi/3
        >>> R(QQbar(-1).sqrt())
        Traceback (most recent call last):
          ...
        FlintDomainError: Root a = 1.00000*I of a^2+1 is not defined in Real field (lazy towers)
        >>> (R.pi().sqrt() - 1).sgn()
        1

    Values are computed through the complex numbers where that is the
    natural way (the representation may then involve nonreal
    generators, as for ``cos(1)`` below); values involving roots of unity
    are rewritten in real terms, and read back only when they are real.
    The real trigonometric constants at rational multiples of pi are
    polynomials in one generator tan(pi/M) per level (the tangent normal
    form; square roots when it is quadratic):

        >>> R = ComplexField_tower().real_field()
        >>> R.cos(R.pi() / 5), R.sin(R.pi() / 12)
        ((-a1^2+7)/8 {a1 = tan(pi/5)}, (a2^2+8*a2+1)/8 {a2 = tan(pi/24)})
        >>> R.tan(R.pi() / 12)
        -a3+2 {a3 = sqrt(3)}
        >>> PolynomialRing(R, "x")([1, 0, 0, 0, 1]).factor()
        (1, [1 + (-a5 {a5 = sqrt(2)})*x + x^2, 1 + (a5 {a5 = sqrt(2)})*x + x^2], [1, 1])
        >>> R("a1 {a1 = exp(2*pi*i/8)}")     # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
          ...
        FlintDomainError
        >>> R("a1^2+a1^6 {a1 = exp(2*pi*i/8)}")
        0

    The trigonometric functions of real arguments (other than rational
    multiples of pi) are expressed through the real generators tan and
    atan, whose relations Richardson's algorithm finds on their angles:

        >>> R.cos(R(1)), R.atan(R(2))
        ((-t1^2+1)/(t1^2+1) {t1 = tan(1/2)}, t2 {t2 = atan(2)})
        >>> 4 * R.atan(R(1)/5) - R.atan(R(1)/239) == R.pi() / 4
        True
        >>> R.atan(R.tan(R(2))), R.tan(R.atan(R(2)))
        (-pi+2, 2)
        >>> R.sin(R(1) + R(2).sqrt()) == R.sin(1) * R.cos(R(2).sqrt()) + R.cos(1) * R.sin(R(2).sqrt())
        True
    """
    _flags = GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_REAL

class ComplexAlgebraicField_tower(gr_tower_lazy_ctx):
    r"""
    The algebraic numbers, represented in towers: exponentials,
    logarithms, trigonometric functions and pi are available only at
    the arguments where their values are algebraic (which is decided by
    the Lindemann-Weierstrass and Gelfond-Schneider theorems).

        >>> A = ComplexAlgebraicField_tower()
        >>> A
        Complex algebraic field (lazy towers)
        >>> A(-2).sqrt(), A(2)**(QQ(1)/3)
        (a1 {a1 = sqrt(-2)}, a2 {a2 = root(2, 3)})
        >>> A.exp(0), A.log(1), A.cos(0), A(1) ** A(2).sqrt()
        (1, 0, 1, 1)
        >>> A.exp(1)
        Traceback (most recent call last):
          ...
        FlintDomainError: exp(x) is not an element of {Complex algebraic field (lazy towers)} for {x = 1}
        >>> A.pi()
        Traceback (most recent call last):
          ...
        FlintDomainError: pi is not an element of {Complex algebraic field (lazy towers)}
        >>> A(2).sqrt() ** A(2).sqrt()
        Traceback (most recent call last):
          ...
        FlintDomainError: x ** y is not an element of {Complex algebraic field (lazy towers)} for {x = a3 {a3 = sqrt(2)}}, {y = a3 {a3 = sqrt(2)}}
        >>> A(QQbar(2).sqrt() + QQbar(3).sqrt()) == A(2).sqrt() + A(3).sqrt()
        True
    """
    _flags = GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_ALGEBRAIC

class RealAlgebraicField_tower(gr_tower_lazy_ctx):
    r"""
    The real algebraic numbers, represented in towers.

        >>> A = RealAlgebraicField_tower()
        >>> A
        Real algebraic field (lazy towers)
        >>> (1 + A(5).sqrt()) / 2
        (a1+1)/2 {a1 = sqrt(5)}
        >>> A(-2).sqrt()
        Traceback (most recent call last):
          ...
        FlintDomainError: sqrt(x) is not an element of {Real algebraic field (lazy towers)} for {x = -2}
        >>> PolynomialRing(A)([-1, -1, 0, 0, 0, 1]).roots()
        ([a2 {a2 = root(-1 - a2 + a2^5, 1.16730)}], [1])
    """
    _flags = GR_TOWER_MERGE_EXPRESS | GR_TOWER_LAZY_REAL | GR_TOWER_LAZY_ALGEBRAIC


class gr_tower:
    r"""
    A fixed tower of fields `\mathbb{Q}(t_1, \ldots)(a_1)(a_2)\cdots`
    (``gr_tower_t``), built step by step by adjoining square roots,
    roots, algebraic numbers, exponentials, logarithms and pi, and the
    field of its top level (:meth:`field`), a ``gr`` context with
    complete zero tests (dynamic evaluation refines the tower when a
    defining polynomial turns out to be reducible). Elements of the field
    are polynomials in the last generator over the field below, printed
    in nested form.

        >>> T = gr_tower()
        >>> a1 = T.adjoin_sqrt(2); a1
        a1
        >>> a2 = T.adjoin_sqrt(3); a2
        a2
        >>> K = T.field(); K
        Tower field Rational field (fmpq)(a1, a2)
        >>> a1, a2 = K.gens()
        >>> (a1 + a2) ** 2
        5 + (2*a1)*a2
        >>> 1 / (a1 + a2)
        -a1 + a2
        >>> a3 = T.adjoin_sqrt(a1 + a2)
        >>> a3 ** 4
        (5 + (2*a1)*a2)
        >>> a3 + a1      # elements of the earlier field convert to the new one
        a1 + a3
        >>> T
        Tower of degree 8 over Rational field (fmpq)
          a1 = 1.414213562  root_2(2)  =  root of  -2 + a1^2  (irreducible)
          a2 = 1.732050808  root_2(3)  =  root of  -3 + a2^2  (irreducible)
          a3 = 1.773771228  root_2(a2+a1)  =  root of  (-a1 - a2) + a3^2  (irreducible)

    A step adjoined with a polynomial which is not irreducible is
    refined when an operation exposes the factorization (the field
    context always represents a field):

        >>> T = gr_tower()
        >>> a1 = T.adjoin_sqrt(2)
        >>> K = T.field()
        >>> a2 = T.adjoin_algebraic(PolynomialRing(K)([-2, 0, 1]), CC(1.4142))    # x^2 - 2 again, dynamic
        >>> T.degree()
        4
        >>> a2 == a1     # the zero test exposes the factorization
        True
        >>> T.degree()   # and the tower was refined
        2
        >>> a2 - a1
        0

    Transcendental generators make the base a rational function field
    (the tower is rebuilt: field objects and elements created before
    become stale):

        >>> T = gr_tower()
        >>> pi = T.adjoin_pi(); pi
        pi
        >>> a1 = T.adjoin_sqrt(2)
        >>> t2 = T.adjoin_exp(a1)
        >>> t2.parent()
        Tower field Fraction field of multivariate polynomials over Integer ring (fmpz) in 2 variables, lex order(a1)
        >>> T
        Tower of degree 2 over Fraction field of multivariate polynomials over Integer ring (fmpz) in 2 variables, lex order
          pi = 3.141592654  pi  (transcendental)
          a1 = 1.414213562  root_2(2)  =  root of  -2 + a1^2  (irreducible)
          t2 = 4.113250379  exp(a1)  (conjecturally transcendental)
        >>> a1.parent()(1)
        Traceback (most recent call last):
          ...
        ValueError: this field object is stale: the tower was rebuilt (transcendental generator adjoined); use tower.field() again
        >>> pi, a1, t2 = T.gens()
        >>> (t2 - pi) * a1
        (-pi+t2)*a1

    Polynomials over the field of a tower factor into irreducibles
    (Trager's method), and their roots in the field are the linear
    factors:

        >>> T = gr_tower()
        >>> a1 = T.adjoin_sqrt(2); a2 = T.adjoin_sqrt(3)
        >>> R = PolynomialRing(T.field(), "x"); x = R.gen()
        >>> (x**4 + 1).factor()
        (1, [1 + a1*x + x^2, 1 - a1*x + x^2], [1, 1])
        >>> (x**4 - 10*x**2 + 1).roots()
        ([-a1 + a2, -a1 - a2, a1 - a2, a1 + a2], [1, 1, 1, 1])
        >>> (x**3 - 2).factor()
        (1, [-2 + x^3], [1])

    Over the lazy fields, polynomials split into linear factors (and
    quadratic factors for the pairs of nonreal roots over the real
    fields):

        >>> x = PolynomialRing(RR_tower, "x").gen()
        >>> ((x**2 - 2) * (x**2 + 1)).factor()
        (1, [(-a1 {a1 = sqrt(2)}) + x, (a1 {a1 = sqrt(2)}) + x, 1 + x^2], [1, 1, 1])

    The tower of an element of a lazy field can be inspected the same way:

        >>> x = CC_tower(2).sqrt() + CC_tower(3).sqrt()
        >>> x.tower()      # doctest: +ELLIPSIS
        Tower of degree 4 over Rational field (fmpq)
          a... = 1.414213562  root_2(2)  =  root of  -2 + a...^2  (irreducible)
          a... = 1.732050808  root_2(3)  =  root of  -3 + a...^2  (irreducible)
    """

    def __init__(self, _ptr=None, _owned=True):
        if _ptr is None:
            self._ptr = libgr.gr_tower_heap_init(QQ._ref)
            self._owned = True
        else:
            self._ptr = _ptr
            self._owned = _owned
        self._fields = []
        # (the field contexts of the tower hold references to it: the
        # tower is freed when the Python object and all of them are gone,
        # whatever the order of finalization at interpreter shutdown)
        self._refcount = 1

    def _decrement_refcount(self):
        self._refcount -= 1
        if not self._refcount and getattr(self, "_owned", False):
            libgr.gr_tower_heap_clear(self._ptr)

    def __del__(self):
        if hasattr(self, "_refcount"):
            self._decrement_refcount()

    def __repr__(self):
        arr = ctypes.c_char_p()
        libgr.gr_tower_get_str(ctypes.byref(arr), self._ptr)
        try:
            s = ctypes.cast(arr, ctypes.c_char_p).value.decode("ascii")
        finally:
            libflint.flint_free(arr)
        return s.rstrip("\n")

    def field(self):
        """
        The field of the top level of the tower, as a context object.
        A field object remains valid when algebraic generators are
        adjoined afterwards (its elements then lie in a subfield of the
        new top field, and convert to it), but not when a transcendental
        generator is adjoined, which rebuilds the tower over a new base
        field: field objects created before are then stale, and their
        elements unusable.
        """
        K = TowerField(self)
        self._fields.append(weakref.ref(K))
        return K

    def _rebuilt(self):
        for r in self._fields:
            K = r()
            if K is not None:
                K._stale = True
        self._fields = []

    def degree(self):
        return libgr.gr_tower_degree(self._ptr)

    def length(self):
        """
        The number of algebraic generators.
        """
        return libgr.gr_tower_length_si(self._ptr)

    def num_gens(self):
        return libgr.gr_tower_num_gens_si(self._ptr)

    def gen_names(self):
        return [ctypes.cast(libgr.gr_tower_gen_name(self._ptr, d), ctypes.c_char_p).value.decode("ascii") for d in range(self.num_gens())]

    def gen(self, d):
        """
        The generator with definition order *d* (0-based, in the order
        of adjunction), as an element of the current top field.
        """
        K = self.field()
        x = K(0)
        status = libgr.gr_tower_gen_get(x._ref, self._ptr, d)
        if status:
            raise FlintUnableError("unable to get the generator")
        return x

    def gens(self):
        """
        All generators (algebraic and transcendental) in the order of
        adjunction, as elements of the current top field.
        """
        return [self.gen(d) for d in range(self.num_gens())]

    def _last_gen(self):
        return self.gen(self.num_gens() - 1)

    def _element(self, x):
        return self.field()(x)

    def _check(self, status, what):
        if status:
            if status & GR_DOMAIN: raise FlintDomainError(f"cannot adjoin {what}")
            raise FlintUnableError(f"unable to adjoin {what}")

    def adjoin_sqrt(self, x, name=None):
        """
        Adjoins the principal square root of *x* (an element of the top
        field, or something convertible to it) and returns its name.
        """
        K = self.field()
        x = K(x)
        self._check(libgr.gr_tower_adjoin_root_ui(self._ptr, x._ref, 2, _name_arg(name)), "sqrt(%s)" % x)
        return self._last_gen()

    def adjoin_root(self, x, n, name=None):
        """
        Adjoins the principal *n*-th root of *x*.
        """
        K = self.field()
        x = K(x)
        self._check(libgr.gr_tower_adjoin_root_ui(self._ptr, x._ref, n, _name_arg(name)), "root(%s, %s)" % (x, n))
        return self._last_gen()

    def adjoin_qqbar(self, x, name=None):
        """
        Adjoins the algebraic number *x* (a ``qqbar``, or something
        convertible to one) through its minimal polynomial over QQ.
        """
        x = QQbar(x)
        self._check(libgr.gr_tower_adjoin_qqbar(self._ptr, x._ref, _name_arg(name)), "%s" % x)
        return self._last_gen()

    def adjoin_root_of_unity(self, n, name=None):
        self._check(libgr.gr_tower_adjoin_root_of_unity(self._ptr, n, _name_arg(name)), "exp(2 pi i / %s)" % n)
        return self._last_gen()

    def adjoin_algebraic(self, poly, z, proven=False, name=None):
        """
        Adjoins the root of the polynomial *poly* (over the top field)
        isolated by the complex enclosure *z* (an ``acb``, or something
        convertible). The polynomial need not be known irreducible
        (``proven=False``): the tower is refined later if it is not.
        """
        K = self.field()
        poly = PolynomialRing(K)(poly)
        z = CC(z)
        status = GR_TOWER_STATUS_PROVEN if proven else GR_TOWER_STATUS_DYNAMIC
        self._check(libgr.gr_tower_adjoin_algebraic(self._ptr, poly._ref, z._ref, status, _name_arg(name)), "root of %s near %s" % (poly, z))
        return self._last_gen()

    def adjoin_pi(self, name=None):
        """
        Adjoins pi as a transcendental generator (this rebuilds the
        tower: see :meth:`field`).
        """
        self._check(libgr.gr_tower_adjoin_pi(self._ptr, _name_arg(name)), "pi")
        self._rebuilt()
        return self._last_gen()

    def adjoin_exp(self, x, name=None):
        """
        Adjoins exp(*x*) as a transcendental generator (conjecturally
        transcendental over the tower unless a relation is found; this
        rebuilds the tower: see :meth:`field`).
        """
        K = self.field()
        x = K(x)
        status = libgr.gr_tower_adjoin_exp(self._ptr, x._ref, _name_arg(name))
        self._rebuilt()
        self._check(status, "exp(...)")
        return self._last_gen()

    def adjoin_log(self, x, name=None):
        """
        Adjoins log(*x*) as a transcendental generator (see :meth:`adjoin_exp`).
        """
        K = self.field()
        x = K(x)
        status = libgr.gr_tower_adjoin_log(self._ptr, x._ref, _name_arg(name))
        self._rebuilt()
        self._check(status, "log(...)")
        return self._last_gen()

GR_TOWER_STATUS_PROVEN = 0
GR_TOWER_STATUS_DYNAMIC = 1

def _name_arg(name):
    if name is None:
        return ctypes.c_char_p(None)
    return ctypes.c_char_p(name.encode("ascii"))

class TowerField(gr_ctx):
    r"""
    The field of the top level of a :class:`gr_tower`
    (``gr_ctx_init_tower_field``). See :class:`gr_tower`.
    """
    def __init__(self, tower):
        gr_ctx.__init__(self)
        self._tower = tower
        self._stale = False
        libgr.gr_ctx_init_tower_field(self._ref, tower._ptr)
        tower._refcount += 1
        self._elem_type = gr_tower_field_elem
        assert libgr.gr_ctx_sizeof_elem(self._ref) <= ctypes.sizeof(gr_tower_field_struct)

    def _decrement_refcount(self):
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
            self._tower._decrement_refcount()

    def __call__(self, *args, **kwargs):
        if self._stale:
            raise ValueError("this field object is stale: the tower was rebuilt (transcendental generator adjoined); use tower.field() again")
        return gr_ctx.__call__(self, *args, **kwargs)

    def tower(self):
        return self._tower




PADIC_RADIX_SIGNED = 1     # allow signed units for exact elements
PADIC_RADIX_DECIMAL = 4    # print the unit as a plain decimal integer
PADIC_RADIX_PREC_INF = WORD_MAX    # context precision sentinel for "infinite"

class Qp_padic_radix(gr_ctx):
    r"""
    The field of p-adic numbers Q_p, implemented with radix arithmetic
    arithmetic (padic_radix). A nonzero element is stored canonically as
    u * p^v + O(p^N): the unit u has p-adic valuation 0, v is the valuation,
    and N is the absolute precision (N == +inf marks an exactly represented
    element, printed with no error term).

        >>> Q7 = Qp_padic_radix(7, rel_prec=30)
        >>> d1 = Mat(Q7, 20, 20)().hilbert().det()
        >>> d2 = Q7(Mat(QQ, 20, 20)().hilbert().det())
        >>> d1
        (5138346895024451929101583) * 7^-19 + O(7^11)
        >>> d2
        (5138346895024451929101583) * 7^-19 + O(7^11)
        >>> d1 - d2
        0 + O(7^11)
        >>> d1.overlaps(d2)
        True

    The relative precision ``rel_prec`` bounds N - v and the absolute precision
    ``abs_prec`` bounds N; either may be left infinite. With both infinite the
    structure holds only exactly representable numbers, so an inexact result
    such as 1/3 is reported as not computable rather than silently truncated.
    ``p`` must be a word-size prime.

        >>> Q7 = Qp_padic_radix(7, rel_prec=None)
        >>> Q7
        Radix 7-adic numbers (rel prec inf, abs prec inf)
        >>> Q7(0)
        0
        >>> Q7(5)
        5
        >>> Q7(14)            # 14 = 2 * 7
        (2) * 7^1
        >>> Q7(98)            # 98 = 2 * 7^2
        (2) * 7^2
        >>> Q7(-1)
        -1

    Arithmetic on exactly representable elements stays exact:

        >>> Q7(2) + Q7(3)
        5
        >>> Q7(3) * Q7(7)
        (3) * 7^1
        >>> Q7(7) * Q7(7)
        (1) * 7^2
        >>> Q7(2) - Q7(2)
        0
        >>> Q7(6) / Q7(2)
        3
        >>> Q7(1) / Q7(7)     # 7^-1 is exact
        (1) * 7^-1

    A non-unit denominator gives an infinite expansion, which an exact ring
    cannot represent:

        >>> Q7(1) / Q7(3)
        Traceback (most recent call last):
          ...
        FlintUnableError: ...

    A finite relative precision truncates such expansions to N - v digits,
    recording the error term O(p^N):

        >>> Q7r = Qp_padic_radix(7, rel_prec=8)
        >>> Q7r
        Radix 7-adic numbers (rel prec 8, abs prec inf)
        >>> Q7r(1) / Q7r(2)               # 2^-1 (mod 7^8)
        2882401 + O(7^8)
        >>> Q7r(1) / Q7r(3)               # 3^-1 (mod 7^8)
        3843201 + O(7^8)
        >>> Q7r(1) / Q7r(14)              # (2*7)^-1: valuation -1, 8 digits
        (2882401) * 7^-1 + O(7^7)

    Exactly representable elements stay exact even in a finite-precision ring,
    and a finite *relative* precision does not bound the valuation:

        >>> Q7r(5)
        5
        >>> Q7r(7**10)
        (1) * 7^10
        >>> Q7r(1) / Q7r(7)
        (1) * 7^-1

    A finite *absolute* precision instead bounds N directly: anything at or
    below the horizon p^abs_prec collapses to zero:

        >>> Q7a = Qp_padic_radix(7, rel_prec=None, abs_prec=5)
        >>> Q7a
        Radix 7-adic numbers (rel prec inf, abs prec 5)
        >>> Q7a(7**3)                     # valuation 3 < 5: still exact
        (1) * 7^3
        >>> Q7a(7**5)                     # at the horizon
        0 + O(7^5)
        >>> Q7a(7**10)
        0 + O(7^5)
        >>> Q7a(1) / Q7a(3)               # 3^-1 (mod 7^5)
        11205 + O(7^5)

    """

    def __init__(self, p, rel_prec=10, abs_prec=None, signed=True, decimal=True):
        gr_ctx.__init__(self)
        prec_rel = PADIC_RADIX_PREC_INF if rel_prec is None else int(rel_prec)
        prec_abs = PADIC_RADIX_PREC_INF if abs_prec is None else int(abs_prec)
        flags = 0
        if signed:
            flags |= PADIC_RADIX_SIGNED
        if decimal:
            flags |= PADIC_RADIX_DECIMAL
        p = int(p)
        if p < 0 or p > UWORD_MAX or not libflint.n_is_prime(p):
            raise FlintUnableError("p must be a word-size prime")
        libgr.gr_ctx_init_padic_radix(self._ref, p, prec_rel, prec_abs, flags)
        self._elem_type = padic_radix

DECIMAL_RND_DOWN = 0
DECIMAL_RND_UP = 1
DECIMAL_RND_FLOOR = 2
DECIMAL_RND_CEIL = 3
DECIMAL_RND_NEAR = 4
DECIMAL_RND_NEAR_AWAY = 5
DECIMAL_RND_NEAR_ZERO = 6

_decimal_rnd_names = ["down", "up", "floor", "ceil", "near", "near_away", "near_zero"]

DECIMAL_ALLOW_INF = 1
DECIMAL_ALLOW_NAN = 2
DECIMAL_ALLOW_UNDERFLOW = 4
DECIMAL_SLOPPY_RADIUS = 8
DECIMAL_WRITE_SCIENTIFIC = 16

DECIMAL_PREC_EXACT = WORD_MAX

def _decimal_rnd(rnd):
    if isinstance(rnd, str):
        return _decimal_rnd_names.index(rnd)
    rnd = int(rnd)
    if not 0 <= rnd < len(_decimal_rnd_names):
        raise ValueError("invalid rounding mode")
    return rnd

class gr_decimal_ctx(gr_ctx):
    r"""
    Base class for the decimal floating-point and decimal ball contexts.
    The keyword arguments common to both constructors are:

    * ``prec``: the working precision in decimal digits, or ``None`` for
      exact arithmetic (an operation whose result is not exactly
      representable then fails with ``FlintUnableError``).
    * ``rnd``: the rounding mode, one of ``"down"`` (toward zero),
      ``"up"`` (away from zero), ``"floor"``, ``"ceil"``, ``"near"``
      (ties to even), ``"near_away"`` or ``"near_zero"``.
    * ``inf``, ``nan``: whether infinities and NaN are representable
      (for balls, infinite midpoints and indeterminate values are
      always representable; the flag only decides whether overflow
      with exponent limits produces an infinity).
    * ``underflow``: whether results below the smallest exponent are
      flushed to zero (otherwise they fail).
    * ``exp_limits``: a pair ``(emin, emax)`` bounding the scientific
      exponent `E` (where `10^E \le |x| < 10^{E+1}`) of nonzero values,
      or ``None`` for unbounded exponents.
    * ``scientific``: always print in scientific notation.
    * ``limb_digits``: the number of digits per internal limb (0 selects
      the largest possible value, 19 on 64-bit machines).
    """

    def _init(self, which, prec, rnd, rad_prec, inf, nan, underflow, sloppy_radius, scientific, exp_limits, limb_digits, rnd_im=None):
        gr_ctx.__init__(self)
        flags = 0
        if inf: flags |= DECIMAL_ALLOW_INF
        if nan: flags |= DECIMAL_ALLOW_NAN
        if underflow: flags |= DECIMAL_ALLOW_UNDERFLOW
        if sloppy_radius: flags |= DECIMAL_SLOPPY_RADIUS
        if scientific: flags |= DECIMAL_WRITE_SCIENTIFIC
        prec = DECIMAL_PREC_EXACT if prec is None else int(prec)
        libgr._gr_ctx_init_decimal(self._ref, which, int(limb_digits), prec, _decimal_rnd(rnd), flags)
        if rnd_im is not None:
            libgr.decimal_ctx_set_rnd_im(self._ref, _decimal_rnd(rnd_im))
        if rad_prec is not None:
            libgr.decimal_ctx_set_rad_prec(self._ref, int(rad_prec))
        if exp_limits is not None:
            self.exp_limits = exp_limits

    @property
    def digits(self):
        """
        The working precision in digits (``None`` for exact arithmetic).
        Unlike ``prec``, which is measured in bits, this can be assigned
        exactly.

            >>> R = RealFloat_decfloat(10)
            >>> R.digits, R.prec
            (10, 32)
            >>> R.digits = 3
            >>> R(1) / 3
            0.333
            >>> R.digits = None
            >>> R(1) / 3
            Traceback (most recent call last):
              ...
            FlintUnableError: ...
        """
        p = libgr.decimal_ctx_get_prec(self._ref)
        return None if p == DECIMAL_PREC_EXACT else p

    @digits.setter
    def digits(self, prec):
        libgr.decimal_ctx_set_prec(self._ref, DECIMAL_PREC_EXACT if prec is None else int(prec))

    @property
    def rnd(self):
        """
        The rounding mode, as a string.

            >>> R = RealFloat_decfloat(3)
            >>> R.rnd
            'near'
            >>> R(2) / 3, R(-2) / 3
            (0.667, -0.667)
            >>> R.rnd = "down"
            >>> R(2) / 3, R(-2) / 3
            (0.666, -0.666)
            >>> R.rnd = "floor"
            >>> R(2) / 3, R(-2) / 3
            (0.666, -0.667)
            >>> R.rnd = "ceil"
            >>> R(2) / 3, R(-2) / 3
            (0.667, -0.666)
            >>> R.rnd = "up"
            >>> R(2) / 3, R(-2) / 3
            (0.667, -0.667)
        """
        return _decimal_rnd_names[libgr.decimal_ctx_get_rnd(self._ref)]

    @rnd.setter
    def rnd(self, rnd):
        libgr.decimal_ctx_set_rnd(self._ref, _decimal_rnd(rnd))

    @property
    def rnd_im(self):
        """
        The rounding mode for imaginary parts (complex contexts), as a
        string. Setting ``rnd`` also sets ``rnd_im``.

            >>> C = ComplexFloat_deccfloat(3)
            >>> C.rnd, C.rnd_im
            ('near', 'near')
            >>> C.rnd_im = "floor"
            >>> C
            Complex decimal floating-point numbers (prec 3, rnd near, rnd im floor)
            >>> C("(2 + 2*I) / 3"), C("(-2 - 2*I) / 3")
            ((0.667 + 0.666*I), (-0.667 - 0.667*I))
            >>> C.rnd = "ceil"
            >>> C.rnd, C.rnd_im
            ('ceil', 'ceil')
        """
        return _decimal_rnd_names[libgr.decimal_ctx_get_rnd_im(self._ref)]

    @rnd_im.setter
    def rnd_im(self, rnd):
        libgr.decimal_ctx_set_rnd_im(self._ref, _decimal_rnd(rnd))

    @property
    def exp_limits(self):
        """
        The exponent limits as a pair ``(emin, emax)``, or ``None``.

            >>> R = RealFloat_decfloat(5, exp_limits=(-3, 3), inf=True, underflow=True)
            >>> R.exp_limits
            (-3, 3)
            >>> R("1234"), R("12345")
            (1234, inf)
            >>> R("0.001"), R("0.0001")
            (0.001, 0)
            >>> R.exp_limits = None
            >>> R("12345"), R("0.0001")
            (12345, 0.0001)
        """
        emin = c_slong(); emax = c_slong()
        libgr.decimal_ctx_get_exp_limits(ctypes.byref(emin), ctypes.byref(emax), self._ref)
        if emin.value == WORD_MIN and emax.value == WORD_MAX:
            return None
        return (emin.value, emax.value)

    @exp_limits.setter
    def exp_limits(self, limits):
        if limits is None:
            emin, emax = WORD_MIN, WORD_MAX
        else:
            emin, emax = limits
            emin = WORD_MIN if emin is None else int(emin)
            emax = WORD_MAX if emax is None else int(emax)
        libgr.decimal_ctx_set_exp_limits(self._ref, emin, emax)

    @property
    def flags(self):
        return libgr.decimal_ctx_get_flags(self._ref)

    @flags.setter
    def flags(self, flags):
        libgr.decimal_ctx_set_flags(self._ref, int(flags))

    @property
    def limb_digits(self):
        """
        The number of decimal digits per internal limb.

            >>> RealFloat_decfloat(limb_digits=4).limb_digits
            4
        """
        return libgr.decimal_ctx_get_limb_digits(self._ref)

    def _float_context(self):
        """A real floating-point context with the same settings."""
        return RealFloat_decfloat(prec=self.digits, rnd=self.rnd, inf=bool(self.flags & DECIMAL_ALLOW_INF),
            nan=bool(self.flags & DECIMAL_ALLOW_NAN), underflow=bool(self.flags & DECIMAL_ALLOW_UNDERFLOW),
            scientific=bool(self.flags & DECIMAL_WRITE_SCIENTIFIC), exp_limits=self.exp_limits, limb_digits=self.limb_digits)

    def _complex_float_context(self):
        """A complex floating-point context with the same settings."""
        return ComplexFloat_deccfloat(prec=self.digits, rnd=self.rnd, rnd_im=self.rnd_im, inf=bool(self.flags & DECIMAL_ALLOW_INF),
            nan=bool(self.flags & DECIMAL_ALLOW_NAN), underflow=bool(self.flags & DECIMAL_ALLOW_UNDERFLOW),
            scientific=bool(self.flags & DECIMAL_WRITE_SCIENTIFIC), exp_limits=self.exp_limits, limb_digits=self.limb_digits)

    def _ball_context(self):
        """A real ball context with the same settings."""
        return RealField_decball(prec=self.digits, rad_prec=libgr.decimal_ctx_get_rad_prec(self._ref), rnd=self.rnd,
            inf=bool(self.flags & DECIMAL_ALLOW_INF), nan=bool(self.flags & DECIMAL_ALLOW_NAN),
            underflow=bool(self.flags & DECIMAL_ALLOW_UNDERFLOW), sloppy_radius=bool(self.flags & DECIMAL_SLOPPY_RADIUS),
            scientific=bool(self.flags & DECIMAL_WRITE_SCIENTIFIC), exp_limits=self.exp_limits, limb_digits=self.limb_digits)


class RealFloat_decfloat(gr_decimal_ctx):
    r"""
    Decimal floating-point numbers with a given number of significant
    digits and correctly rounded arithmetic (decfloat). Values are stored
    with leading and trailing zeros stripped, so exactly representable
    numbers print exactly regardless of the precision:

        >>> R = RealFloat_decfloat(10)
        >>> R
        Decimal floating-point numbers (prec 10, rnd near)
        >>> R(1), R(-2), R("0.5"), R("1e100"), R("123.4500")
        (1, -2, 0.5, 1e100, 123.45)
        >>> R(1) / 3
        0.3333333333
        >>> R(2) / 3
        0.6666666667
        >>> R(10) ** 30
        1e30
        >>> R("1e-7") + 1
        1.0000001
        >>> R("1e-10") + 1
        1
        >>> R("1e-10") + 1 - 1
        0
        >>> R(2).sqrt()
        1.414213562
        >>> R(2).sqrt() ** 2
        1.999999999
        >>> R(2).sqrt() ** 2 - 2
        -1e-9

    Numbers print in positional notation when the scientific exponent is
    between -6 and 20 and in scientific notation otherwise; both forms are
    accepted as input, as are arithmetic expressions:

        >>> R("0.000001"), R("0.0000001")
        (0.000001, 1e-7)
        >>> R("1" + "0" * 20), R("1" + "0" * 21)
        (100000000000000000000, 1e21)
        >>> R("1.5e3 * 2 + 1/4")
        3000.25
        >>> RealFloat_decfloat(10, scientific=True)("123.5")
        1.235e2

    The precision and rounding mode can be given in the constructor or
    changed later:

        >>> RealFloat_decfloat(3, rnd="floor")(1) / 3
        0.333
        >>> RealFloat_decfloat(3, rnd="ceil")(1) / 3
        0.334

    Rounding to a given number of digits is done with :meth:`decfloat.round`;
    with ``prec=None`` the context is exact and inexact operations fail:

        >>> X = RealFloat_decfloat(None)
        >>> X("1e20") + X("1e-20")
        100000000000000000000.00000000000000000001
        >>> X(1) / 4
        0.25
        >>> X(1) / 3
        Traceback (most recent call last):
          ...
        FlintUnableError: ...

    Infinities and NaN are not representable unless enabled:

        >>> R(1) / 0
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> R.inf()
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> Rx = RealFloat_decfloat(10, inf=True, nan=True)
        >>> Rx(1) / 0, Rx(-1) / 0, Rx(0) / 0
        (inf, -inf, nan)
        >>> Rx("1e400") ** 3, Rx("inf") - Rx("inf")
        (1e1200, nan)

    Exponents are unbounded by default, but limits on the scientific
    exponent can be imposed; values beyond the limits overflow to infinity
    or underflow to zero if permitted by the ``inf`` and ``underflow`` flags:

        >>> Rl = RealFloat_decfloat(10, exp_limits=(-99, 99))
        >>> Rl("1e99") * 10
        Traceback (most recent call last):
          ...
        FlintUnableError: ...
        >>> Rl = RealFloat_decfloat(10, exp_limits=(-99, 99), inf=True, underflow=True)
        >>> Rl("1e99") * 10, Rl("1e-99") / 10, Rl("-9.999999999e99") * 2
        (inf, 0, -inf)

    Conversions to and from other types round to the context precision
    (binary fractions are converted exactly in an exact context):

        >>> R(QQ(1)/8), R(ZZ(10)**25), R(2**70), R(0.1)
        (0.125, 1e25, 1.180591621e21, 0.1)
        >>> X = RealFloat_decfloat(None)
        >>> X(2**70), X(0.1)
        (1.180591620717411303424e21, 0.1000000000000000055511151231257827021181583404541015625)
        >>> QQ(R("0.125")), ZZ(R("1e5")), float(R("0.125"))
        (1/8, 100000, 0.125)
        >>> QQ(R("0.1")), int(R("3.7")), R("3.7").floor(), R("-3.7").ceil(), R("3.5").nint(), R("2.5").nint()
        (1/10, 3, 3, -3, 4, 2)
        >>> ZZ(R("0.5"))
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> RR(R("0.1")), RF(R("0.1"))
        ([0.100000000000000 +/- 2.23e-17], 0.1000000000000000)
        >>> R(RR("0.1")), R(RF("0.1")), X(RF("0.1"))
        (0.1, 0.1, 0.1000000000000000055511151231257827021181583404541015625)
        >>> R(RR(1) / 3)             # the midpoint is used
        0.3333333333
        >>> X(RR(1) / 3)             # not allowed in an exact context
        Traceback (most recent call last):
          ...
        FlintUnableError: ...
        >>> CC(R("0.1")), RealFloat_arf(20)(R(1) / 3)
        ([0.100000000000000 +/- 2.23e-17], 0.3333335)

    Elementary functions are computed with correct rounding, via arb:

        >>> R.pi()
        3.141592654
        >>> R(1).exp(), R(2).log(), R(1).sin(), R(1).atan()
        (2.718281828, 0.6931471806, 0.8414709848, 0.7853981634)
        >>> RealFloat_decfloat(50).pi()
        3.1415926535897932384626433832795028841971693993751
        >>> R(2) ** R("0.5")
        1.414213562
        >>> R(1).tan(), R(1).sinh(), R(1).cosh(), R(1).tanh(), R("1e-5").expm1(), R("1e-5").log1p()
        (1.557407725, 1.175201194, 1.543080635, 0.761594156, 0.00001000005, 0.00000999995)

    Correct rounding is guaranteed even for tiny arguments, where the
    result lies extremely close to a representable number (a case where
    Ziv's strategy alone would not terminate), and for exact powers,
    which ``arb`` does not detect:

        >>> Rd = RealFloat_decfloat(10, rnd="down")
        >>> Rd("1e-1000000000").sin(), Rd("-1e-1000000000").sin(), Rd("1e-1000000000").exp(), Rd("-1e-1000000000").exp()
        (9.999999999e-1000000001, -9.999999999e-1000000001, 1, 0.9999999999)
        >>> Rd("1e-1000000000").cos(), Rd("1e-1000000000").sinh(), Rd("1e-1000000000").atan()
        (0.9999999999, 1e-1000000000, 9.999999999e-1000000001)
        >>> Rd("1e100") ** Rd("0.5"), Rd(8) ** Rd("0.125"), Rd("0.0001") ** Rd("0.25"), Rd("1.5") ** 3, Rd("1.5") ** -3
        (1e50, 1.296839554, 0.1, 3.375, 0.2962962962)
        >>> Rd(4) ** Rd("1e-1000000000"), Rd("0.25") ** Rd("1e-1000000000")
        (1, 0.9999999999)

    The full table of elementary and special functions of ``arb`` is
    available with correct rounding, including exact special values
    (which ``arb`` does not detect) and asymptotic cases:

        >>> R.gamma(5), R.gamma(R("0.5")), R.zeta(-3), R.zeta(2), R.log10(1000), R.log2(R("0.125")), R.sin_pi(R("0.25")), R.tan_pi(R("0.25")), R.asin_pi(R("0.5"))
        (24, 1.772453851, 0.008333333333, 1.644934067, 3, -3, 0.7071067812, 1, 0.1666666667)
        >>> Rd.gamma(Rd("1e-1000000000")), Rd.gamma(Rd("-1e-1000000000")), Rd.digamma(Rd("1e-1000000000")), Rd.cot(Rd("1e-1000000000"))
        (9.999999999e999999999, -1e1000000000, -1e1000000000, 9.999999999e999999999)
        >>> Rd.tanh(1000000), Rd.tanh(-1000000), Rd.erf(1000000), Rd.erfc(-1000000), Rd.expm1(-1000000), Rd.zeta(1000000), Rd.acot(Rd("1e1000000000"))
        (0.9999999999, -0.9999999999, 0.9999999999, 1.999999999, -0.9999999999, 1, 9.999999999e-1000000001)
        >>> R.euler(), R.catalan(), R.erf(1), R.erfinv(R("0.5")), R.lambertw(1), R.dilog(1), R.bessel_j(1, 2), R.agm(R("0.5"), 1), R.hurwitz_zeta(2, 3), R.polylog(2, R("0.5")), R.atan2(1, 1), R.atan2(0, -1)
        (0.5772156649, 0.9159655942, 0.8427007929, 0.4769362762, 0.5671432904, 1.644934067, 0.5767248078, 0.7283955155, 0.3949340668, 0.5822405265, 0.7853981634, 3.141592654)
        >>> R.fac(20), R.fac(100), R.gamma(QQ(1)/3), R.rising(R("0.5"), 5), R.lambertw(R("-0.25"), -1), R.fresnel_s(1), R.airy_ai(1)
        (2432902008000000000, 9.332621544e157, 2.678938535, 29.53125, -2.153292364, 0.3102683017, 0.1352924163)
        >>> R.gamma(0)
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> R.zeta(1)
        Traceback (most recent call last):
          ...
        FlintDomainError: ...

    Polynomials and matrices work with decimal coefficients (the
    matrix computations are done with floating-point Gaussian elimination,
    without any error control):

        >>> Mat(R)([[1, 2], [3, 4]]).det()
        -2
        >>> Mat(R, 3, 3)().hilbert().inv()
        [[9.00000017, -36.00000101, 30.00000097],
        [-36.00000101, 192.0000056, -180.0000054],
        [30.00000098, -180.0000054, 180.0000052]]
        >>> Mat(RealFloat_decfloat(30), 3, 3)().hilbert().inv()
        [[9.0000000000000000000000000017, -36.0000000000000000000000000101, 30.0000000000000000000000000097],
        [-36.0000000000000000000000000101, 192.000000000000000000000000056, -180.000000000000000000000000054],
        [30.0000000000000000000000000098, -180.000000000000000000000000054, 180.000000000000000000000000052]]
        >>> Mat(R, 4, 4)().hilbert().det()
        1.6534431e-7
        >>> QQ(1) / 6048000
        1/6048000
        >>> Mat(RealFloat_decfloat(None), 2, 2)([[1, 2], [3, 4]]).inv()
        [[-2, 1],
        [1.5, -0.5]]
        >>> Mat(RealFloat_decfloat(None), 2, 2)([[1, 2], [3, 5]]).inv()
        [[-5, 2],
        [3, -1]]
        >>> Mat(RealFloat_decfloat(None), 2, 2)([[1, 2], [3, 6]]).inv()
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> Mat(RealFloat_decfloat(None), 2, 2)([[1, 2], [3, 9]]).inv()
        Traceback (most recent call last):
          ...
        FlintUnableError: ...
        >>> P = PolynomialRing(R)
        >>> P([1, R("0.5")]) ** 2
        1 + x + 0.25*x^2
        >>> P("(x - 1.5) * (x + 0.25)")
        -0.375 - 1.25*x + x^2
        >>> P("(x - 1.5) * (x + 0.25)")(R("1.5"))
        0
        >>> P("x^3 - 2")(R(2).sqrt())
        0.828427123

    """

    def __init__(self, prec=20, rnd="near", inf=False, nan=False, underflow=False, scientific=False, exp_limits=None, limb_digits=0):
        self._init(0, prec, rnd, None, inf, nan, underflow, False, scientific, exp_limits, limb_digits)
        self._elem_type = decfloat


class RealField_decball(gr_decimal_ctx):
    r"""
    Real numbers represented as decimal balls (decball): a midpoint that is
    a decimal floating-point number with a given number of digits and a
    radius with a small fixed number of digits (``rad_prec``, between 1 and
    9). Operations produce balls that are guaranteed to contain the exact
    result.

        >>> R = RealField_decball(10)
        >>> R
        Decimal balls (prec 10, rad prec 4)
        >>> R(1), R("0.5"), R("1e100")
        (1, 0.5, 1e100)
        >>> R(1) / 3
        [0.3333333333 +/- 3.334e-11]
        >>> R(2).sqrt()
        [1.414213562 +/- 3.731e-10]
        >>> R(2).sqrt() ** 2
        [1.999999999 +/- 1.113e-9]
        >>> R.pi()
        [3.141592654 +/- 4.104e-10]
        >>> R.pi().sin()
        [-4.102067611e-10 +/- 4.106e-10]
        >>> R(10) ** 20 + 1
        [100000000000000000000 +/- 1]
        >>> RealField_decball(21)(10) ** 20 + 1
        100000000000000000001
        >>> (R(1) / 3) * 3
        [0.9999999999 +/- 1.001e-10]
        >>> (R(1) / 3) * 3 == 1
        Traceback (most recent call last):
          ...
        Undecidable: ...
        >>> (R(1) / 3) * 3 == 2
        False
        >>> R(1) / 3 == R(1) / 3
        Traceback (most recent call last):
          ...
        Undecidable: ...
        >>> R(1) / 4 == R(2) / 8
        True
        >>> (R(1) / 3) * 3 != 2, R(1) / 3 < R(1) / 2, R(1) / 3 > R(1) / 2
        (True, True, False)
        >>> R(1) / 3 < R(1) / 3
        Traceback (most recent call last):
          ...
        Undecidable: ...

    Balls can be written as ``mid +/- rad`` or ``[mid +/- rad]``, which
    is also the output format, and as arbitrary arithmetic expressions:

        >>> R("[1.5 +/- 0.01]")
        [1.5 +/- 0.01]
        >>> R("1.5 +/- 0.01")
        [1.5 +/- 0.01]
        >>> R("+/- 1e-30")
        [0 +/- 1e-30]
        >>> R("(1 +/- 0.1) * (2 +/- 0.1)")
        [2 +/- 0.31]
        >>> R("1.23456789012345 +/- 1e-20")
        [1.23456789 +/- 1.236e-10]
        >>> R("1.2345678901 +/- 1e-20")
        [1.23456789 +/- 1.001e-10]
        >>> R("123456789012345 +/- 1e-20")
        [123456789000000 +/- 12360]
        >>> R("1e100 +/- 1e50"), R("[-1e-100 +/- 1e-150]")
        ([1e100 +/- 1e50], [-1e-100 +/- 1e-150])
        >>> R("[1.5 +/- 0.01]") + R("[2.5 +/- 0.02]")
        [4 +/- 0.03]
        >>> R("1 +/- 0.001").contains(R("1.0005 +/- 0.0001")), R("1 +/- 0.001").contains(R("1.0005 +/- 0.001"))
        (True, False)
        >>> R("1 +/- 0.001").overlaps(R("1.0005 +/- 0.001")), R("1 +/- 0.001").overlaps(R("1.01"))
        (True, False)
        >>> R("1 +/- 0.001").contains(1), R("1 +/- 0.001").contains(QQ(1)/2)
        (True, False)

    The radius precision can be chosen between 1 and 9 digits; radii
    are always rounded up. By default the actual rounding error of each
    operation is tracked (rounded up to the radius precision); with
    ``sloppy_radius=True`` it is instead bounded by half an ulp of the
    midpoint (or a full ulp in directed rounding modes), which is
    marginally cheaper:

        >>> RealField_decball(10, rad_prec=1)(1) / 3
        [0.3333333333 +/- 4e-11]
        >>> RealField_decball(10, rad_prec=9)(1) / 3
        [0.3333333333 +/- 3.33333334e-11]
        >>> RealField_decball(10, rad_prec=2)(1) / 3
        [0.3333333333 +/- 3.4e-11]
        >>> RealField_decball(10, rad_prec=9, sloppy_radius=True)(1) / 3
        [0.3333333333 +/- 5e-11]
        >>> RealField_decball(10)("1e-25") + 1
        [1 +/- 1e-25]
        >>> RealField_decball(10, sloppy_radius=True)("1e-25") + 1
        [1 +/- 5e-10]

    The midpoint and radius are available separately, as decimal
    floating-point numbers, and the radius can be enlarged:

        >>> x = R(1) / 3
        >>> x.mid(), x.rad()
        (0.3333333333, 3.334e-11)
        >>> x.mid().parent()
        Decimal floating-point numbers (prec 10, rnd near)
        >>> x.is_exact(), R(1).is_exact()
        (False, True)
        >>> x.add_error(R("0.001"))
        [0.3333333333 +/- 0.001001]
        >>> R(1).add_error_10exp(-5)
        [1 +/- 1e-5]
        >>> R("[1.234567890 +/- 0.01]").trim()
        [1.234568 +/- 0.01001]
        >>> R("[1.234567890 +/- 0.01]").rel_accuracy_digits()
        2

    Conversions to and from ``arb`` balls are rigorous:

        >>> RR(R(1) / 3)
        [0.3333333333 +/- 3.34e-11]
        >>> R(RR(1) / 3)
        [0.3333333333 +/- 3.335e-11]
        >>> R(RR.pi())
        [3.141592654 +/- 4.104e-10]
        >>> R(RR("[1.5 +/- 1e-20]"))
        [1.5 +/- 1.001e-20]
        >>> RealField_decball(30)(RR.pi())
        [3.14159265358979311599796346854 +/- 2.222e-16]
        >>> RealField_decball(30)(RealField_arb(200).pi())
        [3.14159265358979323846264338328 +/- 4.973e-31]
        >>> RealField_arb(200)(RealField_decball(30).pi())
        [3.141592653589793238462643383280 +/- 4.98e-31]
        >>> RR(R("[1e100 +/- 1e90]")), R(RR("[1e100 +/- 1e90]"))
        ([1.00000000e+100 +/- 1.01e+90], [1e100 +/- 1.002e90])
        >>> CC(R("[1 +/- 0.1]")), RF(R("[1 +/- 0.1]")), QQ(R("[0.5 +/- 0]")), ZZ(R("1e10"))
        ([1e+0 +/- 0.101], 1.000000000000000, 1/2, 10000000000)
        >>> QQ(R("[0.5 +/- 0.1]"))
        Traceback (most recent call last):
          ...
        FlintUnableError: ...

    Elementary and special functions are computed via arb:

        >>> R(1).exp(), R(2).log(), R(1).sin(), R(1).atan()
        ([2.718281828 +/- 4.592e-10], [0.6931471806 +/- 4.007e-11], [0.8414709848 +/- 7.898e-12], [0.7853981634 +/- 2.553e-12])
        >>> R("0.5").gamma(), R(2).zeta()
        ([1.772453851 +/- 9.45e-11], [1.644934067 +/- 1.519e-10])
        >>> R(2) ** R("0.5"), R(2) ** (R(1)/2) - R(2).sqrt()
        ([1.414213562 +/- 3.732e-10], [0 +/- 7.463e-10])
        >>> RealField_decball(40).pi()
        [3.141592653589793238462643383279502884197 +/- 1.695e-40]
        >>> R("[1 +/- 0.001]").exp()
        [2.718283188 +/- 0.00272]
        >>> R("[1000 +/- 1]").sin()
        [0.4867696208 +/- 0.5134]

    Polynomials and matrices over decimal balls, including string
    conversions in both directions:

        >>> P = PolynomialRing(R)
        >>> f = P([R(1)/3, R(1)/7, 1])
        >>> f
        [0.3333333333 +/- 3.334e-11] + [0.1428571429 +/- 4.286e-11]*x + x^2
        >>> P(str(f))
        [0.3333333333 +/- 3.334e-11] + [0.1428571429 +/- 4.286e-11]*x + x^2
        >>> P("(x + [1 +/- 0.001])^2")
        [1 +/- 0.002001] + [2 +/- 0.002]*x + x^2
        >>> P("x^2 - 2")(R(2).sqrt())
        [-1e-9 +/- 1.113e-9]
        >>> M = Mat(R, 2, 2)([[R(1)/3, 2], [3, 4]])
        >>> M
        [[[0.3333333333 +/- 3.334e-11], 2],
        [3, 4]]
        >>> str(Mat(R, 2, 2)(str(M))) == str(M)
        True
        >>> M.det()
        [-4.666666667 +/- 3.334e-10]
        >>> Mat(R, 4, 4)().hilbert().det()
        [1.6534431e-7 +/- 5.939e-12]
        >>> Mat(RealField_decball(30), 4, 4)().hilbert().det()
        [1.653439153439153439153439197e-7 +/- 4.824e-32]
        >>> Mat(RealField_decball(30), 4, 4)().hilbert().inv()[3, 3]
        [2800.00000000000000000000001954 +/- 9.595e-23]

    """

    def __init__(self, prec=20, rad_prec=4, rnd="near", inf=False, nan=False, underflow=False, sloppy_radius=False, scientific=False, exp_limits=None, limb_digits=0):
        self._init(1, prec, rnd, rad_prec, inf, nan, underflow, sloppy_radius, scientific, exp_limits, limb_digits)
        self._elem_type = decball

    @property
    def rad_prec(self):
        """
        The number of digits used for radii. Like the precision, it can be
        changed at any time; existing balls keep their radii, which are
        rounded to the new radius precision by subsequent operations.

            >>> R = RealField_decball(10, rad_prec=2)
            >>> R.rad_prec
            2
            >>> x = R("1 +/- 0.123456")
            >>> x
            [1 +/- 0.13]
            >>> R.rad_prec = 4
            >>> y = R("1 +/- 0.123456")
            >>> x, y, x + y, x.union(y), x.contains(y)
            ([1 +/- 0.13], [1 +/- 0.1235], [2 +/- 0.2535], [1 +/- 0.13], True)
            >>> R.rad_prec = 1
            >>> x + y, x.union(y), x.contains(y)
            ([2 +/- 0.4], [1 +/- 0.2], True)
        """
        return libgr.decimal_ctx_get_rad_prec(self._ref)

    @rad_prec.setter
    def rad_prec(self, rad_prec):
        libgr.decimal_ctx_set_rad_prec(self._ref, int(rad_prec))


class ComplexFloat_deccfloat(gr_decimal_ctx):
    r"""
    Complex numbers represented as pairs of decimal floating-point numbers
    (deccfloat), with correctly rounded arithmetic: each operation rounds
    the real and imaginary parts of the exact result. The rounding mode
    for imaginary parts (``rnd_im``) can differ from that of real parts.

        >>> C = ComplexFloat_deccfloat(10)
        >>> C
        Complex decimal floating-point numbers (prec 10, rnd near)
        >>> C(1), C.i(), C("1 + 2*I"), C("(1.5 - 2*I)"), C("-3*I")
        (1, 1*I, (1 + 2*I), (1.5 - 2*I), -3*I)
        >>> C("(1 + 2*I) * (3 - 4*I)"), C("(1 + 2*I) / (3 - 4*I)")
        ((11 + 2*I), (-0.2 + 0.4*I))
        >>> C("1 + I") / 3
        (0.3333333333 + 0.3333333333*I)
        >>> C("(1 + I) / 3") * 3
        (0.9999999999 + 0.9999999999*I)
        >>> C.i() ** 2, C.i() ** 3, C("1+I") ** 100
        (-1, -1*I, -1125899907000000)
        >>> ComplexFloat_deccfloat(None)("(1+I)^100")
        -1125899906842624

    Square roots and integer powers are exact when the result is a
    Gaussian decimal number, and otherwise correctly rounded
    componentwise:

        >>> C("3 + 4*I").sqrt(), C("-4").sqrt(), C("2*I").sqrt(), C("1 + I").sqrt()
        ((2 + 1*I), 2*I, (1 + 1*I), (1.098684113 + 0.4550898606*I))
        >>> C("3 + 4*I").abs(), C("1 + I").abs(), C("3 + 4*I").sgn(), C("3 + 4*I").arg()
        (5, 1.414213562, (0.6 + 0.8*I), 0.927295218)
        >>> C("3 + 4*I") ** -3, C("3 + 4*I") ** C("0.5"), C("-8") ** (QQ(1) / 3)
        ((-0.007488 - 0.002816*I), (2 + 1*I), (1 + 1.732050807*I))
        >>> C(2) ** C("1 + I"), C("1 + I") ** C("1 + I")
        ((1.538477803 + 1.277922553*I), (0.2739572538 + 0.5837007588*I))
        >>> C("1 + I").re(), C("1 + I").im(), C("1 + I").conj(), C("1 + I").real(), C("1 + I").imag()
        (1, 1, (1 - 1*I), 1, 1)
        >>> C("1 + I").real().parent()
        Decimal floating-point numbers (prec 10, rnd near)

    Elementary and special functions are computed with correct rounding
    of both parts, via acb. Real arguments are handled by the real
    functions (so the imaginary part is exactly zero), purely imaginary
    arguments are reduced to real functions where possible, and tiny
    arguments are handled with Taylor expansions:

        >>> C.i().exp(), C(1).exp(), C("1 + I").exp()
        ((0.5403023059 + 0.8414709848*I), 2.718281828, (1.46869394 + 2.287355287*I))
        >>> C(-1).log(), C("-1e-100").log(), C.i().log(), C("1 + I").log()
        (3.141592654*I, (-230.2585093 + 3.141592654*I), 1.570796327*I, (0.3465735903 + 0.7853981634*I))
        >>> C(2).acos(), C(-2).acos(), C("0.5").acosh(), C(-2).acosh(), C(2).asin(), C(2).atanh()
        (1.316957897*I, (3.141592654 - 1.316957897*I), 1.047197551*I, (1.316957897 + 3.141592654*I), (1.570796327 - 1.316957897*I), (0.5493061443 - 1.570796327*I))
        >>> C("3*I").sin(), C("3*I").cos(), C("3*I").tan(), C("3*I").atan(), C("0.5*I").atan()
        (10.01787493*I, 10.067662, 0.9950547537*I, (1.570796327 + 0.3465735903*I), 0.5493061443*I)
        >>> C("1 + I").gamma(), C("1 + I").zeta(), C("1 + I").erf(), C("2*I").erf(), C.lambertw(C("1 + I"))
        ((0.4980156681 - 0.1549498283*I), (0.5821580598 - 0.9268485643*I), (1.316151282 + 0.1904534692*I), 18.56480241*I, (0.6569660692 + 0.3254503394*I))
        >>> C.exp_pi_i(QQ(1)/2), C.exp_pi_i(QQ(-3)/2), C.exp_pi_i(QQ(1)/4), C.exp_pi_i(C("1 + I"))
        (1*I, 1*I, (0.7071067812 + 0.7071067812*I), -0.04321391826)
        >>> C.gamma(5), C.zeta(2), C.pi(), C.bessel_j(1, C.i()), C.agm(C("1 + I"), 2), C.hurwitz_zeta(2, C("1+I"))
        (24, 1.644934067, 3.141592654, 0.565159104*I, (1.527316275 + 0.5710047826*I), (0.4630000966 - 0.7942335428*I))
        >>> C.rising(C("1 + I"), 3), C.fac(10), C.lambertw(C("1 + I"), -1)
        (10*I, 3628800, (-0.9869695732 - 3.663857003*I))

    Tiny nonreal arguments are handled without Ziv's loop (which would
    not terminate in the directed rounding modes):

        >>> Cd = ComplexFloat_deccfloat(10, rnd="down")
        >>> z = Cd("1e-1000000000 + 1e-1000000000*I")
        >>> z.sin(), z.exp(), z.gamma()
        ((1e-1000000000 + 9.999999999e-1000000001*I), (1 + 1e-1000000000*I), (4.999999999e999999999 - 4.999999999e999999999*I))
        >>> z.cos(), z.log1p(), z.tan(), z.atan()
        ((0.9999999999 - 9.999999999e-2000000001*I), (9.999999999e-1000000001 + 9.999999999e-1000000001*I), (9.999999999e-1000000001 + 1e-1000000000*I), (1e-1000000000 + 9.999999999e-1000000001*I))
        >>> Cd.acot(1 / z)
        (1e-1000000000 + 9.999999999e-1000000001*I)

    Comparisons of nonreal numbers fail, except for equality and
    comparisons of absolute values:

        >>> C("1 + I") == C("1 + I"), C("1 + I") != C("1 - I"), C(1) < C(2)
        (True, True, True)
        >>> C("1 + I") < C(2)
        Traceback (most recent call last):
          ...
        ValueError: ...
        >>> abs(C("1 + I")) < abs(C("1.5")), abs(C("3 + 4*I")) == abs(C(-5))
        (True, True)

    Conversions:

        >>> C(CC("1 + 2*I")), C(CF("0.5 - I")), C(RR_decball("1 +/- 0.1")), C(RF_decfloat(1) / 3)
        ((1 + 2*I), (0.5 - 1*I), 1, 0.3333333333)
        >>> x = QQbar(2).sqrt() + QQbar("0.0023") * QQbar.i()
        >>> ComplexFloat_deccfloat(2, rnd="floor")(x), ComplexFloat_deccfloat(3, rnd="ceil")(QQbar("0.125") + x)
        ((1.4 + 0.0023*I), (1.54 + 0.0023*I))
        >>> CC(C("1 + 2*I")), CF(C("1 + 2*I")), ZZ(C(3)), QQ(C("0.5")), RF_decfloat(C(2))
        ((1.000000000000000 + 2.000000000000000*I), (1.000000000000000 + 2.000000000000000*I), 3, 1/2, 2)
        >>> ZZ(C("1 + 2*I"))
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> C(1) / 0
        Traceback (most recent call last):
          ...
        FlintDomainError: ...
        >>> Cx = ComplexFloat_deccfloat(10, inf=True, nan=True)
        >>> Cx(1) / 0, Cx.i() / 0, Cx("inf*I"), Cx("inf") * Cx.i(), Cx(0) / 0
        (inf, inf*I, inf*I, inf*I, (nan + nan*I))

    Polynomials and matrices:

        >>> P = PolynomialRing(C)
        >>> P("(x - I) * (x + I)")
        1 + x^2
        >>> P("x^2 + 1")(C.i())
        0
        >>> Mat(C)([[1, C.i()], [C.i(), 1]]).det()
        2
        >>> Mat(C)([[1, C.i()], [C.i(), 1]]).inv()
        [[0.5, -0.5*I],
        [-0.5*I, 0.5]]

    """

    def __init__(self, prec=20, rnd="near", rnd_im=None, inf=False, nan=False, underflow=False, scientific=False, exp_limits=None, limb_digits=0):
        self._init(2, prec, rnd, None, inf, nan, underflow, False, scientific, exp_limits, limb_digits, rnd_im=rnd_im)
        self._elem_type = deccfloat


class ComplexField_deccball(gr_decimal_ctx):
    r"""
    Complex numbers represented as pairs of decimal balls (deccball):
    rectangular enclosures with a real and an imaginary ball, like ``acb``.

        >>> C = ComplexField_deccball(10)
        >>> C
        Complex decimal balls (prec 10, rad prec 4)
        >>> C(1), C.i(), C("1 + 2*I")
        (1, 1*I, (1 + 2*I))
        >>> C("1 + I") / 3
        ([0.3333333333 +/- 3.334e-11] + [0.3333333333 +/- 3.334e-11]*I)
        >>> C("1 + I") / 3 * 3
        ([0.9999999999 +/- 1.001e-10] + [0.9999999999 +/- 1.001e-10]*I)
        >>> C("1 + I") / 3 * 3 == C("1 + I")
        Traceback (most recent call last):
          ...
        Undecidable: ...
        >>> C("1 + I") / 3 * 3 == 1, C("1 + I") / 3 * 3 != 1
        (False, True)
        >>> C.i() ** 2, C.i() ** 3, C("3 + 4*I").sqrt(), C("-4").sqrt()
        (-1, -1*I, (2 + 1*I), 2*I)
        >>> C("[1 +/- 0.1] + [2 +/- 0.2]*I")
        ([1 +/- 0.1] + [2 +/- 0.2]*I)
        >>> C("[1 + 2*I +/- 0.01]")
        ([1 +/- 0.01] + 2*I)
        >>> C("([1 +/- 0.1] + [2 +/- 0.2]*I) * (3 - I)")
        ([5 +/- 0.5] + [5 +/- 0.7]*I)
        >>> C("1 + I").sqrt(), C("1 + I").abs(), C("1 + I").arg(), C("1 + I").sgn()
        (([1.098684113 +/- 4.68e-10] + [0.4550898606 +/- 3.779e-11]*I), [1.414213562 +/- 3.732e-10], [0.7853981634 +/- 2.553e-12], ([0.7071067812 +/- 1.347e-11] + [0.7071067812 +/- 1.347e-11]*I))

    Elementary and special functions are computed via acb, with real
    arguments handled by the real functions:

        >>> C.i().exp(), C(1).exp(), C(-1).log(), C("1 + I").gamma()
        (([0.5403023059 +/- 3.188e-11] + [0.8414709848 +/- 7.898e-12]*I), [2.718281828 +/- 4.592e-10], [3.141592654 +/- 4.104e-10]*I, ([0.4980156681 +/- 1.837e-11] + [-0.1549498283 +/- 1.812e-12]*I))
        >>> C.gamma(5), C.zeta(C("0.5 + 14.13472514*I"))
        (24, ([2.163160215e-10 +/- 1.056e-17] + [-1.358779596e-9 +/- 1.226e-17]*I))
        >>> ComplexField_deccball(30).zeta(ComplexField_deccball(30)("0.5 + 14.134725141734693790*I"))
        ([5.70192445022535683266896974805e-20 +/- 2.741e-37] + [-3.58163883768763925565723762546e-19 +/- 3.021e-37]*I)

    The parts, midpoints and radii are available separately:

        >>> x = C("1 + I") / 3
        >>> x.mid(), x.real(), x.imag()
        ((0.3333333333 + 0.3333333333*I), [0.3333333333 +/- 3.334e-11], [0.3333333333 +/- 3.334e-11])
        >>> x.mid().parent(), x.real().parent()
        (Complex decimal floating-point numbers (prec 10, rnd near), Decimal balls (prec 10, rad prec 4))
        >>> x.is_exact(), C(1).is_exact(), x.rel_accuracy_digits(), C(1).rel_accuracy_digits() is None
        (False, True, 10, True)
        >>> x.contains(QQ(1)/3 + QQ(1)/3 * C.i()), x.contains(1), x.overlaps(C("0.3333 + 0.3333*I")), x.overlaps(C("[0.3333 +/- 0.0001] + [0.3333 +/- 0.0001]*I"))
        (True, False, False, True)
        >>> C(1).add_error(C("0.001")), C(1).add_error_10exp(-5), C("1 + I").add_error(C("0.1 + 0.2*I"))
        ([1 +/- 0.001], ([1 +/- 1e-5] + [0 +/- 1e-5]*I), ([1 +/- 0.1] + [1 +/- 0.2]*I))
        >>> (C(1) / 3 + C("+/- 0.001")).trim()
        [0.3333333 +/- 0.001002]

    Conversions:

        >>> C(CC("1 + 2*I")), C(CC(1) / 3), C(RR(1) / 3 * CC.i()), C(RF_decfloat(1) / 3), C(RR_decball(1) / 3)
        ((1 + 2*I), [0.3333333333 +/- 3.335e-11], [0.3333333333 +/- 3.335e-11]*I, [0.3333333333 +/- 3.334e-11], [0.3333333333 +/- 3.335e-11])
        >>> x = QQbar(2).sqrt() + QQbar("0.0023") * QQbar.i()
        >>> ComplexField_deccball(1)(x), ComplexField_deccball(2)(x), ComplexField_deccball(2)(QQbar.exp_pi_i(QQ(1) / 3))
        (([1 +/- 0.4144] + [0.002 +/- 3.002e-4]*I), ([1.4 +/- 0.01423] + 0.0023*I), (0.5 + [0.87 +/- 0.003976]*I))
        >>> RealField_decball(5)(QQbar(2).sqrt()), RealFloat_decfloat(5)(QQbar(2).sqrt())
        ([1.4142 +/- 1.358e-5], 1.4142)
        >>> CC(C("1 + I") / 3), CF(C("1 + I") / 3), ZZ(C(3)), QQ(C("0.5"))
        (([0.3333333333 +/- 3.34e-11] + [0.3333333333 +/- 3.34e-11]*I), (0.3333333333000000 + 0.3333333333000000*I), 3, 1/2)
        >>> ZZ(C("3 +/- 0.1"))
        Traceback (most recent call last):
          ...
        FlintUnableError: ...
        >>> ZZ(C("3 + I"))
        Traceback (most recent call last):
          ...
        FlintDomainError: ...

    Polynomials and matrices:

        >>> P = PolynomialRing(C)
        >>> f = P("(x - I) * (x + [1 +/- 0.001])")
        >>> f
        ([-1 +/- 0.001]*I) + ([1 +/- 0.001] - 1*I)*x + x^2
        >>> P(str(f))
        ([-1 +/- 0.001]*I) + ([1 +/- 0.001] - 1*I)*x + x^2
        >>> f(C.i())
        [0 +/- 0.002]*I
        >>> Mat(C)([[1, C.i()], [C.i(), C(1) / 3]]).det()
        [1.333333333 +/- 3.334e-10]

    """

    def __init__(self, prec=20, rad_prec=4, rnd="near", inf=False, nan=False, underflow=False, sloppy_radius=False, scientific=False, exp_limits=None, limb_digits=0):
        self._init(3, prec, rnd, rad_prec, inf, nan, underflow, sloppy_radius, scientific, exp_limits, limb_digits)
        self._elem_type = deccball

    @property
    def rad_prec(self):
        """
        The number of digits used for radii.

            >>> ComplexField_deccball(10, rad_prec=2).rad_prec
            2
        """
        return libgr.decimal_ctx_get_rad_prec(self._ref)

    @rad_prec.setter
    def rad_prec(self, rad_prec):
        libgr.decimal_ctx_set_rad_prec(self._ref, int(rad_prec))


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

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
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

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
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

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
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


class radix_integer(gr_elem):

    _struct_type = radix_integer_struct

    def size_limbs(self):
        """
            >>> Z = IntegerRing_radix_integer(10, 6)
            >>> Z(10**20).size_limbs()
            4
        """
        libflint.radix_integer_size_limbs.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
        libflint.radix_integer_size_limbs.restype = c_slong
        return libflint.radix_integer_size_limbs(self._ref, libgr.gr_ctx_data_as_ptr(self._ctx))

    def size_digits(self):
        """
            >>> Z = IntegerRing_radix_integer(10, 6)
            >>> Z(10**20).size_digits()
            21
            >>> Z(10**20-1).size_digits()
            20

        """
        libflint.radix_integer_size_digits.argtypes = (ctypes.c_void_p, ctypes.c_void_p)
        libflint.radix_integer_size_digits.restype = c_slong
        return libflint.radix_integer_size_digits(self._ref, libgr.gr_ctx_data_as_ptr(self._ctx))

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

class gr_tower_lazy(gr_elem):
    _struct_type = gr_tower_lazy_elem_struct

    @staticmethod
    def _default_context():
        return CC_tower

    def tower(self):
        """
        The tower (a :class:`gr_tower`, owned by the context) in which
        this element is represented. (The tower object keeps the element
        alive: the context collects towers in which no element lives.)
        """
        level = c_slong()
        ptr = libgr.gr_tower_lazy_get_tower(ctypes.byref(level), self._ref, self._ctx)
        t = gr_tower(_ptr=ptr, _owned=False)
        t._keep = self
        return t

class gr_tower_field_elem(gr_elem):
    _struct_type = gr_tower_field_struct

    @staticmethod
    def _default_context():
        return None

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

    def secondary_zeta(self):
        """
        Secondary zeta function (ad-hoc wrapper for testing).

            >>> CC(0.5+10j).secondary_zeta()
            ([0.1725005546943535 +/- 4.73e-17] + [-0.1680692210280708 +/- 4.12e-17]*I)
            >>> CC(2).secondary_zeta()
            [0.02310499311541897 +/- 3.75e-18]
            >>> CC("-20.001").secondary_zeta()
            [1.32e+18 +/- 6.41e+15]
            >>> ComplexField_acb(prec=128)("-20.001").secondary_zeta()
            [1.31697271383159e+18 +/- 8.21e+3]
            >>> [raises(lambda: CC(s).secondary_zeta(), FlintUnableError)
            ...     for s in [1, -1, -3, "1 +/- 0.001", "-5 +/- 0.001"]]
            [True, True, True, True, True]

        Dyadic values:

            >>> for n in range(10):
            ...     print(CC(-2*n).secondary_zeta())
            ... 
            0.8750000000000000
            -0.2812500000000000
            0.02343750000000000
            -0.1347656250000000
            -0.6723632812500000
            -6.168090820312500
            [-82.48159790039063 +/- 5.00e-15]
            [-1521.003639221191 +/- 4.07e-13]
            [-36986.37416267395 +/- 1.96e-13]
            [-1146735.990261555 +/- 2.82e-10]
        """
        C = self.parent()
        prec = C.prec
        res = C()
        libflint.acb_dirichlet_secondary_zeta(res._ref, self._ref, prec)
        if not libflint.acb_is_finite(res._ref):
            raise FlintUnableError("unable")
        return res

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
    """
        >>> RealFloat_nfloat(128)
        Floating-point numbers with prec = 128 (nfloat)
        >>> RealFloat_nfloat(128).pi()
        3.14159265358979323846264338327950288420
        >>> RealFloat_nfloat(10000)
        Traceback (most recent call last):
          ...
        FlintUnableError: precision out of range for nfloat
    """
    def __init__(self, prec=128):
        gr_ctx.__init__(self)
        if libflint.nfloat_ctx_init(self._ref, prec, 0) != GR_SUCCESS:
            raise FlintUnableError("precision out of range for nfloat")
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
    """
        >>> ComplexFloat_nfloat_complex(128)
        Complex floating-point numbers with prec = 128 (nfloat_complex)
        >>> ComplexFloat_nfloat_complex(128).i()
        1.00000000000000000000000000000000000000*I
        >>> ComplexFloat_nfloat_complex(10000)
        Traceback (most recent call last):
          ...
        FlintUnableError: precision out of range for nfloat_complex
    """
    def __init__(self, prec=128):
        gr_ctx.__init__(self)
        if libflint.nfloat_complex_ctx_init(self._ref, prec, 0) != GR_SUCCESS:
            raise FlintUnableError("precision out of range for nfloat_complex")
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
    """

        >>> IntegersMod_mpn_mod(10**20 + 1)
        Integers mod 100000000000000000001 (mpn)
        >>> IntegersMod_mpn_mod(10**1000)
        Traceback (most recent call last):
          ...
        FlintUnableError: n is not in range for the mpn_mod implementation

    """
    def __init__(self, n, n_is_prime=None):
        n = self._as_fmpz(n)
        gr_ctx.__init__(self)
        if libgr.gr_ctx_init_mpn_mod(self._ref, n._ref) != GR_SUCCESS:
            raise FlintUnableError("n is not in range for the mpn_mod implementation")
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
    """
    Finite field (fq representation).

        >>> K = FiniteField_fq(5, 3)
        >>> K
        GF(5^3) (fq)
        >>> (1 + K.gen()) ** 10
        4*a^2+2
    """

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
    """
    Finite field (fq_nmod representation).

        >>> K = FiniteField_fq_nmod(5, 3)
        >>> K
        GF(5^3) (fq_nmod)
        >>> (1 + K.gen()) ** 10
        4*a^2+2
    """

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
    """
    Finite field (Zech logarithm representation).

        >>> K = FiniteField_fq_zech(5, 3)
        >>> K
        GF(5^3) (fq_zech)
        >>> (1 + K.gen()) ** 10
        a^92
    """

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

    def resultant(self, other, algorithm=None):
        """
            >>> I, x = PolynomialRing(ZZi, "x").gens(recursive=True)
            >>> f = ((2+3*I) + x)**3
            >>> g = ((3+4*I) + x)**3
            >>> f.resultant(g)
            (16+16*I)
            >>> g.resultant(f)
            (-16-16*I)
            >>> (f * (x + 1)).resultant(g * (x + 1))
            0
            >>> f.resultant(g, algorithm="subresultant")
            (16+16*I)
            >>> f.resultant(g, algorithm="sylvester")
            (16+16*I)
            >>> f.resultant(g, algorithm="euclidean")
            Traceback (most recent call last):
              ...
            ValueError
            >>> Kx = PolynomialRing(Fraction_gr_fraction(ZZi), "x")
            >>> Kx(f).resultant(g, algorithm="euclidean")
            ((16+16*I)) / (1)

        """
        Rx = self.parent()
        R = Rx._coefficient_ring
        # fixme:
        other = Rx(other)
        res = R()
        if algorithm is None:
            status = libgr.gr_poly_resultant(res._ref, self._ref, other._ref, R._ref)
        elif algorithm == "euclidean":
            status = libgr.gr_poly_resultant_euclidean(res._ref, self._ref, other._ref, R._ref)
        elif algorithm == "subresultant":
            status = libgr.gr_poly_resultant_subresultant(res._ref, self._ref, other._ref, R._ref)
        elif algorithm == "sylvester":
            status = libgr.gr_poly_resultant_sylvester(res._ref, self._ref, other._ref, R._ref)
        else:
            raise ValueError
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return res

    # todo: want gr_xgcd for generic elements
    def xgcd(self, other, algorithm=None):
        """
            >>> x = ZZx.gen()
            >>> f = (x+1)**2; g = (x-1)**3
            >>> G, S, T = f.xgcd(g)
            >>> (G, S, T); G == S*f + T*g
            (64, 12*x^2-40*x+44, -12*x-20)
            True
            >>> f = (x**2 + 2) * (x+1)**3; g = (x**2 + 2) * (x-1)**2
            >>> f.gcd(g)
            x^2+2
            >>> G, S, T = f.xgcd(g)
            >>> G
            64*x^2+128
            >>> S
            -12*x+20
            >>> T
            12*x^2+40*x+44
            >>> G == S*f + T*g
            True
            >>> f.xgcd(g, algorithm="subresultant")
            (64*x^2+128, -12*x+20, 12*x^2+40*x+44)
            >>> f.xgcd(g, algorithm="euclidean")
            Traceback (most recent call last):
              ...
            ValueError
            >>> Kx = PolynomialRing_gr_poly(QQ, "x")
            >>> Kx(f).xgcd(g, algorithm="euclidean")
            (2 + x^2, (5/16) + (-3/16)*x, (11/16) + (5/8)*x + (3/16)*x^2)
        """
        Rx = self.parent()
        R = Rx._coefficient_ring
        # fixme:
        other = Rx(other)
        G = Rx()
        S = Rx()
        T = Rx()
        if algorithm is None:
            status = libgr.gr_poly_xgcd(G._ref, S._ref, T._ref, self._ref, other._ref, R._ref)
        elif algorithm == "euclidean":
            status = libgr.gr_poly_xgcd_euclidean(G._ref, S._ref, T._ref, self._ref, other._ref, R._ref)
        elif algorithm == "subresultant":
            status = libgr.gr_poly_xgcd_subresultant(G._ref, S._ref, T._ref, self._ref, other._ref, R._ref)
        else:
            raise ValueError
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (G, S, T)

    def roots(self, domain=None):
        """
        Computes the roots in the coefficient ring of this polynomial,
        returning a tuple (``roots``, ``multiplicities``).
        If the ring is not algebraically closed, the sum of multiplicities
        can be smaller than the degree of the polynomial.
        If ``domain`` is given, returns roots in that ring instead.

            >>> (ZZx([3,2]) * ZZx([15,1])**2 * ZZx([-10,1])).roots()
            ([-15, 10], [2, 1])
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
            ([[-1.414213562373095 +/- 6.23e-17], [1.414213562373095 +/- 6.23e-17], -1.500000000000000], [1, 1, 2])
            >>> f.roots(domain=CC)      # complex ball roots
            ([[-1.414213562373095 +/- 4.89e-17], [1.414213562373095 +/- 4.89e-17], 1.000000000000000*I, -1.000000000000000*I, -1.500000000000000], [1, 1, 1, 1, 2])
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
            FlintDomainError: f.exp_series(n) is not an element of {Ring of polynomials over Rational field (fmpq)} for {f = 1 + x}, {n = 2}
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
    """

        >>> f = QQser("exp(x)")
        >>> f
        1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 + (1/120)*x^5 + O(x^6)
        >>> f[5]
        1/120
        >>> f[6]
        Traceback (most recent call last):
            ...
        Undecidable: coefficient is not known

    """

    _struct_type = gr_series_struct

    def __init__(self, val=None, error=None, context=None, random=False):
        """
        If error is not None, add O(x^error) to the input val.

            >>> QQser("exp(-x)")
            1 - x + (1/2)*x^2 + (-1/6)*x^3 + (1/24)*x^4 + (-1/120)*x^5 + O(x^6)
            >>> QQser("exp(-x)", error=3)
            1 - x + (1/2)*x^2 + O(x^3)
            >>> QQser("exp(-x)", error=10)
            1 - x + (1/2)*x^2 + (-1/6)*x^3 + (1/24)*x^4 + (-1/120)*x^5 + O(x^6)
        """
        gr_elem.__init__(self, val, context, random)
        if error is not None:
            cur_err = self._error()
            if cur_err is None or cur_err > error:
                libgr._gr_series_set_error(self._ref, error, self.parent()._ref)

    def _error(self):
        """
        Returns the exponent n in the O(x^n) error term of self.
        If self is exact, returns None.

            >>> x = ZZser.gen()
            >>> f = 1+x
            >>> f._error()
            >>> g = (1+x)**10
            >>> g
            1 + 10*x + 45*x^2 + 120*x^3 + 210*x^4 + 252*x^5 + O(x^6)
            >>> g._error()
            6
        """
        if libgr._gr_series_is_exact(self._ref, self.parent()._ref) == T_TRUE:
            return None
        else:
            return libgr._gr_series_get_error(self._ref, self.parent()._ref)

    def _set_error(self, n):
        """
        Set the error of self to O(x^n) in-place, truncating
        any higher terms present. If n is None, makes self exact.

            >>> x = ZZser.gen()
            >>> g = (1+x)**10
            >>> g
            1 + 10*x + 45*x^2 + 120*x^3 + 210*x^4 + 252*x^5 + O(x^6)
            >>> g._set_error(3)
            >>> g
            1 + 10*x + 45*x^2 + O(x^3)
            >>> g._set_error(None)
            >>> g
            1 + 10*x + 45*x^2
            >>> g._set_error(1000)
            >>> g
            1 + 10*x + 45*x^2 + O(x^1000)
            >>> g._set_error(-10)
            >>> g
            0 + O(x^0)
        """
        if n is None:
            libgr._gr_series_make_exact(self._ref, self.parent()._ref)
        else:
            libgr._gr_series_set_error(self._ref, n, self.parent()._ref)

    def __getitem__(self, i):
        assert i >= 0
        error = self._error()
        if error is not None and i >= error:
            raise Undecidable("coefficient is not known")
        R = self.parent()._coefficient_ring
        c = R()
        # XXX: hack; want a gr_series method
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

        >>> DirichletGroup(10**16+61)
        Traceback (most recent call last):
          ...
        NotImplementedError: modulus with prime factor p > 10^16 is not currently supported

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

    Construction from strings:

        >>> Mat(QQ)("[[1, 1/2], [1/3, 1/4]]")
        [[1, 1/2],
        [1/3, 1/4]]
        >>> Mat(QQ)("[[0, 0, 0], [0, 0, 0]]")
        [[0, 0, 0],
        [0, 0, 0]]
        >>> Mat(QQ, 2, 3)("[[0, 0, 0], [0, 0, 0]]")
        [[0, 0, 0],
        [0, 0, 0]]
        >>> Mat(QQ, 2, 3)("[[0, 0], [0, 0]]")  # input does not match shape
        Traceback (most recent call last):
          ...
        ValueError


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

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
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
                elif isinstance(val, str):
                    status = libgr.gr_set_str(self._ref, ctypes.c_char_p(str(val).encode('ascii')), self._ctx)
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

    def permanent(self, algorithm=None):
        """
        Permanent of this matrix.

            >>> MatZZ(3, 3, ZZ.fac_vec(9)).permanent()
            1995840
            >>> Mat(RealField_arb(64))(8, 8, ZZ.fac_vec(64)).permanent()
            [+/- 4.07e+321]
            >>> Mat(RealField_arb(64))(8, 8, ZZ.fac_vec(64)).permanent(algorithm="cofactor")
            [5.5848931822182876e+307 +/- 6.08e+290]
            >>> Mat(RealField_arb(128))(8, 8, ZZ.fac_vec(64)).permanent()
            [5.5849e+307 +/- 2.12e+302]
        """
        element_ring = self.parent()._element_ring
        res = element_ring()
        if algorithm is None:
            status = libgr.gr_mat_permanent(res._ref, self._ref, element_ring._ref)
        elif algorithm == "cofactor":
            status = libgr.gr_mat_permanent_cofactor(res._ref, self._ref, element_ring._ref)
        elif algorithm == "ryser":
            status = libgr.gr_mat_permanent_ryser(res._ref, self._ref, element_ring._ref)
        elif algorithm == "glynn":
            status = libgr.gr_mat_permanent_glynn(res._ref, self._ref, element_ring._ref)
        elif algorithm == "glynn_threaded":
            status = libgr.gr_mat_permanent_glynn_threaded(res._ref, self._ref, element_ring._ref)
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
        r = (c_slong * 1)()
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
            ([-5, 2], [1, 1])
            >>> Mat(ZZ)([[1,2],[3,4]]).eigenvalues(domain=QQbar)
            ([Root a = 5.37228 of a^2-5*a-2, Root a = -0.372281 of a^2-5*a-2], [1, 1])
            >>> Mat(ZZ)([[1,2],[3,4]]).eigenvalues(domain=RR)
            ([[-0.3722813232690143 +/- 3.00e-17], [5.372281323269014 +/- 3.31e-16]], [1, 1])
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
            [[2, 0],
            [0, 3]]
            >>> D
            [2, 3]
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

    def qr(self):
        """
        QR decomposition.

            >>> A = Mat(RF)([[1,2,3],[0,0,1],[0,1,0],[2,2,1]])
            >>> Q, R = A.qr()
            >>> Q
            [[0.4472135954999579, 0.5962847939999439, 0.5716619504750294],
            [0, 0, 0.5144957554275265],
            [0, 0.7453559924999299, -0.5716619504750294],
            [0.8944271909999159, -0.2981423969999719, -0.2858309752375147]]
            >>> R
            [[2.236067977499790, 2.683281572999748, 2.236067977499790],
            [0, 1.341640786499874, 1.490711984999860],
            [0, 0, 1.943650631615100]]
            >>> (Q * R - A).norm_max() < 1e-15
            True

            >>> A = Mat(QQbar)([[1,2,3],[0,0,1],[0,1,0],[2,2,1]])
            >>> Q, R = A.qr()
            >>> Q[1,2]
            Root a = 0.514496 of 34*a^2-9
            >>> Q * R - A
            [[0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]]

            >>> A = Mat(RF, 100, 100)([[i+j+i//(1+j) for i in range(100)] for j in range(100)])
            >>> Q, R = A.qr()
            >>> (Q * R - A).norm_max() < 1e-13
            True

        """
        Cmat = self.parent()
        C = Cmat._element_ring
        m = self.nrows()
        n = self.ncols()
        Q = gr_mat(m, n, context=self.parent())
        R = gr_mat(n, n, context=self.parent())
        status = libgr.gr_mat_qr(Q._ref, R._ref, self._ref, C._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (Q, R)

    def lq(self):
        """
        LQ decomposition.

            >>> A = Mat(RF)([[1,2,3],[0,0,1],[0,1,0],[2,2,1]]).transpose()
            >>> L, Q = A.lq()
            >>> L
            [[2.236067977499790, 0, 0],
            [2.683281572999748, 1.341640786499874, 0],
            [2.236067977499790, 1.490711984999860, 1.943650631615100]]
            >>> Q
            [[0.4472135954999579, 0, 0, 0.8944271909999159],
            [0.5962847939999439, 0, 0.7453559924999299, -0.2981423969999719],
            [0.5716619504750294, 0.5144957554275265, -0.5716619504750294, -0.2858309752375147]]
            >>> (L * Q - A).norm_max() < 1e-15
            True
        """
        Cmat = self.parent()
        C = Cmat._element_ring
        m = self.nrows()
        n = self.ncols()
        L = gr_mat(m, m, context=self.parent())
        Q = gr_mat(m, n, context=self.parent())
        status = libgr.gr_mat_lq(L._ref, Q._ref, self._ref, C._ref)
        if status:
            if status & GR_UNABLE: raise NotImplementedError
            if status & GR_DOMAIN: raise ValueError
        return (L, Q)

    #def __getitem__(self, i):
    #    pass



libgr.gr_mat_entry_ptr.argtypes = (ctypes.c_void_p, c_slong, c_slong, ctypes.POINTER(gr_ctx_struct))
libgr.gr_mat_entry_ptr.restype = ctypes.POINTER(ctypes.c_char)

libgr.gr_vec_entry_ptr.restype = ctypes.POINTER(ctypes.c_char)


# todo singleton/cached domains (also for matrices, etc...)
class Vec(gr_ctx):
    """
    Parent class for vector domains.

        >>> Vec(ZZ)
        Vectors (any length) over Integer ring (fmpz)
        >>> Vec(ZZ, 5)
        Space of length 5 vectors over Integer ring (fmpz)
        >>> VecZZ([0, 5, 10])
        [0, 5, 10]
        >>> VecZZ(range(3, 20, 3))
        [3, 6, 9, 12, 15, 18]

    Construction of vectors from strings:

        >>> Vec(QQ)("[1/3, 1/5]")
        [1/3, 1/5]
        >>> Vec(ZZ, 3)("[1, 2, 3]")
        [1, 2, 3]
        >>> Vec(ZZ, 3)("[1, 2]")     # input does not match size
        Traceback (most recent call last):
          ...
        ValueError
        >>> Vec(RR)("[1, 1/3, 1/3 +/- exp(-10)]")
        [1, [0.3333333333333333 +/- 7.04e-17], [0.3333 +/- 7.88e-5]]
        >>> v = Vec(Mat(QQ, 2, 2))("[[[1,0],[0,1]], [[0,1],[-1,0]]]")
        >>> v[0]
        [[1, 0],
        [0, 1]]
        >>> v[1]
        [[0, 1],
        [-1, 0]]

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

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
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
                elif isinstance(val, str):
                    status = libgr.gr_set_str(self._ref, ctypes.c_char_p(str(val).encode('ascii')), self._ctx)
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


class fmpq_mpoly(gr_elem):
    _struct_type = fmpq_mpoly_struct

class PolynomialRing_fmpq_mpoly(gr_ctx):

    def __init__(self, nvars, vars=None):
        gr_ctx.__init__(self)
        nvars = gr_ctx._as_si(nvars)
        assert nvars >= 0
        libgr.gr_ctx_init_fmpq_mpoly(self._ref, nvars, 0)
        self._elem_type = fmpq_mpoly
        if vars is not None:
            assert len(vars) == nvars
            self._set_gen_names(vars)

    @property
    def _coefficient_ring(self):
        return QQ


class gr_mpoly(gr_elem):
    _struct_type = gr_mpoly_struct

class PolynomialRing_gr_mpoly(gr_ctx):
    def __init__(self, coefficient_ring, nvars, vars=None):
        assert isinstance(coefficient_ring, gr_ctx)
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

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
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



class Fraction_gr_fraction(gr_ctx):
    """
    Fractions with GCD reduction:

        >>> Q = Fraction_gr_fraction(ZZi)
        >>> Q(ZZi.i())
        (I) / (1)
        >>> Q(Fraction_gr_fraction(QQbar).i())
        (I) / (1)
        >>> I, = Q.gens(recursive=True)
        >>> s = sum(1/(1+j*I) for j in range(10))
        >>> s
        ((32798-28121*I)) / ((14651+1807*I))

    Fractions without GCD reduction

        >>> Q2 = Fraction_gr_fraction(ZZi, reduction=False)
        >>> I2, = Q2.gens(recursive=True)
        >>> s2 = sum(1/(1+j*I2) for j in range(10))
        >>> s2
        ((-468880-1874340*I)) / ((365300-549900*I))
        >>> Q2(s)
        ((32798-28121*I)) / ((14651+1807*I))
        >>> Q2(s) == s2
        True
        >>> s - s2
        (0) / (1)
        >>> s2 - s
        (0) / ((6345679600-7396487800*I))
        >>> s2 - s == 0
        True

    """

    def __init__(self, base_ring, reduction=True, strongly_canonical=False):
        assert isinstance(base_ring, gr_ctx)
        gr_ctx.__init__(self)

        flags = 0
        if not reduction:
            flags |= 1
        if strongly_canonical:
            flags |= 2

        libgr.gr_ctx_init_gr_fraction(self._ref, base_ring, flags)

        class _gr_fraction_struct(ctypes.Structure):
            _fields_ = [('data', ctypes.c_ubyte * libgr.gr_ctx_sizeof_elem(self._ref))]

        class gr_fraction(gr_elem):
            _struct_type = _gr_fraction_struct

        self._elem_type = gr_fraction

        base_ring._refcount += 1
        self._base_ring = base_ring
        self._elem_type = gr_fraction

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
            self._base_ring._decrement_refcount()



class Complex_gr_complex(gr_ctx):
    """
        >>> C = Complex_gr_complex(QQ)
        >>> C
        Complex algebra over Rational field (fmpq)
        >>> x = C("2+3*I"); x
        (2) + (3) * I
        >>> x / 5
        (2/5) + (3/5) * I
        >>> 1 / (1 / x)
        (2) + (3) * I
        >>> x.re(); x.im(); x.conj()
        (2) + (0) * I
        (3) + (0) * I
        (2) + (-3) * I

        >>> C = Complex_gr_complex(ZZx)
        >>> x, I = C.gens(recursive=True)
        >>> (2+x*I)**5
        (10*x^4-80*x^2+32) + (x^5-40*x^3+80*x) * I

        >>> A = RealAlgebraicField_qqbar()
        >>> C = Complex_gr_complex(A)
        >>> C("2+3*I")
        (2) + (3) * I
        >>> C(A(2).sqrt())
        (Root a = 1.41421 of a^2-2) + (0) * I
        >>> C(2 + Complex_gr_complex(QQ).i()/5)
        (2) + (1/5) * I
        >>> abs(C("2+3*I"))
        (Root a = 3.60555 of a^2-13) + (0) * I
        >>> C(QQbar(-1) ** (QQ(1) / 5))
        (Root a = 0.809017 of 4*a^2-2*a-1) + (Root a = 0.587785 of 16*a^4-20*a^2+5) * I
        >>> _**5
        (-1) + (0) * I

        >>> C = Complex_gr_complex(RealFloat_nfloat(64))
        >>> (C.pi() + C.i())**2
        (8.8696044010893586177) + (6.2831853071795864766) * I
    """

    def __init__(self, real_ctx):
        assert isinstance(real_ctx, gr_ctx)
        gr_ctx.__init__(self)

        libgr.gr_ctx_init_gr_complex(self._ref, real_ctx)

        class _gr_complex_struct(ctypes.Structure):
            _fields_ = [('data', ctypes.c_ubyte * libgr.gr_ctx_sizeof_elem(self._ref))]

        class gr_complex(gr_elem):
            _struct_type = _gr_complex_struct

        self._elem_type = gr_complex

        real_ctx._refcount += 1
        self._real_ctx = real_ctx
        self._elem_type = gr_complex

    def _decrement_refcount(self):
        # (the base context is released when this context is cleared,
        # after its last element, not when the Python object dies)
        self._refcount -= 1
        if not self._refcount:
            libgr.gr_ctx_clear(self._ref)
            self._real_ctx._decrement_refcount()


class padic_radix(gr_elem):

    _struct_type = padic_radix_struct

    def valuation(self):
        """
        The p-adic valuation v of this element (the exponent of the leading
        power of p). The valuation of zero is reported as ``None``.

            >>> Q7 = Qp_padic_radix(7)
            >>> Q7(14).valuation()
            1
            >>> Q7(98).valuation()
            2
            >>> (Q7(1) / Q7(7)).valuation()
            -1
            >>> Q7(5).valuation()
            0
            >>> Q7(0).valuation() is None
            True
        """
        if self._data.u.size == 0:
            return None
        return int(self._data.v)

    def precision(self):
        """
        The absolute precision N: the element is known modulo p^N. An exactly
        represented element returns ``None`` (infinite precision).

            >>> Q7r = Qp_padic_radix(7, rel_prec=8)
            >>> (Q7r(1) / Q7r(2)).precision()
            8
            >>> (Q7r(1) / Q7r(14)).precision()
            7
            >>> Q7r(5).precision() is None        # exact
            True
            >>> Q7r(1) / Q7r(7)                    # 7^-1 is exact
            (1) * 7^-1
            >>> (Q7r(1) / Q7r(7)).precision() is None
            True
        """
        if int(self._data.N) == PADIC_RADIX_PREC_INF:
            return None
        return int(self._data.N)

class decfloat(gr_elem):
    """
    Element of a :class:`RealFloat_decfloat` context.
    """

    _struct_type = decfloat_struct

    @staticmethod
    def _default_context():
        return RF_decfloat

    def __hash__(self):
        return hash(str(self))

    def round(self, prec, rnd="near"):
        """
        Round to *prec* significant digits with the given rounding mode.

            >>> x = RF_decfloat(1) / 7
            >>> x
            0.14285714285714285714
            >>> x.round(5), x.round(5, "down"), x.round(5, "up"), x.round(1)
            (0.14286, 0.14285, 0.14286, 0.1)
            >>> x.round(50)
            0.14285714285714285714
            >>> RF_decfloat("123456789").round(3), RF_decfloat("-0.5").round(1, "floor")
            (123000000, -0.5)
        """
        res = type(self)(context=self._ctx_python)
        status = libgr.decfloat_set_round(res._ref, self._ref, DECIMAL_PREC_EXACT if prec is None else int(prec), _decimal_rnd(rnd), self._ctx)
        if status:
            _handle_error(self.parent(), status, "round(x)", self)
        return res

    def digits(self):
        """
        Number of digits of the integer mantissa `M` in `x = \\pm M 10^v`
        (`M` not divisible by 10), i.e. the precision needed to represent
        *x* exactly (zero for zero and special values).

            >>> RF_decfloat("123.4500").digits(), RF_decfloat("1e100").digits(), RF_decfloat(0).digits()
            (5, 1, 0)
        """
        return libgr.decfloat_digits(self._ref, self._ctx)

    def limbs(self):
        """
        Number of limbs of the mantissa.

            >>> RealFloat_decfloat(30, limb_digits=3)("123.45").limbs(), RF_decfloat(0).limbs()
            (2, 0)
        """
        return libgr.decfloat_limbs(self._ref, self._ctx)

    def digit(self, k):
        """
        The digit of `|x|` at position `10^k`.

            >>> x = RF_decfloat("-123.45")
            >>> [x.digit(k) for k in range(3, -4, -1)]
            [0, 1, 2, 3, 4, 5, 0]
        """
        return libgr.decfloat_get_digit_si(self._ref, k, self._ctx)

    def set_digit(self, k, d):
        """
        A copy of *x* with the digit of `|x|` at position `10^k` replaced by
        *d* (exactly, without rounding to the context precision).

            >>> x = RF_decfloat("-123.45")
            >>> x.set_digit(1, 9), x.set_digit(-5, 7), x.set_digit(2, 0).set_digit(1, 0).set_digit(0, 0)
            (-193.45, -123.45007, -0.45)
        """
        res = type(self)(context=self._ctx_python)
        status = libgr.decfloat_set_digit_si(res._ref, self._ref, k, d, self._ctx)
        if status:
            _handle_error(self.parent(), status, "set_digit(x)", self)
        return res

    def exponent(self):
        """
        The scientific exponent `E` such that `10^E \\le |x| < 10^{E+1}`,
        or ``None`` for zero and special values.

            >>> RF_decfloat("123.45").exponent(), RF_decfloat("0.001").exponent(), RF_decfloat("1e100").exponent()
            (2, -3, 100)
            >>> RF_decfloat(0).exponent() is None
            True
        """
        if self._data.m.size == 0:
            return None
        E = c_slong()
        if libgr.decfloat_get_sci_exp_si(ctypes.byref(E), self._ref, self._ctx):
            return E.value
        return None

    def str_sci(self):
        """
        String representation in scientific notation.

            >>> RF_decfloat("123.45").str_sci(), RF_decfloat(1).str_sci(), RF_decfloat("0.001").str_sci()
            ('1.2345e2', '1', '1e-3')
        """
        ptr = libgr.decfloat_get_str_sci(self._ref, self._ctx)
        try:
            return ctypes.cast(ptr, ctypes.c_char_p).value.decode("ascii")
        finally:
            libflint.flint_free(ptr)

    def is_finite(self):
        """
            >>> R = RealFloat_decfloat(inf=True, nan=True)
            >>> R(1).is_finite(), R("inf").is_finite(), R("nan").is_finite()
            (True, False, False)
        """
        return bool(libgr._decfloat_is_finite(self._ref))


class decball(gr_elem):
    """
    Element of a :class:`RealField_decball` context.
    """

    _struct_type = decball_struct

    @staticmethod
    def _default_context():
        return RR_decball

    def mid(self):
        """
        The midpoint, as a decimal floating-point number.

            >>> (RR_decball(1) / 3).mid()
            0.33333333333333333333
        """
        ctx = self._ctx_python._float_context()
        res = decfloat(context=ctx)
        status = libgr.decball_get_mid(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "mid(x)", self)
        return res

    def rad(self):
        """
        The radius, as a decimal floating-point number.

            >>> (RR_decball(1) / 3).rad()
            3.334e-21
            >>> RR_decball(1).rad()
            0
        """
        ctx = self._ctx_python._float_context()
        ctx.digits = None
        res = decfloat(context=ctx)
        status = libgr._decmag_get_decfloat(res._ref, ctypes.byref(self._data.rad), self._ctx)
        if status:
            _handle_error(self.parent(), status, "rad(x)", self)
        return res

    def is_exact(self):
        """
            >>> RR_decball(1).is_exact(), (RR_decball(1) / 3).is_exact()
            (True, False)
        """
        return self._data.rad.m == 0

    def contains(self, other):
        """
        Whether this ball contains every point of *other* (converted to a
        ball in the same context).

            >>> R = RealField_decball(10)
            >>> x = R("1 +/- 0.001")
            >>> x.contains(1), x.contains(R("1.001")), x.contains(R("1.0011")), x.contains(QQ(1000)/999)
            (True, True, False, False)
            >>> x.contains(R("[1 +/- 0.001]")), x.contains(R("[1 +/- 0.0011]")), x.contains(R("[1.0005 +/- 0.0005]"))
            (True, False, True)
            >>> R("[1 +/- 1e-5]").contains(R("[1 +/- 1e-5]") ** 2)
            False
            >>> R("[1 +/- 1e-5]").contains(R("[1 +/- 1e-5]").sqrt())
            True
        """
        if not isinstance(other, decball) or other._ctx_python is not self._ctx_python:
            other = type(self)(other, context=self._ctx_python)
        return bool(libgr._decball_contains(self._ref, other._ref, self._ctx))

    def overlaps(self, other):
        """
        Whether this ball and *other* have a common point.

            >>> R = RealField_decball(10)
            >>> R("[1 +/- 0.1]").overlaps(R("[1.2 +/- 0.1]")), R("[1 +/- 0.1]").overlaps(R("[1.2 +/- 0.09]"))
            (True, False)
        """
        if not isinstance(other, decball) or other._ctx_python is not self._ctx_python:
            other = type(self)(other, context=self._ctx_python)
        return bool(libgr._decball_overlaps(self._ref, other._ref, self._ctx))

    def add_error(self, err):
        """
        Return a copy with the radius increased by an upper bound for
        `|err|`.

            >>> RR_decball(1).add_error(RR_decball("0.001")), RR_decball(1).add_error(RR_decball("-1e-30"))
            ([1 +/- 0.001], [1 +/- 1e-30])
            >>> RR_decball("[1 +/- 0.1]").add_error(RR_decball("0.1 +/- 0.01"))
            [1 +/- 0.21]
        """
        if not isinstance(err, decball) or err._ctx_python is not self._ctx_python:
            err = type(self)(err, context=self._ctx_python)
        res = type(self)(context=self._ctx_python)
        status = libgr.decball_set_interval_mid_rad(res._ref, self._ref, err._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "x.add_error(err)", self, err)
        return res

    def add_error_10exp(self, e):
        """
        Return a copy with the radius increased by `10^e`.

            >>> RR_decball(1).add_error_10exp(-5), RR_decball(1).add_error_10exp(2)
            ([1 +/- 1e-5], [1 +/- 100])
        """
        res = type(self)(self, context=self._ctx_python)
        status = libgr.decball_add_error_10exp_si(res._ref, int(e), self._ctx)
        if status:
            _handle_error(self.parent(), status, "x.add_error_10exp(e)", self)
        return res

    def trim(self):
        """
        Round the midpoint to the number of digits justified by the radius.

            >>> R = RealField_decball(20)
            >>> x = R(1) / 3 + R("+/- 0.001")
            >>> x
            [0.33333333333333333333 +/- 0.001001]
            >>> x.trim()
            [0.3333333 +/- 0.001002]
            >>> x.trim().contains(x)
            True
        """
        res = type(self)(context=self._ctx_python)
        status = libgr.decball_trim(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "trim(x)", self)
        return res

    def _unary_ball(self, func, name):
        res = type(self)(context=self._ctx_python)
        status = func(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, name, self)
        return res

    def lower(self):
        """
        Lower endpoint rounded toward `-\\infty` to the context precision,
        as an exact ball.

            >>> R = RealField_decball(5)
            >>> x = R("[1.23456 +/- 0.001]")
            >>> x.lower(), x.upper(), x.abs_lower(), x.abs_upper()
            (1.2335, 1.2357, 1.2335, 1.2357)
            >>> (-x).lower(), (-x).upper(), (-x).abs_lower(), (-x).abs_upper()
            (-1.2357, -1.2335, 1.2335, 1.2357)
            >>> R("[0.5 +/- 1]").abs_lower(), R("[0.5 +/- 1]").abs_upper()
            (0, 1.5)
        """
        return self._unary_ball(libgr.decball_lower, "lower(x)")

    def upper(self):
        """
        Upper endpoint rounded toward `+\\infty` to the context precision,
        as an exact ball. See :meth:`.lower`.
        """
        return self._unary_ball(libgr.decball_upper, "upper(x)")

    def abs_lower(self):
        """
        Lower bound for the absolute value, rounded toward zero to the
        context precision, as an exact ball. See :meth:`.lower`.
        """
        return self._unary_ball(libgr.decball_abs_lower, "abs_lower(x)")

    def abs_upper(self):
        """
        Upper bound for the absolute value, rounded toward `+\\infty` to the
        context precision, as an exact ball. See :meth:`.lower`.
        """
        return self._unary_ball(libgr.decball_abs_upper, "abs_upper(x)")

    def shell(self):
        """
        The ball with the same radius centered at zero.

            >>> RR_decball("[1.5 +/- 0.25]").shell(), RR_decball("[1.5 +/- 0.25]").mid_ball()
            ([0 +/- 0.25], 1.5)
        """
        return self._unary_ball(libgr.decball_shell, "shell(x)")

    def mid_ball(self):
        """
        The midpoint as an exact ball. See :meth:`.shell`.
        """
        return self._unary_ball(libgr.decball_mid, "mid_ball(x)")

    def rad_ball(self):
        """
        The radius as an exact ball.

            >>> RR_decball("[1.5 +/- 0.25]").rad_ball()
            0.25
        """
        return self._unary_ball(libgr.decball_rad, "rad_ball(x)")

    def union(self, other):
        """
        The smallest ball containing both balls.

            >>> RR_decball("[1 +/- 0.5]").union(RR_decball(3))
            [1.75 +/- 1.25]
            >>> RR_decball(3).union(RR_decball("[1 +/- 0.5]"))
            [1.75 +/- 1.25]
        """
        other = self.parent()(other)
        res = type(self)(context=self._ctx_python)
        status = libgr.decball_set_interval(res._ref, self._ref, other._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "union(x, y)", self, other)
        return res

    def add_rad(self, other):
        """
        Adds the absolute value of *other* (a ball) to the radius.

            >>> RR_decball("[1 +/- 0.5]").add_rad(RR_decball("[-1 +/- 0.5]"))
            [1 +/- 2]
        """
        other = self.parent()(other)
        res = type(self)(context=self._ctx_python)
        status = libgr.decball_add_rad(res._ref, self._ref, other._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "add_rad(x, y)", self, other)
        return res

    def round2(self, prec, rad_prec):
        """
        Rounds the midpoint to *prec* digits and then the radius to
        *rad_prec* digits.

            >>> x = RR_decball(1) / 3
            >>> x.round2(5, 1), x.round2(5, 3), x.round2(None, 1), x.round2(1, 3)
            ([0.33333 +/- 4e-6], [0.33333 +/- 3.34e-6], [0.33333333333333333333 +/- 4e-21], [0.3 +/- 0.0334])
        """
        if prec is None:
            prec = libgr.decimal_ctx_get_prec(self._ctx)
        res = type(self)(context=self._ctx_python)
        status = libgr.decball_set_round2(res._ref, self._ref, prec, rad_prec, self._ctx)
        if status:
            _handle_error(self.parent(), status, "round2(x, prec, rad_prec)", self)
        return res

    def rel_accuracy_digits(self):
        """
        Number of correct significant digits of the midpoint given the
        radius: the difference between the scientific exponents of the
        midpoint and the radius. Returns ``None`` for an exact ball.

            >>> R = RealField_decball(20)
            >>> (R(1) / 3).rel_accuracy_digits(), R("[1 +/- 0.01]").rel_accuracy_digits(), R("[123.4 +/- 0.01]").rel_accuracy_digits()
            (20, 2, 4)
            >>> R("[+/- 1]").rel_accuracy_digits(), R("[1 +/- 10]").rel_accuracy_digits(), R("[+/- 0.001]").rel_accuracy_digits()
            (0, -1, 3)
            >>> R(1).rel_accuracy_digits() is None
            True
        """
        acc = libgr.decball_rel_accuracy_digits(self._ref, self._ctx)
        if acc == DECIMAL_PREC_EXACT:
            return None
        return acc


class deccfloat(gr_elem):
    """
    Element of a :class:`ComplexFloat_deccfloat` context.
    """

    _struct_type = deccfloat_struct

    @staticmethod
    def _default_context():
        return CF_deccfloat

    def __hash__(self):
        return hash(str(self))

    def round(self, prec, rnd="near", rnd_im=None):
        """
        Round both parts to *prec* significant digits with the given
        rounding modes (*rnd_im* defaults to *rnd*).

            >>> x = CF_deccfloat("(1 + 2*I) / 7")
            >>> x
            (0.14285714285714285714 + 0.28571428571428571429*I)
            >>> x.round(5), x.round(5, "down"), x.round(5, "up", "floor"), x.round(1)
            ((0.14286 + 0.28571*I), (0.14285 + 0.28571*I), (0.14286 + 0.28571*I), (0.1 + 0.3*I))
        """
        if rnd_im is None:
            rnd_im = rnd
        res = type(self)(context=self._ctx_python)
        status = libgr.deccfloat_set_round(res._ref, self._ref, DECIMAL_PREC_EXACT if prec is None else int(prec), _decimal_rnd(rnd), _decimal_rnd(rnd_im), self._ctx)
        if status:
            _handle_error(self.parent(), status, "round(x)", self)
        return res

    def real(self):
        """
        The real part, as a decimal floating-point number.

            >>> CF_deccfloat("1.5 - 2*I").real()
            1.5
        """
        ctx = self._ctx_python._float_context()
        res = decfloat(context=ctx)
        status = libgr.deccfloat_get_re(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "real(x)", self)
        return res

    def imag(self):
        """
        The imaginary part, as a decimal floating-point number.

            >>> CF_deccfloat("1.5 - 2*I").imag()
            -2
        """
        ctx = self._ctx_python._float_context()
        res = decfloat(context=ctx)
        status = libgr.deccfloat_get_im(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "imag(x)", self)
        return res

    def is_real(self):
        """
            >>> CF_deccfloat("1.5 - 2*I").is_real(), CF_deccfloat("1.5").is_real()
            (False, True)
        """
        return bool(libgr._deccfloat_is_real(self._ref))

    def is_finite(self):
        """
            >>> C = ComplexFloat_deccfloat(inf=True, nan=True)
            >>> C(1).is_finite(), C("inf*I").is_finite(), C("nan").is_finite()
            (True, False, False)
        """
        return bool(libgr._deccfloat_is_finite(self._ref))


class deccball(gr_elem):
    """
    Element of a :class:`ComplexField_deccball` context.
    """

    _struct_type = deccball_struct

    @staticmethod
    def _default_context():
        return CC_deccball

    def mid(self):
        """
        The midpoint, as a complex decimal floating-point number.

            >>> (CC_deccball("1 + I") / 3).mid()
            (0.33333333333333333333 + 0.33333333333333333333*I)
        """
        ctx = self._ctx_python._complex_float_context()
        res = deccfloat(context=ctx)
        status = libgr.deccball_get_mid(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "mid(x)", self)
        return res

    def real(self):
        """
        The real part, as a decimal ball.

            >>> (CC_deccball("1 + I") / 3).real()
            [0.33333333333333333333 +/- 3.334e-21]
        """
        ctx = self._ctx_python._ball_context()
        res = decball(context=ctx)
        status = libgr.decball_set(res._ref, ctypes.byref(self._data.re), self._ctx)
        if status:
            _handle_error(self.parent(), status, "real(x)", self)
        return res

    def imag(self):
        """
        The imaginary part, as a decimal ball.

            >>> (CC_deccball("1 + 2*I") / 3).imag()
            [0.66666666666666666667 +/- 3.334e-21]
        """
        ctx = self._ctx_python._ball_context()
        res = decball(context=ctx)
        status = libgr.decball_set(res._ref, ctypes.byref(self._data.im), self._ctx)
        if status:
            _handle_error(self.parent(), status, "imag(x)", self)
        return res

    def is_exact(self):
        """
            >>> CC_deccball("1 + I").is_exact(), (CC_deccball("1 + I") / 3).is_exact()
            (True, False)
        """
        return self._data.re.rad.m == 0 and self._data.im.rad.m == 0

    def is_real(self):
        """
        Whether the imaginary part is exactly zero.

            >>> CC_deccball("1 + I").is_real(), CC_deccball("[1 +/- 0.1]").is_real(), CC_deccball("[+/- 0.1]*I").is_real()
            (False, True, False)
        """
        return bool(libgr._deccball_is_real(self._ref, self._ctx))

    def contains(self, other):
        """
        Whether this ball contains every point of *other* (converted to a
        ball in the same context).

            >>> C = ComplexField_deccball(10)
            >>> x = C("[1 +/- 0.001] + [2 +/- 0.001]*I")
            >>> x.contains(C("1 + 2*I")), x.contains(C("1.001 + 2*I")), x.contains(C("1.0011 + 2*I")), x.contains(1)
            (True, True, False, False)
        """
        if not isinstance(other, deccball) or other._ctx_python is not self._ctx_python:
            other = type(self)(other, context=self._ctx_python)
        return bool(libgr._deccball_contains(self._ref, other._ref, self._ctx))

    def overlaps(self, other):
        """
        Whether this ball and *other* have a common point.

            >>> C = ComplexField_deccball(10)
            >>> C("[1 +/- 0.1] + I").overlaps(C("[1.2 +/- 0.1] + I")), C("[1 +/- 0.1] + I").overlaps(C("[1.2 +/- 0.1] + 1.01*I"))
            (True, False)
        """
        if not isinstance(other, deccball) or other._ctx_python is not self._ctx_python:
            other = type(self)(other, context=self._ctx_python)
        return bool(libgr._deccball_overlaps(self._ref, other._ref, self._ctx))

    def add_error(self, err):
        """
        Return a copy with the radii increased by upper bounds for the
        absolute values of the real and imaginary parts of *err*.

            >>> CC_deccball(1).add_error(CC_deccball("0.001")), CC_deccball(1).add_error(CC_deccball("-1e-30 + 1e-20*I"))
            ([1 +/- 0.001], ([1 +/- 1e-30] + [0 +/- 1e-20]*I))
        """
        if not isinstance(err, deccball) or err._ctx_python is not self._ctx_python:
            err = type(self)(err, context=self._ctx_python)
        res = type(self)(context=self._ctx_python)
        status = libgr.deccball_set_interval_mid_rad(res._ref, self._ref, err._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "x.add_error(err)", self, err)
        return res

    def add_error_10exp(self, e):
        """
        Return a copy with both radii increased by `10^e`.

            >>> CC_deccball("1 + I").add_error_10exp(-5)
            ([1 +/- 1e-5] + [1 +/- 1e-5]*I)
        """
        res = type(self)(self, context=self._ctx_python)
        status = libgr.decball_add_error_10exp_si(ctypes.byref(res._data.re), int(e), self._ctx)
        status |= libgr.decball_add_error_10exp_si(ctypes.byref(res._data.im), int(e), self._ctx)
        if status:
            _handle_error(self.parent(), status, "x.add_error_10exp(e)", self)
        return res

    def trim(self):
        """
        Round the midpoints to the number of digits justified by the radii.

            >>> C = ComplexField_deccball(20)
            >>> x = C("1 + I") / 3 + C("+/- 0.001")
            >>> x
            ([0.33333333333333333333 +/- 0.001001] + [0.33333333333333333333 +/- 3.334e-21]*I)
            >>> x.trim()
            ([0.3333333 +/- 0.001002] + [0.33333333333333333333 +/- 3.334e-21]*I)
            >>> x.trim().contains(x)
            True
        """
        res = type(self)(context=self._ctx_python)
        status = libgr.deccball_trim(res._ref, self._ref, self._ctx)
        if status:
            _handle_error(self.parent(), status, "trim(x)", self)
        return res

    def rel_accuracy_digits(self):
        """
        Number of correct significant digits of the larger part given the
        larger radius. Returns ``None`` for an exact ball.

            >>> C = ComplexField_deccball(20)
            >>> (C("1 + I") / 3).rel_accuracy_digits(), C("[1 +/- 0.01] + [1000 +/- 0.1]*I").rel_accuracy_digits(), C("[+/- 1]*I").rel_accuracy_digits()
            (20, 4, 0)
            >>> C("1 + I").rel_accuracy_digits() is None
            True
        """
        acc = libgr.deccball_rel_accuracy_digits(self._ref, self._ctx)
        if acc == DECIMAL_PREC_EXACT:
            return None
        return acc


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
libflint.fexpr_set_si.argtypes = fexpr, c_slong
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
CC_tower = ComplexField_tower()
RR_tower = CC_tower.real_field()
QQbar_tower = ComplexAlgebraicField_tower()
AA_tower = QQbar_tower.real_field()

RF = RealFloat_arf()
CF = ComplexFloat_acf()

RF_decfloat = RealFloat_decfloat()
RR_decball = RealField_decball()
CF_deccfloat = ComplexFloat_deccfloat()
CC_deccball = ComplexField_deccball()

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

    M = Mat(CC)
    M2 = Mat(RR)
    A = M([[1,2],[3,4]])
    A2 = M2([[1,2],[3,4]])
    B = M([[1,2,3],[3,4,5]])
    B2 = M2([[1,2,3],[3,4,5]])
    C = M([[2,3],[4,5]])
    C2 = M2([[2,3],[4,5]])
    c = M([[2, 0], [0, 2]])
    cc = M([[0.5, 0], [0, 0.5]])
    for T in [QQ, ZZ, RR, CC]:
        assert A + T(2) == A + c
        assert T(2) + A == c + A
        assert A - T(2) == A - c
        assert T(2) - A == c - A
        assert A * T(2) == A * c
        assert T(2) * A == c * A
        assert A / T(2) == A * cc
    assert A + C2 == A + C
    assert A2 + C == A + C
    assert A - C2 == A - C
    assert A2 - C == A - C
    assert raises(lambda: A + B, ValueError)
    assert raises(lambda: A + B2, ValueError)
    assert raises(lambda: A2 + B, ValueError)
    assert raises(lambda: A - B, ValueError)
    assert raises(lambda: A - B2, ValueError)
    assert raises(lambda: A2 - B, ValueError)
    assert A * B == M([[7,10,13],[15,22,29]])
    assert A * B2 == M([[7,10,13],[15,22,29]])
    assert A2 * B == M([[7,10,13],[15,22,29]])
    assert raises(lambda: B * A, ValueError)
    assert raises(lambda: B2 * A, ValueError)
    assert raises(lambda: B * A2, ValueError)

    M = Mat(ZZ)
    MM = Mat(Mat(ZZ))
    A = M([[1,2],[3,4]])
    B = M([[2,3],[4,5]])
    C = M([[3,4],[5,6]])
    D = M([[4,5],[6,7]])
    assert MM([[A,B],[C,D]]) * MM([[B,C],[D,A]]) == MM([[A*B + B*D, A*C + B*A], [C*B + D**2, C**2 + D*A]])

    A = MatCC([[5,2],[3,4]])
    with optimistic_logic:
        assert (A ** 2) ** (QQ(1) / 2) == A
        assert (A ** 3) ** (RR(1) / 3) == A
        assert (A ** CC.i()) ** (1 / CC.i()) == A
        assert A ** (-5) == A ** ZZ(-5)
        assert A ** QQ(-5) == (A ** ZZ(5)).inv()
        assert A ** (-5) == ((A.log() * 5).exp()).inv()

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


def test_tower():
    """
    The ca test cases, run against the lazy tower field (gamma, erf and
    the extended values are not available there).
    """
    C = ComplexField_tower()
    sqrt = C.sqrt
    exp = C.exp
    log = C.log
    tan = C.tan
    sin = C.sin
    cos = C.cos
    acos = C.acos
    arg = C.arg
    atan = C.atan
    re = C.re
    im = C.im
    floor = C.floor
    ceil = C.ceil
    pi = C.pi()
    i = C.i()
    e = C.exp(1)

    def gd(x):
        return 2*atan(exp(x))-pi/2

    def cosh(x):
        y = exp(x)
        return (y + 1/y)/2

    def sinh(x):
        y = exp(x)
        return (y - 1/y)/2

    def tanh(x):
        return sinh(x)/cosh(x)

    assert floor(sqrt(2)) == 1
    assert ceil(sqrt(2)) == 2
    assert C.nint(C(5)/2) == 2
    assert C.trunc(-sqrt(2)) == -1

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

    assert sin(gd(1)) == tanh(1)
    assert tan(gd(1)) == sinh(1)
    assert sin(gd(sqrt(2))) == tanh(sqrt(2))
    assert tan(gd(1)/2) - tanh(C(1)/2) == 0

    assert C(qqbar(sqrt(2))) == sqrt(2)
    assert C(QQbar(2).sqrt()) == sqrt(2)

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

    # the tower field decides these (ca does not)
    assert acos(cos(1)) == 1
    assert acos(cos(sqrt(2) - 1)) == sqrt(2) - 1
    assert tan(sqrt(pi*2))*tan(sqrt(pi*3)) == \
        (cos(sqrt(pi*(5-2*sqrt(6)))) - cos(sqrt(pi*(5+2*sqrt(6)))))/(cos(sqrt(pi*(5-2*sqrt(6)))) + cos(sqrt(pi*(5+2*sqrt(6)))))
    assert sqrt(exp(2*sqrt(2)) + exp(-2*sqrt(2)) - 2) == (exp(2*sqrt(2))-1)/sqrt(exp(2*sqrt(2)))

    # Some examples from Stoutemyer
    z = pi
    assert i*((3-5*i)*z+1)/(((5+3*i)*z+i)*z) == 1/pi
    assert C(-1)**(C(1)/8) * sqrt(i+1) / C(2)**(C(3)/4) + i*exp(i*pi/2) == (i-1)/2
    assert sin(2*atan(z)) + (15*sqrt(3)+26)**(C(1)/3) == 2*z/(z**2+1) + sqrt(3) + 2

    # roots of unity and radicals
    assert exp(2*pi*i/3) * exp(2*pi*i/5) == exp(16*pi*i/15)
    assert C(6)**(C(1)/6) == sqrt(2) / C(2)**(C(1)/3) * C(3)**(C(1)/6)
    assert sqrt(-2) == i*sqrt(2)
    assert log(1+i) == log(2)/2 + pi*i/4
    assert log(-1) == pi*i
    assert log(3-4*i) == log(5) - 2*i*atan(C(1)/2)

    # polynomial roots: symmetric functions of the roots of an irreducible
    # polynomial reduce to the coefficients
    x = PolynomialRing(C).gen()
    f = x**5 - x - 1
    roots = f.roots()
    assert len(roots[0]) == 5
    assert sum(r for r in roots[0]) == 0
    prod = C(1)
    for r in roots[0]:
        prod = prod * r
    assert prod == 1
    assert sum(r**2 for r in roots[0]) == 0
    assert sum(1/r for r in roots[0]) == -1
    assert (x**2 - 2).roots()[0] == [sqrt(2), -sqrt(2)]

    # radicals of primes over fields unramified at them are new steps
    # (no search); over Q(zeta_12), sqrt(3) is found
    z16 = exp(2*pi*i/16)
    u = z16 + sqrt(2) + sqrt(3) + sqrt(5) + sqrt(7)
    v = sqrt(11) * u
    assert v**2 == 11 * u**2
    assert v != u * sqrt(13)
    z12 = exp(2*pi*i/12)
    assert sqrt(3) == z12 + 1/z12
    assert sqrt(3) * z12 == z12**2 + 1

    # square roots of written-out squares, without searching the tower
    assert sqrt((1-pi)**2) == pi-1
    assert sqrt(-(pi+1)**2) == i*(pi+1)
    assert sqrt((i*pi - C.exp(1))**2) == C.exp(1) - i*pi

    # roots over the field of the coefficients only, even when the towers
    # of the coefficients are shared with unrelated generators (the roots
    # of other polynomials, other radicals): a splitting tower of degree
    # 240 for x^5 - sqrt(2) x - 1, rather than a tower of degree 5760
    s = sum(sqrt(p) for p in [3,5,7,11,13]) + C(2)**(C(1)/3)
    c = (s + sqrt(2)) - s
    for f, n in [(x**4 - c*x - 1, 4), (x**3 - (1+c)*x**2 + 1, 3), (x**5 - c*x - 1, 5)]:
        rts = f.roots()[0]
        assert len(rts) == n
        assert all(f(r) == 0 for r in rts)
        assert sum(rts) == -f[n-1]
        # the roots live in a splitting tower over the small field of the
        # generators of c (of degree at most n! over it): Q(sqrt(2)), or
        # Q(zeta_16) when sqrt(2) is written through the roots of unity
        # of the shared tower, not over the tower of s (the last root, by
        # Vieta's formula, may involve the latter)
        assert sum(r.tower().degree() <= 8 * math.factorial(n) for r in rts) >= n - 1

def test_tower_trigonometric():
    C = ComplexField_tower()
    sqrt = C.sqrt
    pi = C.pi()
    xsin = C.sin
    xcos = C.cos
    xtan = C.tan

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
        assert sum(xcos(n*a) for n in range(1,N+1)) == xsin((N+C(1)/2)*a)/(2*xsin(a/2)) - C(1)/2

    assert xcos(a) == -sqrt(1 - xsin(a)**2)
    assert xsin(a/2) == sqrt((1-xcos(a))/2)

    assert xsin(3*a) == 4*xsin(a)*xsin(pi/3-a)*xsin(pi/3+a)
    assert xtan((a+b)/2) == (xsin(a) + xsin(b)) / (xcos(a) + xcos(b))
    assert xtan(a)*xtan(b) == ((xcos(a-b)-xcos(a+b))/(xcos(a-b)+xcos(a+b)))


def test_tower_trigonometric_real():
    # the same identities in the real field, where the trigonometric
    # functions of real arguments are tangents and arctangents

    C = ComplexField_tower().real_field()
    sqrt = C.sqrt
    pi = C.pi()
    xsin = C.sin
    xcos = C.cos
    xtan = C.tan

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
        assert sum(xcos(n*a) for n in range(1,N+1)) == xsin((N+C(1)/2)*a)/(2*xsin(a/2)) - C(1)/2

    assert xcos(a) == -sqrt(1 - xsin(a)**2)
    assert xsin(a/2) == sqrt((1-xcos(a))/2)

    assert xsin(3*a) == 4*xsin(a)*xsin(pi/3-a)*xsin(pi/3+a)
    assert xtan((a+b)/2) == (xsin(a) + xsin(b)) / (xcos(a) + xcos(b))
    assert xtan(a)*xtan(b) == ((xcos(a-b)-xcos(a+b))/(xcos(a-b)+xcos(a+b)))

def test_tower_calcium_issues():
    """
    Open issues from the Calcium issue tracker
    (https://github.com/flintlib/calcium/issues), on the lazy tower field.
    """
    C = ComplexField_tower()
    sqrt = C.sqrt
    exp = C.exp
    log = C.log
    cos = C.cos
    acos = C.acos
    pi = C.pi()
    i = C.i()

    # #38: cancellation of powers with a large denominator
    a = 417 / (962 * pi + 80808)
    assert a**50 - a**51 * a**-1 == 0

    # #24: exp((log(2i) - pi i/2)/2) = sqrt(2)
    assert exp((log(2*i) - pi*i/2)/2) == sqrt(2)

    # #33: cos(acos(sqrt2-sqrt3)/3) is an algebraic number of degree 12
    v = cos(acos(sqrt(2) - sqrt(3))/3)
    assert 1 - 90*v**2 + 321*v**4 - 592*v**6 + 864*v**8 - 768*v**10 + 256*v**12 == 0
    assert str(QQbar(v)).endswith("256*a^12-768*a^10+864*a^8-592*a^6+321*a^4-90*a^2+1")

    # #23: identities that ca cannot decide
    a = exp(2*sqrt(2)); b = exp(-2*sqrt(2))
    assert sqrt(a + b - 2) - (a-1)/sqrt(a) == 0
    assert (pi-1) / (sqrt(pi) - 1) == sqrt(pi) + 1
    M = Mat(C)([[1,-1],[-1,-1]])
    assert M.exp().log() == M

    # #25: normalization of sqrt(-3) and roots of unity
    z = (-1 + sqrt(-3))/2
    assert z**3 == 1 and z != 1 and sqrt(-3) == i*sqrt(3)

def test_tower_sage_examples():
    """
    Examples from the documentation of Sage's QQbar (sage.rings.qqbar)
    and related tickets, on the lazy tower field and on qqbar.
    """
    for R in [ComplexField_tower(), ComplexAlgebraicField_qqbar()]:
        sqrt = R.sqrt
        i = R.i()
        Rx = PolynomialRing(R)
        x = Rx.gen()

        # golden ratio
        s5 = sqrt(5); phi = (1 + s5)/2; tau = (1 - s5)/2
        assert phi**2 == phi + 1 and tau**2 == tau + 1 and phi + tau == 1

        # nested radicals
        assert (sqrt(5 + 2*sqrt(6)) - sqrt(3))**2 == 2
        assert sqrt(R(2)/3) * sqrt(R(3)/5) == sqrt(R(2)/5)

        # (-8)^(1/3): principal branch, absolute value, norm
        r = R(-8)**(R(1)/3)
        assert r**3 == -8 and abs(r) == 2 and r * r.conj() == 4
        assert r != -2

        # a cube root of unity
        z = -R(1)/2 + i*sqrt(3)/2
        assert z**3 == 1 and z**2 != 1

        # the roots of x^5 - x - 1
        f = x**5 - x - 1
        rts = f.roots()[0]
        assert len(rts) == 5
        assert all(r**5 - r - 1 == 0 for r in rts)
        assert sum(rts, R(0)) == 0

        # (sqrt2 + sqrt3)^5
        assert (sqrt(2) + sqrt(3))**5 == 109*sqrt(2) + 89*sqrt(3)

        # |3/5 + 4i/5| = 1
        assert abs(R(3)/5 + 4*i/5) == 1

        # a symbolic identity: (2/(3 sqrt3) + 10/27)^(1/3) - 2/(9 t) + 1/3 = 1
        s3 = sqrt(3)
        t = (2/(3*s3) + R(10)/27)**(R(1)/3)
        assert t - 2/(9*t) + R(1)/3 == 1

        # discriminant of a cubic: formula versus roots
        def disc1(b, c, d): return b**2*c**2 - 4*b**3*d - 4*c**3 + 18*b*c*d - 27*d**2
        def disc2(s1, s2, s3): return ((s1-s2)*(s1-s3)*(s2-s3))**2
        polys = [x*(x-2)*(x-4), x*(x-2)*(x-4) + 1]
        if not isinstance(R, ComplexAlgebraicField_qqbar):   # takes ~40 seconds with qqbar
            polys.append((x - sqrt(2))*(x - R(2)**(R(1)/3))*(x - sqrt(3)))
        for p in polys:
            rts = p.roots()[0]
            assert disc1(p[2], p[1], p[0]) == disc2(rts[0], rts[1], rts[2])

        # the regular 34-gon: rotating (1,0) 34 times by the angle 2pi/34,
        # expressed in radicals, gives (1,0) back
        rt17 = sqrt(17); rt2 = sqrt(2)
        eps = sqrt(17 + rt17); epss = sqrt(17 - rt17)
        delta = rt17 - 1
        alpha = sqrt(34 + 6*rt17 + rt2*delta*epss - 8*rt2*eps)
        cx = rt2*sqrt(15 + rt17 + rt2*(alpha + epss))/8
        cy = rt2*sqrt(epss**2 - rt2*(alpha + epss))/8
        px, py = R(1), R(0)
        for n in range(34):
            px, py = cx*px - cy*py, cx*py + cy*px
        assert px == 1 and py == 0

        # the same cos(2pi/34) as a root of a polynomial
        p = 256*x**8 - 128*x**7 - 448*x**6 + 192*x**5 + 240*x**4 - 80*x**3 - 40*x**2 + 8*x + 1
        cx2 = [r for r in p.roots()[0] if abs(complex(r) - 0.98297) < 1e-3][0]
        assert cx == cx2 and cy == sqrt(1 - cx2**2)

    # an identity from the ARPREC documentation: alpha^630 - 1 as a product
    # of cyclotomic-like factors, alpha the largest real root of a
    # degree-10 polynomial (Lehmer's polynomial)
    C = ComplexField_tower()
    x = PolynomialRing(C).gen()
    p = x**10 + x**9 - x**7 - x**6 - x**5 - x**4 - x**3 + x + 1
    a = [r for r in p.roots()[0] if abs(complex(r) - 1.17628) < 1e-3][0]
    lhs = a**630 - 1
    num = (a**315 - 1) * (a**210 - 1) * (a**126 - 1)**2 * (a**90 - 1) * (a**3 - 1)**3 * (a**2 - 1)**5 * (a - 1)**3
    den = (a**35 - 1) * (a**15 - 1)**2 * (a**14 - 1)**2 * (a**5 - 1)**6 * a**68
    assert lhs == num / den

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

def test_series():

    x = ZZser.gen()
    assert (x**2 / x) == x
    assert ((x**3 - x**4) / x**2) == x - x**2
    assert (3 * x) / (-3) == -x
    assert (3 * x**2) / (-3 * x) == -x
    assert (3 * x**2) / (-3 * x**2) == -1
    assert raises((lambda: (3 + 3*x**6) / 3), FlintUnableError)
    assert str((3 + 3*x**6) / (-1)) == "-3 + O(x^6)"

    assert (10 * x) / (5 * x) == 2
    assert raises(lambda: (10 * x) / (3 * x), FlintDomainError)

    assert (4 * x**0).sqrt() == 2
    assert raises(lambda: (4 + x**5).sqrt(), FlintUnableError)

    x = QQser.gen()
    assert str((3 + 3*x**6) / 3) == "1 + O(x^6)"
    assert str((4 + x**5).sqrt()) == "2 + (1/4)*x^5 + O(x^6)"
    assert str((4 + x**6).sqrt()) == "2 + O(x^6)"
    assert str((4 + 3*x).sqrt() * (4 + 3*x).rsqrt()) == "1 + O(x^6)"

    x = PowerSeriesRing(IntegersMod_nmod(17)).gen()
    assert raises(lambda: x.sqrt(), FlintUnableError)
    assert str((2 + x).sqrt()) == "6 + 10*x + 3*x^2 + 12*x^3 + 9*x^4 + 13*x^5 + O(x^6)"
    assert str((1 + 3*x).sqrt() * (1 + 3*x).rsqrt()) == "1 + O(x^6)"

    x = PowerSeriesRing(IntegersMod_nmod(16)).gen()
    assert raises(lambda: (1 + x).sqrt(), FlintUnableError)
    assert raises(lambda: (1 + x).rsqrt(), FlintUnableError)

    assert CCser(1+ZZser.gen()) == 1 + RRser.gen()

    # Test that trigonometric and hyperbolic functions use numerically
    # stable formulas for large arguments
    S = PowerSeriesRing(RR,  4)
    x = S.gen()
    c = RR(10)
    v = c + x + x**2
    d = 14
    assert S.tanh(v)[3].nstr(d) == '-1.0992819274357e-8'
    assert S.tanh(-v)[3].nstr(d) == '1.0992819274357e-8'
    assert S.coth(v)[3].nstr(d) == '1.0992819364988e-8'
    assert S.coth(-v)[3].nstr(d) == '-1.0992819364988e-8'
    S = PowerSeriesRing(CC,  3)
    x = S.gen()
    c = CC(0.25+10j)
    v = c + x + x**2
    d = 14
    assert S.tan(v)[2].nstr(d) == '(3.2826512022173e-9 + 1.1188008582726e-8*I)'
    assert S.tan(-v)[2].nstr(d) == '(-3.2826512022173e-9 - 1.1188008582726e-8*I)'
    assert S.cot(v)[2].nstr(d) == '(3.2826511245479e-9 + 1.1188008713377e-8*I)'
    assert S.cot(-v)[2].nstr(d) == '(-3.2826511245479e-9 - 1.1188008713377e-8*I)'
    assert S.tan_pi(v)[2].nstr(d) == '(-2.0362573263061e-26 + 6.4816083777740e-27*I)'
    assert S.tan_pi(-v)[2].nstr(d) == '(2.0362573263061e-26 - 6.4816083777740e-27*I)'
    assert S.cot_pi(v)[2].nstr(d) == '(-2.0362573263061e-26 + 6.4816083777740e-27*I)'
    assert S.cot_pi(-v)[2].nstr(d) == '(2.0362573263061e-26 - 6.4816083777740e-27*I)'
    v *= CC(1j)
    assert S.tanh(v)[2].nstr(d) == '(-1.1188008582726e-8 + 3.2826512022173e-9*I)'
    assert S.tanh(-v)[2].nstr(d) == '(1.1188008582726e-8 - 3.2826512022173e-9*I)'
    assert S.coth(v)[2].nstr(d) == '(1.1188008713377e-8 - 3.2826511245479e-9*I)'
    assert S.coth(-v)[2].nstr(d) == '(-1.1188008713377e-8 + 3.2826511245479e-9*I)'


def test_float():
    assert RF(5).mul_2exp(-1) == RF(2.5)
    assert CF(2+3j).mul_2exp(-1) == CF(1+1.5j)


def test_decimal_extreme():
    # found by out-of-tree fuzzing: no crashes, hangs or blowups
    import time
    def unable(f):
        try:
            f()
        except (FlintUnableError, FlintDomainError):
            return True
        return False
    t0 = time.time()
    R = RealField_decball(10)
    assert str(R("+/- inf").ceil()) == "[0 +/- inf]"
    assert str(ComplexField_deccball(10)(R("+/- inf")).floor()) == "[0 +/- inf]"
    assert str(RealFloat_decfloat(10)(3).mul_2exp(10**15)) == "4.702567019e301029995663981"
    assert str(R(3).mul_2exp(-10**15)) == "[1.913848322e-301029995663981 +/- 4.879e-301029995663991]"
    assert unable(lambda: RealFloat_decfloat(None)(3).mul_2exp(10**15))
    assert unable(lambda: RealFloat_decfloat(None)(3) ** (10**15))
    assert unable(lambda: ComplexFloat_deccfloat(None)(3) ** (10**15))
    C = ComplexFloat_deccfloat(20)
    assert unable(lambda: C("(0.5 + 1e20*I)").zeta())
    assert unable(lambda: C.riemann_xi(C("(0.5 + 1e7*I)")))
    assert unable(lambda: C("(1e100 + 1*I)").zeta())
    assert unable(lambda: C("1e10000000").gamma())
    RR.polylog(RR("1e25"), RR("1e-100"))
    assert time.time() - t0 < 10

def test_deccomplex():
    C = ComplexFloat_deccfloat(10)
    X = ComplexFloat_deccfloat(None)
    B = ComplexField_deccball(10)
    I = C.i()

    # parsing and printing
    for s, t in [("0", "0"), ("1", "1"), ("I", "1*I"), ("-I", "-1*I"), ("2*I", "2*I"), ("1+2*I", "(1 + 2*I)"),
                 ("1 - 2*I", "(1 - 2*I)"), ("(1 - 2*I)", "(1 - 2*I)"), ("-1.5e-7 + 2.5e3*I", "(-1.5e-7 + 2500*I)"),
                 ("(1+2*I)*(3-4*I)", "(11 + 2*I)"), ("(1+2*I)/(3-4*I)", "(-0.2 + 0.4*I)"), ("I^2", "-1"), ("I^3", "-1*I"),
                 ("(1+I)^2", "2*I"), ("(1+I)^-2", "-0.5*I"), ("sqrt(-4)", "2*I"), ("sqrt(2*I)", "(1 + 1*I)"),
                 ("(3+4*I)^(1/2)", "(2 + 1*I)"), ("(-3+4*I)^(1/2)", "(1 + 2*I)"), ("(-3-4*I)^(1/2)", "(1 - 2*I)"),
                 ("1/(1+I)", "(0.5 - 0.5*I)"), ("(1+I)/I", "(1 - 1*I)"), ("I*I*I*I", "1"), ("1e100*I", "1e100*I")]:
        assert str(X(s)) == t, (s, str(X(s)), t)
        assert str(X(str(X(s)))) == t
        assert str(C(s)) == t
        assert str(B(s)) == t
    for s in ["", "I I", "1 +", "(1 + I", "J"]:
        assert raises(lambda: C(s), (FlintUnableError, FlintDomainError, ValueError)), s
    assert str(X("1/(1+2*I)")) == "(0.2 - 0.4*I)" and str(X("1/(1+3*I)")) == "(0.1 - 0.3*I)"
    assert raises(lambda: X("1/(1+4*I)"), FlintUnableError)
    assert str(C("1/(1+4*I)")) == "(0.05882352941 - 0.2352941176*I)"
    Ci = ComplexFloat_deccfloat(10, inf=True, nan=True)
    assert str(Ci("inf*I")) == "inf*I" and str(Ci("-inf")) == "-inf" and str(Ci("(1 + nan*I)")) == "(nan + nan*I)"
    assert str(Ci("inf") * Ci.i()) == "inf*I" and str(Ci("inf*I") * Ci("inf*I")) == "-inf"
    assert not Ci("inf*I").is_finite() and Ci("1+I").is_finite()
    assert str(Ci(1) / 0) == "inf" and str(Ci.i() / 0) == "inf*I" and str(Ci(0) / 0) == "(nan + nan*I)"
    assert raises(lambda: C(1) / 0, FlintDomainError)

    # rounding modes, separately for the parts
    for rnd, rnd_im, expected in [("near", None, "(0.667 + 0.667*I)"), ("down", None, "(0.666 + 0.666*I)"),
                                  ("floor", "ceil", "(0.666 + 0.667*I)"), ("ceil", "floor", "(0.667 + 0.666*I)"),
                                  ("up", "down", "(0.667 + 0.666*I)")]:
        C3 = ComplexFloat_deccfloat(3, rnd=rnd, rnd_im=rnd_im)
        assert str(C3("(2 + 2*I) / 3")) == expected, (rnd, rnd_im)
        assert C3.rnd == rnd and C3.rnd_im == (rnd if rnd_im is None else rnd_im)
        assert str(C3("2 + 2*I") / 3) == expected
    C3 = ComplexFloat_deccfloat(3)
    C3.rnd_im = "floor"
    assert str(C3("(-2 - 2*I) / 3")) == "(-0.667 - 0.667*I)"
    assert str(C3("(2 + 2*I) / 3")) == "(0.667 + 0.666*I)"
    C3.rnd = "ceil"
    assert C3.rnd_im == "ceil"
    x = C("(1 + 2*I) / 7")
    assert str(x.round(3)) == "(0.143 + 0.286*I)"
    assert str(x.round(3, "down")) == "(0.142 + 0.285*I)"
    assert str(x.round(3, "down", "up")) == "(0.142 + 0.286*I)"
    assert str(x.round(None)) == str(x)

    # exact arithmetic: correct rounding against rationals
    Q = QQ
    for a, b, c, d, exact in [(1, 2, 3, 4, True), (Q(1)/3, Q(2)/7, Q(-5)/11, Q(3)/13, False), (10**20, 1, -1, Q(10)**-20, True), (0, 1, 0, 1, True), (7, 0, 0, 3, True), (Q(1)/4, Q(3)/8, Q(-5)/16, Q(1)/1000, True)]:
        if exact:
            Xa = X(a) + X(b) * X.i()
            Xc = X(c) + X(d) * X.i()
            prod = Xa * Xc
            assert str(prod) == str(X(Q(a)*Q(c) - Q(b)*Q(d)) + X(Q(a)*Q(d) + Q(b)*Q(c)) * X.i())
        Ca = C(Q(a)) + C(Q(b)) * I
        Cc = C(Q(c)) + C(Q(d)) * I
        # the inputs are rounded to 10 digits; the products and quotients
        # of the rounded inputs must be correctly rounded
        a, b, c, d = QQ(Ca.real()), QQ(Ca.imag()), QQ(Cc.real()), QQ(Cc.imag())
        pr = Q(a)*Q(c) - Q(b)*Q(d)
        pi = Q(a)*Q(d) + Q(b)*Q(c)
        assert str((Ca * Cc).real()) == str(RealFloat_decfloat(10)(pr))
        assert str((Ca * Cc).imag()) == str(RealFloat_decfloat(10)(pi))
        if c != 0 or d != 0:
            den = Q(c)**2 + Q(d)**2
            qr = (Q(a)*Q(c) + Q(b)*Q(d)) / den
            qi = (Q(b)*Q(c) - Q(a)*Q(d)) / den
            assert str((Ca / Cc).real()) == str(RealFloat_decfloat(10)(qr))
            assert str((Ca / Cc).imag()) == str(RealFloat_decfloat(10)(qi))
    assert str(C("3+4*I") * C("3+4*I")) == "(-7 + 24*I)" and str(C("3+4*I") ** 2) == "(-7 + 24*I)"
    assert str(C("3+4*I") ** 3) == "(-117 + 44*I)" and str(C("3+4*I") ** -3) == "(-0.007488 - 0.002816*I)"
    assert str(C("3+4*I") ** 0) == "1" and str(C("3+4*I") ** 1) == "(3 + 4*I)"
    assert str(C(0) ** 0) == "1" and str(C(0) ** C("1+I")) == "0"
    assert raises(lambda: C(0) ** C("-1+I"), FlintDomainError)
    assert str(C("2*I") ** 3) == "-8*I" and str(C("2*I") ** 4) == "16" and str(C("-2*I") ** -1) == "0.5*I"
    assert str(C("1+I") ** 100) == "-1125899907000000" and str(X("(1+I)^100")) == "-1125899906842624"
    assert str(C("1+I") ** 1000000) == "1.024e150515" or str(C("1+I") ** 1000000).endswith("e150514")
    assert str(C(2) ** C("0.5")) == "1.414213562" and str(C(-2) ** C("0.5")) == "1.414213562*I"
    assert str(C(-4) ** C("1.5")) == "-8*I" and str(C(-4) ** C("-0.5")) == "-0.5*I" and str(C(-8) ** C("2.5")) == "181.019336*I"
    assert str(C("3+4*I") ** C("-0.5")) == "(0.4 - 0.2*I)"

    # exact and rounded square roots, abs, sgn, arg
    assert str(C("3+4*I").sqrt()) == "(2 + 1*I)" and str(C("-3+4*I").sqrt()) == "(1 + 2*I)"
    assert str(C("3-4*I").sqrt()) == "(2 - 1*I)" and str(C("-3-4*I").sqrt()) == "(1 - 2*I)"
    assert str(C("0.0625*I").sqrt()) == "(0.1767766953 + 0.1767766953*I)"
    assert str(C("-0.0625*I").sqrt()) == "(0.1767766953 - 0.1767766953*I)"
    assert str(C("0.5*I").sqrt()) == "(0.5 + 0.5*I)" and str(C("-4").sqrt()) == "2*I" and str(C("-2").sqrt()) == "1.414213562*I"
    assert str(C("3+4*I").rsqrt()) == "(0.4 - 0.2*I)" and str(C(-4).rsqrt()) == "-0.5*I" and str(C(4).rsqrt()) == "0.5"
    assert str(C("3+4*I").abs()) == "5" and str(C("-4*I").abs()) == "4" and str(C("1+I").abs()) == "1.414213562"
    assert str(abs(C("1e100+1e100*I"))) == "1.414213562e100" and str(abs(C("1e-100+1e-100*I"))) == "1.414213562e-100"
    assert str(C("3+4*I").sgn()) == "(0.6 + 0.8*I)" and str(C("-2*I").sgn()) == "-1*I" and str(C(-3).sgn()) == "-1" and str(C(0).sgn()) == "0"
    assert str(C("1+I").sgn()) == "(0.7071067812 + 0.7071067812*I)"
    assert str(C("3+4*I").csgn()) == "1" and str(C("-3+4*I").csgn()) == "-1" and str(C("4*I").csgn()) == "1" and str(C("-4*I").csgn()) == "-1"
    assert str(C(1).arg()) == "0" and str(C(-1).arg()) == "3.141592654" and str(C.i().arg()) == "1.570796327" and str(C("-1-I").arg()) == "-2.35619449"
    assert raises(lambda: C(0).arg(), FlintDomainError)
    assert str(C("1+2*I").conj()) == "(1 - 2*I)" and str(C("1+2*I").re()) == "1" and str(C("1+2*I").im()) == "2"
    assert str(C("1+2*I").real()) == "1" and str(C("1+2*I").imag()) == "2" and str(C("1+2*I").real().parent()) == str(RealFloat_decfloat(10))
    assert str(-C("1+2*I")) == "(-1 - 2*I)" and str(C("1+2*I") - C("1+2*I")) == "0"
    assert C("1+2*I") == C("1+2*I") and C("1+2*I") != C("1-2*I") and C(1) == 1 and C.i() != 1
    assert raises(lambda: C("1+I") < C(1), ValueError)
    assert C(1) < C(2) and C(-1) <= C(-1)
    assert abs(C("1+I")) < abs(C("1.5")) and abs(C("3+4*I")) == abs(C(-5))
    assert str(C("3.7 + 0*I").floor()) == "3" and raises(lambda: C("3.7 + I").floor(), FlintDomainError)
    for K in [C, B]:
        assert [str(getattr(K(s), f)()) for s, f in [("3.7", "floor"), ("3.2", "ceil"), ("-3.7", "trunc"), ("2.5", "nint")]] == ["3", "4", "-3", "2"]
        assert str(K("-3+4*I").csgn()) == "-1" and not K.is_exact()
    assert str(Ci.neg_inf()) == "-inf" and raises(B.neg_inf, FlintDomainError) and raises(C.neg_inf, FlintDomainError)
    assert str(Ci.undefined()) == "(nan + nan*I)" and raises(B.undefined, FlintDomainError) and raises(B.unknown, FlintDomainError)
    # complex balls represent complex numbers; an infinite radius is the whole plane
    assert str(B("([+/- inf] + [0 +/- inf]*I)")) == "([0 +/- inf] + [0 +/- inf]*I)" and str(B("[1 +/- inf]*I")) == "[1 +/- inf]*I"
    assert str(B("([+/- inf] + 2*I)") * B.i()) == "(-2 + [0 +/- inf]*I)" and raises(lambda: B("(inf + 2*I)"), FlintUnableError)
    assert raises(lambda: B("[3.7 +/- 0.1] + [+/- 0.1]*I").floor(), FlintUnableError) and str(RR_decball("[3.7 +/- 0.1]").nint()) == "4"

    # conversions
    assert str(C(CC("1+2*I"))) == "(1 + 2*I)" and str(C(CF("0.5-I"))) == "(0.5 - 1*I)"
    assert str(C(RR_decball("1 +/- 0.1"))) == "1" and str(C(RF_decfloat(1) / 3)) == "0.3333333333"
    assert str(CC(C("1+2*I"))) == "(1.000000000000000 + 2.000000000000000*I)"
    assert str(CF(C("1+2*I"))) == "(1.000000000000000 + 2.000000000000000*I)"
    assert ZZ(C(3)) == 3 and QQ(C("0.5")) == QQ(1)/2 and str(RF_decfloat(C(2))) == "2" and str(RR_decball(C(2))) == "2"
    assert raises(lambda: ZZ(C("1+2*I")), FlintDomainError) and raises(lambda: RF_decfloat(C("1+2*I")), FlintDomainError)
    assert raises(lambda: RR(C("1+2*I")), (FlintDomainError, FlintUnableError))
    assert str(X(CF("0.1"))) == "0.1000000000000000055511151231257827021181583404541015625"
    assert str(ComplexFloat_deccfloat(5)(C("1.23456789 + 9.87654321*I"))) == "(1.2346 + 9.8765*I)"
    assert str(ComplexFloat_deccfloat(5, limb_digits=2)(C("1.23456789 + 9.87654321*I"))) == "(1.2346 + 9.8765*I)"
    assert str(C(QQbar(-1) ** (QQ(1)/2))) == "1*I"
    assert str(C(B("1 + I") / 3)) == "(0.3333333333 + 0.3333333333*I)"
    assert raises(lambda: X(B("1 + I") / 3), FlintUnableError)
    assert str(B(C("1 + I") / 3)) == "(0.3333333333 + 0.3333333333*I)"

    # functions: real and imaginary arguments, tiny arguments
    assert str(C(1).exp()) == "2.718281828" and str(C.i().exp()) == "(0.5403023059 + 0.8414709848*I)"
    assert str(C("3*I").sin()) == "10.01787493*I" and str(C("3*I").cos()) == "10.067662" and str(C("3*I").sinh()) == "0.1411200081*I"
    assert str(C("3*I").cosh()) == "-0.9899924966" and str(C("0.5*I").tan()) == "0.4621171573*I" and str(C("0.5*I").tanh()) == "0.5463024898*I"
    assert str(C("0.5*I").asin()) == "0.4812118251*I" and str(C("0.5*I").atan()) == "0.5493061443*I" and str(C("2*I").atan()) == "(1.570796327 + 0.5493061443*I)"
    assert str(C("0.5*I").asinh()) == "0.5235987756*I" and str(C("2*I").asinh()) == "(1.316957897 + 1.570796327*I)" and str(C("0.5*I").atanh()) == "0.463647609*I"
    assert str(C("2*I").erf()) == "18.56480241*I" and str(C("2*I").erfi()) == "0.995322265*I"
    assert str(C(-1).log()) == "3.141592654*I" and str(C(-2).log()) == "(0.6931471806 + 3.141592654*I)" and str(C.i().log()) == "1.570796327*I"
    assert str(C(2).acos()) == "1.316957897*I" and str(C(-2).acos()) == "(3.141592654 - 1.316957897*I)"
    assert str(C("0.5").acosh()) == "1.047197551*I" and str(C(-2).acosh()) == "(1.316957897 + 3.141592654*I)" and str(C(-1).acosh()) == "3.141592654*I"
    assert str(C(4).atanh()) == "(0.2554128119 - 1.570796327*I)" and str(C(4).asin()) == "(1.570796327 - 2.063437069*I)"
    assert str(C(-4).sqrt()) == "2*I" and str(C(-1).sqrt()) == "1*I"
    assert str(C.gamma(5)) == "24" and str(C("0.5").gamma()) == "1.772453851" and raises(lambda: C(0).gamma(), (FlintDomainError, FlintUnableError))
    assert str(C.zeta(-3)) == "0.008333333333" and str(C.log10(1000)) == "3" and str(C.sin_pi(QQ(1)/2)) == "1" and str(C.sin_pi(QQ(1)/6)) == "0.5000000001"
    assert str(C.exp_pi_i(QQ(1)/2)) == "1*I" and str(C.exp_pi_i(1)) == "-1" and str(C.exp_pi_i(QQ(1)/4)) == "(0.7071067812 + 0.7071067812*I)"
    assert str(C.exp_pi_i(C.i())) == "0.04321391826" and str(C.exp_pi_i(C("1+I"))) == "-0.04321391826"
    assert str(C.fac(20)) == "2432902008000000000" and str(C.rising(C("1+I"), 3)) == "10*I" and str(C.rising(C("1+I"), 0)) == "1"
    assert str(C.rising(C("0.5"), 5)) == "29.53125" and str(C.rising(C("1+I"), C("0.5"))) == "(1.003009581 + 0.4891951308*I)"
    assert str(C.lambertw(C("1+I"))) == "(0.6569660692 + 0.3254503394*I)" and str(C.lambertw(1)) == "0.5671432904"
    assert str(C.pi()) == "3.141592654" and str(C.euler()) == "0.5772156649"
    assert str(C.bessel_j(0, C.i())) == "1.266065878" and str(C.hurwitz_zeta(2, C("1+I"))) == "(0.4630000966 - 0.7942335428*I)"
    assert str(C.polylog(2, C("0.5"))) == "0.5822405265" and str(C.dilog(C.i())) == "(-0.2056167584 + 0.9159655942*I)"
    assert str(C.agm(C("1+I"), 2)) == "(1.527316275 + 0.5710047826*I)" and str(C.agm(1, C.i())) == "(0.5990701174 + 0.5990701174*I)"
    assert str(C("1+I").gamma()) == "(0.4980156681 - 0.1549498283*I)" and str(C("1+I").zeta()) == "(0.5821580598 - 0.9268485643*I)"
    assert str(C("1+I").erf()) == "(1.316151282 + 0.1904534692*I)" and str(C("1+I").lgamma()) == "(-0.6509231993 - 0.3016403205*I)"
    assert str(C("1+I").digamma()) == "(0.09465032062 + 1.076674047*I)" and str(C("1+I").exp()) == "(1.46869394 + 2.287355287*I)"
    assert str(C("1+I").sin()) == "(1.298457581 + 0.6349639148*I)" and str(C("1+I").cos()) == "(0.8337300251 - 0.9888977058*I)"
    assert str(C("1+I").tan()) == "(0.2717525853 + 1.083923327*I)" and str(C("1+I").atan()) == "(1.017221968 + 0.4023594781*I)"
    assert str(C("1+I").log()) == "(0.3465735903 + 0.7853981634*I)" and str(C("1+I").sqrt()) == "(1.098684113 + 0.4550898606*I)"
    assert str(C.sec(C("2+I"))) == "(-0.4131493443 + 0.6875274387*I)" and str(C.acot(C("2+I"))) == "(0.3926990817 - 0.1732867951*I)"
    assert str(C.csc(C("2+I"))) == "(0.6354937993 + 0.2215009309*I)" and str(C.sech(C("2+I"))) == "(0.1511762983 - 0.2269736754*I)"
    assert str(C.csch(C("2+I"))) == "(0.1413630216 - 0.2283750656*I)" and str(C.cot(C("2+I"))) == "(-0.1713836129 - 0.8213297975*I)"
    Cd = ComplexFloat_deccfloat(10, rnd="down")
    Cu = ComplexFloat_deccfloat(10, rnd="up")
    Cc = ComplexFloat_deccfloat(10, rnd="ceil")
    z = Cd("1e-1000000000 + 1e-1000000000*I")
    assert str(z.sin()) == "(1e-1000000000 + 9.999999999e-1000000001*I)"
    assert str(Cu(z).sin()) == "(1.000000001e-1000000000 + 1e-1000000000*I)"
    assert str(Cc(z).sin()) == "(1.000000001e-1000000000 + 1e-1000000000*I)"
    assert str(Cc(-z).sin()) == "(-1e-1000000000 - 9.999999999e-1000000001*I)"
    assert str(z.exp()) == "(1 + 1e-1000000000*I)" and str(Cu(z).exp()) == "(1.000000001 + 1.000000001e-1000000000*I)"
    assert str(z.cos()) == "(0.9999999999 - 9.999999999e-2000000001*I)" and str(Cu(z).cos()) == "(1 - 1e-2000000000*I)"
    assert str(z.gamma()) == "(4.999999999e999999999 - 4.999999999e999999999*I)"
    assert str(Cu(z).gamma()) == "(5e999999999 - 5e999999999*I)"
    assert str(z.log1p()) == "(9.999999999e-1000000001 + 9.999999999e-1000000001*I)"
    assert str(z.tan()) == "(9.999999999e-1000000001 + 1e-1000000000*I)" and str(z.atan()) == "(1e-1000000000 + 9.999999999e-1000000001*I)"
    assert str(Cd.zeta(z + 1)) == "(0.5772156649 - 9.999999999e999999999*I)"    # z + 1 rounds to 1 + 1e-1000000000*I
    assert str(Cd.acot(1 / z)) == "(1e-1000000000 + 9.999999999e-1000000001*I)"
    assert str(Cd.acsc(1 / z)) == "(9.999999999e-1000000001 + 1e-1000000000*I)"
    assert str(Cd.lambertw(z)) == "(9.999999999e-1000000001 + 9.999999999e-1000000001*I)"
    assert str(Cd.rgamma(z)) == "(1e-1000000000 + 1e-1000000000*I)"
    assert str(Cd.sinc(z)) == "(0.9999999999 - 3.333333333e-2000000001*I)"
    w = Cd("1e-1000000000 + 1e-2000000000*I")
    assert str(w.sin()) == "(9.999999999e-1000000001 + 9.999999999e-2000000001*I)"
    assert str(w.exp()) == "(1 + 1e-2000000000*I)" and str(w.cos()) == "(0.9999999999 - 9.999999999e-3000000001*I)"
    assert str(Cu(w).exp()) == "(1.000000001 + 1.000000001e-2000000000*I)"
    assert str(Cd("1e-2000000000 + 1e-1000000000*I").sin()) == "(1e-2000000000 + 1e-1000000000*I)"

    # balls
    x = B("1 + I") / 3
    assert str(x) == "([0.3333333333 +/- 3.334e-11] + [0.3333333333 +/- 3.334e-11]*I)"
    assert str(x * 3) == "([0.9999999999 +/- 1.001e-10] + [0.9999999999 +/- 1.001e-10]*I)"
    assert raises(lambda: x * 3 == B("1+I"), Undecidable) and x * 3 != 1 and not (x * 3 == 1)
    assert str(x.mid()) == "(0.3333333333 + 0.3333333333*I)" and str(x.mid().parent()) == str(ComplexFloat_deccfloat(10))
    assert str(x.real()) == "[0.3333333333 +/- 3.334e-11]" and str(x.imag()) == "[0.3333333333 +/- 3.334e-11]"
    assert str(x.real().parent()) == str(RealField_decball(10))
    assert not x.is_exact() and B("1+I").is_exact() and x.rel_accuracy_digits() == 10 and B("1+I").rel_accuracy_digits() is None
    assert x.contains(QQ(1)/3 + QQ(1)/3 * B.i()) and not x.contains(1) and x.overlaps(x) and not x.overlaps(B("0.3333 + 0.3333*I"))
    assert B("1+I").is_real() == False and B("1").is_real() and B("[1 +/- 0.1]").is_real() and not B("[+/- 0.1]*I").is_real()
    assert str(B("[1 +/- 0.1] + [2 +/- 0.2]*I")) == "([1 +/- 0.1] + [2 +/- 0.2]*I)"
    assert str(B("[1 + 2*I +/- 0.01]")) == "([1 +/- 0.01] + 2*I)" and str(B("[1 + 2*I +/- (0.01 + 0.02*I)]")) == "([1 +/- 0.01] + [2 +/- 0.02]*I)"
    assert str(B("[1 +/- 0.1]*I")) == "[1 +/- 0.1]*I" and str(B("[1 +/- 0.1] - 2*I")) == "([1 +/- 0.1] - 2*I)"
    assert str(B("([1 +/- 0.1] + [2 +/- 0.2]*I) * (3 - I)")) == "([5 +/- 0.5] + [5 +/- 0.7]*I)"
    assert str(B("([1 +/- 0.1] + [2 +/- 0.2]*I) / (3 - I)")) == "([0.1 +/- 0.05] + [0.7 +/- 0.07]*I)"
    assert str(B("(1+I)^2")) == "2*I" and str(B("(3+4*I)^(1/2)")) == "(2 + 1*I)" and str(B("I^3")) == "-1*I"
    assert str(B.i() ** 100) == "1" and str(B("(1+I)") ** 100) == "[-1125899907000000 +/- 157400]"
    assert str(B("1+I").sqrt()) == "([1.098684113 +/- 4.68e-10] + [0.4550898606 +/- 3.779e-11]*I)"
    assert str(B("1+I").abs()) == "[1.414213562 +/- 3.732e-10]" and str(B("3+4*I").abs()) == "5" and str(B("-4*I").abs()) == "4"
    assert str(B(-4).sqrt()) == "2*I" and str(B("-2").sqrt()) == "[1.414213562 +/- 3.731e-10]*I"
    assert str(B(-1).log()) == "[3.141592654 +/- 4.104e-10]*I" and str(B.i().exp()) == "([0.5403023059 +/- 3.188e-11] + [0.8414709848 +/- 7.898e-12]*I)"
    assert str(B.gamma(5)) == "24" and str(B("0.5").gamma()) == "[1.772453851 +/- 9.45e-11]"
    assert str(B("1+I").gamma()) == "([0.4980156681 +/- 1.837e-11] + [-0.1549498283 +/- 1.812e-12]*I)"
    assert str(B(2).acos()) == "[1.316957897 +/- 7.52e-11]*I"
    assert B("[1 +/- 0.001] + [2 +/- 0.001]*I").contains(B("1.001 + 2*I")) and not B("[1 +/- 0.001] + [2 +/- 0.001]*I").contains(B("1.0011 + 2*I"))
    assert str(B(1).add_error(B("0.001"))) == "[1 +/- 0.001]" and str(B("1+I").add_error(B("0.1 + 0.2*I"))) == "([1 +/- 0.1] + [1 +/- 0.2]*I)"
    assert str(B("1+I").add_error_10exp(-5)) == "([1 +/- 1e-5] + [1 +/- 1e-5]*I)"
    assert str((B(1) / 3 + B("+/- 0.001")).trim()) == "[0.3333333 +/- 0.001002]"
    assert str(B(CC("1+2*I"))) == "(1 + 2*I)" and str(B(CC(1) / 3)) == "[0.3333333333 +/- 3.335e-11]"
    assert str(CC(B("1+I") / 3)) == "([0.3333333333 +/- 3.34e-11] + [0.3333333333 +/- 3.34e-11]*I)"
    assert str(CF(B("1+I") / 3)) == "(0.3333333333000000 + 0.3333333333000000*I)"
    assert ZZ(B(3)) == 3 and QQ(B("0.5")) == QQ(1)/2 and raises(lambda: ZZ(B("3 +/- 0.1")), FlintUnableError) and raises(lambda: ZZ(B("3+I")), FlintDomainError)
    assert str(RR_decball(B(2))) == "2" and raises(lambda: RR_decball(B("2+I")), FlintDomainError) and raises(lambda: RR_decball(B("2 + [+/- 1]*I")), FlintUnableError)
    assert str(ComplexField_deccball(10, rad_prec=2)("1+I") / 3) == "([0.3333333333 +/- 3.4e-11] + [0.3333333333 +/- 3.4e-11]*I)"
    assert str(ComplexField_deccball(10, sloppy_radius=True)("1+I") / 3) == "([0.3333333333 +/- 5e-11] + [0.3333333333 +/- 5e-11]*I)"
    P = PolynomialRing(B)
    f = P("(x - I) * (x + [1 +/- 0.001])")
    assert str(f) == "([-1 +/- 0.001]*I) + ([1 +/- 0.001] - 1*I)*x + x^2" and str(P(str(f))) == str(f)
    assert str(f(B.i())) == "[0 +/- 0.002]*I"
    assert str(Mat(B)([[1, B.i()], [B.i(), B(1) / 3]]).det()) == "[1.333333333 +/- 3.334e-10]"
    P = PolynomialRing(C)
    assert str(P("(x - I) * (x + I)")) == "1 + x^2" and str(P("x^2 + 1")(C.i())) == "0"
    assert str(Mat(C)([[1, C.i()], [C.i(), 1]]).inv()) == "[[0.5, -0.5*I],\n[-0.5*I, 0.5]]"
    assert str(Mat(C)([[1, C.i()], [C.i(), 1]]).det()) == "2"

def test_decimal_functions():
    # Every elementary and special function on real and complex decimal
    # floats and balls, compared with 200-bit arb/acb: floats must be
    # correctly rounded (within half an ulp in each component, tested with
    # one ulp) and balls must contain the value. The complex contexts are
    # also given a non-real last argument.
    D = 15
    R2, C2 = RealField_arb(200), ComplexField_acb(200)
    x, y, z = "0.375", "1.25", "2.5"
    real = [("pi", ()), ("euler", ()), ("catalan", ()), ("khinchin", ()), ("glaisher", ()),
        ("exp", (x,)), ("expm1", (x,)), ("exp2", (x,)), ("exp10", (x,)), ("log", (z,)), ("log1p", (x,)),
        ("log2", (z,)), ("log10", (z,)), ("sin", (x,)), ("cos", (x,)), ("tan", (x,)), ("cot", (x,)),
        ("sec", (x,)), ("csc", (x,)), ("sinc", (x,)), ("sin_pi", (x,)), ("cos_pi", (x,)), ("tan_pi", (x,)),
        ("cot_pi", (x,)), ("sec_pi", (x,)), ("csc_pi", (x,)), ("sinc_pi", (x,)), ("sin_cos", (x,)),
        ("sin_cos_pi", (x,)), ("asin", (x,)), ("acos", (x,)), ("atan", (x,)), ("acot", (z,)), ("asec", (z,)),
        ("acsc", (z,)), ("asin_pi", (x,)), ("acos_pi", (x,)), ("atan_pi", (y,)), ("acot_pi", (y,)),
        ("asec_pi", (z,)), ("acsc_pi", (z,)), ("sinh", (x,)), ("cosh", (x,)), ("sinh_cosh", (x,)),
        ("tanh", (x,)), ("coth", (x,)), ("sech", (x,)), ("csch", (x,)), ("asinh", (x,)), ("acosh", (z,)),
        ("atanh", (x,)), ("acoth", (z,)), ("asech", (x,)), ("acsch", (x,)), ("lambertw", (x,)),
        ("gamma", (z,)), ("rgamma", (z,)), ("lgamma", (z,)), ("digamma", (z,)), ("barnes_g", (z,)),
        ("log_barnes_g", (z,)), ("zeta", (z,)), ("erf", (x,)), ("erfc", (x,)), ("erfi", (x,)),
        ("fresnel", (x,)), ("fresnel_s", (x,)), ("fresnel_c", (x,)), ("exp_integral_ei", (z,)),
        ("sin_integral", (z,)), ("cos_integral", (z,)), ("sinh_integral", (z,)), ("cosh_integral", (z,)),
        ("log_integral", (z,)), ("dilog", (x,)), ("agm", (y, z)), ("airy", (x,)), ("airy_ai", (x,)),
        ("airy_bi", (x,)), ("airy_ai_prime", (x,)), ("airy_bi_prime", (x,)), ("rising", (x, z)),
        ("bessel_j", (y, z)), ("bessel_y", (y, z)), ("bessel_i", (y, z)), ("bessel_k", (y, z)),
        ("polylog", (z, x)), ("hurwitz_zeta", (z, y)), ("exp_integral", (y, z)), ("gamma_upper", (y, z)),
        ("gamma_lower", (y, z)), ("beta_lower", (y, z, x)), ("chebyshev_t", (z, x)), ("chebyshev_u", (z, x)),
        ("hermite_h", (z, x)), ("gegenbauer_c", (z, y, x)), ("laguerre_l", (z, y, x)), ("jacobi_p", (z, x, y, x)),
        ("legendre_p", (z, y, x)), ("legendre_q", (z, y, x)), ("hypgeom_0f1", (y, x)), ("hypgeom_1f1", (y, z, x)),
        ("hypgeom_u", (y, z, x)), ("hypgeom_2f1", (y, x, z, x)), ("coulomb_f", (y, x, z)), ("coulomb_g", (y, x, z))]
    real = [(f, a, {}) for (f, a) in real] + [("bessel_i", (y, z), {"scaled": True}), ("bessel_k", (y, z), {"scaled": True}),
        ("fresnel", (x,), {"normalized": True}), ("log_integral", (z,), {"offset": True}),
        ("gamma_upper", (y, z), {"regularized": 1}), ("hypgeom_2f1", (y, x, z, x), {"regularized": True}),
        ("lambertw", ("-0.25",), {"k": -1}), ("erfinv", (x,), {}), ("erfcinv", (y,), {}), ("atan2", (x, "-1.25"), {})]
    tau = "0.25 + 1.5*I"
    cplx = [("modular_j", (tau,)), ("modular_lambda", (tau,)), ("modular_delta", (tau,)), ("dedekind_eta", (tau,)),
        ("exp_pi_i", (x,)), ("log_pi_i", (z,)), ("dirichlet_eta", (z,)), ("riemann_xi", (z,)), ("polygamma", (y, z)),
        ("lerch_phi", (x, z, y)), ("elliptic_k", (x,)), ("elliptic_e", (x,)), ("elliptic_pi", (x, z)),
        ("elliptic_f", (x, z)), ("elliptic_e_inc", (x, z))]
    cplx = [(f, a, {}) for (f, a) in cplx]
    # functions with only real arguments supported in the complex types
    # (and atan2, which is not provided for complex numbers)
    only_real = ["erfinv", "erfcinv"]

    def check(ctx, name, args, kw):
        ref_ctx = C2 if ctx.is_complex_vector_space() else R2
        try:
            ref = getattr(ref_ctx, name)(*[ref_ctx(a) for a in args], **kw)
        except (FlintUnableError, FlintDomainError):
            ref_ctx = R2
            ref = getattr(ref_ctx, name)(*[ref_ctx(a) for a in args], **kw)
        val = getattr(ctx, name)(*[ctx(a) for a in args], **kw)
        if not isinstance(ref, tuple):
            ref, val = (ref,), (val,)
        for r, v in zip(ref, val):
            if ctx.is_canonical():
                size = R2(abs(r.re()) + abs(r.im())) if ref_ctx is C2 else abs(r)
                assert R2(abs(ref_ctx(v) - r)) < size * R2(10) ** (1 - D) + R2(10) ** -100, (ctx, name, args, v, r)
            else:
                assert ref_ctx(v).overlaps(r), (ctx, name, args, v, r)

    for ctx in [RealFloat_decfloat(D), RealField_decball(D), ComplexFloat_deccfloat(D), ComplexField_deccball(D)]:
        is_complex = ctx.is_complex_vector_space()
        for (name, args, kw) in real + (cplx if is_complex else []):
            if is_complex and name == "atan2":
                continue
            check(ctx, name, args, kw)
            if is_complex and args and name not in only_real:
                check(ctx, name, args[:-1] + (args[-1] + " + 0.125*I",), kw)
    # exact and special-argument cases through the same interface
    F, B, C, CB = RealFloat_decfloat(D), RealField_decball(D), ComplexFloat_deccfloat(D), ComplexField_deccball(D)
    assert F.fac(20) == 2432902008176640000 and C.fac(ZZ(20)) == 2432902008176640000 and str(B.rising(B(3), 4)) == "360"
    assert B.fac(ZZ(20)) == 2432902008176640000 and CB.fac(ZZ(20)) == 2432902008176640000 and str(CB.rising(CB(3), 4)) == "360"
    assert str(B.gamma(QQ(1)/2)) == "[1.77245385090552 +/- 3.974e-15]" and str(CB.gamma(QQ(-1)/2)) == "[-3.54490770181103 +/- 2.056e-15]"
    assert str(F.gamma(QQ(1)/2)) == "1.77245385090552" and str(C.gamma(QQ(-1)/2)) == "-3.54490770181103"
    assert str(F.sec_pi(QQ(1)/3)) == "2" and str(C.asin_pi(QQ(1)/2)) == "0.166666666666667" and str(F.acos_pi(-1)) == "1"
    assert str(C.hermite_h(3, 4)) == "464" and str(CB.hermite_h(3, CB("I"))) == "-20*I"
    assert raises(lambda: C.erfinv(C("0.5 + 0.5*I")), FlintUnableError) and raises(lambda: C.atan2(C.i(), 1), FlintUnableError)
    assert str(C.lambertw(C("0.5 + I"), k=-1)) == "(-1.1303865182239 - 3.27266733025251*I)"

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

    assert str(sum(PolynomialRing_gr_mpoly(ZZi, 20).gens())) == "x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + x19 + x20"

    RA = PolynomialRing_gr_mpoly(ZZi, 2)
    RB = PolynomialRing_gr_mpoly(QQbar, 2)
    IA, xA, yA = RA.gens(recursive=True)
    xB, yB = RB.gens()
    IB = QQbar.i()
    cA = 2 - 3*IA
    cB = 2 - 3*IB
    assert xA == xB
    assert yA == yB
    assert xA != yB
    assert xB != yA
    assert RA(cB) == RB(cA)
    assert xA + cB == cA + xB
    assert RA(cB*yB + xB) == RB(cA*yA + xA)

    RA = PolynomialRing_gr_mpoly(ZZ, 2, ["x", "y"])
    RB = PolynomialRing_gr_mpoly(ZZmod(5), 2, ["x", "y"])
    xA, yA = RA.gens()
    xB, yB = RB.gens()
    assert RB((xA+yA)**10) == xB**10 + 2*(xB*yB)**5 + yB**10

    RA = PolynomialRing_gr_mpoly(RR, 2, ["x", "y"])
    RB = PolynomialRing_gr_mpoly(CC, 2, ["x", "y"])
    xA, yA = RA.gens()
    xB, yB = RB.gens()
    c = RR("0 +/- 1e-10")
    v = (xA + yA + c)**3 - (xB + yB)**3
    assert str(v) == "[+/- 3.01e-10]*x^2 + [+/- 6.01e-10]*x*y + [+/- 3.01e-20]*x + [+/- 3.01e-10]*y^2 + [+/- 3.01e-20]*y + [+/- 1.01e-30]"

    RA = PolynomialRing_gr_mpoly(ZZi, 2, ["x", "y"])
    RB = PolynomialRing_gr_mpoly(ZZi, 3, ["z", "y", "x"])
    xA, yA = RA.gens()
    zB, yB, xB = RB.gens()
    assert RA(xB) == xA
    assert RA(yB) == yA
    assert RB(xA) == xB
    assert RB(yA) == yB
    assert RA((-3+xB+2*yB)**3) == (-3+xA+2*yA)**3
    assert RB((-3+xA+2*yA)**3) == (-3+xB+2*yB)**3
    assert RA(xB * 0) == 0
    assert RB(xA * 0) == 0
    assert RA(xB ** 0) == 1
    assert RB(xA ** 0) == 1
    assert raises(lambda: RA(zB), NotImplementedError)   # todo: domain error

    RA2 = PolynomialRing_gr_mpoly(FiniteField_fq(2, 3), 2, ["x", "y"])
    assert raises(lambda: RA(RA2(0)), NotImplementedError)

    RA = PolynomialRing_gr_mpoly(ZZi, 2, ["x", "y"])
    RB = PolynomialRing_fmpz_mpoly(3, ["z", "y", "x"])
    xA, yA = RA.gens()
    zB, yB, xB = RB.gens()
    assert RA(xB) == xA
    assert RA(yB) == yA
    assert raises(lambda: RA(zB), NotImplementedError)   # todo: domain error
    assert RA((-3+xB+2*yB)**3) == (-3+xA+2*yA)**3

    # todo: match index when variables are not named ?
    RA = PolynomialRing_gr_mpoly(ZZi, 2)
    RB = PolynomialRing_gr_mpoly(ZZi, 3)
    assert raises(lambda: RA(RB.gens()[0]), NotImplementedError)
    assert raises(lambda: RB(RA.gens()[0]), NotImplementedError)
    RA = PolynomialRing_gr_mpoly(ZZi, 2)
    RB = PolynomialRing_gr_mpoly(ZZi, 2, ["x", "y"])
    assert raises(lambda: RA(RB.gens()[0]), NotImplementedError)
    assert raises(lambda: RB(RA.gens()[0]), NotImplementedError)

    x, y, z = PolynomialRing_gr_mpoly(QQbar, 3, ["x", "y", "z"]).gens()
    I = QQbar.i()
    assert ((x+I)*(x-I)*(y+I)*(y-I)) / ((x+I)*(y-I)) == (x-I)*(y+I)
    assert raises(lambda: ((x+I)*(x-I)*(y+I)*(y-I)) / ((x+I)*(z+I)), FlintDomainError)

    x, y, z = PolynomialRing_gr_mpoly(RR, 3, ["x", "y", "z"]).gens()
    assert raises(lambda: x / (x**0 - y**0), FlintDomainError)
    assert raises(lambda: (x**4 * y**3) / (RR("0 +/- 0.1")*(x**3 * y**2)), FlintUnableError)

    R = PolynomialRing_gr_mpoly(QQx, 1, "y");
    f = R("(x/3 + x^2/5)*y^2")
    assert f == R(str(f))


def test_fmpq_mpoly():
    QQxyz = PolynomialRing_fmpq_mpoly(3)
    x, y, z = QQxyz.gens()
    f = (-72 * (1+x)**2 * (y+z+1)) / 5
    assert f.numerator() == (-72 * (1+x)**2 * (y+z+1))
    assert f.denominator() == 5
    assert x.denominator() == 1
    assert (x * 0).numerator() == 0
    assert (x * 0).denominator() == 1
    c, fac, exp = f.factor()
    assert c == QQ(-72) / 5
    assert ((fac, exp) == ([1+x, y+z+1], [2, 1])) or ((fac, exp) == ([y+z+1, 1+x], [1, 2]))
    assert f.gcd((-100-100*x) / 17) == (1+1*x)
    assert str(QQxyz.gens()) == '[x1, x2, x3]'
    assert str(PolynomialRing_fmpq_mpoly(2, ["a", "b"]).gens()[1]) == "b"

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
    # "+/- inf" is parsed without evaluating inf as an element of the ring
    assert str(RR("1 +/- inf")) == "[+/- inf]" and str(RR("[+/- inf]")) == "[+/- inf]" and str(RR("+/- Inf + 2")) == "[+/- inf]"
    assert str(CC("1 + 2*i +/- inf")) == "([+/- inf] + 2.000000000000000*I)" and str(CC("(1+i) * [0 +/- inf]*i")) == "([+/- inf] + [+/- inf]*I)"
    assert str(RRx("x^2 + 1 +/- inf")) == "[+/- inf] + x^2" and str(RRx("[1 +/- inf]*x")) == "[+/- inf] + [+/- inf]*x"
    assert raises(lambda: RR("1 +/- infx"), FlintUnableError) and raises(lambda: QQ("1 +/- inf"), FlintUnableError)
    assert raises(lambda: ZZ("+/- inf"), FlintUnableError) and raises(lambda: QQ("2 +/- 1"), FlintUnableError)
    assert str(RealField_decball(10)("3 +/- inf")) == "[3 +/- inf]" and str(ComplexField_deccball(10)("(3 + 4*I) +/- inf")) == "([3 +/- inf] + 4*I)"

    assert raises(lambda: RR("foo"), FlintUnableError)
    assert raises(lambda: RR("expexp2"), FlintUnableError)
    assert raises(lambda: RR("sqrt(1"), FlintUnableError)
    assert raises(lambda: RR("sqrt1)"), FlintUnableError)
    assert raises(lambda: RR("foo(3)"), FlintUnableError)
    assert ZZ("sqrt 1") == 1
    assert CC("sqrt -1") == 1j
    assert ZZ("sqrt(5 * 5)") == 5
    assert raises(lambda: ZZ("sqrt(5 * 5 + 1)"), FlintUnableError)
    assert ZZ("fac(10) / fac(9)") == 10
    assert CC_ca("cos(1)^2 + sin(1)^2") == 1
    assert QQ("abs(floor(-11/2))") == 6
    assert QQ("ceil(-11/2)") == -5
    assert QQ("rsqrt(16)") == 0.25
    assert raises(lambda: QQ("rsqrt(0)"), FlintUnableError)
    assert QQbar("sinpi(1/4)/2 + 2*cospi(1/4) + tanpi(-1/3)^2/2") == QQbar("3/2 + 5*sqrt(2)/4")
    assert raises(lambda: QQbar("tanpi(1/2)"), FlintUnableError)
    assert RR("gamma(5)") == 24
    assert QQbar("re(2-7*i) * im(2-7*I)") == -14
    assert QQbar("conj(3+4*I)") == QQbar(3+4j).conj()

    with optimistic_logic:
        assert RR("exp(log(10))") == 10
        assert RR("tan(atan(1))") == 1
        assert RR("sin(asin(0.5))") == 0.5
        assert RR("cos(acos(0.5))") == 0.5
        assert CC("arg(sgn(1+I)) - pi/4") == 0
        assert RR("sqrt2 + sqrt3") == RR(2).sqrt() + RR(3).sqrt()
        assert RR("sqrt 2 + sqrt 3") == RR(2).sqrt() + RR(3).sqrt()
        assert RR("log log log (10^100)") == (RR(10)**100).log().log().log()
        assert RR("zeta(2)") == RR("pi^2/6")

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

    assert str(RR("(+/- 1e-3) + 0.1")) == "[0.10 +/- 1.01e-3]"
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

    assert str(CCx("(pi+I*x)^2")) == "[9.86960440108936 +/- 6.96e-15] + ([6.283185307179586 +/- 6.77e-16]*I)*x - x^2"

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

    # regression test
    assert CCser(1).cos_pi() == -1
    assert str(CCser(1).cos()) == "[0.540302305868140 +/- 4.59e-16]"
    assert CCser(1).sin_pi() == 0
    assert str(CCser(1).sin()) == "[0.841470984807897 +/- 6.08e-16]"


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

def test_is_vector_space():
    for R in [ZZ, ZZi, ZZx, PolynomialRing(ZZi), PolynomialRing_gr_mpoly(ZZi, 2), \
                Mat(ZZ), Mat(QQ), IntegersMod_nmod(5), Mat(IntegersMod_nmod(5), 3, 3), \
                PowerSeriesRing(ZZ), PowerSeriesModRing(ZZ, 3)]:
        assert not R.is_rational_vector_space()
        assert not R.is_real_vector_space()
        assert not R.is_complex_vector_space()
    for R in [QQ, QQx, AA, QQbar, Mat(QQ, 2, 3), NumberField(QQx.gen()**2 + 1), FractionField_fmpz_mpoly_q(2), PowerSeriesRing(QQ), PowerSeriesModRing(QQ, 3)]:
        assert R.is_rational_vector_space()
        assert not R.is_real_vector_space()
        assert not R.is_complex_vector_space()
    for R in [RR, RRx, RR_ca, Mat(RR, 2, 3), Vec(RR, 2)]:
        assert R.is_rational_vector_space()
        assert R.is_real_vector_space()
        assert not R.is_complex_vector_space()
    for R in [CC, CCx, CC_ca, Mat(CC, 2, 3), Vec(CC, 1)]:
        assert R.is_rational_vector_space()
        assert R.is_real_vector_space()
        assert R.is_complex_vector_space()
    for R in [PowerSeriesModRing(ZZ, 0), Mat(ZZ, 0, 0), Vec(ZZ, 0),
                PowerSeriesModRing(IntegersMod_nmod(5), 0), Mat(IntegersMod_nmod(5), 0, 0), Vec(IntegersMod_nmod(5), 0),
            IntegersMod_nmod(1)]:
        assert R.is_rational_vector_space()
        assert R.is_real_vector_space()
        assert R.is_complex_vector_space()

def test_padic():
    Q7 = Qp_padic_radix(7, rel_prec=10)
    Q2 = Qp_padic_radix(2, rel_prec=10)

    assert str(Q7(7) * Q7(7)) == "(1) * 7^2"
    assert Q7(0).sqrt() == 0
    assert Q7(1).sqrt() == 1
    assert Q7(4).sqrt() == 2
    assert str(Q7(2).sqrt()) == "266983762 + O(7^10)"
    assert str((1/(1/Q7(4))).sqrt()) == "2 + O(7^10)"
    assert raises(Q7(3).sqrt, FlintDomainError)
    assert raises(Q7(5).sqrt, FlintDomainError)
    assert raises(Q7(6).sqrt, FlintDomainError)

    assert Q7(0).exp() == 1
    assert str(Q7(7).exp()) == "182289612 + O(7^10)"
    assert str(Q7("7 + O(7^11)").exp()) == "182289612 + O(7^10)"
    assert str(Q7("7 + O(7^9)").exp()) == "20875184 + O(7^9)"
    assert str(Q7("0 + O(7^2)").exp()) == "1 + O(7^2)"
    assert str(Q7("0 + O(7^1)").exp()) == "1 + O(7^1)"
    assert raises(Q7("0 + O(7^0)").exp, FlintUnableError)
    assert raises(Q7("0 + O(7^-1)").exp, FlintUnableError)
    assert raises((Q7(7)**(-5) + Q7("O(7^-2)")).exp, FlintDomainError)

    assert Q7(1).log() == 0
    assert str(Q7("1 + O(7^5)").log()) == "0 + O(7^5)"
    assert str(Q7(1000 * 7).exp().log()) == "(1000) * 7^1 + O(7^10)"
    assert raises(Q7(2).log, FlintDomainError)
    assert raises(Q7(0).log, FlintDomainError)
    assert str(Q7("1 + O(7)").log()) == "0 + O(7^1)"
    assert raises(Q7("O(7)").log, FlintDomainError)
    assert raises(Q7("0 + O(7^0)").log, FlintUnableError)
    assert raises(Q7("0 + O(7^-1)").log, FlintUnableError)
    assert raises(Q7("7^-5 + O(7^-2)").log, FlintDomainError)

    assert Q2(0).exp() == 1
    assert str(Q2(4).exp()) == "333 + O(2^10)"
    assert str(Q2("4 + O(2^12)").exp()) == "333 + O(2^10)"
    assert str(Q2("4 + O(2^8)").exp()) == "77 + O(2^8)"
    assert str(Q2("0 + O(2^2)").exp()) == "1 + O(2^2)"
    assert raises(Q2("0 + O(2^1)").exp, FlintUnableError)
    assert raises(Q2("0 + O(2^0)").exp, FlintUnableError)
    assert raises(Q2("0 + O(2^1)").exp, FlintUnableError)

    assert Q2(1).log() == 0
    assert str(Q2("1 + O(2^5)").log()) == "0 + O(2^5)"
    assert str(Q2("1 + 2 * 3").log()) == "(59) * 2^3 + O(2^10)"
    assert str(Q2(123 * 4).exp().log()) == "(123) * 2^2 + O(2^10)"
    assert raises(Q2(2).log, FlintDomainError)
    assert raises(Q2(0).log, FlintDomainError)
    assert str(Q2("1 + O(2)").log()) == "0 + O(2^1)"
    assert raises(Q2("O(2)").log, FlintDomainError)
    assert raises(Q2("0 + O(2^0)").log, FlintUnableError)
    assert raises(Q2("0 + O(2^-1)").log, FlintUnableError)
    assert raises(Q2("2^-5 + O(2^-2)").log, FlintDomainError)

    Q5 = Qp_padic_radix(5, rel_prec=10)
    f = PowerSeriesRing(Q5)("5 * 1234 + x")
    ef = f.exp()
    stref = "(1779246 + O(5^10)) + (1779246 + O(5^10))*x + (889623 + O(5^10))*x^2 + (296541 + O(5^10))*x^3 + (7398354 + O(5^10))*x^4 + ((7398354) * 5^-1 + O(5^9))*x^5 + O(x^6)"
    assert str(ef) == stref
    assert str((PowerSeriesRing(Q5))(stref)) == stref
    assert str(f.exp() * (-f).exp()) == "(1 + O(5^10)) + (0 + O(5^10))*x + (0 + O(5^10))*x^2 + (0 + O(5^10))*x^3 + (0 + O(5^10))*x^4 + (0 + O(5^9))*x^5 + O(x^6)"


def test_decimal():
    R = RealFloat_decfloat(10)
    R40 = RealFloat_decfloat(40)
    X = RealFloat_decfloat(None)
    B = RealField_decball(10)

    # parsing and printing
    for s, t in [("0", "0"), ("-0", "0"), ("1", "1"), ("-1", "-1"), ("1.0", "1"), ("0.50", "0.5"),
                 ("000123.4500", "123.45"), ("1e3", "1000"), ("1.5E3", "1500"), ("-2.5e-3", "-0.0025"),
                 ("1e20", "100000000000000000000"), ("1e21", "1e21"), ("123e19", "1.23e21"),
                 ("0.000001", "0.000001"), ("0.0000001", "1e-7"), ("123.456e-9", "1.23456e-7"),
                 ("1e-100", "1e-100"), ("1e+100", "1e100"), ("+7", "7"), ("  3.25  ", "3.25"),
                 ("1/8", "0.125"), ("2^10", "1024"), ("(1+2)*3", "9"), ("10^-3", "0.001"),
                 ("1e1000000000000000000000", "1e1000000000000000000000"),
                 ("12345678901234567890123456789", "1.2345678901234567890123456789e28"),
                 ("1234567890123456789012", "1.234567890123456789012e21"),
                 ("123456789012345678901", "123456789012345678901")]:
        assert str(X(s)) == t, (s, str(X(s)), t)
        assert str(X(str(X(s)))) == t
    assert str(X("123.5").str_sci()) == "1.235e2"
    assert str(RealFloat_decfloat(None, scientific=True)("123.5")) == "1.235e2"
    assert str(RealFloat_decfloat(None, scientific=True)("0")) == "0"
    for s in ["", "abc", "1..2", "1e", "e5", "[1", "1 +/- 2", "1 +", "*1"]:
        assert raises(lambda: X(s), (FlintUnableError, FlintDomainError, ValueError)), s
    assert raises(lambda: R("inf"), (FlintUnableError, FlintDomainError, ValueError))
    assert raises(lambda: R("nan"), (FlintUnableError, FlintDomainError, ValueError))
    # exact values and domain errors detected from the digits, for any size
    X = RealFloat_decfloat(None); F = RealFloat_decfloat(20); Fd = RealFloat_decfloat(20, rnd="down"); Fu = RealFloat_decfloat(20, rnd="up")
    assert F.zeta(X("-1.5e1000000")) == 0 and str(F.zeta(-3)) == "0.0083333333333333333333" and raises(lambda: F.zeta(1), FlintDomainError)
    assert F.log2(X(2) ** -1000) == -1000 and F.log2(X(2) ** 1000) == 1000 and F.log10(X("1e-1000000")) == -1000000
    assert F.sin_pi(X("1234567890123.5")) == -1 and F.cos_pi(X("-98765432109876543211")) == -1 and F.tan_pi(X("12345678901.25")) == 1
    assert raises(lambda: F.cot_pi(X("-5e999")), FlintDomainError) and F.sinc_pi(X("7e1000000")) == 0 and raises(lambda: F.gamma(X("-7e1000000")), FlintDomainError)
    for f, t in [("asin_pi", "1.5"), ("acos_pi", "-7e1000000"), ("acosh", "0.5"), ("log", "-2"), ("atanh", "-3"), ("asec", "0.5"), ("erfinv", "2")]:
        assert raises(lambda: getattr(F, f)(X(t)), FlintDomainError), f
    # directed rounding with exactly representable leading terms
    t = X("3.25e-1000000")
    assert str(Fd.barnes_g(t)) == "3.25e-1000000" and str(Fu.barnes_g(t)) == "3.2500000000000000001e-1000000"
    assert str(Fd.zeta(t)) == "-0.5" and str(Fu.zeta(t)) == "-0.50000000000000000001" and str(Fu.sec_pi(t)) == "1.0000000000000000001"
    assert str(Fd.acos_pi(t)) == "0.49999999999999999999" and str(Fu.acot_pi(-t)) == "-0.5" and str(Fd.atan_pi(X("4e1000000"))) == "0.49999999999999999999"
    assert str(Fu.asec_pi(X("-4e1000000"))) == "0.50000000000000000001"
    # digits and limbs
    y = RealFloat_decfloat(30, limb_digits=3)("-123.45")
    assert y.digits() == 5 and y.limbs() == 2 and [y.digit(k) for k in range(3, -4, -1)] == [0, 1, 2, 3, 4, 5, 0]
    assert str(y.set_digit(3, 9)) == "-9123.45" and str(y.set_digit(-1, 0)) == "-123.05" and str(y.set_digit(-2, 0).set_digit(-1, 0)) == "-123"
    Ri = RealFloat_decfloat(10, inf=True, nan=True)
    assert str(Ri("inf")) == "inf" and str(Ri("-inf")) == "-inf" and str(Ri("nan")) == "nan"
    assert str(Ri("Infinity")) == "inf"
    assert not Ri("inf").is_finite() and Ri(1).is_finite()
    assert [str(Ri(t).rsqrt()) for t in ["inf", "0", "nan"]] == ["0", "inf", "nan"] and raises(lambda: Ri("-inf").rsqrt(), FlintDomainError)

    # rounding
    cases = [("1/3", {"near": "0.3333333333", "down": "0.3333333333", "up": "0.3333333334", "floor": "0.3333333333", "ceil": "0.3333333334"}),
             ("-1/3", {"near": "-0.3333333333", "down": "-0.3333333333", "up": "-0.3333333334", "floor": "-0.3333333334", "ceil": "-0.3333333333"}),
             ("2/3", {"near": "0.6666666667", "down": "0.6666666666", "up": "0.6666666667", "floor": "0.6666666666", "ceil": "0.6666666667"}),
             ("12345678905", {"near": "12345678900", "near_away": "12345678910", "near_zero": "12345678900", "down": "12345678900", "up": "12345678910"}),
             ("12345678915", {"near": "12345678920", "near_away": "12345678920", "near_zero": "12345678910"}),
             ("-12345678915", {"near": "-12345678920", "near_away": "-12345678920", "near_zero": "-12345678910", "floor": "-12345678920", "ceil": "-12345678910"}),
             ("12345678905e20", {"near": "1.23456789e30", "near_away": "1.234567891e30", "near_zero": "1.23456789e30", "down": "1.23456789e30", "up": "1.234567891e30"}),
             ("-12345678915e20", {"near": "-1.234567892e30", "near_away": "-1.234567892e30", "near_zero": "-1.234567891e30", "floor": "-1.234567892e30", "ceil": "-1.234567891e30"}),
             ("9999999999.5", {"near": "10000000000", "down": "9999999999", "up": "10000000000", "near_zero": "9999999999"}),
             ("9999999999.5e20", {"near": "1e30", "down": "9.999999999e29", "up": "1e30", "near_zero": "9.999999999e29"}),
             ("1.00000000001", {"near": "1", "down": "1", "up": "1.000000001", "floor": "1", "ceil": "1.000000001"}),
             ("-1.00000000001", {"near": "-1", "down": "-1", "up": "-1.000000001", "floor": "-1.000000001", "ceil": "-1"}),
             ("123", {"near": "123", "up": "123", "down": "123"})]
    for s, d in cases:
        for rnd, t in d.items():
            Rr = RealFloat_decfloat(10, rnd=rnd)
            assert str(Rr(s)) == t, (s, rnd, str(Rr(s)), t)
            assert str(R40(s).round(10, rnd)) == t, (s, rnd)
            assert str(R40(s).round(10, rnd)) == str(Rr(R40(s)))
    assert str(R40("1/3").round(3)) == "0.333"
    assert str(R40("1/3").round(1)) == "0.3"
    assert str(X("0.5").round(1, "up")) == "0.5"
    assert str(X("0.95").round(1, "near")) == "1"
    assert str(X("0.95").round(1, "down")) == "0.9"
    assert str(X("999").round(2, "up")) == "1000"
    assert str(X("999").round(2, "down")) == "990"
    assert str(X("-999").round(2, "floor")) == "-1000"
    assert str(X("-999").round(2, "ceil")) == "-990"
    assert R40("1/3").round(None) == R40("1/3")
    assert raises(lambda: X("1").round(1, "sideways"), ValueError)

    # exact arithmetic
    assert str(X("0.1") + X("0.2")) == "0.3"
    assert str(X("0.1") * X("0.2")) == "0.02"
    assert str(X("1e20") + X("1e-20")) == "100000000000000000000.00000000000000000001"
    assert str(X("1e20") - X("1e-20")) == "99999999999999999999.99999999999999999999"
    assert str(X("1e-20") - X("1e20")) == "-99999999999999999999.99999999999999999999"
    assert str(X("1e20") * X("1e-20")) == "1"
    assert str(X("1e20") / X("1e-20")) == "1e40"
    assert str(X("1") / X("1e-20")) == "100000000000000000000"
    assert str(X("1") / X("1e-21")) == "1e21"
    assert str(X("125") / X("1000")) == "0.125"
    assert str(X("1") / X("64")) == "0.015625"
    assert str(X(0) / 3) == "0" and str(X(0) * X("1e100")) == "0"
    assert raises(lambda: X(1) / 3, FlintUnableError)
    assert raises(lambda: X(1) / 0, FlintDomainError)
    assert str(X("1e100") ** 20) == "1e2000"
    assert str(X("1.5") ** 40) == "11057332.3209400121422731899656355381011962890625"
    assert QQ(X("1.5") ** 40) == QQ(3) ** 40 / 2 ** 40
    assert str(X("0.0001").sqrt()) == "0.01"
    assert str(X("1e-40").sqrt()) == "1e-20"
    assert str(X("225").sqrt()) == "15"
    assert raises(lambda: X(2).sqrt(), FlintUnableError)
    assert raises(lambda: X(-1).sqrt(), FlintDomainError)
    assert str(X("12345.678") - X("12345.678")) == "0"
    assert str(X("1e300") * X("1e-300") - 1) == "0"

    # rounded arithmetic and correct rounding
    assert str(R(2).sqrt()) == "1.414213562"
    assert str(RealFloat_decfloat(10, rnd="up")(2).sqrt()) == "1.414213563"
    assert str(RealFloat_decfloat(50)(2).sqrt()) == "1.4142135623730950488016887242096980785696718753769"
    assert str(RealFloat_decfloat(50, rnd="down")(2).sqrt()) == "1.4142135623730950488016887242096980785696718753769"
    assert str(RealFloat_decfloat(50, rnd="up")(2).sqrt()) == "1.414213562373095048801688724209698078569671875377"
    assert str(R("1e10") + 1) == "10000000000"
    assert str(RealFloat_decfloat(10, rnd="up")("1e10") + 1) == "10000000010"
    assert str(RealFloat_decfloat(10, rnd="floor")("-1e10") - 1) == "-10000000010"
    assert str(R("1e10") + R("0.5")) == "10000000000"
    assert str(R("1e10") + R("5")) == "10000000000"
    assert str(R("1e10") + R("5.000000001")) == "10000000010"
    assert str(R("1e10") + R("15")) == "10000000020"
    assert str(R("9999999999") + 1) == "10000000000"
    assert str(R("9999999999") + R("0.5")) == "10000000000"
    assert str(R("1e30") + 1) == "1e30"
    assert str(RealFloat_decfloat(10, rnd="up")("1e30") + 1) == "1.000000001e30"
    assert str(RealFloat_decfloat(10, rnd="floor")("-1e30") - 1) == "-1.000000001e30"
    assert str(R("1e30") + R("5e20")) == "1e30"
    assert str(R("1e30") + R("5.000000001e20")) == "1.000000001e30"
    assert str(R("1e30") + R("1.5e21")) == "1.000000002e30"
    assert str(R("9.999999999e29") + R("1e20")) == "1e30"
    assert str(R("9999999999") + R("0.4")) == "9999999999"
    assert str(R("1.5") * R("1.5")) == "2.25"
    assert str(R("1234567890") * R("9876543210")) == "12193263110000000000"
    assert str(RealFloat_decfloat(10, rnd="down")("1234567890") * R("9876543210")) == "12193263110000000000"
    assert str(RealFloat_decfloat(10, rnd="up")("1234567890") * R("9876543210")) == "12193263120000000000"
    assert str(R("1234567890e10") * R("9876543210")) == "1.219326311e29"
    assert str(RealFloat_decfloat(10, rnd="up")("1234567890e10") * R("9876543210")) == "1.219326312e29"
    assert str(X("1234567890") * X("9876543210")) == "12193263111263526900"
    assert str(R(1) / 7) == "0.1428571429"
    assert str(R(22) / 7) == "3.142857143"
    assert str(R(1) / R("1e-30")) == "1e30"
    assert str(R("1e30").inv()) == "1e-30"
    assert str(R(3).inv()) == "0.3333333333"
    assert str(-R(3).inv()) == "-0.3333333333"
    assert str(abs(-R(3).inv())) == "0.3333333333"
    assert str(R(3).inv() * 3) == "0.9999999999"
    assert str(R("1e100") * R("1e100")) == "1e200"
    assert str(R(1) / 3 - R(1) / 3) == "0"
    assert str(R("1e-100") + R("1e100") - R("1e100")) == "0"

    # properties and comparisons
    assert R("123.4500").digits() == 5
    assert R("123.45").exponent() == 2 and R("0.00123").exponent() == -3 and R("1e100").exponent() == 100
    assert R(0).exponent() is None and R(0).digits() == 0
    assert R(1) < R(2) and R(-1) < R(1) and not (R(1) < R(1)) and R(1) <= R(1)
    assert R("1e100") > R("9.99e99") and R("-1e100") < R("-9.99e99")
    assert R("0.1") == R("0.10") and R("1e2") == 100 and R("0.5") == QQ(1) / 2
    assert R(1) != R(2) and R("1e-30") != 0 and R("1e-30") > 0
    assert R("1e-30").sgn() == 1 and R("-1e-30").sgn() == -1 and R(0).sgn() == 0
    assert Ri("inf") > Ri("1e1000") and Ri("-inf") < Ri("-1e1000")
    assert hash(R(1)) == hash(R("1.0"))

    # conversions
    assert int(R("123.9")) == 123 and int(R("-123.9")) == -123
    assert ZZ(R("1e10")) == 10 ** 10 and ZZ(X("1e100")) == 10 ** 100
    assert QQ(R("0.1")) == QQ(1) / 10 and QQ(X("1e-100")) == QQ(1) / 10 ** 100
    assert float(R("0.5")) == 0.5 and float(R("0.1")) == 0.1 and float(R("1e300")) == 1e300
    assert float(R(1) / 3) == 0.3333333333
    assert str(X(0.1)) == "0.1000000000000000055511151231257827021181583404541015625"
    assert str(X(2.0 ** 100)) == "1.267650600228229401496703205376e30"
    assert str(X(2 ** 100)) == "1.267650600228229401496703205376e30"
    assert str(X(-2 ** 100)) == "-1.267650600228229401496703205376e30"
    assert str(X(QQ(1) / 1024)) == "0.0009765625"
    assert raises(lambda: X(QQ(1) / 3), FlintUnableError)
    assert str(R(QQ(1) / 3)) == "0.3333333333"
    assert str(R(2 ** 100)) == "1.2676506e30"
    assert str(RealFloat_decfloat(10, rnd="up")(2 ** 100)) == "1.267650601e30"
    assert str(R(RF("0.1"))) == "0.1"
    assert str(X(RF("0.1"))) == "0.1000000000000000055511151231257827021181583404541015625"
    assert str(RF(R("0.1"))) == "0.1000000000000000"
    assert RealFloat_arf(200)(X("1e-100")) < RF("1e-100") and RealFloat_arf(200)(X("1e-100")) > RF("0.999e-100")
    assert RR(X("0.125")) == QQ(1) / 8 and str(RR(X("0.125"))) == "0.1250000000000000" and RR(X(2**100)) == 2**100
    assert str(RR(R("0.1"))) == "[0.100000000000000 +/- 2.23e-17]" and RR(R("0.1")).overlaps(RR(1) / 10)
    assert str(R(RR("0.1"))) == "0.1"
    assert raises(lambda: X(RR(1) / 3), FlintUnableError)
    assert raises(lambda: ZZ(R("0.5")), FlintDomainError)
    assert raises(lambda: ZZ(X("1e1000000000000")), FlintUnableError)
    assert str(R(1) + 1) == "2" and str(1 + R(1)) == "2" and str(R(1) + QQ(1) / 2) == "1.5"
    assert str(R("0.5") * ZZ(3)) == "1.5" and str(ZZ(3) * R("0.5")) == "1.5"
    assert str(R(1) / ZZ(3)) == "0.3333333333"
    assert str(R(2) ** 3) == "8" and str(R(2) ** -3) == "0.125" and str(R(2) ** QQ(1)/2) == "1"
    assert str(R(2) ** (QQ(1) / 2)) == "1.414213562"
    assert str(R(9) ** (QQ(1) / 2)) == "3"
    assert str(R(2) ** R("0.5")) == "1.414213562"
    assert str(R(2) ** R(10)) == "1024"
    assert str(R(3).floor()) == "3" and str(R("3.7").floor()) == "3" and str(R("-3.7").floor()) == "-4"
    assert str(R("3.2").ceil()) == "4" and str(R("-3.2").ceil()) == "-3"
    assert str(R("3.5").nint()) == "4" and str(R("2.5").nint()) == "2" and str(R("-2.5").nint()) == "-2"
    assert str(R("3.7").trunc()) == "3" and str(R("-3.7").trunc()) == "-3"
    assert str(R("1e100").floor()) == "1e100"
    assert str(R("123456789012345").floor()) == "123456789000000"

    # exponent limits
    Rl = RealFloat_decfloat(5, exp_limits=(-10, 10))
    assert str(Rl("1e10")) == "10000000000" and str(Rl("9.9999e10")) == "99999000000" and str(Rl("1e-10")) == "1e-10"
    assert raises(lambda: Rl("1e11"), FlintUnableError)
    assert raises(lambda: Rl("1e-11"), FlintUnableError)
    assert raises(lambda: Rl("1e10") * 10, FlintUnableError)
    assert raises(lambda: Rl("9.9999e10") + Rl("1e6"), FlintUnableError)
    assert str(Rl("9.9999e10") + Rl("1e5")) == "99999000000"
    assert raises(lambda: Rl("1e-10") / 10, FlintUnableError)
    assert raises(lambda: Rl("1e-10") / 2, FlintUnableError)
    assert str(Rl("1e-10") * 2) == "2e-10"
    assert str(Rl.exp_limits) == "(-10, 10)"
    Rl = RealFloat_decfloat(5, exp_limits=(-10, 10), inf=True, underflow=True)
    assert str(Rl("1e11")) == "inf" and str(Rl("-1e11")) == "-inf" and str(Rl("1e-11")) == "0"
    assert str(Rl("1e10") * 10) == "inf" and str(Rl("1e-10") / 10) == "0" and str(Rl("1e-10") / 2) == "0"
    assert str(Rl("99999e6") + Rl("1e6")) == "inf"
    assert str(Rl(2) ** 100) == "inf" and str(Rl(2) ** -100) == "0"
    Rl.exp_limits = (None, 3)
    assert str(Rl("1e-100")) == "1e-100" and str(Rl("1e4")) == "inf"
    Rl.exp_limits = None
    assert str(Rl("1e4")) == "10000"
    Ru = RealFloat_decfloat(5, exp_limits=(-10, 10), underflow=True)
    assert str(Ru("1e-11")) == "0" and raises(lambda: Ru("1e11"), FlintUnableError)

    # special values
    inf = Ri("inf"); nan = Ri("nan")
    assert str(inf + 1) == "inf" and str(inf - 1) == "inf" and str(1 - inf) == "-inf" and str(-inf) == "-inf"
    assert str(inf - inf) == "nan" and str(inf * 0) == "nan" and str(inf * -2) == "-inf" and str(inf / inf) == "nan"
    assert str(1 / inf) == "0" and str(Ri(1) / 0) == "inf" and str(Ri(-1) / 0) == "-inf" and str(Ri(0) / 0) == "nan"
    assert str(inf.sqrt()) == "inf" and str(nan + 1) == "nan" and str(nan * 0) == "nan"
    assert raises(lambda: nan == nan, Undecidable)
    assert inf == inf and inf != -inf and inf > Ri(1) and -inf < Ri(1)
    assert str(Ri("1e400") * Ri("1e400")) == "1e800"
    assert str(Ri(0).inv()) == "inf"
    assert raises(lambda: R(0).inv(), FlintDomainError)

    # elementary functions with correct rounding
    assert str(R.pi()) == "3.141592654"
    assert str(RealFloat_decfloat(10, rnd="down").pi()) == "3.141592653"
    assert str(RealFloat_decfloat(10, rnd="up").pi()) == "3.141592654"
    assert str(RealFloat_decfloat(30).pi()) == "3.14159265358979323846264338328"
    assert str(RealFloat_decfloat(30, rnd="down").pi()) == "3.14159265358979323846264338327"
    assert str(R(1).exp()) == "2.718281828"
    assert str(RealFloat_decfloat(10, rnd="up")(1).exp()) == "2.718281829"
    assert str(R(0).exp()) == "1" and str(R(0).sin()) == "0" and str(R(1).log()) == "0" and str(R(0).atan()) == "0"
    assert str(R(2).log()) == "0.6931471806" and str(R(10).log()) == "2.302585093"
    assert str(R("1e100").log()) == "230.2585093"
    assert str(R("1e-100").exp()) == "1"
    assert str(RealFloat_decfloat(10, rnd="up")("1e-100").exp()) == "1.000000001"
    assert str(RealFloat_decfloat(10, rnd="down")("-1e-100").exp()) == "0.9999999999"
    assert str(R(100).exp()) == "2.688117142e43" and str(R(-100).exp()) == "3.720075976e-44"
    assert str(R(1).sin()) == "0.8414709848" and str(R(1).cos()) == "0.5403023059"
    assert str(R("1e20").sin()) == "-0.6452512853"
    assert str(R(1).atan()) == "0.7853981634" and str(R("1e30").atan()) == "1.570796327"
    assert raises(lambda: R(0).log(), FlintDomainError)
    assert raises(lambda: R(-1).log(), FlintDomainError)
    # tiny arguments (Ziv's strategy alone cannot decide these in directed modes)
    for rnd, sinp, sinn, expp, expn, cosp, atanp, atann in [
            ("near", "1e-1000000000", "-1e-1000000000", "1", "1", "1", "1e-1000000000", "-1e-1000000000"),
            ("down", "9.999999999e-1000000001", "-9.999999999e-1000000001", "1", "0.9999999999", "0.9999999999", "9.999999999e-1000000001", "-9.999999999e-1000000001"),
            ("up", "1e-1000000000", "-1e-1000000000", "1.000000001", "1", "1", "1e-1000000000", "-1e-1000000000"),
            ("floor", "9.999999999e-1000000001", "-1e-1000000000", "1", "0.9999999999", "0.9999999999", "9.999999999e-1000000001", "-1e-1000000000"),
            ("ceil", "1e-1000000000", "-9.999999999e-1000000001", "1.000000001", "1", "1", "1e-1000000000", "-9.999999999e-1000000001")]:
        Rt = RealFloat_decfloat(10, rnd=rnd)
        x = Rt("1e-1000000000"); xn = -x
        assert str(x.sin()) == sinp and str(xn.sin()) == sinn and str(x.exp()) == expp and str(xn.exp()) == expn
        assert str(x.cos()) == cosp and str(xn.cos()) == cosp and str(x.atan()) == atanp and str(xn.atan()) == atann
        assert str(x.tan()) == str(x.sinh()) == str(x.expm1()) and str(x.tanh()) == str(x.atan()) == str(x.log1p()) == str(x.sin())
        assert str(x.cosh()) == ("1.000000001" if rnd in ("up", "ceil") else "1")
        assert str(Rt(2) ** x) == expp and str(Rt("0.5") ** x) == expn and str(Rt(2) ** xn) == expn and str(Rt("0.5") ** xn) == expp
        assert str(Rt("1e-100000000000000000000").sin()) == sinp.replace("1000000000", "100000000000000000000").replace("1000000001", "100000000000000000001")
        assert str(Rt("123456789e-1000000000").sin()) in ("1.23456789e-999999992", "1.234567889e-999999992", "1.23456789e-999999992")
        assert str(Rt("-123456789e-1000000000").cos()) == cosp
    x = RealFloat_decfloat(3, rnd="down")("1e-100")
    assert str(x.sin()) == "9.99e-101" and str((x + 1).log()) == "0" and str(x.log1p()) == "9.99e-101" and str(x.expm1()) == "1e-100"
    assert str(RealFloat_decfloat(3, rnd="up")("1e-100").sin()) == "1e-100" and str(RealFloat_decfloat(3, rnd="up")("1e-100").expm1()) == "1.01e-100"
    assert str(R("1e-30").sin()) == "1e-30" and str(R("1e-30").exp()) == "1" and str(RealFloat_decfloat(10, rnd="up")("1e-30").exp()) == "1.000000001"
    assert str(RealFloat_decfloat(10, rnd="up")("1e-9").sin()) == "1e-9" and str(RealFloat_decfloat(10, rnd="down")("1e-9").sin()) == "9.999999999e-10"
    assert str(RealFloat_decfloat(10, rnd="down")("1e-6").sin()) == "9.999999999e-7" and str(RealFloat_decfloat(10, rnd="down")("1e-6").tan()) == "0.000001"
    # exact powers
    for rnd in _decimal_rnd_names:
        Rt = RealFloat_decfloat(10, rnd=rnd)
        assert str(Rt(8) ** Rt("0.125")) == ("1.296839555" if rnd in ("up", "ceil", "near", "near_away", "near_zero") else "1.296839554")
        assert str(Rt(256) ** Rt("0.125")) == "2" and str(Rt(256) ** Rt("-0.125")) == "0.5" and str(Rt("0.0001") ** Rt("0.25")) == "0.1"
        assert str(Rt("1e100") ** Rt("0.5")) == "1e50" and str(Rt("1e100") ** Rt("-0.5")) == "1e-50" and str(Rt("1e-100") ** Rt("0.5")) == "1e-50"
        assert str(Rt("1.5") ** 3) == "3.375" and str(Rt("1.5") ** 4) == "5.0625" and str(Rt("6.25") ** Rt("1.5")) == "15.625"
        assert str(Rt("1.5") ** -3) == ("0.2962962963" if rnd in ("up", "ceil", "near", "near_away", "near_zero") else "0.2962962962")
        assert str(Rt("-1.5") ** -3) == ("-0.2962962963" if rnd in ("up", "floor", "near", "near_away", "near_zero") else "-0.2962962962")
        assert str(Rt("1e-5") ** -(10**15)) == "1e5000000000000000" and str(Rt("1e-10") ** 10**8) == "1e-1000000000"
        assert str(Rt(3) ** 10**9) == ("5.243997032e477121254" if rnd in ("down", "floor") else "5.243997033e477121254")
        assert str(Rt(3) ** -10**9) == ("1.906942346e-477121255" if rnd in ("up", "ceil") else "1.906942345e-477121255")
        assert str(Rt(4) ** Rt("1e-1000000000")) == ("1.000000001" if rnd in ("up", "ceil") else "1")
        assert str(RealFloat_decfloat(40, rnd=rnd)("1.0000000000000000000000000000001") ** RealFloat_decfloat(40, rnd=rnd)("1e10")) == ("1.000000000000000000001000000000000000001" if rnd in ("up", "ceil") else "1.000000000000000000001")
        assert str(Rt("1.000000001") ** Rt("1e100"))[:11] == ("5.601255459" if rnd in ("up", "ceil") else "5.601255458")
        assert str(Rt("0.9999999999") ** Rt("1e100"))[:11] == ("2.049262995" if rnd in ("down", "floor") else "2.049262996")
        assert str(Rt("0.9999999999") ** Rt("1e100"))[11:] == "e-434294481924966551747739158572280941286215521827993337769370509468333092239211207443855568"
        assert str(Rt(1) ** Rt("1e100")) == "1" and str(Rt(1) ** Rt("-1e100")) == "1" and str(Rt(10) ** 1000000000) == "1e1000000000"
        assert str(Rt(2) ** 100) == ("1.267650601e30" if rnd in ("up", "ceil") else "1.2676506e30")
        assert str(Rt(2) ** -100) == ("7.888609053e-31" if rnd in ("up", "ceil") else "7.888609052e-31")
    assert str(X(2) ** 100) == "1.267650600228229401496703205376e30" and raises(lambda: X("1.5") ** -4, FlintUnableError)
    assert str(X("0.5") ** -4) == "16" and str(X(2) ** -4) == "0.0625" and str(X(256) ** X("0.125")) == "2" and str(X("0.0001") ** X("0.25")) == "0.1"
    assert raises(lambda: X(2) ** X("0.5"), FlintUnableError) and raises(lambda: X(3) ** -1, FlintUnableError)
    # special functions and special points
    assert str(R.gamma(5)) == "24" and str(R.gamma(3000)) == "1.383119868e9127" and str(R.rgamma(5)) == "0.04166666667" and str(R.rgamma(-4)) == "0"
    assert str(R.zeta(-1)) == "-0.08333333333" and str(R.zeta(-3)) == "0.008333333333" and str(R.zeta(-4)) == "0" and str(R.zeta(0)) == "-0.5" and str(R.zeta(-1001)) == "-1.348590824e1771"
    assert str(R.log10(R("1e100"))) == "100" and str(R.log10(R("1e-100"))) == "-100" and str(R.log2(2 ** 100)) == "100" and str(R.log2(R("0.5") ** 100)) == "-100"
    assert str(R.lgamma(1)) == "0" and str(R.lgamma(2)) == "0" and str(R.barnes_g(1)) == "1" and str(R.barnes_g(4)) == "2" and str(R.barnes_g(6)) == "288"
    assert str(R.sin_pi(R("0.5"))) == "1" and str(R.sin_pi(1)) == "0" and str(R.sin_pi(R("1.5"))) == "-1" and str(R.cos_pi(R("0.5"))) == "0" and str(R.cos_pi(1)) == "-1"
    assert str(R.tan_pi(R("0.25"))) == "1" and str(R.tan_pi(R("0.75"))) == "-1" and str(R.cot_pi(R("0.25"))) == "1" and str(R.csc_pi(R("0.5"))) == "1" and str(R.sec_pi(1)) == "-1"
    assert str(R.asin_pi(1)) == "0.5" and str(R.asin_pi(R("-0.5"))) == "-0.1666666667" and str(R.acos_pi(-1)) == "1" and str(R.acos_pi(0)) == "0.5" and str(R.atan_pi(-1)) == "-0.25"
    assert str(R.acos(1)) == "0" and str(R.acosh(1)) == "0" and str(R.sinc(0)) == "1" and str(R.sinc_pi(3)) == "0" and str(R.lambertw(0)) == "0" and str(R.erfinv(0)) == "0"
    assert raises(lambda: R.tan_pi(R("0.5")), FlintDomainError) and raises(lambda: R.cot_pi(1), FlintDomainError) and raises(lambda: R.gamma(-2), FlintDomainError)
    assert raises(lambda: R.zeta(1), FlintDomainError) and raises(lambda: R.lgamma(0), FlintDomainError) and raises(lambda: R.log2(0), FlintDomainError)
    for rnd in ["down", "up", "floor", "ceil", "near"]:
        Rt = RealFloat_decfloat(10, rnd=rnd)
        # S + tiny rounds to S + ulp only when rounding away; S - tiny rounds to S - ulp only when rounding toward zero
        up = rnd in ("up", "ceil")            # S + tiny, S > 0
        upn = rnd in ("up", "ceil", "near")   # S - tiny, S > 0: rounds to S unless toward zero
        neg_away = rnd in ("up", "floor")     # for S < 0: away from zero
        neg_stay = rnd in ("up", "floor", "near")   # S < 0, S + tiny (toward zero): stays at S unless toward zero
        x = Rt("1e-1000000000")
        assert str(Rt.gamma(x)) == ("1e1000000000" if upn else "9.999999999e999999999")
        assert str(Rt.gamma(-x)) == ("-1.000000001e1000000000" if neg_away else "-1e1000000000")
        assert str(Rt.digamma(x)) == ("-1.000000001e1000000000" if neg_away else "-1e1000000000")
        assert str(Rt.cot(x)) == ("1e1000000000" if upn else "9.999999999e999999999")
        assert str(Rt.csc(x)) == ("1.000000001e1000000000" if up else "1e1000000000")
        assert str(Rt.coth(x)) == ("1.000000001e1000000000" if up else "1e1000000000")
        assert str(Rt.csch(x)) == ("1e1000000000" if upn else "9.999999999e999999999")
        assert str(Rt.rgamma(x)) == ("1.000000001e-1000000000" if up else "1e-1000000000")
        assert str(Rt.sec(x)) == str(Rt.cosh(x)) == ("1.000000001" if up else "1")
        assert str(Rt.sech(x)) == str(Rt.cos_pi(x)) == str(Rt.sinc(x)) == str(Rt.sinc_pi(x)) == ("1" if upn else "0.9999999999")
        assert str(Rt.lambertw(x)) == str(Rt.asinh(x)) == str(Rt.sin_integral(x)) == ("1e-1000000000" if upn else "9.999999999e-1000000001")
        assert str(Rt.asin(x)) == str(Rt.atanh(x)) == str(Rt.dilog(x)) == str(Rt.sinh_integral(x)) == ("1.000000001e-1000000000" if up else "1e-1000000000")
        assert str(Rt.tanh(1000000)) == ("1" if upn else "0.9999999999") and str(Rt.tanh(-1000000)) == ("-1" if neg_stay else "-0.9999999999")
        assert str(Rt.coth(1000000)) == ("1.000000001" if up else "1") and str(Rt.erf(-1000000)) == ("-1" if neg_stay else "-0.9999999999")
        assert str(Rt.erfc(-1000000)) == ("2" if upn else "1.999999999") and str(Rt.expm1(-1000000)) == ("-1" if neg_stay else "-0.9999999999")
        assert str(Rt.zeta(10 ** 9)) == ("1.000000001" if up else "1") and str(Rt.acot(Rt("1e100"))) == ("1e-100" if upn else "9.999999999e-101")
        assert str(Rt.acoth(Rt("1e100"))) == str(Rt.acsc(Rt("1e100"))) == ("1.000000001e-100" if up else "1e-100") and str(Rt.acsch(Rt("-1e100"))) == ("-1e-100" if neg_stay else "-9.999999999e-101")
    R110 = RealFloat_decfloat(110, rnd="down")
    assert str(R110.zeta(R110("1e-100") + 1)) == "1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000577215664e100"
    assert str(R110.gamma(R110("1e-100"))) == "9.999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999422784335e99"
    assert str(Ri(0).log()) == "-inf" and str(Ri("inf").exp()) == "inf" and str(Ri("-inf").exp()) == "0"
    assert str(Ri("inf").log()) == "inf" and str(Ri("inf").atan()) == "1.570796327"
    assert raises(lambda: Ri("inf").sin(), FlintUnableError)

    # polynomials and matrices
    P = PolynomialRing(R)
    f = P("x^2 - 3*x + 1.25")
    assert str(f) == "1.25 - 3*x + x^2"
    assert str(P(str(f))) == str(f)
    assert str(f(R("0.5"))) == "0" and str(f(R(2.5))) == "0"
    assert str(f * f) == "1.5625 - 7.5*x + 11.5*x^2 - 6*x^3 + x^4"
    assert str(P("(x + 1/3)^2")) == "0.1111111111 + 0.6666666666*x + x^2"
    assert str(P("(x + 1/3) * (x - 1/3)")) == "-0.1111111111 + x^2"
    assert str(P("x^2 - 2") % P("x - 1.5")) == "0.25"
    assert str(P("x^3 + 2*x") / P("x")) == "2 + x^2"
    assert str(P("x^2-2").derivative()) == "2*x"
    M = Mat(R, 2, 2)([[1, 2], [3, 4]])
    assert str(M.det()) == "-2" and str(M.inv()) == "[[-1.999999998, 0.9999999992],\n[1.499999999, -0.4999999996]]"
    assert str(Mat(R, 2, 2)(str(M))) == str(M)
    assert str(Mat(R, 2, 2)("[[1/3, 2], [3, 4]]")) == "[[0.3333333333, 2],\n[3, 4]]"
    assert str(M * M.inv()) == "[[1, 0],\n[2e-9, 1]]"
    assert str(Mat(R, 3, 3)().hilbert().det()) == "0.000462962953"
    assert str(Mat(RealFloat_decfloat(20), 3, 3)().hilbert().det()) == "0.000462962962962962962"
    assert str(Mat(X, 2, 2)([[1, 2], [3, 5]]).inv()) == "[[-5, 2],\n[3, -1]]"
    assert str(Mat(X, 2, 2)([[1, 2], [3, 4]]).inv()) == "[[-2, 1],\n[1.5, -0.5]]"
    assert raises(lambda: Mat(X, 2, 2)([[1, 2], [3, 9]]).inv(), FlintUnableError)
    assert raises(lambda: Mat(X, 2, 2)([[1, 2], [2, 4]]).inv(), FlintDomainError)
    S = PowerSeriesRing(R, 6)
    assert str(S("(1 + x/3)^2")) == "1 + 0.6666666666*x + 0.1111111111*x^2"
    assert str(S("1/(1 - x/3)")) == "1 + 0.3333333333*x + 0.1111111111*x^2 + 0.03703703703*x^3 + 0.01234567901*x^4 + 0.004115226336*x^5 + O(x^6)"
    assert str(S("1 + x").exp()) == "2.718281828 + 2.718281828*x + 1.359140914*x^2 + 0.4530469713*x^3 + 0.1132617428*x^4 + 0.02265234856*x^5 + O(x^6)"
    assert str(S("1 + x").log()) == "x - 0.5*x^2 + 0.3333333333*x^3 - 0.25*x^4 + 0.2*x^5 + O(x^6)"
    assert str(S("(1 + x)^10")) == "1 + 10*x + 45*x^2 + 120*x^3 + 210*x^4 + 252*x^5 + O(x^6)"

    # context settings
    Rc = RealFloat_decfloat(7, rnd="floor", limb_digits=3)
    assert str(Rc) == "Decimal floating-point numbers (prec 7, rnd floor, limb 10^3)"
    assert Rc.digits == 7 and Rc.rnd == "floor" and Rc.limb_digits == 3 and Rc.exp_limits is None
    assert str(Rc(1) / 3) == "0.3333333" and str(Rc(-1) / 3) == "-0.3333334"
    assert str(Rc("1e100") + Rc("1e-100")) == "1e100"
    Rc.digits = 3; Rc.rnd = "up"
    assert str(Rc(1) / 3) == "0.334" and str(Rc(2).sqrt()) == "1.42" and str(Rc(2).sqrt() ** 2) == "2.02"
    assert Rc.prec == 8
    Rc.prec = 100
    assert Rc.digits == 31
    Rc.digits = None
    assert Rc.digits is None
    assert raises(lambda: Rc(1) / 3, FlintUnableError)
    for e in [1, 2, 3, 5, 9, 18]:
        Re = RealFloat_decfloat(25, limb_digits=e)
        assert str(Re("1e100") + Re("1e-100")) == "1e100"
        assert str(Re(1) / 7) == "0.1428571428571428571428571"
        assert str(Re(2).sqrt()) == "1.414213562373095048801689"
        assert str(Re("1e-100").exp()) == "1"
        assert str(Re("123456789.123456789") * Re("987654321.987654321")) == "121932631356500531.3472032"
        assert str(RealFloat_decfloat(None, limb_digits=e)("123456789.123456789") * Re("987654321.987654321")) == "121932631356500531.347203169112635269"
        assert str(Re("123456789.123456789") - Re("123456789.123456788")) == "1e-9"
        assert str(Re.pi()) == "3.141592653589793238462643"
        assert str(PolynomialRing(Re)("(x - 1/3) * (x + 1/3)")) == "-0.1111111111111111111111111 + x^2"

    # balls
    assert str(B) == "Decimal balls (prec 10, rad prec 4)"
    x = B(1) / 3
    assert str(x) == "[0.3333333333 +/- 3.334e-11]"
    assert str(x.mid()) == "0.3333333333" and str(x.rad()) == "3.334e-11"
    assert x.mid().parent().digits == 10 and x.rad().parent().digits is None
    assert not x.is_exact() and B(1).is_exact() and B("1e100").is_exact()
    assert x.contains(QQ(1) / 3) and not x.contains(QQ(1) / 2) and x.contains(x)
    assert (x * 3).contains(1) and (x * 3).overlaps(1) and not (x * 3).contains(2)
    assert x.overlaps(x) and not x.overlaps(x + 1)
    assert raises(lambda: x == x, Undecidable) and raises(lambda: x != x, Undecidable)
    assert x != x + 1 and x < x + 1 and x > x - 1 and raises(lambda: x < x, Undecidable)
    assert B(1) / 3 * 3 != 2 and not (B(1) / 3 * 3 == 2)
    assert B(1) == B(1) and B("0.5") == QQ(1) / 2
    assert B("[1 +/- 0.1]").contains(B("[1 +/- 0.1]")) and not B("[1 +/- 0.1]").contains(B("[1 +/- 0.11]"))
    assert B("[1 +/- 0.1]").contains(B("[1.05 +/- 0.05]")) and not B("[1 +/- 0.1]").contains(B("[1.05 +/- 0.06]"))
    assert B("[1 +/- 0.1]").overlaps(B("[1.2 +/- 0.1]")) and not B("[1 +/- 0.1]").overlaps(B("[1.2 +/- 0.09]"))
    assert B("[1 +/- 0.1]").contains(B("[1 +/- 0.1]").mid())
    assert str(B("[1 +/- 0.1]").mid()) == "1" and str(B("[1 +/- 0.1]").rad()) == "0.1"
    assert str(B("[1 +/- 0.1]").add_error(B("0.05"))) == "[1 +/- 0.15]"
    assert str(B("[1 +/- 0.1]").add_error(B("[0.05 +/- 0.01]"))) == "[1 +/- 0.16]"
    assert str(B("[1 +/- 0.1]").add_error(-B("0.05"))) == "[1 +/- 0.15]"
    assert str(B(1).add_error_10exp(-3)) == "[1 +/- 0.001]" and str(B(1).add_error_10exp(-30)) == "[1 +/- 1e-30]"
    assert str(B(0).add_error_10exp(5)) == "[0 +/- 100000]"
    assert str(B("[1234.5678 +/- 0.5]").trim()) == "[1234.5678 +/- 0.5]"
    assert str(B("[1234.56789012 +/- 0.5]").trim()) == "[1234.56789 +/- 0.5001]"
    assert B("[1234.56789012 +/- 0.5]").trim().contains(B("[1234.56789012 +/- 0.5]"))
    assert B("[1234.5678 +/- 0.5]").rel_accuracy_digits() == 4
    assert B("[1234.5678 +/- 0.001]").rel_accuracy_digits() == 6
    assert B("[1e100 +/- 1e90]").rel_accuracy_digits() == 10
    assert B("[1e100 +/- 1e120]").rel_accuracy_digits() == -20
    assert B("[0 +/- 1e-20]").rel_accuracy_digits() == 20
    assert B(1).rel_accuracy_digits() is None
    assert (B(1) / 3).rel_accuracy_digits() == 10
    # string roundtrips
    for s, t in [("[1 +/- 0.1]", "[1 +/- 0.1]"), ("1 +/- 0.1", "[1 +/- 0.1]"), ("[1+/-0.1]", "[1 +/- 0.1]"),
                 ("+/- 0.1", "[0 +/- 0.1]"), ("[+/- 0.1]", "[0 +/- 0.1]"), ("[-1.5e3 +/- 1e-3]", "[-1500 +/- 0.001]"),
                 ("[1 +/- 0.123456]", "[1 +/- 0.1235]"), ("[1 +/- 1e-100]", "[1 +/- 1e-100]"),
                 ("[1 +/- 0]", "1"), ("[1]", "1"), ("(1)", "1"), ("1e100 +/- 1e50", "[1e100 +/- 1e50]"),
                 ("1/3 +/- 1e-30", "[0.3333333333 +/- 3.335e-11]"), ("(1/3) +/- 1e-30", "[0.3333333333 +/- 3.335e-11]"),
                 ("1/3", "[0.3333333333 +/- 3.334e-11]"), ("1/4", "0.25"), ("1e-30 + 1", "[1 +/- 1e-30]"),
                 ("(1 +/- 0.1)^2", "[1 +/- 0.21]"), ("(1 +/- 0.1)*(1 +/- 0.1)", "[1 +/- 0.21]"),
                 ("[1 +/- [0.1 +/- 0.01]]", "[1 +/- 0.11]"), ("[[1 +/- 0.1] +/- 0.1]", "[1 +/- 0.2]"),
                 ("-[1 +/- 0.1]", "[-1 +/- 0.1]"), ("2 - [1 +/- 0.1]", "[1 +/- 0.1]"),
                 ("123456789012 +/- 1", "[123456789000 +/- 13]"), ("12345678901234567890", "[12345678900000000000 +/- 1.235e9]"),
                 ("123456789012345678901", "[123456789000000000000 +/- 1.235e10]"), ("-1234567890123456789012 +/- 1e5", "[-1.23456789e21 +/- 1.236e11]"),
                 ("[1 +/- 0.001]", "[1 +/- 0.001]"), ("[1 +/- 0.0001]", "[1 +/- 1e-4]"), ("[1 +/- 12345]", "[1 +/- 12350]"), ("[1 +/- 123456]", "[1 +/- 123500]"), ("[1 +/- 1234567]", "[1 +/- 1.235e6]")]:
        assert str(B(s)) == t, (s, str(B(s)), t)
        assert B(str(B(s))).contains(B(s)), s
    for s in ["[1 +/- x]", "[1 +/-]", "[1 +/- 1", "+/-", "[]", "[1 +/- 1e]"]:
        assert raises(lambda: B(s), (FlintUnableError, FlintDomainError, ValueError)), s
    # balls represent real numbers: no infinite or undefined midpoints (an
    # infinite radius denotes the whole real line)
    assert str(B("[1 +/- inf]")) == "[1 +/- inf]" and str(B("[+/- inf]")) == "[0 +/- inf]" and str(B("+/- inf")) == "[0 +/- inf]"
    assert str(B("[+/- inf]") + 1) == "[1 +/- inf]" and str(B("[+/- inf]") * 0) == "0" and str(B("[1 +/- inf]") + 1) == "[2 +/- inf]"
    assert B("[1 +/- inf]").contains(10 ** 100) and B("[1 +/- inf]").contains(B("[+/- inf]")) and not B("[1 +/- inf]").is_exact()
    assert str(B(B("[+/- inf]"))) == "[0 +/- inf]" and str(RR(B("[+/- inf]"))) == "[+/- inf]" and str(B(RR("[+/- inf]"))) == "[0 +/- inf]"
    for s in ["inf", "-inf", "nan", "[inf +/- 1]", "1/0", "0/0", "[+/- inf] * 0 + inf"]:
        assert raises(lambda: B(s), (FlintUnableError, FlintDomainError)), s
    assert raises(B.neg_inf, FlintDomainError) and raises(B.undefined, FlintDomainError) and raises(B.unknown, FlintDomainError)
    assert raises(lambda: B("[+/- inf]").lower(), FlintDomainError) and raises(lambda: B("[+/- inf]").abs_upper(), FlintDomainError)
    assert str(B("[+/- inf]").abs_lower()) == "0" and raises(lambda: B("[+/- inf]").rad_ball(), FlintDomainError)
    assert raises(lambda: B(0).log(), FlintDomainError) and raises(lambda: B(0).rsqrt(), FlintDomainError) and raises(lambda: B.gamma(B(0)), FlintDomainError)
    assert raises(lambda: B("[0 +/- 0.1]").log(), FlintUnableError) and raises(lambda: B(1) / B("[0 +/- 0.1]"), FlintUnableError)
    Bi = RealField_decball(10, inf=True)   # the inf and nan flags have no effect on balls
    assert raises(lambda: Bi("1/0"), (FlintUnableError, FlintDomainError)) and raises(lambda: Bi("inf"), (FlintUnableError, FlintDomainError))
    assert str(Bi("1e400") ** 3) == "1e1200" and str(Bi) == "Decimal balls (prec 10, rad prec 4, inf)"
    # containment through arithmetic
    third = QQ(1) / 3
    for a, b in [(B(1) / 3, third), (B(1) / 3 + B(1) / 7, third + QQ(1) / 7), (B(1) / 3 * (B(1) / 7), third / 7),
                 ((B(1) / 3) / (B(1) / 7), third * 7), (B(2).sqrt() ** 2, 2), (B(1) / 3 - B(1) / 3, 0),
                 (B("1e-30") + 1, QQ(10) ** -30 + 1), (B("1e30") + 1, QQ(10) ** 30 + 1), ((B(1) / 3) ** 10, third ** 10),
                 (B(2).sqrt().inv(), None), ((B(1) / 3).abs(), third), (-(B(1) / 3), -third), ((B(1) / 3) * 3 - 1, 0),
                 (B("[1 +/- 0.5]") * B("[1 +/- 0.5]"), QQ(9) / 4), (B("[1 +/- 0.5]") - B("[1 +/- 0.5]"), 1),
                 (B("[1 +/- 0.5]") / B("[1 +/- 0.5]"), 3), (B("[1 +/- 0.5]").sqrt(), QQ(3) / 4),
                 (B("[1 +/- 0.5]").sqrt(), 1), (B(10) ** -20, QQ(10) ** -20), (B(3).inv(), third)]:
        if b is not None:
            assert a.contains(b), (a, b)
    assert not (B(1) / 3 + B(1) / 7).contains(third + QQ(1) / 7 + QQ(1) / 10 ** 9)
    assert B(2).sqrt().inv().overlaps(B(2).sqrt() / 2)
    assert raises(lambda: B("[1 +/- 2]").sqrt(), FlintUnableError)
    assert raises(lambda: B(0) / B("[1 +/- 2]"), FlintUnableError)
    assert raises(lambda: B(1) / B("[1 +/- 1]"), FlintUnableError)
    assert raises(lambda: B("[-1 +/- 0.5]").sqrt(), FlintDomainError)
    assert raises(lambda: B(-1).sqrt(), FlintDomainError)
    assert raises(lambda: B(1) / B(0), FlintDomainError)
    assert str(B(0).sqrt()) == "0" and raises(lambda: B("[0 +/- 1e-20]").sqrt(), FlintUnableError)
    assert str(B("[4 +/- 1e-20]").sqrt()) == "[2 +/- 2.502e-21]"
    assert str(B(4).sqrt()) == "2" and str(B("1e-40").sqrt()) == "1e-20" and str(B(2).sqrt()) == "[1.414213562 +/- 3.731e-10]"
    assert str(B(1) / 4) == "0.25" and str(B(1) / 8) == "0.125" and str(B("1e20") + B("1e-20")) == "[100000000000000000000 +/- 1e-20]"
    assert str(B("1e30") + B("1e-20")) == "[1e30 +/- 1e-20]"
    assert str(B(1) / 3 - B(1) / 3) == "[0 +/- 6.668e-11]"
    assert str(B("[1 +/- 0.5]") * 0) == "0" and str(B("[1 +/- 0.5]") - B("[1 +/- 0.5]")) == "[0 +/- 1]"
    assert str(B("[0 +/- 1]") ** 2) == "[0 +/- 1]" and str(B("[0 +/- 1]") ** 3) == "[0 +/- 1]"
    assert str(B("[0 +/- 1]") * B("[0 +/- 2]")) == "[0 +/- 2]"
    assert str(B("[1 +/- 1]") ** 2) == "[1 +/- 3]"
    assert str(B("[-1 +/- 1]").abs()) == "[1 +/- 1]" and str(B("[-1 +/- 0.5]").abs()) == "[1 +/- 0.5]"
    assert str(B("[10 +/- 1]").floor()) == "[10 +/- 2]" and str(B("[10.5 +/- 0.2]").floor()) == "10"
    assert str(B("[10.5 +/- 0.2]").ceil()) == "11" and str(B("[10.5 +/- 0.2]").nint()) == "[10 +/- 1.2]"
    assert str(B("[10.4 +/- 0.05]").nint()) == "10" and str(B("[10.4 +/- 0.05]").trunc()) == "10" and str(B("[-10.4 +/- 0.05]").trunc()) == "-10"
    assert B("[10 +/- 1]").floor().contains(9) and B("[10 +/- 1]").floor().contains(11)
    # radius precision and precise mode
    for rp in range(1, 10):
        Bp = RealField_decball(10, rad_prec=rp)
        assert Bp.rad_prec == rp
        assert str(Bp(1) / 3) == "[0.3333333333 +/- %se-11]" % ("4" if rp == 1 else "3." + "3" * (rp - 2) + "4")
        assert str(Bp("[1 +/- 0.123456789]")) == "[1 +/- %s]" % ("0.123456789" if rp == 9 else "0." + "123456789"[:rp - 1] + str(int("123456789"[rp - 1]) + 1))
        assert (Bp(1) / 3).contains(third) and (Bp(2).sqrt() ** 2).contains(2)
        assert (Bp(1) / 3 + Bp(1) / 7).contains(third + QQ(1) / 7)
        assert (Bp("[1 +/- 0.5]") ** 5).contains(QQ(3) ** 5 / 2 ** 5)
    assert str(RealField_decball(10, rad_prec=0)) == "Decimal balls (prec 10, rad prec 1)"
    assert str(RealField_decball(10, rad_prec=100)) == "Decimal balls (prec 10, rad prec 9)"
    Bq = RealField_decball(10, rad_prec=9, sloppy_radius=False)
    assert str(Bq(1) / 3) == "[0.3333333333 +/- 3.33333334e-11]"
    assert str(Bq(2) / 3) == "[0.6666666667 +/- 3.33333334e-11]"
    assert str(Bq(1) / 4) == "0.25" and str(Bq(2).sqrt()) == "[1.414213562 +/- 3.73095049e-10]"
    assert str(Bq("1e-30") + 1) == "[1 +/- 1e-30]" and str(Bq("1e30") + 1) == "[1e30 +/- 1]"
    assert str(Bq("12345678901") * 1) == "[12345678900 +/- 1]" and str(Bq("12345678905") * 1) == "[12345678900 +/- 5]"
    assert str(Bq("12345678905.0001") * 1) == "[12345678910 +/- 4.9999]"
    assert str(Bq(1) / 3 * 3) == "[0.9999999999 +/- 1.00000001e-10]"
    assert (Bq(1) / 3 * 3).contains(1)
    Bd = RealField_decball(10, rnd="down")
    assert str(Bd(1) / 3) == "[0.3333333333 +/- 3.334e-11]" and str(Bd(2) / 3) == "[0.6666666666 +/- 6.667e-11]"
    assert (Bd(2) / 3).contains(QQ(2) / 3)
    Bd = RealField_decball(10, rnd="up", sloppy_radius=False)
    assert str(Bd(2) / 3) == "[0.6666666667 +/- 3.334e-11]" and (Bd(2) / 3).contains(QQ(2) / 3)
    # exponent limits for balls
    Bl = RealField_decball(10, exp_limits=(-20, 20), inf=True, underflow=True)
    assert str(Bl("1e-30")) == "[0 +/- 1e-20]" and raises(lambda: Bl("1e30"), FlintUnableError) and raises(lambda: Bl("-1e30"), FlintUnableError)
    assert str(Bl("1e-15") * Bl("1e-15")) == "[0 +/- 1e-20]" and raises(lambda: Bl("1e-15") * Bl("1e-15") != 0, Undecidable)
    assert (Bl("1e-15") * Bl("1e-15")).contains(RealField_decball(10)("1e-30")) and Bl("1e-30").contains(RealField_decball(10)("1e-30"))
    assert raises(lambda: Bl("1e15") * Bl("1e15"), FlintUnableError)
    assert str(Bl("[1 +/- 1e-30]")) == "[1 +/- 1e-20]" and Bl("[1 +/- 1e-30]").contains(1) and Bl("[1 +/- 1e-30]").overlaps(RealField_decball(40)("1 + 1e-30"))
    assert raises(lambda: Bl(2) ** 100, FlintUnableError) and str(Bl(2) ** -100) == "[0 +/- 1e-20]"
    Bl = RealField_decball(10, exp_limits=(-20, 20))
    assert raises(lambda: Bl("1e-30"), FlintUnableError) and raises(lambda: Bl("1e30"), FlintUnableError)
    assert raises(lambda: Bl("[1 +/- 1e-30]"), FlintUnableError)
    # conversions
    assert str(RR(B(1) / 3)) == "[0.3333333333 +/- 3.34e-11]"
    assert RR(B(1) / 3).overlaps(RR(1) / 3) and RR(B("[1 +/- 0.1]")).overlaps(RR("[1 +/- 0.1]"))
    assert RR(B("[1e100 +/- 1e90]")).overlaps(RR("1e100") + RR("1e90")) and not RR(B("[1e100 +/- 1e90]")).overlaps(RR("1e100") + RR("1.1e90"))
    assert B(RR(1) / 3).contains(third) and B(RR.pi()).contains(B.pi()) and B.pi().contains(RealField_decball(50).pi())
    assert RealField_decball(50)(RealField_arb(300).pi()).contains(RealField_decball(60).pi())
    assert str(B(RR("[1 +/- 1e-30]"))) == "[1 +/- 1.001e-30]" and str(B(RR("[1e100 +/- 1e50]"))) == "[1e100 +/- 3.205e83]"
    assert str(B(RR("[1e-100 +/- 1e-150]"))) == "[1e-100 +/- 4.021e-117]"
    assert B(RR("[1e-100 +/- 1e-150]")).contains(QQ(10) ** -100 + QQ(10) ** -150)
    assert str(RealField_decball(40)(RR("1e-100"))) == "[1.000000000000000019991899802602883619648e-100 +/- 2.022e-117]"
    assert str(RealField_decball(40)(RealField_arb(300)("1e-100"))) == "[1e-100 +/- 8.991e-192]"
    assert str(B(2 ** 100)) == "[1.2676506e30 +/- 2.283e20]" and str(RealField_decball(31)(2 ** 100)) == "1.267650600228229401496703205376e30"
    assert str(B(QQ(1) / 3)) == "[0.3333333333 +/- 3.334e-11]" and str(B(QQ(1) / 4)) == "0.25"
    assert str(B(0.1)) == "[0.1 +/- 5.552e-18]" and str(RealField_decball(60)(0.1)) == "0.1000000000000000055511151231257827021181583404541015625"
    assert str(B(RF("0.1"))) == "[0.1 +/- 5.552e-18]"
    assert QQ(B("0.125")) == QQ(1) / 8 and ZZ(B("1e5")) == 100000 and float(B("0.5")) == 0.5
    assert raises(lambda: QQ(B(1) / 3), FlintUnableError) and raises(lambda: ZZ(B("0.5")), FlintDomainError)
    assert raises(lambda: ZZ(B("[1 +/- 0.1]")), FlintUnableError) and raises(lambda: ZZ(B("[1.5 +/- 0.1]")), FlintDomainError)
    assert int(B("3.7")) == 3
    assert str(B(R(1) / 3)) == "0.3333333333" and str(R(B(1) / 3)) == "0.3333333333"
    assert str(X(B("0.125"))) == "0.125" and raises(lambda: X(B(1) / 3), FlintUnableError)
    assert str(RF(B(1) / 3)) == "0.3333333333000000" and str(CC(B(1) / 3)) == "[0.3333333333 +/- 3.34e-11]"
    assert str(RealFloat_decfloat(5)(B(1) / 3)) == "0.33333"
    assert str(RealField_decball(5)(B(1) / 3)) == "[0.33333 +/- 3.335e-6]"
    assert str(RealField_decball(5)(R(1) / 3)) == "[0.33333 +/- 3.334e-6]"
    assert str(RealField_decball(5, limb_digits=2)(B(1) / 3)) == "[0.33333 +/- 3.341e-6]"
    assert str(RealField_decball(5, limb_digits=2, rad_prec=2)(B(1) / 3)) == "[0.33333 +/- 3.5e-6]"
    assert str(RealField_decball(12, limb_digits=7)(B(1) / 3)) == "[0.3333333333 +/- 3.334e-11]"
    assert str(B(1) + 1) == "2" and str(1 + B(1)) == "2" and str(B(1) + QQ(1) / 2) == "1.5" and str(B(1) + R("0.5")) == "1.5"
    assert str(B(1) / 3 + R(1) / 3) == "[0.6666666666 +/- 3.334e-11]" and str(R(1) / 3 + B(1) / 3) == "[0.6666666666 +/- 3.334e-11]"
    # functions
    assert str(B.pi()) == "[3.141592654 +/- 4.104e-10]"
    assert str(RealField_decball(10, sloppy_radius=False, rad_prec=9).pi()) == "[3.141592654 +/- 4.10206763e-10]"
    assert str(RealField_decball(30).pi()) == "[3.14159265358979323846264338328 +/- 4.973e-31]"
    assert RealField_decball(10, rad_prec=9, sloppy_radius=False).pi().contains(RealField_decball(50).pi())
    assert str(B(1).exp()) == "[2.718281828 +/- 4.592e-10]" and str(B(0).exp()) == "1" and str(B(1).log()) == "0"
    assert str(B(2).log()) == "[0.6931471806 +/- 4.007e-11]" and str(B(1).sin()) == "[0.8414709848 +/- 7.898e-12]"
    assert str(B(0).sin()) == "0" and str(B(0).cos()) == "1" and str(B(0).atan()) == "0" and str(B(0).sinh()) == "0"
    assert str(B("[0 +/- 1e-20]").sin()) == "[0 +/- 1.001e-20]" and str(B("[0 +/- 1e-20]").exp()) == "[1 +/- 1.001e-20]"
    assert str(B("[1e30 +/- 1]").sin()) == "[0 +/- 1]" and str(B("1e30").sin()) == "[0 +/- 1]"
    assert str(RealField_decball(40)("1e30").sin()) == "[-0.0901169019121380580303864289529873302744 +/- 3.669e-42]"
    assert str(B(100).exp()) == "[2.688117142e43 +/- 1.84e33]" and str(B(-100).exp()) == "[3.720075976e-44 +/- 2.085e-55]"
    assert str(B("[1 +/- 1e-5]").exp()) == "[2.718281828 +/- 2.72e-5]" or B("[1 +/- 1e-5]").exp().contains(B(1).exp())
    assert B("[1 +/- 1e-5]").exp().contains(B(1 + QQ(1) / 10 ** 5).exp())
    assert B(3).gamma().contains(2) and B(5).gamma().contains(24) and B(QQ(1) / 2).gamma().overlaps(B.pi().sqrt())
    assert B(2).zeta().overlaps(B.pi() ** 2 / 6) and B(-1).zeta().contains(QQ(-1) / 12) and B(0).zeta().contains(QQ(-1) / 2)
    assert str(B(1).tan()) == "[1.557407725 +/- 3.452e-10]" and str(B(1).tanh()) == "[0.761594156 +/- 4.425e-11]"
    assert str(B(1).cosh()) == "[1.543080635 +/- 1.849e-10]" and str(B(1).sinh()) == "[1.175201194 +/- 3.563e-10]"
    assert str(B("1e-20").expm1()) == "[1e-20 +/- 3.339e-39]" and str(B("1e-20").log1p()) == "[1e-20 +/- 4.802e-39]"
    assert str(B(2) ** B("0.5")) == "[1.414213562 +/- 3.732e-10]" and str(B(2) ** 10) == "1024" and str(B(2) ** -1) == "0.5"
    assert str(B(2) ** (QQ(1) / 2)) == "[1.414213562 +/- 3.731e-10]" and str(B(4) ** (QQ(1) / 2)) == "2"
    assert str(B(9) ** (B(1) / 2)) == "3" and raises(lambda: B(-8) ** (QQ(1) / 3), FlintDomainError) and str(B(-8) ** 3) == "-512"
    assert (B(2) ** B(3)).contains(8) and str(B(2) ** B(3)) == "8"
    assert raises(lambda: B(0).log(), FlintDomainError) and raises(lambda: B("[0 +/- 1]").log(), (FlintDomainError, FlintUnableError))
    assert raises(lambda: B(-1).log(), FlintDomainError) and raises(lambda: B("[1 +/- 2]").log(), (FlintDomainError, FlintUnableError))
    assert raises(lambda: B(-2) ** B("0.5"), FlintDomainError)
    # polynomials, matrices, series over balls
    P = PolynomialRing(B)
    f = P("x^2 - [3 +/- 0.001]*x + 1.25")
    assert str(f) == "1.25 + [-3 +/- 0.001]*x + x^2"
    assert str(P(str(f))) == str(f)
    assert str(P(str(f * f))) == str(f * f)
    assert str(P("(x - 1/3)^3")) == "[-0.03703703703 +/- 1.523e-11] + [0.3333333333 +/- 1.001e-10]*x + [-0.9999999999 +/- 1.001e-10]*x^2 + x^3"
    assert str(P("(x - 1/3)^3")(B(1) / 3)) == "[0 +/- 8.238e-11]"
    assert str(P("x^2 - 2") % P("x - [1.5 +/- 0.1]")) == "[0.25 +/- 0.31]"
    assert str(P("x^2 - 2")(B(2).sqrt())) == "[-1e-9 +/- 1.113e-9]"
    assert P("x^2 - 2")(B(2).sqrt()).contains(0)
    assert P("x^3 - 2")(B(2) ** (QQ(1) / 3)).contains(0)
    assert str(P("x^2 - 2").derivative()) == "2*x"
    M = Mat(B, 2, 2)([[B(1) / 3, 2], [3, 4]])
    assert str(M) == "[[[0.3333333333 +/- 3.334e-11], 2],\n[3, 4]]"
    assert str(Mat(B, 2, 2)(str(M))) == str(M)
    assert str(M.det()) == "[-4.666666667 +/- 3.334e-10]" and M.det().contains(QQ(-14) / 3)
    assert (M * M.inv() - 1)[1, 1].contains(0) and (M * M.inv())[0, 1].contains(0)
    assert Mat(B, 3, 3)().hilbert().det().contains(QQ(1) / 2160)
    assert Mat(B, 4, 4)().hilbert().det().contains(QQ(1) / 6048000)
    assert Mat(RealField_decball(40), 6, 6)().hilbert().det().contains(QQ(1) / 186313420339200000)
    assert str(Mat(B, 2, 2)("[[1 +/- 0.1, 2], [3, 4 +/- 0.1]]")) == "[[[1 +/- 0.1], 2],\n[3, [4 +/- 0.1]]]"
    assert str(Mat(B, 2, 2)("[[[1 +/- 0.1], 2], [3, [4 +/- 0.1]]]")) == "[[[1 +/- 0.1], 2],\n[3, [4 +/- 0.1]]]"
    assert str(Mat(B, 2, 2)("[[1 +/- 0.1, 2], [3, 4 +/- 0.1]]").det()) == "[-2 +/- 0.51]"
    assert str(Mat(B, 1, 1)("[[1/3]]")) == "[[[0.3333333333 +/- 3.334e-11]]]"
    assert str(Mat(B, 1, 1)("[[[1/3]]]")) == "[[[0.3333333333 +/- 3.334e-11]]]"
    assert str(Mat(B, 1, 1)("[[[1/3 +/- 1e-20]]]")) == "[[[0.3333333333 +/- 3.335e-11]]]"
    S = PowerSeriesRing(B, 4)
    assert str(S("1 + x").exp()) == "[2.718281828 +/- 4.592e-10] + [2.718281828 +/- 4.592e-10]*x + [1.359140914 +/- 2.296e-10]*x^2 + [0.4530469713 +/- 1.099e-10]*x^3 + O(x^4)"
    assert str(S("1 + x").log()) == "x - 0.5*x^2 + [0.3333333333 +/- 3.334e-11]*x^3 + O(x^4)"
    ef = S("[1 +/- 1e-20] + x").exp()
    assert str(S(str(ef))) == str(ef)
    assert str(S("1/(1 - x)")) == "1 + x + x^2 + x^3 + O(x^4)"
    assert str(S("1/(1 - x/3)")) == "1 + [0.3333333333 +/- 3.334e-11]*x + [0.1111111111 +/- 3.337e-11]*x^2 + [0.03703703703 +/- 1.523e-11]*x^3 + O(x^4)"
    PP = PolynomialRing(PolynomialRing(B, "x"), "y")
    g = PP("(x + [1 +/- 0.1]*y)^2 + [1/3]")
    assert str(PP(str(g))) == str(g)
    assert PP(str(g)) - g != 1
    # random string roundtrips through polynomials
    for i in range(30):
        h = P(random=True)
        assert str(P(str(h))) == str(h), str(h)
        h = P("1 + x") ** 3 * P(random=True) + P(random=True)
        assert str(P(str(h))) == str(h), str(h)
    for i in range(30):
        h = PolynomialRing(R)(random=True)
        assert str(PolynomialRing(R)(str(h))) == str(h), str(h)
    for i in range(30):
        Bx = RealField_decball(random.randint(9, 30), rad_prec=random.randint(1, 9), sloppy_radius=random.randint(0, 1))
        Px = PolynomialRing(Bx)
        h = Px(random=True) * Px(random=True) + Px(random=True)
        assert str(Px(str(h))) == str(h), str(h)
        Mx = Mat(Bx, 2, 2)(random=True) * Mat(Bx, 2, 2)(random=True)
        assert str(Mat(Bx, 2, 2)(str(Mx))) == str(Mx), str(Mx)
        Rx = RealFloat_decfloat(random.randint(1, 30), rnd=random.choice(_decimal_rnd_names))
        Mx = Mat(Rx, 2, 2)(random=True) * Mat(Rx, 2, 2)(random=True)
        assert str(Mat(Rx, 2, 2)(str(Mx))) == str(Mx), str(Mx)


def test_big_o():
    Q = Qp_padic_radix(7, rel_prec=3)
    assert raises(lambda: Q("O(5^3)"), FlintUnableError)
    assert str(Q("O(7^-3)")) == "0 + O(7^-3)"
    # should saturate the precision limit
    assert str(Q("O(7^100000000000000000000000000000000000000)")) == str(Q(7) ** (10**100))
    assert raises(lambda: Q("O(7^-100000000000000000000000000000000000000)"), FlintUnableError)
    Qx = PolynomialRing(Q, "x")
    Qxy = PowerSeriesRing(Qx, 2, "y")
    Qxyz = PolynomialRing(Qxy, "z")
    Qxyzt = PowerSeriesRing(Qxyz, 3, "t")
    f = Qxyzt("(3 + O(7^2))*t + (2 + O(y^2))*t + (5 + O(7^1))*t^2")
    assert str(f) == "((5 + O(7^2)) + O(y^2))*t + (5 + O(7^1))*t^2"
    for e in range(8):
        f = Qxyzt("3+x+y+z+t")**e
        fs = str(f)
        assert str(Qxyzt(fs)) == fs
    assert raises(lambda: Q("O(5^3)"), FlintUnableError)
    assert raises(lambda: Qxy("O(y^-2)"), FlintUnableError)
    assert str(Qxy("O(y^0)")) == "0 + O(y^0)"
    assert str(Qxy("O(y)")) == "0 + O(y^1)"

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
