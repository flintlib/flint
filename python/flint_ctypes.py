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
                ('den', ctypes.c_long),
                ('alloc', ctypes.c_long),
                ('length', ctypes.c_long)]

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

libgr.gr_set_str.argtypes = (ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_get_str.argtypes = (ctypes.POINTER(ctypes.c_char_p), ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmp.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))
libgr.gr_cmpabs.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

libgr.gr_heap_clear.argtypes = (ctypes.c_void_p, ctypes.POINTER(gr_ctx_struct))

class gr_ctx:

    def __init__(self, which, **kwargs):
        self._data = gr_ctx_struct()
        self._ref = ctypes.byref(self._data)
        if which == "ZZ":
            libgr.gr_ctx_init_fmpz(self._ref)
            self._elem_type = fmpz
        elif which == "QQ":
            libgr.gr_ctx_init_fmpq(self._ref)
            self._elem_type = fmpq
        elif which == "QQbar":
            libgr.gr_ctx_init_complex_qqbar(self._ref)
            self._elem_type = qqbar
        self._str = self._repr()

    def _repr(self):
        arr = ctypes.c_char_p()
        if libgr.gr_ctx_get_str(ctypes.byref(arr), self._ref) != GR_SUCCESS:
            raise NotImplementedError
        try:
            return ctypes.cast(arr, ctypes.c_char_p).value.decode("ascii")
        finally:
            libflint.flint_free(arr)

    def __call__(self, value=None):
        return self._elem_type(value, self)

    def __repr__(self):
        return self._str

    def __del__(self):
        # todo: refcounting
        # libgr.gr_ctx_clear(self._ref)
        pass

class gr_elem:

    def __init__(self, val=None, context=None):
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
                #c = libgr.gr_ctx_cmp_coercion(self._ctx, other._ctx)
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
            elif self._ctx_python is not other._ctx_python:
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
    def _binary_predicate(self, other, op, rstr):
        self, other = gr_elem._binary_coercion(self, other)
        truth = op(self._ref, other._ref, self._ctx)
        if truth == T_TRUE: return True
        if truth == T_FALSE: return False
        raise Undecidable(f"{rstr} cannot be decided for x = {self}, y = {other} over {self.parent()}")

    @staticmethod
    def _unary_op(self, op, rstr):
        elem_type = type(self)
        res = elem_type(None, self._ctx_python)
        status = op(res._ref, self._ref, self._ctx)
        if status:
            if status & GR_UNABLE: raise NotImplementedError(f"{rstr} is not implemented for x = {self} over {self.parent()}")
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
        return self._binary_op(self, other, libgr.gr_pow, "x ** y")

    def __rpow__(self, other):
        return self._binary_op(other, self, libgr.gr_pow, "x ** y")


class fmpz(gr_elem):
    _struct_type = fmpz_struct

    def __int__(self):
        return fmpz_to_python_int(self._ref)

class fmpq(gr_elem):
    _struct_type = fmpq_struct

class qqbar(gr_elem):
    _struct_type = qqbar_struct

ZZ = gr_ctx("ZZ")
QQ = gr_ctx("QQ")
QQbar = gr_ctx("QQbar")

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
    assert str(xy) == "Root a = 0.629961 + 1.09112i of a^3+2"
    i = QQbar(-1) ** (QQ(1)/2)
    assert str(i) == 'Root a = 1.00000i of a^2+1'
    assert str(-i) == 'Root a = -1.00000i of a^2+1'
    assert str(1-i) == 'Root a = 1.00000 - 1.00000i of a^2-2*a+2'
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

