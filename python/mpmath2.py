# just for fun

from mpmath.ctx_base import StandardBaseContext
import flint

import mpmath.libmp
from mpmath.libmp import prec_to_dps, dps_to_prec

class RFContext(StandardBaseContext):

    def __init__(ctx, RR, QQ):
        ctx.R = RR
        ctx.RR = RR
        ctx.QQ = QQ
        ctx._type = ctx.R()
        ctx.R.prec = 53
        ctx._dps = prec_to_dps(53)
        StandardBaseContext.__init__(ctx)
        #ctx.pretty = False
        #ctx._init_aliases()

        ctx.zero = ctx.R()
        ctx.one = ctx.R(1)
        ctx.inf = ctx.R("inf")
        ctx.ninf = -ctx.inf
        ctx.nan = ctx.R("nan")

        # .j = ...

    @property
    def pi(ctx):
        return ctx.zero.pi()

    @property
    def eps(ctx):
        return ctx.one * ctx.R(2)**(1 - ctx.prec)

    def mpf(ctx, val):
        return ctx.R(val)

    def mpc(ctx, val):
        return ctx.R(val)

    def _mpq(ctx, x):
        a, b = x
        return ctx.QQ(a) / b

    def convert(ctx, x):
        # todo: avoid copying
        return ctx.R(x)

    NoConvergence = mpmath.libmp.NoConvergence

    def _get_prec(ctx):
        return ctx.R.prec

    def _set_prec(ctx, p):
        p = max(1, int(p))
        ctx.R.prec = p
        ctx._dps = prec_to_dps(p)

    def _get_dps(ctx):
        return ctx._dps

    def _set_dps(ctx, d):
        d = max(1, int(d))
        ctx.R.prec = dps_to_prec(d)
        ctx._dps = d

    _fixed_precision = False

    prec = property(_get_prec, _set_prec)
    dps = property(_get_dps, _set_dps)

    def is_special(ctx, x):
        return x - x != 0.0

    def isnan(ctx, x):
        return x != x

    def isinf(ctx, x):
        return abs(x) == ctx.inf

    def isnormal(ctx, x):
        if x:
            return x - x == 0.0
        return False

    # todo
    def _is_real_type(ctx, x):
        return True

    # todo:
    def ldexp(ctx, x, n):
        return x * ctx.R(2) ** n

    def sqrt(ctx, x, prec=0):
        if not (type(x) is ctx._type and x.parent() is ctx.R):
            x = ctx.R(x)
        return x.sqrt()

    def log(ctx, x, prec=0):
        if not (type(x) is ctx._type and x.parent() is ctx.R):
            x = ctx.R(x)
        return x.log()

    ln = log

    def exp(ctx, x, prec=0):
        if not (type(x) is ctx._type and x.parent() is ctx.R):
            x = ctx.R(x)
        return x.exp()

    def sin(ctx, x, prec=0):
        if not (type(x) is ctx._type and x.parent() is ctx.R):
            x = ctx.R(x)
        return x.sin()

    def gamma(ctx, x, prec=0):
        if not (type(x) is ctx._type and x.parent() is ctx.R):
            x = ctx.R(x)
        return x.gamma()


    '''
    # Called by SpecialFunctions.__init__()
    @classmethod
    def _wrap_specfun(cls, name, f, wrap):
        if wrap:
            def f_wrapped(ctx, *args, **kwargs):
                convert = ctx.convert
                args = [convert(a) for a in args]
                return f(ctx, *args, **kwargs)
        else:
            f_wrapped = f
        f_wrapped.__doc__ = function_docs.__dict__.get(name, f.__doc__)
        setattr(cls, name, f_wrapped)

    def bernoulli(ctx, n):
        cache = ctx._bernoulli_cache
        if n in cache:
            return cache[n]
        cache[n] = to_float(mpf_bernoulli(n, 53, 'n'), strict=True)
        return cache[n]

    pi = math2.pi
    e = math2.e
    euler = math2.euler
    sqrt2 = 1.4142135623730950488
    sqrt5 = 2.2360679774997896964
    phi = 1.6180339887498948482
    ln2 = 0.69314718055994530942
    ln10 = 2.302585092994045684
    euler = 0.57721566490153286061
    catalan = 0.91596559417721901505
    khinchin = 2.6854520010653064453
    apery = 1.2020569031595942854
    glaisher = 1.2824271291006226369

    absmin = absmax = abs


    def isnpint(ctx, x):
        if type(x) is complex:
            if x.imag:
                return False
            x = x.real
        return x <= 0.0 and round(x) == x

    power = staticmethod(math2.pow)
    sqrt = staticmethod(math2.sqrt)
    exp = staticmethod(math2.exp)
    ln = log = staticmethod(math2.log)
    cos = staticmethod(math2.cos)
    sin = staticmethod(math2.sin)
    tan = staticmethod(math2.tan)
    cos_sin = staticmethod(math2.cos_sin)
    acos = staticmethod(math2.acos)
    asin = staticmethod(math2.asin)
    atan = staticmethod(math2.atan)
    cosh = staticmethod(math2.cosh)
    sinh = staticmethod(math2.sinh)
    tanh = staticmethod(math2.tanh)
    gamma = staticmethod(math2.gamma)
    rgamma = staticmethod(math2.rgamma)
    fac = factorial = staticmethod(math2.factorial)
    floor = staticmethod(math2.floor)
    ceil = staticmethod(math2.ceil)
    cospi = staticmethod(math2.cospi)
    sinpi = staticmethod(math2.sinpi)
    cbrt = staticmethod(math2.cbrt)
    _nthroot = staticmethod(math2.nthroot)
    _ei = staticmethod(math2.ei)
    _e1 = staticmethod(math2.e1)
    _zeta = _zeta_int = staticmethod(math2.zeta)

    # XXX: math2
    def arg(ctx, z):
        z = complex(z)
        return math.atan2(z.imag, z.real)

    def expj(ctx, x):
        return ctx.exp(ctx.j*x)

    def expjpi(ctx, x):
        return ctx.exp(ctx.j*ctx.pi*x)

    ldexp = math.ldexp
    frexp = math.frexp

    def mag(ctx, z):
        if z:
            return ctx.frexp(abs(z))[1]
        return ctx.ninf

    def isint(ctx, z):
        if hasattr(z, "imag"):   # float/int don't have .real/.imag in py2.5
            if z.imag:
                return False
            z = z.real
        try:
            return z == int(z)
        except:
            return False

    def nint_distance(ctx, z):
        if hasattr(z, "imag"):   # float/int don't have .real/.imag in py2.5
            n = round(z.real)
        else:
            n = round(z)
        if n == z:
            return n, ctx.ninf
        return n, ctx.mag(abs(z-n))

    def _convert_param(ctx, z):
        if type(z) is tuple:
            p, q = z
            return ctx.mpf(p) / q, 'R'
        if hasattr(z, "imag"):    # float/int don't have .real/.imag in py2.5
            intz = int(z.real)
        else:
            intz = int(z)
        if z == intz:
            return intz, 'Z'
        return z, 'R'

    def _is_real_type(ctx, z):
        return isinstance(z, float) or isinstance(z, int_types)

    def _is_complex_type(ctx, z):
        return isinstance(z, complex)

    def hypsum(ctx, p, q, types, coeffs, z, maxterms=6000, **kwargs):
        coeffs = list(coeffs)
        num = range(p)
        den = range(p,p+q)
        tol = ctx.eps
        s = t = 1.0
        k = 0
        while 1:
            for i in num: t *= (coeffs[i]+k)
            for i in den: t /= (coeffs[i]+k)
            k += 1; t /= k; t *= z; s += t
            if abs(t) < tol:
                return s
            if k > maxterms:
                raise ctx.NoConvergence

    def atan2(ctx, x, y):
        return math.atan2(x, y)

    def psi(ctx, m, z):
        m = int(m)
        if m == 0:
            return ctx.digamma(z)
        return (-1)**(m+1) * ctx.fac(m) * ctx.zeta(m+1, z)

    digamma = staticmethod(math2.digamma)

    def harmonic(ctx, x):
        x = ctx.convert(x)
        if x == 0 or x == 1:
            return x
        return ctx.digamma(x+1) + ctx.euler

    nstr = str

    def to_fixed(ctx, x, prec):
        return int(math.ldexp(x, prec))

    def rand(ctx):
        import random
        return random.random()

    _erf = staticmethod(math2.erf)
    _erfc = staticmethod(math2.erfc)

    def sum_accurately(ctx, terms, check_step=1):
        s = ctx.zero
        k = 0
        for term in terms():
            s += term
            if (not k % check_step) and term:
                if abs(term) <= 1e-18*abs(s):
                    break
            k += 1
        return s
    '''


mp = RFContext(flint.RF, flint.QQ)

sqrt = mp.sqrt
log = mp.log
ln = mp.ln
sin = mp.sin
exp = mp.exp
gamma = mp.gamma

inf = mp.inf
quad = mp.quad
quadosc = mp.quadosc
diff = mp.diff
nsum = mp.nsum

"""
from mpmath2 import *
mp.dps = 100
nsum(lambda n: 1/n**2, [1, inf])
diff(lambda x: gamma(x), 1)
quadosc(lambda x: sin(x)/x, [0, inf], omega=1)

"""
