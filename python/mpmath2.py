# just for fun

from mpmath.ctx_base import StandardBaseContext
import flint2

import mpmath.libmp
from mpmath.libmp import prec_to_dps, dps_to_prec

class RFContext(StandardBaseContext):

    def __init__(ctx, RR, CC, QQ, ZZ):
        ctx.R = RR
        ctx.RR = RR
        ctx.CC = CC
        ctx.QQ = QQ
        ctx.ZZ = ZZ
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
        return ctx.R.pi()

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

    def ldexp(ctx, x, n):
       return ctx.R.mul_2exp(x, n)

    def sqrt(ctx, x, prec=0):
       return ctx.R.sqrt(x)

    def log(ctx, x, prec=0):
       return ctx.R.log(x)

    ln = log

    # todo
    def nthroot(ctx, x, n):
        return ctx.R(x) ** (1 / ctx.R(n))

    def exp(ctx, x, prec=0):
        return ctx.R.exp(x)

    def sin(ctx, x, prec=0):
        return ctx.R.sin(x)

    def atan(ctx, x, prec=0):
        return ctx.R.atan(x)

    def gamma(ctx, x, prec=0):
        return ctx.R.gamma(x)

    def factorial(ctx, x, prec=0):
        return ctx.R.fac(x)

    fac = factorial

    def bernoulli(ctx, n):
        return ctx.R.bernoulli(n)

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


mp = RFContext(flint2.RF, flint2.CF, flint2.QQ, flint2.ZZ)

sqrt = mp.sqrt
log = mp.log
ln = mp.ln
sin = mp.sin
exp = mp.exp
atan = mp.atan
gamma = mp.gamma
zeta = mp.zeta
fac = mp.fac
root = mp.root

chop = mp.chop
inf = mp.inf
quad = mp.quad
quadosc = mp.quadosc
diff = mp.diff
nsum = mp.nsum
nprod = mp.nprod
limit = mp.limit

"""
from mpmath2 import *
mp.dps = 50
nsum(lambda n: 1/n**2, [1, inf])
nsum(lambda n: 1/n**2, [1, inf], method='e')
diff(lambda x: gamma(x), 1)
quadosc(lambda x: sin(x)/x, [0, inf], omega=1)

from mpmath2 import *
mp.dps = 50; mp.pretty = True

+pi
180*degree
4*atan(1)
16*acot(5)-4*acot(239)
48*acot(49)+128*acot(57)-20*acot(239)+48*acot(110443)
chop(2*j*log((1-j)/(1+j)))
chop(-2j*asinh(1j))
chop(ci(-inf)/1j)
gamma(0.5)**2
beta(0.5,0.5)
(2/diff(erf, 0))**2
findroot(sin, 3)
findroot(cos, 1)*2
chop(-2j*lambertw(-pi/2))
besseljzero(0.5,1)
3*sqrt(3)/2/hyp2f1((-1,3),(1,3),1,1)
8/(hyp2f1(0.5,0.5,1,0.5)*gamma(0.75)/gamma(1.25))**2
4*(hyp1f2(1,1.5,1,1) / struvel(-0.5, 2))**2
1/meijerg([[],[]], [[0],[0.5]], 0)**2
(meijerg([[],[2]], [[1,1.5],[]], 1, 0.5) / erfc(1))**2
(1-e) / meijerg([[1],[0.5]], [[1],[0.5,0]], 1)
sqrt(psi(1,0.25)-8*catalan)
elliprc(1,2)*4
elliprg(0,1,1)*4
2*agm(1,0.5)*ellipk(0.75)
(gamma(0.75)*jtheta(3,0,exp(-pi)))**4
cbrt(gamma(0.25)**4*agm(1,sqrt(2))**2/8)
sqrt(6*zeta(2))
sqrt(6*(zeta(2,3)+5./4))
sqrt(zeta(2,(3,4))+8*catalan)
exp(-2*zeta(0,1,1))/2
sqrt(12*altzeta(2))
4*dirichlet(1,[0,1,0,-1])
2*catalan/dirichlet(-1,[0,1,0,-1],1)
exp(-dirichlet(0,[0,1,0,-1],1))*gamma(0.25)**2/(2*sqrt(2))
sqrt(7*zeta(3)/(4*diff(lerchphi, (-1,-2,1), (0,1,0))))
sqrt(-12*polylog(2,-1))
sqrt(6*log(2)**2+12*polylog(2,0.5))
chop(root(-81j*(polylog(3,root(1,3,1))+4*zeta(3)/9)/2,3))
2*clsin(1,1)+1
(3+sqrt(3)*sqrt(1+8*clcos(2,1)))/2
root(2,6)*sqrt(e)/(glaisher**6*barnesg(0.5)**4)
nsum(lambda k: 4*(-1)**(k+1)/(2*k-1), [1,inf])
nsum(lambda k: (3**k-1)/4**k*zeta(k+1), [1,inf])
nsum(lambda k: 8/(2*k-1)**2, [1,inf])**0.5
nsum(lambda k: 2*fac(k)/fac2(2*k+1), [0,inf])
nsum(lambda k: fac(k)**2/fac(2*k+1), [0,inf])*3*sqrt(3)/2
nsum(lambda k: fac(k)**2/(phi**(2*k+1)*fac(2*k+1)), [0,inf])*(5*sqrt(phi+2))/2
nsum(lambda k: (4/(8*k+1)-2/(8*k+4)-1/(8*k+5)-1/(8*k+6))/16**k, [0,inf])
2/nsum(lambda k: (-1)**k*(4*k+1)*(fac2(2*k-1)/fac2(2*k))**3, [0,inf])
nsum(lambda k: 72/(k*expm1(k*pi))-96/(k*expm1(2*pi*k))+24/(k*expm1(4*pi*k)), [1,inf])
1/nsum(lambda k: binomial(2*k,k)**3*(42*k+5)/2**(12*k+4), [0,inf])
4/nsum(lambda k: (-1)**k*(1123+21460*k)*fac2(2*k-1)*fac2(4*k-1)/(882**(2*k+1)*32**k*fac(k)**3), [0,inf])
9801/sqrt(8)/nsum(lambda k: fac(4*k)*(1103+26390*k)/(fac(k)**4*396**(4*k)), [0,inf])
426880*sqrt(10005)/nsum(lambda k: (-1)**k*fac(6*k)*(13591409+545140134*k)/(fac(k)**3*fac(3*k)*(640320**3)**k), [0,inf])
4/nsum(lambda k: (6*k+1)*rf(0.5,k)**3/(4**k*fac(k)**3), [0,inf])
(ln(8)+sqrt(48*nsum(lambda m,n: (-1)**(m+n)/(m**2+n**2), [1,inf],[1,inf]) + 9*log(2)**2))/2
-nsum(lambda x,y: (-1)**(x+y)/(x**2+y**2), [-inf,inf], [-inf,inf], ignore=True)/ln2
2*nsum(lambda k: sin(k)/k, [1,inf])+1
quad(lambda x: 2/(x**2+1), [0,inf])
quad(lambda x: exp(-x**2), [-inf,inf])**2
2*quad(lambda x: sqrt(1-x**2), [-1,1])
chop(quad(lambda z: 1/(2j*z), [1,j,-1,-j,1]))
3*(4*log(2+sqrt(3))-quad(lambda x,y: 1/sqrt(1+x**2+y**2), [-1,1],[-1,1]))/2
sqrt(8*quad(lambda x,y: 1/(1-(x*y)**2), [0,1],[0,1]))
sqrt(6*quad(lambda x,y: 1/(1-x*y), [0,1],[0,1]))
sqrt(6*quad(lambda x: x/expm1(x), [0,inf]))
quad(lambda x: (16*x-16)/(x**4-2*x**3+4*x-4), [0,1])
quad(lambda x: sqrt(x-x**2), [0,0.25])*24+3*sqrt(3)/4
mpf(22)/7 - quad(lambda x: x**4*(1-x)**4/(1+x**2), [0,1])
mpf(355)/113 - quad(lambda x: x**8*(1-x)**8*(25+816*x**2)/(1+x**2), [0,1])/3164
2*quadosc(lambda x: sin(x)/x, [0,inf], omega=1)
40*quadosc(lambda x: sin(x)**6/x**6, [0,inf], omega=1)/11
e*quadosc(lambda x: cos(x)/(1+x**2), [-inf,inf], omega=1)
8*quadosc(lambda x: cos(x**2), [0,inf], zeros=lambda n: sqrt(n))**2
2*quadosc(lambda x: sin(exp(x)), [1,inf], zeros=ln)+2*si(e)
exp(2*quad(loggamma, [0,1]))/2
2*nprod(lambda k: sec(pi/2**k), [2,inf])
s=lambda k: sqrt(0.5+s(k-1)/2) if k else 0; 2/nprod(s, [1,inf])
s=lambda k: sqrt(2+s(k-1)) if k else 0; limit(lambda k: sqrt(2-s(k))*2**(k+1), inf)
2*nprod(lambda k: (2*k)**2/((2*k-1)*(2*k+1)), [1,inf])
2*nprod(lambda k: (4*k**2)/(4*k**2-1), [1, inf])
sqrt(6*ln(nprod(lambda k: exp(1/k**2), [1,inf])))
nprod(lambda k: (k**2-1)/(k**2+1), [2,inf])/csch(pi)
nprod(lambda k: (k**2-1)/(k**2+1), [2,inf])*sinh(pi)
nprod(lambda k: (k**4-1)/(k**4+1), [2, inf])*(cosh(sqrt(2)*pi)-cos(sqrt(2)*pi))/sinh(pi)
sinh(pi)/nprod(lambda k: (1-1/k**4), [2, inf])/4
sinh(pi)/nprod(lambda k: (1+1/k**2), [2, inf])/2
(exp(1+euler/2)/nprod(lambda n: (1+1/n)**n * exp(1/(2*n)-1), [1, inf]))**2/2
3*sqrt(2)*cosh(pi*sqrt(3)/2)**2*csch(pi*sqrt(2))/nprod(lambda k: (1+1/k+1/k**2)**2/(1+2/k+3/k**2), [1, inf])
2/e*nprod(lambda k: (1+2/k)**((-1)**(k+1)*k), [1,inf])
limit(lambda k: 16**k/(k*binomial(2*k,k)**2), inf)
limit(lambda x: 4*x*hyp1f2(0.5,1.5,1.5,-x**2), inf)
1/log(limit(lambda n: nprod(lambda k: pi/(2*atan(k)), [n,2*n]), inf),4)
limit(lambda k: 2**(4*k+1)*fac(k)**4/(2*k+1)/fac(2*k)**2, inf)
limit(lambda k: fac(k) / (sqrt(k)*(k/e)**k), inf)**2/2
limit(lambda k: (-(-1)**k*bernoulli(2*k)*2**(2*k-1)/fac(2*k))**(-1/(2*k)), inf)
limit(lambda k: besseljzero(1,k)/k, inf)
1/limit(lambda x: airyai(x)*2*x**0.25*exp(2*x**1.5/3), inf, exp=True)**2
1/limit(lambda x: airybi(x)*x**0.25*exp(-2*x**1.5/3), inf, exp=True)**2

"""
