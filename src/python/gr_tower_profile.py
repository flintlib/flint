"""
Profiling and correctness catalog for exact real and complex numbers:
the lazy tower field (ComplexField_tower) against Calcium
(ComplexField_ca), qqbar (ComplexAlgebraicField_qqbar, on the algebraic
cases) and SymPy.

    python3 gr_tower_profile.py                       # all cases, tower and ca
    python3 gr_tower_profile.py nm.                   # cases whose id contains "nm."
    python3 gr_tower_profile.py --category special    # one category
    python3 gr_tower_profile.py --engines tower,ca,qqbar,sympy --timeout 20
    python3 gr_tower_profile.py --json out.json --html report.html
    python3 gr_tower_profile.py --from-json out.json --html report.html
    python3 gr_tower_profile.py --list                # list the cases

Each case runs once per engine, in a fresh field, in a subprocess with a
time limit (--jobs N runs N subprocesses at a time). The result of a run
is one of

    ok          the engine gives the right answer
    wrong       the engine gives the wrong answer (never acceptable: a near
                miss decided zero, a zero decided nonzero, a false inequality)
    undecided   the engine cannot decide (Undecidable; SymPy: None)
    unable      the engine cannot compute an operation (FlintUnableError)
    error       any other exception, or a crash
    timeout     the time limit was exceeded
    n/a         the engine does not apply (qqbar: a transcendental case;
                SymPy: a program case without a SymPy version)

with the wall time in seconds. A text table is printed as the cases
finish; --json writes the results and --html a report (formulas rendered
with KaTeX from fexpr LaTeX, summary statistics per engine and category,
speed ratios, filters).

Cases are of two kinds. Expression cases are tuples (id, category,
expect, expression, approx, source, what it stresses). Expressions use
SymPy/Python syntax: `**` for powers, `1/3` the rational number, `I` the
imaginary unit, principal branches throughout; real_root(a, n) is the real
n-th root, RootOf(p, k) the k-th root of p in SymPy's CRootOf order (real
roots ascending, then the others by real part, then imaginary part; the
roots are sorted by exact comparisons in the field), zeta(s, a) the Hurwitz
zeta function, Eq/Ne/Lt/Le/Gt/Ge relations. They are evaluated in the
field through Python's ast, operation by operation (expand() is ignored:
the value is the same), and parsed by SymPy for the SymPy engine. The
expectations:

    zero        the expression is zero
    nonzero     the expression is not zero (approx: its value)
    minpoly:P   P(x) = 0 for the expression x
    true        the relation holds
    false       the relation does not hold

Program cases (PROGRAMS) are Python functions of the field returning True
for the right answer: matrices (Cayley-Hamilton, a nullspace from Sage
trac #37927, exp(log(M))) and polynomial roots (Sage's QQbar examples, and
roots of polynomials with transcendental-algebraic coefficients). All
comparisons, including the selection of roots, are done in the exact field.

The families (gauss_sum_cases, sine_product_cases, dft_cases,
swinnerton_dyer_cases) scale: Gauss sums sum (k|p) cos(2 pi k/p) = sqrt(p),
products of sines prod 2 sin(k pi/n) = n, the exact DFT benchmark
x = IDFT(DFT(x)) for six kinds of inputs, and Swinnerton-Dyer polynomials
S_n(sqrt 2 + ... + sqrt p_n) = 0.

The big cases (left out with --no-big): the expressions of
examples/huge_expr.c (read from that file), and cos(pi/257) by Gauss's
construction (gauss_period_cos: written out, 643 kB with 15712 square
roots, generated outside the timings; and with the periods bound by
name, 127 square roots).

The prefix of an id says where the case comes from: ca (FLINT/Calcium
examples and documentation), sp (SymPy tests and issues), mp (SymPy
minimal polynomials), dn (denesting: SymPy's sqrtdenest, Ramanujan,
Shanks, Cavallo, Landau), w (Wester's problems, via SymPy's
test_wester.py), tr (trigonometric constants), el (elementary
functions), nm (near misses), sage (Sage's QQbar), pf (the Calcium issue
tracker, the Calcium notebook and Sage's QQbar documentation), big
(stress), fam (families), prog (programs); and by topic: fl (integer
parts), cmp (comparisons), pt (real and imaginary parts, absolute values),
sg (signs), br (branch cuts), rich (Richardson-type exp-log problems,
algebraic arguments), mix (algebraic and transcendental extensions
together), mac (Machin-like formulas in non-Gaussian number fields), sf
(special function identities), asy (inequalities in asymptotic regimes),
rad (Cardano and Ferrari root formulas of cubics and quartics with
rational, algebraic and transcendental coefficients: in the polynomial,
Vieta's formulas, against closed forms; prog.rad_* against root finding),
nt (algebraic number theory: class number formulas as quotients of
sines, cubic Gauss and Jacobi sums, cubic period polynomials, quadratic
Gauss sums, norms of cyclotomic units), ra (real algebraic geometry:
Mignotte root separation), fs (translation surfaces: the sage-flatsurf
17-gon), prog.mat_* (matrix exponentials, logarithms, roots and powers
through the Jordan form), prog.veech* (Veech group relations).
Programs may carry a SymPy version (root finding with CRootOf or roots,
Matrix.exp, log and powers).
The C test src/gr_tower/test/t-catalog.c runs most of the zero tests.
"""

import sys, os, time, ast, math, json, html, functools, threading
import multiprocessing as mp
from fractions import Fraction

# ---------------------------------------------------------------------------
# expression cases
# ---------------------------------------------------------------------------

CASES = [
    ('ca.euler', 'exp-log', 'zero', 'exp(pi*I) + 1', None,
     'FLINT examples/elementary.c', 'exp at i*pi*Q is algebraic'),
    ('ca.log_m1', 'exp-log', 'zero', 'log(-1)/(pi*I) - 1', None,
     'FLINT examples/elementary.c', 'principal log'),
    ('ca.log_mi', 'exp-log', 'zero', 'log(-I)/(pi*I) + 1/2', None,
     'FLINT examples/elementary.c', 'principal log, negative argument'),
    ('ca.log_pow10', 'exp-log', 'zero', 'log(1/10**123)/log(100) + 123/2', None,
     'FLINT examples/elementary.c', 'log of integers -> prime logs'),
    ('ca.log_unit', 'exp-log', 'zero', 'log(1 + sqrt(2))/log(3 + 2*sqrt(2)) - 1/2', None,
     'FLINT examples/elementary.c', 'multiplicative relation between units (1+sqrt2)^2 = 3+2sqrt2'),
    ('ca.sqrt6', 'algebraic', 'zero', 'sqrt(2)*sqrt(3) - sqrt(6)', None,
     'FLINT examples/elementary.c', 'relations between radical generators'),
    ('ca.exp_sum', 'exp-log', 'zero', 'exp(1 + sqrt(2))*exp(1 - sqrt(2))/exp(1)**2 - 1', None,
     'FLINT examples/elementary.c', 'additive relations among exp arguments'),
    ('ca.i_to_i', 'exp-log', 'zero', 'I**I - exp(-pi/2)', None,
     'FLINT examples/elementary.c; SymPy evalf docs (nsimplify(I**I,[pi]))', 'complex power via exp(I*log(I))'),
    ('ca.exp_sqrt12', 'exp-log', 'zero', 'exp(sqrt(3))**2 - exp(sqrt(12))', None,
     'FLINT examples/elementary.c', 'exp arguments Q-dependent after algebraic reduction'),
    ('ca.log_pi_i', 'exp-log', 'zero', '2*log(pi*I) - 4*log(sqrt(pi)) - pi*I', None,
     'FLINT examples/elementary.c', 'log of transcendental*I; branch bookkeeping'),
    ('ca.bbk_ex1', 'exp-log', 'zero', '-I*pi/8*log(2/3 - 2*I/3)**2 + I*pi/8*log(2/3 + 2*I/3)**2 + pi**2/12*log(-1 - I) + pi**2/12*log(-1 + I) + pi**2/12*log(1/3 - I/3) + pi**2/12*log(1/3 + I/3) + pi**2/48*log(18)', None,
     "FLINT examples/elementary.c (NB: the printed label there says '- Pi^2/48*Log(18)', but the code computes '+'; only '+' is an identity); Bailey-Borwein-Kaiser, JSC 60 (2014), Example 1", 'integer relations among logs of Gaussian rationals, including squared logs'),
    ('ca.denest5_2_6', 'denest', 'zero', 'sqrt(5 + 2*sqrt(6)) - sqrt(2) - sqrt(3)', None,
     'FLINT examples/elementary.c', 'classic denesting'),
    ('ca.sqrt_i', 'branch', 'zero', 'sqrt(I) - (1 + I)/sqrt(2)', None,
     'FLINT examples/elementary.c', 'complex principal sqrt'),
    ('ca.erf_logs', 'special', 'zero', 'erf(2*log(sqrt(1/2 - sqrt(2)/4)) + log(4)) - erf(log(2 - sqrt(2)))', None,
     'FLINT examples/elementary.c (from Mathematica PossibleZeroQ docs)', 'must simplify the log arguments under an opaque function'),
    ('ca.iter_pow', 'branch', 'zero', 'I**I - exp(pi/((sqrt(-2)**sqrt(2))**sqrt(2)))', None,
     'Calcium introduction notebook', 'iterated complex powers; (a^b)^c != a^(bc)'),
    ('ca.trig_tr1', 'trig', 'zero', 'sin(sqrt(2)/2)**2 + cos(1/sqrt(2))**2 - 1', None,
     'Calcium docs/blog', 'Pythagoras with transcendental argument in disguise'),
    ('ca.trig_tr2', 'trig', 'zero', 'sin(3 + pi) + sin(3)', None,
     'Calcium docs/blog', 'shift by pi'),
    ('ca.trig_tr3', 'trig', 'zero', 'tan(1 + pi) - tan(1)', None,
     'Calcium docs/blog', 'period of tan'),
    ('ca.gd', 'trig', 'zero', 'sin(2*atan(tanh(1/2))) - tanh(1)', None,
     'Calcium docs (sin(gd(1)) = tanh(1))', 'gudermannian: mixes circular and hyperbolic'),
    ('ca.atan_alg', 'trig', 'zero', 'atan(1 - sqrt(2)) + pi/8', None,
     'Calcium dev update 2021', 'atan of algebraic giving rational multiple of pi'),
    ('ca.atan_tan', 'branch', 'zero', 'atan(tan(23*pi/27)) + 4*pi/27', None,
     'Calcium dev update 2021', 'branch reduction of atan(tan)'),
    ('ca.asin_sin', 'branch', 'zero', 'asin(sin(sqrt(2) - 1)) - (sqrt(2) - 1)', None,
     'Calcium dev update 2021', 'asin(sin(x)) = x only on principal strip'),
    ('ca.gamma_fe', 'special', 'zero', 'gamma(pi + 1)/gamma(pi) - pi', None,
     'Johansson, NUS 2021 slides', 'functional equation'),
    ('ca.erf_erfc', 'special', 'zero', 'erf(exp(pi*I/3)) - erfc(exp(-2*pi*I/3)) + 1', None,
     'Johansson, NUS 2021 slides', 'erf(-z) = -erf(z) hidden behind root-of-unity argument'),
    ('ca.mixed1', 'denest', 'zero', '(pi + sqrt(2) + sqrt(3))/(pi + sqrt(5 + 2*sqrt(6))) - 1', None,
     'Calcium dev update 2021', 'denesting inside a transcendental field'),
    ('ca.mixed2', 'exp-log', 'zero', 'log(1/exp(sqrt(2) + 1)) + sqrt(2) + 1', None,
     'Calcium dev update 2021', 'log(exp(x)) = x only if |Im x| < pi'),
    ('ca.mixed3', 'branch', 'zero', 'arg(sqrt(-pi*I)) + pi/4', None,
     'Calcium dev update 2021', 'arg of principal sqrt'),
    ('ca.mixed4', 'trig', 'zero', 'sin((1 + sqrt(2))/2) - sqrt((1 - cos(1 + sqrt(2)))/2)', None,
     'Calcium dev update 2021', 'half-angle with sign; needs sin((1+sqrt2)/2) > 0'),
    ('ca.gosper', 'algebraic', 'zero', 'sqrt(36 + 3*(-54 + 35*I*sqrt(3))**(1/3)*3**(1/3) + 117/(-162 + 105*I*sqrt(3))**(1/3))/3 + sqrt(5)*(1296*I + 840*sqrt(3) - 35*3**(5/6)*(-54 + 35*I*sqrt(3))**(1/3) - 54*I*(-162 + 105*I*sqrt(3))**(1/3) + 13*I*(-162 + 105*I*sqrt(3))**(2/3))/(5*(162*I + 105*sqrt(3))) - sqrt(5) - sqrt(7)', None,
     'Calcium introduction notebook (Gosper)', 'casus irreducibilis; complex cube roots must cancel'),
    ('ca.ramanujan_163', 'near-miss', 'nonzero', 'exp(pi*sqrt(163)) - (640320**3 + 744)', '-7.4993e-13',
     'FLINT examples/elementary.c; SymPy test_evalf.py:test_evalf_ramanujan', "Ramanujan's constant"),
    ('ca.exp_tiny', 'near-miss', 'nonzero', 'exp(1/10**10000) - 1', '1.0e-10000',
     "Johansson, 'Computing with metavalues' (blog, 2021)", 'Lindemann says nonzero; Calcium returns Unknown. Structural proof vs ~33000-bit evaluation'),
    ('ca.machin', 'machin', 'zero', '4*atan(1/5) - atan(1/239) - pi/4', None,
     'FLINT examples/machin.c', 'atan -> log of Gaussian integers -> integer relation'),
    ('ca.machin2', 'machin', 'zero', 'atan(1/2) + atan(1/3) - pi/4', None,
     'FLINT examples/machin.c', 'Euler'),
    ('ca.machin3', 'machin', 'zero', '2*atan(1/2) - atan(1/7) - pi/4', None,
     'FLINT examples/machin.c', 'Hermann'),
    ('ca.machin4', 'machin', 'zero', '2*atan(1/3) + atan(1/7) - pi/4', None,
     'FLINT examples/machin.c', 'Hutton'),
    ('ca.machin5', 'machin', 'zero', 'atan(1/2) + atan(1/5) + atan(1/8) - pi/4', None,
     'FLINT examples/machin.c', '3 terms'),
    ('ca.machin6', 'machin', 'zero', 'atan(1/3) + atan(1/4) + atan(1/7) + atan(1/13) - pi/4', None,
     'FLINT examples/machin.c', '4 terms'),
    ('ca.machin7', 'machin', 'zero', '12*atan(1/49) + 32*atan(1/57) - 5*atan(1/239) + 12*atan(1/110443) - pi/4', None,
     'FLINT examples/machin.c (Takano)', 'large arguments'),
    ('ca.hmachin2', 'machin', 'zero', '14*atanh(1/31) + 10*atanh(1/49) + 6*atanh(1/161) - log(2)', None,
     'FLINT examples/machin.c', 'hyperbolic Machin'),
    ('ca.hmachin3', 'machin', 'zero', '22*atanh(1/31) + 16*atanh(1/49) + 10*atanh(1/161) - log(3)', None,
     'FLINT examples/machin.c', 'hyperbolic Machin'),
    ('ca.hmachin5', 'machin', 'zero', '32*atanh(1/31) + 24*atanh(1/49) + 14*atanh(1/161) - log(5)', None,
     'FLINT examples/machin.c', 'hyperbolic Machin'),
    ('ca.hmachin7', 'machin', 'zero', '404*atanh(1/251) + 152*atanh(1/449) - 106*atanh(1/4801) + 174*atanh(1/8749) - log(7)', None,
     'FLINT examples/machin.c', '4-term hyperbolic set; larger sets in FLINT src/mp_real/machin_tab.c (n up to 48)'),
    ('ca.hmachin2b', 'machin', 'zero', '144*atanh(1/251) + 54*atanh(1/449) - 38*atanh(1/4801) + 62*atanh(1/8749) - log(2)', None,
     'FLINT examples/machin.c', '4-term hyperbolic set'),
    ('sp.equals1', 'algebraic', 'zero', '-3 - sqrt(5) + (-sqrt(10)/2 - sqrt(2)/2)**2', None,
     'sympy/core/tests/test_expr.py:test_equals; test_assumptions.py:test_special_assumptions', 'simplify(e<0) and simplify(e>0) both False; needs exact zero'),
    ('sp.equals2', 'branch', 'zero', '-(-1)**(3/4)*6**(1/4) + (-6)**(1/4)*I', None,
     'sympy/core/tests/test_expr.py:test_equals', 'principal roots of negative numbers'),
    ('sp.equals3', 'denest', 'zero', 'sqrt(1 + sqrt(3)) + sqrt(3 + 3*sqrt(3)) - sqrt(10 + 6*sqrt(3))', None,
     'sympy/core/tests/test_expr.py:test_equals; test_arit.py:test_issue_8247_8354', 'hidden zero among 3 nested radicals'),
    ('sp.equals4', 'denest', 'zero', 'expand((3**(1/3) + 3)**3)**(1/3) - (3**(1/3) + 3)', None,
     'sympy/core/tests/test_expr.py:test_equals', 'cube root of an expanded cube'),
    ('sp.equals_branch_zero', 'branch', 'zero', '(2*sqrt(2)*x**(5/2)*(1 + 1/(2*x))**(5/2)/5 + 2*sqrt(2)*x**(3/2)*(1 + 1/(2*x))**(5/2)/(-6 - 3/x) - sqrt(2*x + 1)*(6*x**2 + x - 1)/15).subs(x, -1)', None,
     'sympy/core/tests/test_expr.py:test_equals (from integrate(x*sqrt(1+2*x)))', 'antiderivative difference is 0 at x=-1 ...'),
    ('sp.equals_branch_nonzero', 'branch', 'zero', '(2*sqrt(2)*x**(5/2)*(1 + 1/(2*x))**(5/2)/5 + 2*sqrt(2)*x**(3/2)*(1 + 1/(2*x))**(5/2)/(-6 - 3/x) - sqrt(2*x + 1)*(6*x**2 + x - 1)/15).subs(x, -1/4) - 7*sqrt(2)/120', None,
     'sympy/core/tests/test_expr.py:test_equals', '... but equals 7*sqrt(2)/120 (not 0) at x=-1/4: same formula, branch-dependent'),
    ('sp.cardano93a', 'algebraic', 'zero', '-2**(1/3)*(3*sqrt(93) + 29)**2 - 4*(3*sqrt(93) + 29)**(4/3) + 12*sqrt(93)*(3*sqrt(93) + 29)**(1/3) + 116*(3*sqrt(93) + 29)**(1/3) + 174*2**(1/3)*sqrt(93) + 1678*2**(1/3)', None,
     'sympy/core/tests/test_arit.py:test_issue_8247_8354', 'Cardano-type hidden zero'),
    ('sp.cardano93b', 'algebraic', 'zero', '9*(3*sqrt(93) + 29)**(2/3)*((3*sqrt(93) + 29)**(1/3)*(-2**(2/3)*(3*sqrt(93) + 29)**(1/3) - 2) - 2*2**(1/3))**3 + 72*(3*sqrt(93) + 29)**(2/3)*(81*sqrt(93) + 783) + (162*sqrt(93) + 1566)*((3*sqrt(93) + 29)**(1/3)*(-2**(2/3)*(3*sqrt(93) + 29)**(1/3) - 2) - 2*2**(1/3))**2', None,
     'sympy/core/tests/test_arit.py:test_issue_8247_8354', "'a single _mexpand isn't enough'"),
    ('sp.trig90', 'trig', 'zero', '2*(-3*tan(19*pi/90) + sqrt(3))*cos(11*pi/90)*cos(19*pi/90) - sqrt(3)*(-3 + 4*cos(19*pi/90)**2)', None,
     'sympy/core/tests/test_arit.py:test_issue_8247_8354', "'it's zero and it shouldn't hang' (degree-24 cyclotomic)"),
    ('sp.issue4956_num', 'algebraic', 'zero', '-27*12**(1/3)*sqrt(31)*I + 27*2**(2/3)*3**(1/3)*sqrt(31)*I', None,
     'sympy/core/tests/test_evalf.py:test_issue_4956_5204', 'numerator of an expression evalf could not resolve'),
    ('sp.issue4956_den', 'algebraic', 'nonzero', '-2511*2**(2/3)*3**(1/3) + (29*18**(1/3) + 9*2**(1/3)*3**(2/3)*sqrt(31)*I + 87*2**(1/3)*3**(1/6)*I)**2', '-6.9122e4 + 3.9971e4*I',
     'sympy/core/tests/test_evalf.py:test_issue_4956_5204', 'denominator must be proven nonzero before simplifying the quotient to 0'),
    ('sp.hyperbolic_nullspace', 'exp-log', 'zero', '(-exp(1) - 2*cosh(1/3))*(-2*cosh(1/3) - exp(-1)) - (4*cosh(1/3)**2 - 1)**2', None,
     'SymPy docs, tutorials/intro-tutorial/matrices.rst (zero testing in nullspace), with q=1', "needs exp(1) = exp(1/3)^3 inside a tower; SymPy's default iszerofunc returned None"),
    ('sp.binet_exact', 'algebraic', 'zero', 'fibonacci(1000) - (GoldenRatio**1000 - (GoldenRatio - 1)**1000)/sqrt(5)', None,
     'SymPy docs, modules/evalf.rst', 'evalf gives 0.e-1336 even with maxn=1000: needs exact arithmetic'),
    ('sp.binet_near', 'near-miss', 'nonzero', 'fibonacci(1000) - GoldenRatio**1000/sqrt(5)', '-4.6012e-210',
     'SymPy docs, modules/evalf.rst', 'drop the conjugate term: tiny but nonzero'),
    ('sp.binet5000', 'near-miss', 'nonzero', '((1 + sqrt(5))**5000)/(2**5000*sqrt(5)) - fibonacci(5000)', '5.156e-1046',
     'sympy/core/tests/test_evalf.py:test_evalf_near_integers', '1046-digit cancellation'),
    ('sp.j_series', 'near-miss', 'nonzero', '1 - 262537412640768744*exp(-pi*sqrt(163)) - 196884*exp(-2*pi*sqrt(163)) + 103378831900730205293632*exp(-3*pi*sqrt(163))', '1.6137e-59',
     'sympy/core/tests/test_evalf.py:test_evalf_ramanujan', 'j-invariant q-expansion near-miss'),
    ('sp.sin_163', 'near-miss', 'nonzero', 'sin(pi*exp(pi*sqrt(163)))', '-2.356e-12',
     'sympy/core/tests/test_evalf.py:test_evalf_trig', 'sin near a root; argument ~2.6e17'),
    ('sp.e_rational', 'near-miss', 'nonzero', '45 - 613*E/37 + 35/991', '6.0376e-11',
     'sympy/core/tests/test_evalf.py:test_evalf_near_integers', 'e close to a rational'),
    ('sp.sin2017', 'near-miss', 'nonzero', '1 + sin(2017*2**(1/5))', '2.1432e-17',
     'sympy/core/tests/test_evalf.py; MathWorld Almost Integer', 'sin of algebraic close to -1'),
    ('sp.atan_tan_pole', 'branch', 'zero', 'atan(tan(260515)) - 260515 + 82924*pi', None,
     'sympy/core/tests/test_evalf.py (issue 20076)', '260515/pi = 82924.4999991..., so tan(260515) ~ 3.8e5 sits near a pole; branch index must be exact'),
    ('sp.floor_50e', 'floor', 'zero', 'floor(factorial(50)/E) - 11188719610782480504630258070757734324011354208865721592720336800', None,
     'sympy/core/tests/test_evalf.py', 'certified floor of a 65-digit transcendental'),
    ('sp.floor_binet', 'floor', 'zero', 'floor(GoldenRatio**1000/sqrt(5) + 1/2) - fibonacci(1000)', None,
     'sympy/core/tests/test_evalf.py', 'floor needs 210 digits of cancellation'),
    ('sp.ceiling_pyth', 'floor', 'zero', 'ceiling(10*(sin(1)**2 + cos(1)**2)) - 10', None,
     'sympy/core/tests/test_evalf.py', 'floor/ceiling at an exact integer: needs a zero test, not just precision'),
    ('sp.pi_109', 'near-miss', 'nonzero', '2*2**(22/109)*3**(42/109)*5**(90/109)*7**(71/109)/15 - pi', '4.854e-11',
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_issue_7113 (nsimplify(pi, tolerance=1e-9))', 'degree-109 algebraic number close to pi'),
    ('sp.nsimplify_cosatan', 'trig', 'zero', 'cos(atan(1/3)) - 3*sqrt(10)/10', None,
     'SymPy docs, modules/evalf.rst (nsimplify)', 'trig of inverse trig of rational is algebraic'),
    ('sp.nsimplify_exp_atan', 'trig', 'zero', '2 + exp(2*atan(1/4)*I) - 49/17 - 8*I/17', None,
     'SymPy docs, modules/evalf.rst; simplify/tests/test_simplify.py:test_nsimplify', 'exp(i*atan(q)) algebraic'),
    ('sp.nsimplify_root10', 'cyclotomic', 'zero', '1/(exp(3*pi*I/5) + 1) - (1/2 - I*sqrt(sqrt(5)/10 + 1/4))', None,
     'SymPy docs, modules/evalf.rst', '1/(zeta+1) in Q(zeta_10)'),
    ('sp.gamma_reflect', 'special', 'zero', 'gamma(1/4)*gamma(3/4) - sqrt(2)*pi', None,
     "SymPy docs, modules/evalf.rst (nsimplify(gamma('1/4')*gamma('3/4'), [pi]))", 'reflection formula'),
    ('sp.golden', 'algebraic', 'zero', '4/(1 + sqrt(5)) - (-2 + 2*GoldenRatio)', None,
     'SymPy docs, modules/evalf.rst', 'named algebraic constant'),
    ('mp.issue26903', 'algebraic', 'zero', 'sqrt(10000000000000061**2*10000000000000069) - 10000000000000061*sqrt(10000000000000069)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_issue_26903', 'square factor of a 48-digit radicand not extracted (primes > 10^15)'),
    ('mp.issue14831', 'denest', 'zero', '-3*sqrt(12*sqrt(2) + 17) + 12*sqrt(2) + 17 - 2*sqrt(2)*sqrt(12*sqrt(2) + 17)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_issue_14831', 'sqrt(17+12sqrt2) = 3+2sqrt2 hidden'),
    ('mp.issue19760', 'algebraic', 'minpoly:x**4 - 4*x**3 + 4*x**2 - 2', '1/(sqrt(1 + sqrt(2)) - sqrt(2)*sqrt(1 + sqrt(2))) + 1', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_issue_19760', 'compose=True and False gave different answers'),
    ('mp.issue6868', 'algebraic', 'minpoly:8000*x**2 - 48000*x + 71999', '-1/(800*sqrt(-1/240 + 1/(18000*(-1/17280000 + sqrt(15)*I/28800000)**(1/3)) + 2*(-1/17280000 + sqrt(15)*I/28800000)**(1/3))) + 3', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_compose (issue 6868)', 'casus irreducibilis collapsing to a quadratic'),
    ('mp.issue5934_den', 'algebraic', 'zero', '-36000 - 7200*sqrt(5) + (12*sqrt(10)*sqrt(sqrt(5) + 5) + 24*sqrt(10)*sqrt(-sqrt(5) + 5))**2', None,
     'sympy/polys/numberfields/tests/test_minpoly.py (issue 5934, ZeroDivisionError)', 'hidden division by zero in 1/den + 1'),
    ('mp.sin7_sqrt2', 'cyclotomic', 'minpoly:4096*x**12 - 63488*x**10 + 351488*x**8 - 826496*x**6 + 770912*x**4 - 268432*x**2 + 28561', 'sin(pi/7) + sqrt(2)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_compose', 'compositum trig + radical'),
    ('mp.exp7_sqrt2', 'cyclotomic', 'minpoly:x**12 - 2*x**11 - 9*x**10 + 16*x**9 + 43*x**8 - 70*x**7 - 97*x**6 + 126*x**5 + 211*x**4 - 212*x**3 - 37*x**2 + 142*x + 127', 'exp(I*pi/7) + sqrt(2)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_compose', 'root of unity + radical'),
    ('mp.cos7_ratio', 'cyclotomic', 'zero', '(5*cos(2*pi/7) - 7)/(9*cos(pi/7) - 5*cos(3*pi/7)) + 1/(2*cos(pi/7))', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_compose (both sides have minpoly x^3+2x^2-x-1)', 'same minimal polynomial; must also identify the same conjugate'),
    ('mp.cube_sq', 'denest', 'minpoly:x**8 - 8*x**7 - 56*x**6 + 448*x**5 + 480*x**4 - 5056*x**3 + 1984*x**2 + 7424*x - 3008', 'expand((1 + sqrt(2) - 2*sqrt(3) + sqrt(7))**3)**(1/3)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minimal_polynomial_sq', 'cube root of an expanded cube (radicand > 0); degree-8 minpoly'),
    ('mp.cube_sq2', 'denest', 'zero', 'expand((1 + 5*sqrt(2) + 2*sqrt(3))**3)**(1/3) - (1 + 5*sqrt(2) + 2*sqrt(3))', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minimal_polynomial_sq', 'cube root of an expanded cube'),
    ('mp.mixed48', 'algebraic', 'nonzero', 'sqrt(1 + 2**(1/3)) + sqrt(1 + 2**(1/4)) + sqrt(2)', '4.3971',
     'sympy/polys/numberfields/tests/test_minpoly.py (minpoly degree 48, constant term -16630256576)', 'minpoly has degree 48; use for minpoly/degree tests'),
    ('mp.hi_prec', 'algebraic', 'nonzero', '1/sqrt(1 - 9*sqrt(2) + 7*sqrt(3) + 1/10**30)', '1.588',
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minimal_polynomial_hi_prec', 'minpoly coefficients ~10^120 (x^6 coeff given in test)'),
    ('mp.cos15', 'cyclotomic', 'minpoly:16*x**4 + 8*x**3 - 16*x**2 - 8*x + 1', 'cos(pi/15)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py', 'trig -> qqbar'),
    ('mp.sin11', 'cyclotomic', 'minpoly:1024*x**10 - 2816*x**8 + 2816*x**6 - 1232*x**4 + 220*x**2 - 11', 'sin(pi/11)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py', 'degree 10'),
    ('mp.sin21', 'cyclotomic', 'minpoly:4096*x**12 - 11264*x**10 + 11264*x**8 - 4992*x**6 + 960*x**4 - 64*x**2 + 1', 'sin(pi/21)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py', 'degree 12'),
    ('mp.tan5', 'cyclotomic', 'minpoly:x**4 - 10*x**2 + 5', 'tan(pi/5)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py', 'tan -> qqbar'),
    ('mp.tan10', 'cyclotomic', 'minpoly:5*x**4 - 10*x**2 + 1', 'tan(pi/10)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py', 'tan -> qqbar'),
    ('mp.root_rot', 'cyclotomic', 'minpoly:x**3 - 2', '2**(1/3)*exp(2*I*pi/3)', None,
     'sympy/polys/numberfields/tests/test_minpoly.py', 'non-real conjugate of cbrt(2)'),
    ('mp.cbrt_m1', 'branch', 'zero', '-(-1)**(1/3) + (-1)**(2/3) + 1', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_issue_7574', 'principal roots of -1'),
    ('mp.exp3ipi', 'exp-log', 'zero', 'exp(3*I*pi) + 1', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_issue_8353', 'unevaluated exp(3*pi*I)'),
    ('mp.not_alg1', 'transcendental', 'nonzero', 'cos(pi*sqrt(2))', '-0.2663',
     'sympy/polys/numberfields/tests/test_minpoly.py (NotAlgebraic)', 'Gelfond-Schneider: must not be treated as a root of unity'),
    ('mp.not_alg2', 'transcendental', 'nonzero', 'exp(I*pi*sqrt(2))', '-0.26626 - 0.96390*I',
     'sympy/polys/numberfields/tests/test_minpoly.py (NotAlgebraic)', '(-1)^sqrt(2) is transcendental'),
    ('mp.issue23677', 'algebraic', 'minpoly:x**3 + 59426520028417434406408556687919*x**2 + 1161475464966574421163316896737773190861975156439163671112508400*x + 7467465541178623874454517208254940823818304424383315270991298807299003671748074773558707779600', '7680000000000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**4*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**4 - 614323200000000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**4*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**3 + 18458112576000000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**4*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**2 - 246896663036160000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**4*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1) + 1240473830323209600000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**4 - 614323200000000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**3*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**4 - 1476464424954240000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**3*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**2 - 99225501687553535904000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**3 + 18458112576000000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**2*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**4 - 1476464424954240000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**2*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**3 - 593391458458356671712000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**2*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1) + 2981354896834339226880720000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)**2 - 246896663036160000000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**4 - 593391458458356671712000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0)*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**2 - 39878756418031796275267195200*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 0) + 1240473830323209600000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**4 - 99225501687553535904000000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**3 + 2981354896834339226880720000*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1)**2 - 39878756418031796275267195200*RootOf(4000000*x**3 - 239960000*x**2 + 4782399900*x - 31663998001, 1) + 200361370275616536577343808012', None,
     'sympy/polys/numberfields/tests/test_minpoly.py:test_minpoly_issue_23677 (@slow)', 'clustered roots (tiny discriminant), huge cancellation; minpoly coefficients ~10^94'),
    ('dn.shanks29a', 'denest', 'zero', 'sqrt(16 - 2*sqrt(29) + 2*sqrt(55 - 10*sqrt(29))) - sqrt(5) - sqrt(11 - 2*sqrt(29))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'Shanks-type'),
    ('dn.shanks29b', 'denest', 'zero', 'sqrt(-sqrt(5) + sqrt(-2*sqrt(29) + 2*sqrt(-10*sqrt(29) + 55) + 16)) - (11 - 2*sqrt(29))**(1/4)', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'denests to a 4th root'),
    ('dn.jr43', 'denest', 'zero', 'sqrt(5*sqrt(3) + 6*sqrt(2)) - sqrt(2)*3**(1/4) - 3**(3/4)', None,
     'sympy/simplify/tests/test_sqrtdenest.py; Jeffrey & Rich (1999) eq. 4.3', 'denesting introduces 4th roots'),
    ('dn.mq1', 'denest', 'zero', 'sqrt(-4*sqrt(14) - 2*sqrt(6) + 4*sqrt(21) + 33) - (-sqrt(2) + sqrt(3) + 2*sqrt(7))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'multiquadratic, sign selection'),
    ('dn.mq1_ctrl', 'denest', 'nonzero', 'sqrt(-4*sqrt(14) - 2*sqrt(6) + 4*sqrt(21) + 34) - (-sqrt(2) + sqrt(3) + 2*sqrt(7))', '0.088440',
     'sympy/simplify/tests/test_sqrtdenest.py (non-denestable control, constant +1)', 'one unit away from a denestable radicand'),
    ('dn.mq2', 'denest', 'zero', 'sqrt(-28*sqrt(7) - 14*sqrt(5) + 4*sqrt(35) + 82) - (-7 + sqrt(5) + 2*sqrt(7))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'multiquadratic, negative rational part'),
    ('dn.mq3', 'denest', 'zero', 'sqrt(468*sqrt(3) + 3024*sqrt(2) + 2912*sqrt(6) + 19735) - (9*sqrt(3) + 26 + 56*sqrt(6))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'large coefficients'),
    ('dn.mq4', 'denest', 'zero', 'sqrt(-490*sqrt(3) - 98*sqrt(115) - 98*sqrt(345) - 2107) - I*(7*sqrt(5) + 7*sqrt(15) + 7*sqrt(23))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'negative radicand -> imaginary denesting'),
    ('dn.mq5', 'denest', 'zero', 'sqrt(4*sqrt(15) + 8*sqrt(5) + 12*sqrt(3) + 24) - (1 + sqrt(3) + sqrt(5) + sqrt(15))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', '3-generator multiquadratic'),
    ('dn.mq6', 'denest', 'zero', 'sqrt(sqrt(2*sqrt(6) + 5) + sqrt(2*sqrt(7) + 8)) - sqrt(1 + sqrt(2) + sqrt(3) + sqrt(7))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'inner denesting enables outer equality'),
    ('dn.nonga1', 'denest', 'zero', 'sqrt(13 - 2*sqrt(10) + 2*sqrt(2)*sqrt(-2*sqrt(10) + 11)) - (-1 + sqrt(2) + sqrt(10))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'nested radicand'),
    ('dn.nonga2', 'denest', 'zero', 'sqrt((112 + 70*sqrt(2)) + (46 + 34*sqrt(2))*sqrt(5)) - (sqrt(10) + 5 + 4*sqrt(2) + 3*sqrt(5))', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'biquadratic coefficient field'),
    ('dn.nonga3', 'denest', 'zero', 'sqrt(2*sqrt(2)*sqrt(sqrt(2) + 2) + 5*sqrt(2) + 4*sqrt(sqrt(2) + 2) + 8) - (sqrt(2) + sqrt(sqrt(2) + 2) + 2)', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'non-Galois tower Q(sqrt2, sqrt(2+sqrt2))'),
    ('dn.c55', 'denest', 'zero', 'sqrt(8 - sqrt(2)*sqrt(5 - sqrt(5)) - sqrt(3)*(1 + sqrt(5))) - (-sqrt(15)*sqrt(5 - sqrt(5)) - sqrt(3)*sqrt(5 - sqrt(5)) + sqrt(5 - sqrt(5)) + sqrt(5)*sqrt(5 - sqrt(5)) - sqrt(6) - sqrt(2) + sqrt(10) + sqrt(30))/4', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'denesting over Q(sqrt(5-sqrt5)), degree 16 compositum'),
    ('dn.complex1', 'denest', 'zero', '(3 - sqrt(2)*sqrt(4 + 3*I) + 3*I)/2 - I', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'complex radicand'),
    ('dn.complex2', 'denest', 'zero', '-sqrt(-2 + 2*sqrt(3)*I) - (-1 - sqrt(3)*I)', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'complex radicand'),
    ('dn.complex3', 'denest', 'zero', 'sqrt(-8 - sqrt(63)) - I*(sqrt(14) + 3*sqrt(2))/2', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'negative radicand'),
    ('dn.recip', 'denest', 'zero', 'sqrt(1/(4*sqrt(3) + 7) + 1) - (sqrt(2) + sqrt(6))/(sqrt(3) + 2)', None,
     'sympy/simplify/tests/test_sqrtdenest.py', 'denest under a reciprocal'),
    ('dn.ctrl1', 'denest', 'nonzero', 'sqrt(15 - 2*sqrt(31) + 2*sqrt(55 - 10*sqrt(29)))', '2.4511',
     'sympy/simplify/tests/test_sqrtdenest.py (returned unchanged)', 'does not denest'),
    ('dn.ctrl2', 'denest', 'minpoly:x**8 - 8*x**6 + 20*x**4 - 16*x**2 + 2', 'sqrt(2 + sqrt(2 + sqrt(2)))', None,
     'sympy/simplify/tests/test_sqrtdenest.py (returned unchanged); = 2cos(pi/16)', 'genuinely nested, degree 8'),
    ('dn.ram1', 'denest', 'zero', '(2**(1/3) - 1)**(1/3) - ((1/9)**(1/3) - (2/9)**(1/3) + (4/9)**(1/3))', None,
     "Ramanujan; Landau, SIAM J. Comput. 21 (1992); Wikipedia 'Nested radical'", 'cube-root denesting'),
    ('dn.ram2', 'denest', 'zero', 'sqrt(5**(1/3) - 4**(1/3)) - (2**(1/3) + 20**(1/3) - 25**(1/3))/3', None,
     "Ramanujan; Wikipedia 'Nested radical'", 'cube roots under sqrt'),
    ('dn.ram3', 'denest', 'zero', 'sqrt(28**(1/3) - 27**(1/3)) - (98**(1/3) - 28**(1/3) - 1)/3', None,
     "Ramanujan; Wikipedia 'Nested radical'", 'cube roots under sqrt'),
    ('dn.ram4', 'denest', 'zero', '((3 + 2*5**(1/4))/(3 - 2*5**(1/4)))**(1/4) - (5**(1/4) + 1)/(5**(1/4) - 1)', None,
     "Ramanujan; Wikipedia 'Nested radical'", '4th roots'),
    ('dn.ram6', 'denest', 'zero', '((32/5)**(1/5) - (27/5)**(1/5))**(1/3) - ((1/25)**(1/5) + (3/25)**(1/5) - (9/25)**(1/5))', None,
     "Ramanujan; Wikipedia 'Nested radical'", '5th roots inside a cube root'),
    ('dn.ram7', 'denest', 'zero', '(49 + 20*sqrt(6))**(1/4) + (49 - 20*sqrt(6))**(1/4) - 2*sqrt(3)', None,
     "Wikipedia 'Nested radical'", '4th-root sum'),
    ('dn.ram_cos9', 'denest', 'zero', 'real_root(cos(2*pi/9), 3) + real_root(cos(4*pi/9), 3) + real_root(cos(8*pi/9), 3) - real_root((3*9**(1/3) - 6)/2, 3)', None,
     "Ramanujan (Berndt, Ramanujan's Notebooks IV; cf. Wikipedia 'Nested radical')", 'real cube roots of cyclotomic numbers; one is negative'),
    ('dn.ram_cos7', 'denest', 'zero', 'real_root(cos(2*pi/7), 3) + real_root(cos(4*pi/7), 3) + real_root(cos(8*pi/7), 3) - real_root((5 - 3*7**(1/3))/2, 3)', None,
     "Ramanujan; MathWorld 'Trigonometry Angles--Pi/7'", 'real cube roots of cyclotomic numbers'),
    ('dn.cavallo', 'branch', 'zero', 'real_root(7 - 5*sqrt(2), 3) - (1 - sqrt(2))', None,
     "Cavallo, 'Denesting cubic radicals', arXiv:2403.04776", 'REAL cube root; the principal cube root of 7-5sqrt2 < 0 is not 1-sqrt2'),
    ('dn.cavallo_principal', 'branch', 'nonzero', '(7 - 5*sqrt(2))**(1/3) - (1 - sqrt(2))', '0.6213 + 0.3587*I',
     '(companion to dn.cavallo)', 'principal branch: must NOT be proven zero'),
    ('dn.shanks', 'denest', 'zero', 'sqrt(5) + sqrt(22 + 2*sqrt(5)) - sqrt(11 + 2*sqrt(29)) - sqrt(16 - 2*sqrt(29) + 2*sqrt(55 - 10*sqrt(29)))', None,
     "D. Shanks, 'Incredible identities', Fibonacci Quarterly 12 (1974)", 'two different degree-8 towers, same number'),
    ('dn.jr44', 'denest', 'zero', 'sqrt(12 + 2*sqrt(6) + 2*sqrt(14) + 2*sqrt(21)) - sqrt(2) - sqrt(3) - sqrt(7)', None,
     'Jeffrey & Rich (1999) eq. 4.4', 'multiquadratic'),
    ('dn.bombelli', 'algebraic', 'zero', '(2 + 11*I)**(1/3) + (2 - 11*I)**(1/3) - 4', None,
     "Bombelli, L'Algebra (1572); casus irreducibilis of x^3 = 15x + 4", 'principal complex cube roots summing to an integer'),
    ('w.C14', 'denest', 'zero', 'sqrt(2*sqrt(3) + 4) - 1 - sqrt(3)', None,
     'sympy/utilities/tests/test_wester.py:test_C14', 'basic denest'),
    ('w.C15', 'denest', 'zero', 'sqrt(14 + 3*sqrt(3 + 2*sqrt(5 - 12*sqrt(3 - 2*sqrt(2))))) - 3 - sqrt(2)', None,
     'sympy/utilities/tests/test_wester.py:test_C15', '4 levels collapse completely'),
    ('w.C16', 'denest', 'zero', 'sqrt(10 + 2*sqrt(6) + 2*sqrt(10) + 2*sqrt(15)) - sqrt(2) - sqrt(3) - sqrt(5)', None,
     'sympy/utilities/tests/test_wester.py:test_C16', 'multiquadratic'),
    ('w.C18', 'branch', 'zero', 'sqrt(-2 + sqrt(-5))*sqrt(-2 - sqrt(-5)) - 3', None,
     'sympy/utilities/tests/test_wester.py:test_C18', 'product of principal complex sqrts'),
    ('w.C19', 'denest', 'zero', '(90 + 34*sqrt(7))**(1/3) - 3 - sqrt(7)', None,
     'sympy/utilities/tests/test_wester.py:test_C19 (XFAIL in SymPy)', 'cube-root denesting'),
    ('w.C20', 'denest', 'zero', '((135 + 78*sqrt(3))**(2/3) + 3)*sqrt(3)/(135 + 78*sqrt(3))**(1/3) - 12', None,
     'sympy/utilities/tests/test_wester.py:test_C20', 'cube roots in Q(sqrt3)'),
    ('w.C21', 'denest', 'zero', '(41 + 29*sqrt(2))**(1/5) - 1 - sqrt(2)', None,
     'sympy/utilities/tests/test_wester.py:test_C21', '5th root of a unit power'),
    ('w.C22', 'exp-log', 'zero', '((6 - 4*sqrt(2))*log(3 - 2*sqrt(2)) + (3 - 2*sqrt(2))*log(17 - 12*sqrt(2)) + 32 - 24*sqrt(2))/(48*sqrt(2) - 72) - (sqrt(2)/3 - log(sqrt(2) - 1)/3)', None,
     'sympy/utilities/tests/test_wester.py:test_C22 (XFAIL in SymPy)', 'log(17-12sqrt2) = 4 log(sqrt2-1): logs of units'),
    ('w.C13', 'algebraic', 'zero', '10*(1 + 29/1000)**(1/3)/7 - 3**(1/3)', None,
     'sympy/utilities/tests/test_wester.py:test_C13', 'rational times radical recognition'),
    ('w.K2', 'algebraic', 'zero', 'abs(3 - sqrt(7) + I*sqrt(6*sqrt(7) - 15)) - 1', None,
     'sympy/utilities/tests/test_wester.py:test_K2', 'modulus of nested complex number'),
    ('w.K4', 'exp-log', 'zero', 'log(3 + 4*I) - log(5) - I*atan(4/3)', None,
     'sympy/utilities/tests/test_wester.py:test_K4', 'complex log split'),
    ('w.L1', 'algebraic', 'zero', 'sqrt(997) - (997**3)**(1/6)', None,
     'sympy/utilities/tests/test_wester.py:test_L1', 'rational powers'),
    ('w.L2', 'algebraic', 'zero', 'sqrt(999983) - (999983**3)**(1/6)', None,
     'sympy/utilities/tests/test_wester.py:test_L2', 'rational powers, large prime'),
    ('w.L3', 'algebraic', 'zero', '(2**(1/3) + 4**(1/3))**3 - 6*(2**(1/3) + 4**(1/3)) - 6', None,
     'sympy/utilities/tests/test_wester.py:test_L3', 'cube-root field'),
    ('w.I1', 'cyclotomic', 'zero', 'tan(7*pi/10) + sqrt(1 + 2/sqrt(5))', None,
     'sympy/utilities/tests/test_wester.py:test_I1', 'trig at rational pi'),
    ('w.I2', 'trig', 'zero', 'sqrt((1 + cos(6))/2) + cos(3)', None,
     'sympy/utilities/tests/test_wester.py:test_I2 (XFAIL in SymPy)', 'half-angle; sign requires cos(3) < 0'),
    ('w.K8', 'branch', 'nonzero', 'sqrt(1/(-1)) - 1/sqrt(-1)', '2*I',
     'sympy/utilities/tests/test_wester.py:test_K8-K10', 'sqrt(1/z) != 1/sqrt(z) for z<0'),
    ('tr.gauss11', 'cyclotomic', 'zero', 'tan(3*pi/11) + 4*sin(2*pi/11) - sqrt(11)', None,
     "MathWorld 'Trigonometry Angles--Pi/11'", 'Gauss-sum-like, degree reduction in Q(zeta_44)'),
    ('tr.prod_tan11', 'cyclotomic', 'zero', 'tan(pi/11)*tan(2*pi/11)*tan(3*pi/11)*tan(4*pi/11)*tan(5*pi/11) - sqrt(11)', None,
     "MathWorld 'Trigonometry Angles--Pi/11'", 'product over conjugates'),
    ('tr.prod_sin7', 'cyclotomic', 'zero', 'sin(pi/7)*sin(2*pi/7)*sin(3*pi/7) - sqrt(7)/8', None,
     "MathWorld 'Trigonometry Angles--Pi/7'", 'product over conjugates'),
    ('tr.sum_sin7', 'cyclotomic', 'zero', 'sin(pi/7) - sin(2*pi/7) - sin(4*pi/7) + sqrt(7)/2', None,
     "MathWorld 'Trigonometry Angles--Pi/7'", 'quadratic Gauss sum in sine form'),
    ('tr.alt_cos7', 'cyclotomic', 'zero', 'cos(pi/7) - cos(2*pi/7) + cos(3*pi/7) - 1/2', None,
     'classical (IMO 1963 Problem 5)', 'cyclotomic cancellation'),
    ('tr.morrie', 'cyclotomic', 'zero', 'cos(pi/9)*cos(2*pi/9)*cos(4*pi/9) - 1/8', None,
     "Morrie's law", 'telescoping product'),
    ('tr.gauss17', 'cyclotomic', 'zero', 'cos(2*pi/17) - (-1 + sqrt(17) + sqrt(34 - 2*sqrt(17)) + 2*sqrt(17 + 3*sqrt(17) - sqrt(34 - 2*sqrt(17)) - 2*sqrt(34 + 2*sqrt(17))))/16', None,
     "Gauss (1796); e.g. Wikipedia 'Heptadecagon'", 'radical <-> cyclotomic, degree 8'),
    ('tr.cos153', 'cyclotomic', 'zero', 'cos(pi/153) - sin(pi/9)*sin(2*pi/17) - cos(pi/9)*cos(2*pi/17)', None,
     'sympy/functions/elementary/tests/test_trigonometric.py (cos(pi/9/17).rewrite(sqrt))', 'angle subtraction across coprime conductors'),
    ('tr.tan2pi5', 'cyclotomic', 'zero', 'tan(2*pi/5) - sqrt(sqrt(5)/8 + 5/8)/(-1/4 + sqrt(5)/4)', None,
     'sympy/functions/elementary/tests/test_trigonometric.py:test_tan_rewrite', 'tan -> radicals'),
    ('tr.zeta10', 'cyclotomic', 'zero', 'exp(I*pi/5) - (sqrt(5) + 1)/4 - I*sqrt((5 - sqrt(5))/8)', None,
     'classical (cos 36 deg = phi/2)', 'root of unity -> radicals'),
    ('tr.cos85', 'cyclotomic', 'zero', '((-1/4 + sqrt(5)/4)*(-3*sqrt(17 - sqrt(17))*sqrt(sqrt(17) + 17)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/16 - 3*sqrt(34)*sqrt(sqrt(17) + 17)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/32 - 3*sqrt(17 - sqrt(17))*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 - sqrt(sqrt(17) + 17)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/16 - 9*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/8 - sqrt(34)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 - sqrt(34)*sqrt(17 - sqrt(17))*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/32 + 7*sqrt(2)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + sqrt(17)*sqrt(17 - sqrt(17))*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 11*sqrt(2)*sqrt(sqrt(17) + 17)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/32 + 5*sqrt(17)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/8 + 19*sqrt(2)*sqrt(17 - sqrt(17))*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/32) + sqrt(sqrt(5)/8 + 5/8)*sqrt(-9*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/16 + sqrt(17)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/16 + sqrt(2)*sqrt(17 - sqrt(17))*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)/16 + sqrt(2)*sqrt(sqrt(17)/32 + sqrt(2)*sqrt(17 - sqrt(17))/32 + sqrt(2)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/32 + 15/32)*sqrt(-8*sqrt(2)*sqrt(sqrt(17) + 17) - sqrt(2)*sqrt(17 - sqrt(17)) + sqrt(34)*sqrt(17 - sqrt(17)) + 6*sqrt(17) + 34)/16 + 1/2)) - cos(pi/85)', None,
     'SymPy cos(pi/85).rewrite(sqrt) (added with gr_tower)', 'nested square roots of nominal degree beyond 4096 against Q(zeta_85): split, minimal polynomial'),
    ('tr.cos257', 'cyclotomic', 'zero', 'BIG:COS_PI_257_SQRT - cos(pi/257)', None,
     "Gauss's construction written out (cf. sympy test_sincos_rewrite_sqrt_257, @slow: 55 kB, 2956 square roots)", 'STRESS: 643 kB radical expression with 15712 square roots, degree 128 (Fermat prime 257)'),
    ('tr.cheb5', 'trig', 'zero', 'cos(5*acos(1/3)) - 241/243', None,
     'classical (Chebyshev T5)', 'cos(n*acos(q)) is rational: needs multiple-angle reasoning with transcendental angle'),
    ('el.asinh', 'exp-log', 'zero', 'asinh(1) - log(1 + sqrt(2))', None,
     'classical', 'inverse hyperbolic -> log'),
    ('el.acosh3', 'exp-log', 'zero', 'acosh(3) - 2*log(1 + sqrt(2))', None,
     'classical', 'unit relation (1+sqrt2)^2 = 3+sqrt8'),
    ('el.gs_tower', 'exp-log', 'zero', '(sqrt(2)**sqrt(2))**sqrt(2) - 2', None,
     'classical (Gelfond-Schneider tower)', 'exp(sqrt2*log(exp(sqrt2*log sqrt2))): real, no branch issue'),
    ('el.pow_exp', 'exp-log', 'zero', 'exp(sqrt(2)*log(2)) - 2**sqrt(2)', None,
     'definition', 'power = exp-log'),
    ('el.log_fac', 'exp-log', 'zero', 'log(sqrt(sqrt(2)/3)/2) - log(sqrt(sqrt(2)/3)) - log(1/2)', None,
     "FLINT src/ca/test/t-log.c ('Test a bug')", 'log(xyz) = log x + log y + log z (real positive)'),
    ('el.lambertw1', 'special', 'zero', 'LambertW(-log(2)/2) + log(2)', None,
     'classical: W(-ln2/2) = -ln2 since -ln2*e^(-ln2) = -ln2/2 and -ln2 > -1', 'implicitly defined function (Richardson-type extension)'),
    ('el.lambertw2', 'special', 'zero', 'LambertW(2*log(2)) - log(2)', None,
     'classical', 'W(x*e^x) = x'),
    ('el.omega', 'special', 'zero', 'LambertW(1)*exp(LambertW(1)) - 1', None,
     'definition (omega constant)', 'defining equation'),
    ('nm.h163_3', 'near-miss', 'nonzero', 'exp(pi*sqrt(163)/3) - 640320', '6.05e-10',
     "cube root of Ramanujan's constant; cf. van der Hoeven, 'Zero-testing, witness conjectures...' (2001)", 'Heegner'),
    ('nm.h67', 'near-miss', 'nonzero', 'exp(pi*sqrt(67)) - (5280**3 + 744)', '-1.33e-6',
     "Wikipedia 'Almost integer'", 'Heegner 67'),
    ('nm.h43', 'near-miss', 'nonzero', 'exp(pi*sqrt(43)) - (960**3 + 744)', '-2.22e-4',
     "Wikipedia 'Almost integer'", 'Heegner 43'),
    ('nm.r58', 'near-miss', 'nonzero', 'exp(pi*sqrt(58)) - (396**4 - 104)', '-1.78e-7',
     "Ramanujan; MathWorld 'Almost Integer'", 'class number 2'),
    ('nm.epi', 'near-miss', 'nonzero', 'exp(pi) - pi - 20', '-9.0e-4',
     "MathWorld 'Almost Integer'", 'coincidence'),
    ('nm.pi4', 'near-miss', 'nonzero', '22*pi**4 - 2143', '2.75e-6',
     "MathWorld 'Almost Integer' (Ramanujan)", 'coincidence'),
    ('nm.pi45', 'near-miss', 'nonzero', 'pi**4 + pi**5 - exp(6)', '-1.767e-5',
     "MathWorld 'Pi Approximations' (Castellanos)", 'mixed pi/e coincidence'),
    ('nm.log163', 'near-miss', 'nonzero', '163/log(163) - 32', '-1.26e-6',
     "MathWorld 'Almost Integer'", 'log coincidence'),
    ('nm.phi17', 'near-miss', 'nonzero', 'GoldenRatio**17 - 3571', '2.8e-4',
     "MathWorld 'Almost Integer'", 'Pisot power'),
    ('nm.rootred', 'near-miss', 'nonzero', 'sqrt(2) + sqrt(3) - sqrt(5 + 2*sqrt(6)) + 1/10**20000', '1e-20000',
     'Mathematica PossibleZeroQ docs (returns True, wrongly)', 'hidden zero plus tiny rational'),
    ('nm.sq1', 'sqrt-sum', 'nonzero', 'sqrt(10) + sqrt(11) - sqrt(5) - sqrt(18)', '1.94e-4',
     'TOPP Problem 33 (r(20,2))', 'sum-of-square-roots near tie'),
    ('nm.sq2', 'sqrt-sum', 'nonzero', 'sqrt(5) + sqrt(6) + sqrt(18) - sqrt(4) - 2*sqrt(12)', '-4.8e-6',
     'TOPP Problem 33 (r(20,3))', 'sum-of-square-roots near tie'),
    ('nm.sq3', 'sqrt-sum', 'nonzero', 'sqrt(29) + sqrt(1097) + sqrt(3153) - sqrt(226) - sqrt(2324) - sqrt(987)', '2.84e-20',
     "Steinerberger, 'Sums of square roots that are close to an integer', arXiv:2401.10152", '6 roots within 3e-20'),
    ('nm.sq4', 'sqrt-sum', 'nonzero', 'sqrt(11075) + sqrt(27187) + sqrt(68057) - 531', '-1.26e-15',
     'T. Oliveira e Silva, via Steinerberger arXiv:2401.10152', '3 roots near an integer'),
    ('nm.sq5', 'sqrt-sum', 'nonzero', 'sqrt(3) + sqrt(20) + sqrt(23) - 11', '1.83e-5',
     'Steinerberger arXiv:2401.10152', '3 roots near an integer'),
    ('sage.34gon', 'cyclotomic', 'zero', 'sqrt(2)*sqrt(15 + sqrt(17) + sqrt(2)*(sqrt(34 + 6*sqrt(17) + sqrt(2)*(sqrt(17) - 1)*sqrt(17 - sqrt(17)) - 8*sqrt(2)*sqrt(17 + sqrt(17))) + sqrt(17 - sqrt(17))))/8 - cos(pi/17)', None,
     "Sage src/sage/rings/qqbar.py doctest (34-gon; 'x == x2 ... currently infinitely long'); sage-devel QQbar benchmarks", 'nested radical form of cos(pi/17); Sage QQbar could not decide this equality'),
    ('sage.cardano', 'algebraic', 'zero', '(2/(3*sqrt(3)) + 10/27)**(1/3) - 2/(9*(2/(3*sqrt(3)) + 10/27)**(1/3)) + 1/3 - 1', None,
     'Sage src/sage/rings/qqbar.py doctest (AA)', 'Cardano output of a cubic with a rational root'),
    ('big.huge_expr', 'stress', 'zero', 'BIG:HUGE_N + (1 - abs(BIG:HUGE_M)**2)**2', None,
     'FLINT examples/huge_expr.c; ask.sagemath.org/question/52653', 'STRESS: ~7000 ops, 2388+125 sqrts. Sage QQbar: hours; Calcium: ~8 s'),

]

_G = "gr_tower profile"
CASES += [
# --- integer parts --------------------------------------------------------
("fl.heegner_floor", "floor", "true", "Eq(floor(exp(pi*sqrt(163))), 262537412640768743)", None, _G, "near-integer from below, 7.5e-13 away"),
("fl.heegner_ceil", "floor", "true", "Eq(ceiling(exp(pi*sqrt(163))), 262537412640768744)", None, _G, "ceiling of the same"),
("fl.pell40", "floor", "true", "Eq(floor((1 + sqrt(2))**40), 2046573816377473)", None, _G, "(1+sqrt2)^40 = integer - 1/(1+sqrt2)^40"),
("fl.sqrt2_30", "floor", "true", "Eq(floor(sqrt(2)*10**30), 1414213562373095048801688724209)", None, _G, "31 digits"),
("fl.pi_40", "floor", "true", "Eq(floor(pi*10**40), 31415926535897932384626433832795028841971)", None, _G, "41 digits"),
("fl.log_exact", "floor", "true", "Eq(floor(log(10**100)/log(10)), 100)", None, _G, "exact integer: needs a proof, not precision"),
("fl.log8", "floor", "true", "Eq(floor(log(8)/log(2)), 3)", None, _G, "exact integer"),
("fl.atan_exact", "floor", "true", "Eq(floor(4*atan(1)/pi), 1)", None, _G, "exact integer"),
("fl.gamma_exact", "floor", "true", "Eq(floor(gamma(1/2)**2/pi), 1)", None, _G, "exact integer via Gamma(1/2) = sqrt(pi)"),
("fl.exp_pi_minus_pi", "floor", "true", "Eq(floor(exp(pi) - pi), 19)", None, _G, "19.9990999..."),
("fl.neg_sqrt2", "floor", "true", "Eq(floor(-sqrt(2)) + ceiling(-sqrt(2)), -3)", None, _G, "negative arguments"),
("fl.golden", "floor", "true", "Eq(floor(1000*(1 + sqrt(5))/2), 1618)", None, _G, ""),
("fl.half", "floor", "true", "Eq(floor(sqrt(9/4) + 1/2), 2)", None, _G, "exact half-integer plus 1/2"),
("fl.hidden_zero", "floor", "true", "Eq(floor(sqrt(2) + sqrt(3) - sqrt(5 + 2*sqrt(6))), 0)", None, _G, "floor of a hidden zero"),
("fl.hidden_zero_minus", "floor", "true", "Eq(floor(sqrt(2) + sqrt(3) - sqrt(5 + 2*sqrt(6)) - 1/10**30), -1)", None, _G, "a hidden zero minus 10^-30"),
("fl.ceil_trig", "floor", "true", "Eq(ceiling(10*(sin(1)**2 + cos(1)**2)) - floor(10*(sin(1)**2 + cos(1)**2)), 0)", None, _G, "exact integer 10 from both sides"),

# --- comparisons ----------------------------------------------------------
("cmp.epi_pie", "comparison", "true", "Gt(exp(pi), pi**E)", None, _G, "e^pi > pi^e"),
("cmp.sqrt_tie", "comparison", "true", "Gt(sqrt(10) + sqrt(11), sqrt(5) + sqrt(18))", None, _G, "difference 1.9e-4"),
("cmp.machin_eq_le", "comparison", "true", "Le(atan(1/2) + atan(1/3), pi/4)", None, _G, "equality case of <="),
("cmp.machin_eq_lt", "comparison", "false", "Lt(atan(1/2) + atan(1/3), pi/4)", None, _G, "equality: strict inequality false"),
("cmp.machin_above", "comparison", "true", "Lt(atan(1/2) + atan(1/3), pi/4 + 1/10**50)", None, _G, "10^-50 above an exact equality"),
("cmp.machin_below", "comparison", "true", "Gt(atan(1/2) + atan(1/3), pi/4 - 1/10**50)", None, _G, "10^-50 below an exact equality"),
("cmp.sqrt6_ge", "comparison", "true", "Ge(sqrt(2)*sqrt(3), sqrt(6))", None, _G, "equality case of >="),
("cmp.exp_tiny", "comparison", "true", "Gt(exp(1/10**20), 1)", None, _G, "exp of a tiny argument"),
("cmp.log1p", "comparison", "true", "Lt(log(1 + 1/10**30), 1/10**30)", None, _G, "log(1+x) < x, difference 5e-61"),
("cmp.sin_huge", "comparison", "true", "Lt(sin(10**22), 0)", None, _G, "argument reduction for 10^22"),
("cmp.cos355", "comparison", "true", "Gt(cos(355), -1)", None, _G, "355 = 113 pi + 3e-5"),
("cmp.tan_pole", "comparison", "true", "Lt(tan(pi/2 - 1/10**30), 10**31)", None, _G, "next to a pole of tan"),
("cmp.root_x5", "comparison", "true", "Gt(RootOf(x**5 - x - 1, 0), 11673/10000)", None, _G, "real root of x^5 - x - 1 = 1.16730..."),
("cmp.heegner_ratio", "comparison", "true", "Lt(10**-30, (640320**3 + 744)/exp(pi*sqrt(163)) - 1)", None, "Calcium notebook", "lower bound"),
("cmp.heegner_ratio2", "comparison", "true", "Lt((640320**3 + 744)/exp(pi*sqrt(163)) - 1, 10**-29)", None, "Calcium notebook", "upper bound"),
("cmp.exp_log_2_3", "comparison", "true", "Gt(exp(log(2)*log(3)), 2)", None, "Calcium notebook", ""),
("cmp.qqbar_roots", "comparison", "true", "Gt((1 + 100/2**1000)**(1/100), (1 + 101/2**1000)**(1/101))", None, "Sage QQbar documentation", "difference 4.4e-603"),

# --- real and imaginary parts, absolute values ----------------------------
("pt.re_sqrt", "parts", "zero", "re(sqrt(-3 + 4*I)) - 1", None, _G, "sqrt(-3+4i) = 1 + 2i"),
("pt.im_sqrt_cut", "parts", "zero", "im(sqrt(-3 - 4*I)) + 2", None, _G, "sqrt(-3-4i) = 1 - 2i"),
("pt.abs_exp_irr", "parts", "zero", "Abs(exp(I*sqrt(2)*pi)) - 1", None, _G, ""),
("pt.abs_exp", "parts", "zero", "Abs(exp(1 + I*pi/3)) - E", None, _G, ""),
("pt.re_log", "parts", "zero", "re(log(-1 - I)) - log(2)/2", None, _G, ""),
("pt.im_log", "parts", "zero", "im(log(-1 - I)) + 3*pi/4", None, _G, "principal argument in (-pi, pi]"),
("pt.i_to_i", "parts", "zero", "re(I**I) - exp(-pi/2) + im(I**I)", None, _G, ""),
("pt.abs_alg", "parts", "zero", "Abs(sqrt(2) + sqrt(3)*I) - sqrt(5)", None, _G, ""),
("pt.conj_exp", "parts", "zero", "conjugate(exp(I*sqrt(2))) - exp(-I*sqrt(2))", None, _G, ""),
("pt.abs_chord", "parts", "zero", "Abs(1 - exp(2*pi*I/5)) - sqrt((5 - sqrt(5))/2)", None, _G, "chord of the pentagon"),
("pt.abs_gamma_i", "parts", "zero", "Abs(gamma(I))**2 - pi/sinh(pi)", None, _G, "|Gamma(i)|^2 = pi/sinh(pi)"),
("pt.abs_gamma_half", "parts", "zero", "Abs(gamma(1/2 + I))**2 - pi/cosh(pi)", None, _G, "|Gamma(1/2+i)|^2 = pi/cosh(pi)"),
("pt.re_exp_sum", "parts", "zero", "re(exp(I*sqrt(2)) + exp(-I*sqrt(2))) - 2*cos(sqrt(2))", None, _G, ""),
("pt.arg_root", "parts", "zero", "arg((1 + I*sqrt(3))**5) + pi/3", None, _G, "5 pi/3 reduced to -pi/3"),

# --- signs -----------------------------------------------------------------
("sg.pi_22_7", "sign", "true", "Eq(sign(pi - 22/7), -1)", None, _G, ""),
("sg.heegner", "sign", "true", "Eq(sign(exp(pi*sqrt(163)) - 262537412640768744), -1)", None, _G, "difference -7.5e-13"),
("sg.hidden_zero", "sign", "true", "Eq(sign(sqrt(2) + sqrt(3) - sqrt(5 + 2*sqrt(6))), 0)", None, _G, "sign of a hidden zero"),
("sg.complex", "sign", "zero", "sign(1 + I) - (1 + I)/sqrt(2)", None, _G, "complex sign z/|z|"),
("sg.imag", "sign", "zero", "sign(-3*I*pi) + I", None, _G, ""),
("sg.sin_kpi", "sign", "true", "Eq(sign(sin(10**30*pi)), 0)", None, _G, "sin(k pi) = 0 for k = 10^30"),
("sg.cos_kpi", "sign", "true", "Eq(sign(cos(10**30*pi/3)), -1)", None, _G, "10^30 = 4 mod 6: cos = -1/2"),
("tr.exp_big", "cyclotomic", "zero", "exp(10**40*pi*I/7) - exp(2*pi*I/7)**2", None, _G, "reduction of a huge rational multiple of pi i"),
("sg.cos_1e10", "sign", "true", "Eq(sign(cos(10**10)), 1)", None, _G, "cos(10^10) = 0.873..."),
("sg.re_log", "sign", "true", "Eq(sign(re(log(1/2 + I*sqrt(3)/2))), 0)", None, _G, "|z| = 1: real part of log is zero"),

# --- branch cuts -----------------------------------------------------------
("br.log_neg", "branch", "zero", "log(-2) - log(2) - I*pi", None, _G, ""),
("br.sqrt_neg", "branch", "zero", "sqrt(-4) - 2*I", None, _G, ""),
("br.cbrt_neg1", "branch", "zero", "(-1)**(1/3) - (1 + sqrt(3)*I)/2", None, _G, "principal cube root"),
("br.pow_neg8", "branch", "zero", "(-8)**(2/3) - 4*exp(2*pi*I/3)", None, _G, ""),
("br.log_exp3", "branch", "zero", "log(exp(3*I)) - 3*I", None, _G, "3 in (-pi, pi]"),
("br.log_exp4", "branch", "zero", "log(exp(4*I)) - (4 - 2*pi)*I", None, _G, "4 reduced by 2 pi"),
("br.log_exp10", "branch", "zero", "log(exp(10*I)) - (10 - 4*pi)*I", None, _G, ""),
("br.sqrt_exp4", "branch", "zero", "sqrt(exp(4*I)) + exp(2*I)", None, _G, "principal sqrt of exp(4i) is -exp(2i)"),
("br.atan_tan2", "branch", "zero", "atan(tan(2)) - 2 + pi", None, _G, ""),
("br.asin_sin3", "branch", "zero", "asin(sin(3)) - (pi - 3)", None, _G, ""),
("br.acos_cos4", "branch", "zero", "acos(cos(4)) - (2*pi - 4)", None, _G, ""),
("br.log_neg_e", "branch", "zero", "log(-exp(1)) - 1 - I*pi", None, _G, ""),
("br.near_cut_above", "branch", "zero", "im(log(-1 + I/10**30)) - (pi - atan(1/10**30))", None, _G, "10^-30 above the cut"),
("br.near_cut_below", "branch", "zero", "im(log(-1 - I/10**30)) + (pi - atan(1/10**30))", None, _G, "10^-30 below the cut"),
("br.sqrt_near_cut", "branch", "true", "Lt(im(sqrt(-1 - I/10**30)), 0)", None, _G, "sqrt just below the cut"),
("br.pow_irrational", "branch", "zero", "(-1)**sqrt(2) - exp(I*pi*sqrt(2))", None, _G, "principal power"),
("br.log_exp_irr", "branch", "zero", "log(exp(1 + sqrt(2)*I*pi)) - 1 - (sqrt(2) - 2)*I*pi", None, _G, "sqrt(2) pi reduced by 2 pi"),

# --- Richardson: exponentials and logarithms ------------------------------
("rich.exp_log_sum", "exp-log", "zero", "exp(log(2) + log(3)) - 6", None, _G, ""),
("rich.log_exp_prod", "exp-log", "zero", "log(exp(2)*exp(3)) - 5", None, _G, ""),
("rich.alg_exponent", "exp-log", "zero", "exp(sqrt(2)*log(3))*exp(-sqrt(2)*log(3)) - 1", None, _G, ""),
("rich.pow_pow", "exp-log", "zero", "(2**sqrt(3))**sqrt(3) - 8", None, _G, ""),
("rich.log_pow", "exp-log", "zero", "log(2**sqrt(2)) - sqrt(2)*log(2)", None, _G, ""),
("rich.exp_pi_e", "exp-log", "zero", "log(exp(pi)*exp(E)) - pi - E", None, _G, ""),
("rich.sqrt8", "exp-log", "zero", "exp(pi*sqrt(2))*exp(pi*sqrt(8)) - exp(3*pi*sqrt(2))", None, _G, "algebraic simplification of the arguments"),
("rich.units", "exp-log", "zero", "log(1 + sqrt(2)) + log(sqrt(2) - 1)", None, _G, "logarithms of units"),
("rich.unit_sqrt", "exp-log", "zero", "log(2 + sqrt(3)) - 2*log((sqrt(6) + sqrt(2))/2)", None, _G, ""),
("rich.nested", "exp-log", "zero", "log(log(exp(exp(sqrt(2))))) - sqrt(2)", None, _G, ""),
("rich.exp_unit", "exp-log", "zero", "exp(2*log(sqrt(2) + 1)) - 3 - 2*sqrt(2)", None, _G, ""),
("rich.many_logs", "exp-log", "zero", "log(15) - log(3) - log(5) + log(2)*log(3) - log(3)*log(2)", None, _G, ""),
("rich.exp_i_irr", "exp-log", "zero", "exp(I*pi*sqrt(2))*exp(I*pi*(2 - sqrt(2))) - 1", None, _G, ""),
("rich.exp_sum_alg", "exp-log", "zero", "exp(sqrt(2) + sqrt(3)) - exp(sqrt(2))*exp(sqrt(3))", None, _G, ""),
("rich.exp_denest", "exp-log", "zero", "exp(sqrt(5 + 2*sqrt(6))) - exp(sqrt(2))*exp(sqrt(3))", None, _G, "a nested radical in the argument"),
("rich.log_denest", "exp-log", "zero", "log(exp(sqrt(5 + 2*sqrt(6)))) - sqrt(2) - sqrt(3)", None, _G, ""),
("rich.tanh_log", "exp-log", "zero", "tanh(log(3)) - 4/5", None, _G, ""),
("rich.exp_log_chain", "exp-log", "zero", "exp(log(exp(log(2) + 1)) - 1) - 2", None, _G, ""),
("rich.cube_root_exp", "exp-log", "zero", "exp(1/3)**3 - E", None, _G, "exp(1) = exp(1/3)^3"),
("rich.exp_alg_log", "exp-log", "zero", "exp(log(1 + sqrt(2))/2)**4 - 3 - 2*sqrt(2)", None, _G, "exp of half a logarithm of a unit"),
("rich.log_alg_power", "exp-log", "zero", "log((1 + sqrt(2))**10) - 10*log(1 + sqrt(2))", None, _G, ""),
("rich.log_root3", "exp-log", "zero", "log(2**(1/3)*3**(1/2)) - log(2)/3 - log(3)/2", None, _G, ""),
("rich.exp_mixed", "exp-log", "zero", "exp(pi + sqrt(2)*log(2)) - exp(pi)*2**sqrt(2)", None, _G, ""),
("rich.nonzero_alg", "exp-log", "nonzero", "exp(sqrt(2)*log(3)) - 3**(14142/10000)", "7.0458e-5", _G, "an algebraic exponent against a rational approximation"),
("rich.stress_nested", "exp-log", "zero", "exp(exp(exp(log(log(log(sqrt(2) + 20)))))) - sqrt(2) - 20", None, _G, "three levels each way"),

# --- algebraic and transcendental extensions together --------------------
("mix.denest_pi", "mixed", "zero", "sqrt(pi + 2*sqrt(2)*sqrt(pi) + 2) - sqrt(pi) - sqrt(2)", None, _G, "denesting over Q(pi)"),
("mix.denest_pi_sqrt3", "mixed", "zero", "sqrt(pi**2 + 2*pi*sqrt(3) + 3) - pi - sqrt(3)", None, _G, ""),
("mix.denest_e", "mixed", "zero", "sqrt(exp(2) - 2*exp(1)*sqrt(2) + 2) - exp(1) + sqrt(2)", None, _G, "e - sqrt(2) > 0"),
("mix.rationalize", "mixed", "zero", "1/(pi + sqrt(2)) - (pi - sqrt(2))/(pi**2 - 2)", None, _G, ""),
("mix.cancel", "mixed", "zero", "(pi**2 - 2)/(pi - sqrt(2)) - pi - sqrt(2)", None, _G, ""),
("mix.cbrt_pi", "mixed", "zero", "(pi**3 + 3*pi**2*2**(1/3) + 3*pi*2**(2/3) + 2)**(1/3) - pi - 2**(1/3)", None, _G, "cube root over Q(pi)"),
("mix.exp_alg_coeffs", "mixed", "zero", "(exp(1) + sqrt(2))*(exp(1) - sqrt(2)) - exp(2) + 2", None, _G, ""),
("mix.log_sqrt_pi", "mixed", "zero", "log(sqrt(pi)) - log(pi)/2", None, _G, ""),
("mix.sqrt_pi_i", "mixed", "zero", "sqrt(-pi)**2 + pi", None, "Calcium (profile)", ""),
("mix.root_sum", "mixed", "zero", "RootOf(x**3 - x - 1, 0) + RootOf(x**3 - x - 1, 1) + RootOf(x**3 - x - 1, 2)", None, _G, "sum of the roots"),
("mix.root_exp", "mixed", "zero", "log(exp(RootOf(x**3 - x - 1, 0))) - RootOf(x**3 - x - 1, 0)", None, _G, "an algebraic number of degree 3 in exp/log"),
("mix.root_pi", "mixed", "zero", "(pi + RootOf(x**3 - x - 1, 0))**3 - (pi + RootOf(x**3 - x - 1, 0))**2*(pi + RootOf(x**3 - x - 1, 0))", None, _G, ""),

# --- Machin-like formulas in non-Gaussian fields -------------------------
("mac.sqrtm2_a", "machin", "zero", "atan(2*sqrt(2)) - 2*atan(sqrt(2)/3) + 2*atan(4*sqrt(2)) - pi", None, _G, "arguments in Z[sqrt(-2)] (pslq, verified exactly)"),
("mac.sqrtm2_b", "machin", "zero", "atan(2*sqrt(2)) + 2*atan(2*sqrt(2)/3) + 2*atan(sqrt(2)/7) - pi", None, _G, "Z[sqrt(-2)]"),
("mac.sqrtm2_c", "machin", "zero", "3*atan(sqrt(2)) + atan(sqrt(2)/5) - pi", None, _G, "Z[sqrt(-2)]: (1 + sqrt(-2))^3 (5 + sqrt(-2))"),
("mac.sqrtm3_a", "machin", "zero", "2*atan(2*sqrt(3)) - 2*atan(2*sqrt(3)/3) + 2*atan(5*sqrt(3)/4) - pi", None, _G, "Z[sqrt(-3)]"),
("mac.sqrtm3_b", "machin", "zero", "3*atan(2*sqrt(3)) + 3*atan(3*sqrt(3)) - 3*atan(3*sqrt(3)/8) - 2*pi", None, _G, "Z[sqrt(-3)]"),
("mac.sqrtm7_a", "machin", "zero", "atan(sqrt(7)/3) - 2*atan(sqrt(7)/4) - 2*atan(5*sqrt(7)/3) + pi", None, _G, "Z[sqrt(-7)]"),
("mac.sqrtm7_b", "machin", "zero", "atan(sqrt(7)/3) + 2*atan(2*sqrt(7)) - 2*atan(sqrt(7)/15) - pi", None, _G, "Z[sqrt(-7)]"),
("mac.hyp_sqrt2_a", "machin", "zero", "atanh(sqrt(2)/3) + atanh(2*sqrt(2)/5) + atanh(5*sqrt(2)/13) - 2*log(1 + sqrt(2))", None, _G, "Z[sqrt(2)]: the unit 1 + sqrt(2)"),
("mac.hyp_sqrt2_b", "machin", "zero", "atanh(sqrt(2)/3) + atanh(4*sqrt(2)/7) + atanh(sqrt(2)/11) - 2*log(1 + sqrt(2))", None, _G, "Z[sqrt(2)]"),
("mac.hyp_sqrt5", "machin", "zero", "atanh(sqrt(5)/4) + atanh(3*sqrt(5)/8) + atanh(8*sqrt(5)/23) - 6*log((1 + sqrt(5))/2)", None, _G, "Z[(1+sqrt(5))/2]: the golden ratio"),
("mac.sqrtm2_near", "machin", "nonzero", "atan(2*sqrt(2)) - 2*atan(sqrt(2)/3) + 2*atan(4*sqrt(2)) - pi + 1/10**40", "1e-40", _G, "near miss of mac.sqrtm2_a"),

# --- special functions ----------------------------------------------------
("sf.li2_half", "special", "zero", "polylog(2, 1/2) - pi**2/12 + log(2)**2/2", None, _G, ""),
("sf.li2_golden", "special", "zero", "polylog(2, (sqrt(5) - 1)/2) - pi**2/10 + log((sqrt(5) - 1)/2)**2", None, _G, "Landen"),
("sf.li2_golden2", "special", "zero", "polylog(2, (3 - sqrt(5))/2) - pi**2/15 + log((sqrt(5) - 1)/2)**2", None, _G, ""),
("sf.li2_reflect", "special", "zero", "polylog(2, 1/3) + polylog(2, 2/3) - pi**2/6 + log(1/3)*log(2/3)", None, _G, "reflection at 1/3"),
("sf.li2_reflect_alg", "special", "zero", "polylog(2, sqrt(2) - 1) + polylog(2, 2 - sqrt(2)) - pi**2/6 + log(sqrt(2) - 1)*log(2 - sqrt(2))", None, _G, "reflection at an algebraic point"),
("sf.li2_inversion", "special", "zero", "polylog(2, -3) + polylog(2, -1/3) + pi**2/6 + log(3)**2/2", None, _G, "inversion"),
("sf.li2_third", "special", "zero", "polylog(2, 1/3) - polylog(2, 1/9)/6 - pi**2/18 + log(3)**2/6", None, _G, ""),
("sf.hurwitz_quarter", "special", "zero", "zeta(2, 1/4) - pi**2 - 8*Catalan", None, _G, ""),
("sf.hurwitz_thirds", "special", "zero", "zeta(2, 1/3) + zeta(2, 2/3) - 4*pi**2/3", None, _G, "distribution"),
("sf.hurwitz_half", "special", "zero", "zeta(3, 1/2) - 7*zeta(3)", None, _G, ""),
("sf.hurwitz_sixths", "special", "zero", "zeta(2, 1/6) + zeta(2, 5/6) - 4*pi**2", None, _G, "reflection"),
("sf.gamma_thirds", "special", "zero", "gamma(1/3)*gamma(2/3) - 2*pi/sqrt(3)", None, _G, ""),
("sf.gamma_quarters", "special", "zero", "gamma(1/4)*gamma(3/4) - pi*sqrt(2)", None, _G, ""),
("sf.gamma_fifths", "special", "zero", "gamma(1/5)*gamma(2/5)*gamma(3/5)*gamma(4/5) - 4*pi**2/sqrt(5)", None, _G, "multiplication formula"),
("sf.gamma_sixth", "special", "zero", "gamma(1/6) - 2**(-1/3)*sqrt(3/pi)*gamma(1/3)**2", None, _G, "duplication and reflection combined"),
("sf.gamma_line_shift", "special", "zero", "gamma(sqrt(2) + 1) - sqrt(2)*gamma(sqrt(2))", None, _G, "Gamma on the line sqrt(2) + Z"),
("sf.gamma_line_reflect", "special", "zero", "gamma(sqrt(2))*gamma(1 - sqrt(2)) - pi/sin(pi*sqrt(2))", None, _G, "reflection at an algebraic point"),
("sf.gamma_line_dup", "special", "zero", "gamma(sqrt(2))*gamma(sqrt(2) + 1/2) - 2**(1 - 2*sqrt(2))*sqrt(pi)*gamma(2*sqrt(2))", None, _G, "duplication at an algebraic point"),
("sf.gamma_pi", "special", "zero", "gamma(pi + 2) - (pi + 1)*pi*gamma(pi)", None, _G, "Gamma on the line pi + Z"),
("sf.gamma_complex", "special", "zero", "gamma(1/2 + I)*gamma(1/2 - I) - pi/cosh(pi)", None, _G, ""),
("sf.digamma_shift", "special", "zero", "digamma(sqrt(2) + 1) - digamma(sqrt(2)) - 1/sqrt(2)", None, _G, ""),
("sf.digamma_reflect", "special", "zero", "digamma(1 - sqrt(2)) - digamma(sqrt(2)) - pi*cos(pi*sqrt(2))/sin(pi*sqrt(2))", None, _G, "reflection at an algebraic point"),
("sf.trigamma_quarter", "special", "zero", "polygamma(1, 1/4) - pi**2 - 8*Catalan", None, _G, ""),
("sf.trigamma_reflect", "special", "zero", "polygamma(1, sqrt(2)) + polygamma(1, 1 - sqrt(2)) - pi**2/sin(pi*sqrt(2))**2", None, _G, ""),
("sf.digamma_third", "special", "zero", "digamma(1/3) + EulerGamma + pi/(2*sqrt(3)) + 3*log(3)/2", None, _G, "Gauss's digamma theorem"),
("sf.digamma_quarter", "special", "zero", "digamma(1/4) + EulerGamma + pi/2 + 3*log(2)", None, _G, ""),
("sf.gamma_near", "special", "nonzero", "gamma(1/3)*gamma(2/3) - 2*pi/sqrt(3) + 1/10**30", "1e-30", _G, "near miss of sf.gamma_thirds"),

# --- inequalities in asymptotic regimes -----------------------------------
("asy.tan_100i_upper", "asymptotic", "true", "Lt(Abs(tan(1 + 100*I) - I), 3*exp(-200))", None, _G, "tan(x + iy) - i ~ 2 exp(-2y)"),
("asy.tan_100i_lower", "asymptotic", "true", "Gt(Abs(tan(1 + 100*I) - I), exp(-200))", None, _G, ""),
("asy.tan_1000i", "asymptotic", "true", "Lt(Abs(tan(1 + 1000*I) - I), 3*exp(-2000))", None, _G, "about 3000 bits of cancellation"),
("asy.erfc10_upper", "asymptotic", "true", "Lt(erfc(10), exp(-100)/(10*sqrt(pi)))", None, _G, "erfc(x) < exp(-x^2)/(x sqrt(pi))"),
("asy.erfc10_lower", "asymptotic", "true", "Gt(erfc(10), exp(-100)/(10*sqrt(pi))*(1 - 1/200))", None, _G, "next term of the expansion"),
("asy.erfc30", "asymptotic", "true", "Gt(erfc(30), exp(-900)/(30*sqrt(pi))*(1 - 1/1800))", None, _G, ""),
("asy.erfc100", "asymptotic", "true", "Lt(erfc(100), exp(-10000)/(100*sqrt(pi)))", None, _G, "erfc(100) = 6.4e-4346"),
("asy.stirling_lower", "asymptotic", "true", "Gt(gamma(100), sqrt(2*pi)*100**(199/2)*exp(-100))", None, _G, "Stirling"),
("asy.stirling_upper", "asymptotic", "true", "Lt(gamma(100), sqrt(2*pi)*100**(199/2)*exp(-100)*exp(1/1200))", None, _G, ""),
("asy.stirling_1e6", "asymptotic", "true", "Lt(gamma(10**6), sqrt(2*pi)*(10**6)**(10**6 - 1/2)*exp(-10**6)*exp(1/(12*10**6)))", None, _G, "Gamma(10^6) ~ 10^5565703"),
("asy.exp_exp", "asymptotic", "true", "Lt(exp(-exp(10)), 10**-9565)", None, _G, "exp(-exp(10)) = 1.4e-9566"),
("asy.cos_tiny", "asymptotic", "true", "Gt(1 - cos(1/10**20), 0)", None, _G, "cancellation: 5e-41"),
("asy.sinh_tiny", "asymptotic", "true", "Lt(sinh(1/10**30) - 1/10**30, 10**-89)", None, _G, "x^3/6 = 1.7e-91"),
("asy.zeta_pole", "asymptotic", "true", "Gt(zeta(1 + 1/10**10), 10**10)", None, _G, "zeta(1 + eps) = 1/eps + gamma + ..."),
("asy.atan_inf", "asymptotic", "true", "Gt(atan(10**30), pi/2 - 1/10**29)", None, _G, ""),
("asy.lambertw_large", "asymptotic", "true", "Gt(LambertW(10**100), log(10**100) - log(log(10**100)))", None, _G, "W(x) > log x - log log x"),
]

# ---------------------------------------------------------------------------
# the Calcium issue tracker, the Calcium notebook, Sage's QQbar documentation
# ---------------------------------------------------------------------------

_CI = "Calcium issue tracker (github.com/flintlib/calcium/issues) #"
_NB = "Calcium introduction notebook"
_SQ = "Sage QQbar documentation (src/sage/rings/qqbar.py)"
_FT = "flint_ctypes tests"

def _dirichlet_kernel(n):
    s = " + ".join("cos(%d*a)" % k for k in range(1, n + 1))
    return "(%s - sin(%d*a/2)/(2*sin(a/2)) + 1/2).subs(a, 1 + sqrt(2))" % (s, 2 * n + 1)

CASES += [
    ('pf.floor_sqrt2', 'floor', 'true', 'Eq(floor(sqrt(2)), 1)', None, _FT, ''),
    ('pf.gs_neg', 'branch', 'zero', '(sqrt(-2)**sqrt(2))**sqrt(2) + 2', None, _FT, 'complex Gelfond-Schneider tower'),
    ('pf.gs3', 'exp-log', 'zero', '(sqrt(3)**sqrt(3))**sqrt(3) - 3*sqrt(3)', None, _FT, ''),
    ('pf.log_pi', 'exp-log', 'zero', 'log(1 + pi) - log(pi) - log(1 + 1/pi)', None, _FT, 'logarithms over Q(pi)'),
    ('pf.loglog', 'exp-log', 'zero', 'log(log(-log(log(exp(exp(-exp(exp(3)))))))) - 3', None, _FT, 'exp(-exp(exp(3))) = 1e-230'),
    ('pf.e2', 'exp-log', 'zero', 'E**2 - exp(2)', None, _FT, ''),
    ('pf.gd_sin', 'trig', 'zero', 'sin(2*atan(exp(1)) - pi/2) - tanh(1)', None, _FT, 'sin(gd(1)) = tanh(1)'),
    ('pf.gd_tan', 'trig', 'zero', 'tan(2*atan(exp(1)) - pi/2) - sinh(1)', None, _FT, 'tan(gd(1)) = sinh(1)'),
    ('pf.gd_sqrt2', 'trig', 'zero', 'sin(2*atan(exp(sqrt(2))) - pi/2) - tanh(sqrt(2))', None, _FT, 'sin(gd(sqrt 2)) = tanh(sqrt 2)'),
    ('pf.gd_half', 'trig', 'zero', 'tan((2*atan(exp(1)) - pi/2)/2) - tanh(1/2)', None, _FT, 'tan(gd(1)/2) = tanh(1/2)'),
    ('pf.root8_pi', 'cyclotomic', 'zero', 'sqrt(2)*(1 + I)/2*pi - exp(pi*I/4)*pi', None, _FT, 'root of unity times a transcendental'),
    ('pf.root12_pi', 'cyclotomic', 'zero', '(sqrt(3) + I)/2*pi - exp(pi*I/6)*pi', None, _FT, ''),
    ('pf.abs_exp_sqrt', 'parts', 'zero', 'Abs(exp(sqrt(1 + I))) - exp(re(sqrt(1 + I)))', None, _FT, ''),
    ('pf.tan_pi_sqrt', 'trig', 'zero', 'tan(pi*sqrt(2))*tan(pi*sqrt(3)) - (cos(pi*sqrt(5 - 2*sqrt(6))) - cos(pi*sqrt(5 + 2*sqrt(6))))/(cos(pi*sqrt(5 - 2*sqrt(6))) + cos(pi*sqrt(5 + 2*sqrt(6))))', None, _FT, 'product-to-sum with algebraic angles'),
    ('pf.tan_sqrt_pi', 'trig', 'zero', 'tan(sqrt(2*pi))*tan(sqrt(3*pi)) - (cos(sqrt(pi*(5 - 2*sqrt(6)))) - cos(sqrt(pi*(5 + 2*sqrt(6)))))/(cos(sqrt(pi*(5 - 2*sqrt(6)))) + cos(sqrt(pi*(5 + 2*sqrt(6)))))', None, _FT, 'denesting sqrt(pi (5 + 2 sqrt 6)) = sqrt(2 pi) + sqrt(3 pi)'),
    ('pf.log_exp_i', 'branch', 'zero', 'log(exp(I)/exp(-I)) - 2*I', None, _FT, ''),
    ('pf.acos_cos1', 'branch', 'zero', 'acos(cos(1)) - 1', None, _FT, ''),
    ('pf.acos_cos_alg', 'branch', 'zero', 'acos(cos(sqrt(2) - 1)) - sqrt(2) + 1', None, _FT, ''),
    ('pf.dirichlet16', 'trig', 'zero', _dirichlet_kernel(16), None, _FT, 'Dirichlet kernel: 16 cosines of multiples of 1 + sqrt(2)'),
    ('pf.sin3a', 'trig', 'zero', '(sin(3*a) - 4*sin(a)*sin(pi/3 - a)*sin(pi/3 + a)).subs(a, 1 + sqrt(2))', None, _FT, ''),
    ('pf.sin_sqrt', 'trig', 'zero', '(sin(a) - sqrt(1 - cos(a)**2)).subs(a, 1 + sqrt(2))', None, _FT, 'needs sin(1 + sqrt 2) > 0'),
    ('pf.qqbar_sqrt_sq', 'branch', 'zero', 'sqrt((-23 + 5*I)**2) - (23 - 5*I)', None, _FT, 'sqrt(a^2) = -a'),
    ('pf.stoutemyer1', 'algebraic', 'zero', 'I*((3 - 5*I)*pi + 1)/(((5 + 3*I)*pi + I)*pi) - 1/pi', None, _CI + '1 (Stoutemyer)', 'rational function of pi with Gaussian coefficients'),
    ('pf.stoutemyer2', 'branch', 'zero', '(-1)**(1/8)*sqrt(I + 1)/2**(3/4) + I*exp(I*pi/2) - (I - 1)/2', None, _CI + '1 (Stoutemyer)', ''),
    ('pf.stoutemyer3', 'trig', 'zero', 'sin(2*atan(pi)) + (15*sqrt(3) + 26)**(1/3) - 2*pi/(pi**2 + 1) - sqrt(3) - 2', None, _CI + '1 (Stoutemyer)', 'sin(2 atan(x)) and a cube root denesting'),
    ('pf.log_ratio', 'exp-log', 'zero', 'log(sqrt(2) + sqrt(3))/log(5 + 2*sqrt(6)) - 1/2', None, _NB, ''),
    ('pf.issue38', 'algebraic', 'zero', '(a**50 - a**51*a**(-1)).subs(a, 417/(962*pi + 80808))', None, _CI + '38', 'large powers of a rational function of pi'),
    ('pf.issue24', 'exp-log', 'zero', 'exp((log(2*I) - pi*I/2)/2) - sqrt(2)', None, _CI + '24', ''),
    ('pf.issue33', 'trig', 'minpoly:256*x**12 - 768*x**10 + 864*x**8 - 592*x**6 + 321*x**4 - 90*x**2 + 1', 'cos(acos(sqrt(2) - sqrt(3))/3)', None, _CI + '33', 'cos(acos(a)/3) is algebraic'),
    ('pf.issue23a', 'exp-log', 'zero', 'sqrt(exp(2*sqrt(2)) + exp(-2*sqrt(2)) - 2) - (exp(2*sqrt(2)) - 1)/sqrt(exp(2*sqrt(2)))', None, _CI + '23', 'square root of a square over Q(exp(2 sqrt 2))'),
    ('pf.issue23b', 'algebraic', 'zero', '(pi - 1)/(sqrt(pi) - 1) - sqrt(pi) - 1', None, _CI + '23', ''),
    ('pf.issue25a', 'cyclotomic', 'zero', '((-1 + sqrt(-3))/2)**3 - 1', None, _CI + '25', ''),
    ('pf.issue25b', 'branch', 'zero', 'sqrt(-3) - I*sqrt(3)', None, _CI + '25', ''),
    ('pf.sage_golden', 'algebraic', 'zero', '((1 + sqrt(5))/2)**2 - (1 + sqrt(5))/2 - 1', None, _SQ, ''),
    ('pf.sage_nested', 'denest', 'zero', '(sqrt(5 + 2*sqrt(6)) - sqrt(3))**2 - 2', None, _SQ, ''),
    ('pf.sage_sqrt_frac', 'algebraic', 'zero', 'sqrt(2/3)*sqrt(3/5) - sqrt(2/5)', None, _SQ, ''),
    ('pf.sage_cuberoot', 'branch', 'zero', '((-8)**(1/3))**3 + 8 + Abs((-8)**(1/3)) - 2', None, _SQ, ''),
    ('pf.sage_unity', 'cyclotomic', 'zero', '(-1/2 + I*sqrt(3)/2)**3 - 1', None, _SQ, ''),
    ('pf.sage_pow5', 'algebraic', 'zero', '(sqrt(2) + sqrt(3))**5 - 109*sqrt(2) - 89*sqrt(3)', None, _SQ, ''),
    ('pf.sage_pythag', 'parts', 'zero', 'Abs(3/5 + 4*I/5) - 1', None, _SQ, ''),
]


# ---------------------------------------------------------------------------
# families
# ---------------------------------------------------------------------------

def _primes(n):
    ps, k = [], 2
    while len(ps) < n:
        if all(k % p for p in ps):
            ps.append(k)
        k += 1
    return ps

def gauss_sum_cases(ps=(5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 1009)):
    """sum_{k=1}^{p-1} (k|p) cos(2 pi k/p) = sqrt(p) for a prime p = 1 mod 4"""
    out = []
    for p in ps:
        qr = {(k * k) % p for k in range(1, p)}
        terms = " + ".join(("" if k in qr else "-") + "cos(2*pi*%d/%d)" % (k, p) for k in range(1, p))
        out.append(("fam.gauss%d" % p, "family-cyclotomic", "zero", "(%s) - sqrt(%d)" % (terms, p), None,
                    "Gauss (1805)", "p-1 cosines collapse to sqrt(p)"))
    return out

def sine_product_cases(ns=(7, 12, 30, 60, 105, 210)):
    """prod_{k=1}^{n-1} 2 sin(k pi/n) = n"""
    return [("fam.sinprod%d" % n, "family-cyclotomic", "zero",
             "*".join("2*sin(%d*pi/%d)" % (k, n) for k in range(1, n)) + " - %d" % n, None,
             "classical", "product over all conjugates") for n in ns]

def dft_cases(Ns=(2, 3, 4, 5, 6, 8, 10, 12, 16, 20), kind="log"):
    """component 0 of x - IDFT(DFT(x)) for the inputs of the exact DFT
    benchmark (https://fredrikj.net/blog/2020/09/benchmarking-exact-dft-computation/)"""
    xs = {"int": "%d", "sqrt": "sqrt(%d)", "log": "log(%d)", "root": "exp(2*pi*I/%d)",
          "pi": "1/(1+%d*pi)", "sqrtpi": "1/(1+sqrt(%d)*pi)"}[kind]
    out = []
    for N in Ns:
        w = "exp(-2*pi*I/%d)" % N
        X = ["(" + " + ".join("(%s)*%s**%d" % (xs % (n + 2), w, (n * k) % N) for n in range(N)) + ")" for k in range(N)]
        out.append(("fam.dft_%s_%d" % (kind, N), "family-dft", "zero",
                    "%s - (%s)/%d" % (xs % 2, " + ".join(X), N), None,
                    "Johansson, DFT benchmark (2020)", "cyclotomic field of degree phi(N) with %s" % kind))
    return out

def _swinnerton_dyer(n):
    """coefficients (ascending) of S_n = prod (x +- sqrt 2 +- ... +- sqrt p_n)"""
    P = [0, 1]
    for p in _primes(n):
        # P(x - y) P(x + y) with y^2 = p: P(x + y) = A(x) + y B(x)
        d = len(P) - 1
        A, B = [0] * (d + 1), [0] * (d + 1)
        for k, c in enumerate(P):
            for j in range(k + 1):   # c (x + y)^k
                t = c * math.comb(k, j) * p ** (j // 2)
                if j % 2 == 0:
                    A[k - j] += t
                else:
                    B[k - j] += t
        # A^2 - p B^2
        R = [0] * (2 * d + 1)
        for i in range(d + 1):
            for j in range(d + 1):
                R[i + j] += A[i] * A[j] - p * B[i] * B[j]
        while len(R) > 1 and R[-1] == 0:
            R.pop()
        P = R
    return P

def swinnerton_dyer_cases(ns=(2, 3, 4, 5, 6)):
    """S_n(sqrt 2 + ... + sqrt p_n) = 0 (degree 2^n)"""
    out = []
    for n in ns:
        alpha = "(" + " + ".join("sqrt(%d)" % p for p in _primes(n)) + ")"
        P = _swinnerton_dyer(n)
        expr = " + ".join("%d*%s**%d" % (c, alpha, k) for k, c in enumerate(P) if c != 0)
        out.append(("fam.sd%d" % n, "family-multiquadratic", "zero", expr, None,
                    "Swinnerton-Dyer polynomials", "degree-%d integer polynomial at a multiquadratic number" % 2 ** n))
    return out

def family_cases():
    cs = gauss_sum_cases() + sine_product_cases() + swinnerton_dyer_cases()
    for kind in ("int", "sqrt", "log", "root", "pi", "sqrtpi"):
        cs += dft_cases(kind=kind)
    return cs


# ---------------------------------------------------------------------------
# big expressions
# ---------------------------------------------------------------------------

def _huge_expr(name):
    import os, re
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "examples", "huge_expr.c")
    src = open(path).read()
    m = re.search(r"const char \* EXPR_%s =\s*((?:\"[^\"]*\"\s*)+);" % name, src)
    return "".join(re.findall(r'"([^"]*)"', m.group(1))).replace("^", "**")

def gauss_period_cos(p, m, shared=True):
    """cos(2 pi m/p) for a Fermat prime p as nested square roots, by Gauss's
    construction: the periods of 2^L terms (sums of zeta^(g^k) over k in a
    residue class mod 2^L, g a primitive root) split in pairs, the roots of
    x^2 - s x + P with s a period of the previous level and P an integer
    combination of periods of the previous level. With shared, each period
    and each square root is a name bound once by .subs (one square root
    per pair: (p - 1)/2 - 1 of them); otherwise everything is written out
    (exponential growth: 643 kB with 15712 square roots for p = 257)."""
    import cmath, re
    n = p - 1
    assert n & (n - 1) == 0, "a Fermat prime"
    g = next(a for a in range(2, p) if pow(a, n // 2, p) != 1)
    ind, x = {}, 1
    for k in range(n):
        ind[x] = k
        x = x * g % p
    zeta = lambda a: cmath.exp(2j * cmath.pi * a / p)
    def val(L, i):
        return sum(zeta(pow(g, k, p)) for k in range(i, n, 2 ** L)).real
    Lmax = (n // 2).bit_length() - 1        # periods of two terms: 2 cos
    defs = {}
    for L in range(1, Lmax + 1):
        h, size = 2 ** (L - 1), n // 2 ** (L - 1)
        for i in range(h):
            A = [pow(g, k, p) for k in range(i, n, 2 ** L)]
            B = [pow(g, k, p) for k in range(i + h, n, 2 ** L)]
            coef, const = {}, 0
            for a in A:
                for b in B:
                    t = (a + b) % p
                    if t == 0:
                        const += 1
                    else:
                        c = ind[t] % h
                        coef[c] = coef.get(c, 0) + 1
            terms = ["%d" % const] if const else []
            for c, k in sorted(coef.items()):
                assert k % size == 0
                terms.append("%d*e%d_%d" % (k // size, L - 1, c))
            P = "(" + " + ".join(terms) + ")"
            s = "e%d_%d" % (L - 1, i)
            v1, v2 = val(L, i), val(L, i + h)
            sigma = "+" if v1 > v2 else "-"         # (s + sigma sqrt(D))/2 = v1
            other = "-" if sigma == "+" else "+"
            D = "(%s**2 - 4*%s)" % (s, P)
            if shared:
                defs[(L, i)] = "((%s %s q%d_%d)/2)" % (s, sigma, L, i)
                defs[(L, i + h)] = "((%s %s q%d_%d)/2)" % (s, other, L, i)
                defs[("q", L, i)] = "sqrt(%s)" % D
            else:
                defs[(L, i)] = "((%s %s sqrt(%s))/2)" % (s, sigma, D)
                defs[(L, i + h)] = "((%s %s sqrt(%s))/2)" % (s, other, D)
    c = ind[m % p] % (2 ** Lmax)
    if shared:
        e = "(e%d_%d/2)" % (Lmax, c)
        for L in range(Lmax, 0, -1):
            for i in range(2 ** L):
                e += ".subs(e%d_%d, %s)" % (L, i, defs[(L, i)])
            for i in range(2 ** (L - 1)):
                e += ".subs(q%d_%d, %s)" % (L, i, defs[("q", L, i)])
        return e + ".subs(e0_0, -1)"
    memo = {(0, 0): "(-1)"}
    def expand(L, i):
        if (L, i) not in memo:
            memo[(L, i)] = re.sub(r"e(\d+)_(\d+)", lambda mm: expand(int(mm.group(1)), int(mm.group(2))), defs[(L, i)])
        return memo[(L, i)]
    for L in range(1, Lmax + 1):        # (bottom up: no deep recursion)
        for i in range(2 ** L):
            expand(L, i)
    return "(%s/2)" % expand(Lmax, c)

_cos257_cache = []

def _cos257():
    # cos(pi/257) = -cos(2 pi 128/257)
    if not _cos257_cache:
        _cos257_cache.append("(-%s)" % gauss_period_cos(257, 128, shared=False))
    return _cos257_cache[0]

def big_cases():
    """Gauss's square root forms built here (no SymPy needed)"""
    return [
        ("tr.cos17_gauss", "cyclotomic", "zero", gauss_period_cos(17, 1) + " - cos(2*pi/17)", None,
         "Gauss (1796), periods bound by name", "7 square roots, each period named once"),
        ("tr.cos257_shared", "cyclotomic", "zero", "-" + gauss_period_cos(257, 128) + " - cos(pi/257)", None,
         "Gauss's construction for the Fermat prime 257, periods bound by name",
         "STRESS: 127 square roots in 7 levels (15 kB), degree 128"),
    ]

def _expand_big(expr):
    if "BIG:HUGE_N" in expr or "BIG:HUGE_M" in expr:
        expr = expr.replace("BIG:HUGE_N", "(" + _huge_expr("N") + ")")
        expr = expr.replace("BIG:HUGE_M", "(" + _huge_expr("M") + ")")
    if "BIG:COS_PI_257_SQRT" in expr:
        expr = expr.replace("BIG:COS_PI_257_SQRT", "(" + _cos257() + ")")
    return expr


# ---------------------------------------------------------------------------
# evaluation in a field
# ---------------------------------------------------------------------------

_FUNCS = {"exp": "exp", "log": "log", "sin": "sin", "cos": "cos", "tan": "tan", "atan": "atan",
          "asin": "asin", "acos": "acos", "sinh": "sinh", "cosh": "cosh", "tanh": "tanh",
          "asinh": "asinh", "acosh": "acosh", "atanh": "atanh", "abs": "abs", "Abs": "abs",
          "arg": "arg", "re": "re", "im": "im", "conjugate": "conj", "floor": "floor",
          "ceiling": "ceil", "erf": "erf", "erfc": "erfc", "erfi": "erfi", "gamma": "gamma",
          "sqrt": "sqrt", "sign": "sgn", "digamma": "digamma"}

_RELATIONS = {"Eq": lambda a, b: a == b, "Ne": lambda a, b: a != b, "Lt": lambda a, b: a < b,
              "Le": lambda a, b: a <= b, "Gt": lambda a, b: a > b, "Ge": lambda a, b: a >= b}

# functions whose value at an algebraic number is (in general) transcendental
_TRANSCENDENTAL = {"exp", "log", "sin", "cos", "tan", "atan", "asin", "acos", "sinh", "cosh",
                   "tanh", "asinh", "acosh", "atanh", "arg", "erf", "erfc", "erfi", "gamma",
                   "digamma", "polygamma", "polylog", "zeta", "LambertW"}

class Unsupported(Exception):
    pass

def _poly_coeffs(tree, var="x"):
    """ascending rational coefficients of a polynomial expression in var"""
    def add(a, b, s):
        r = dict(a)
        for k, v in b.items():
            r[k] = r.get(k, 0) + s * v
        return r
    def mul(a, b):
        r = {}
        for i, u in a.items():
            for j, v in b.items():
                r[i + j] = r.get(i + j, 0) + u * v
        return r
    def ev(n):
        if isinstance(n, ast.Constant):
            return {0: Fraction(n.value)}
        if isinstance(n, ast.Name) and n.id == var:
            return {1: Fraction(1)}
        if isinstance(n, ast.UnaryOp):
            return {k: (-v if isinstance(n.op, ast.USub) else v) for k, v in ev(n.operand).items()}
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                a, e = ev(n.left), ev(n.right)
                r = {0: Fraction(1)}
                for _ in range(int(e[0])):
                    r = mul(r, a)
                return r
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): return add(a, b, 1)
            if isinstance(n.op, ast.Sub): return add(a, b, -1)
            if isinstance(n.op, ast.Mult): return mul(a, b)
            if isinstance(n.op, ast.Div): return {k: v / b[0] for k, v in a.items()}
        raise Unsupported("polynomial")
    d = ev(tree)
    return [d.get(k, Fraction(0)) for k in range(max(d) + 1)]

def _rational(node):
    """the rational value of a constant subexpression built from integers, or None"""
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return Fraction(node.value)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        v = _rational(node.operand)
        return None if v is None else (-v if isinstance(node.op, ast.USub) else v)
    if isinstance(node, ast.BinOp):
        a, b = _rational(node.left), _rational(node.right)
        if a is None or b is None:
            return None
        if isinstance(node.op, ast.Add): return a + b
        if isinstance(node.op, ast.Sub): return a - b
        if isinstance(node.op, ast.Mult): return a * b
        if isinstance(node.op, ast.Div): return a / b if b != 0 else None
        if isinstance(node.op, ast.Pow) and b.denominator == 1 and abs(b) < 100000:
            return a ** int(b) if (a != 0 or b >= 0) else None
    return None

def _pi_multiple(node):
    """(q, imaginary) if node is syntactically q*pi or q*pi*I for a rational q"""
    def mono(n):
        # {(power of pi, power of I): rational coefficient}
        r = _rational(n)
        if r is not None:
            return {(0, 0): r}
        if isinstance(n, ast.Name) and n.id in ("pi", "I"):
            return {(1, 0) if n.id == "pi" else (0, 1): Fraction(1)}
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub):
            m = mono(n.operand)
            return None if m is None else {k: -v for k, v in m.items()}
        if isinstance(n, ast.BinOp):
            a, b = mono(n.left), mono(n.right)
            if a is None or b is None:
                return None
            if isinstance(n.op, (ast.Add, ast.Sub)):
                s = 1 if isinstance(n.op, ast.Add) else -1
                r = dict(a)
                for k, v in b.items():
                    r[k] = r.get(k, 0) + s * v
                return {k: v for k, v in r.items() if v != 0}
            if isinstance(n.op, ast.Mult):
                r = {}
                for (p1, i1), u in a.items():
                    for (p2, i2), v in b.items():
                        k = (p1 + p2, (i1 + i2) % 4)
                        r[k] = r.get(k, 0) + u * v
                out = {}
                for (p, i), v in r.items():   # I^2 = -1
                    k, s = (p, i % 2), (-1 if i >= 2 else 1)
                    out[k] = out.get(k, 0) + s * v
                return {k: v for k, v in out.items() if v != 0}
            if isinstance(n.op, ast.Div) and list(b) == [(0, 0)] and b[(0, 0)] != 0:
                return {k: v / b[(0, 0)] for k, v in a.items()}
        return None
    m = mono(node)
    if m is not None and len(m) == 1:
        (p, i), q = next(iter(m.items()))
        if p == 1:
            return q, bool(i)
    return None

def is_algebraic(s):
    """syntactic test: no transcendental constants or functions (sin, cos, tan
    at rational multiples of pi and exp at rational multiples of pi*I allowed)"""
    stack = [ast.parse(_expand_big(s), mode="eval").body]
    while stack:
        n = stack.pop()
        if isinstance(n, ast.Name) and n.id in ("pi", "E", "Catalan", "EulerGamma"):
            return False
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Pow) and _rational(n.right) is None:
            return False
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name):
            fn = n.func.id
            if fn in ("sin", "cos", "tan", "exp") and len(n.args) == 1:
                pm = _pi_multiple(n.args[0])
                if pm is not None and pm[1] == (fn == "exp"):
                    continue
            if fn in _TRANSCENDENTAL:
                return False
        stack.extend(ast.iter_child_nodes(n))
    return True

def _exact_sort(roots):
    """roots in CRootOf order, using exact comparisons in the field: the real
    roots ascending, then the others by real part, then imaginary part"""
    real = [z == z.conj() for z in roots]
    def cmp(a, b):
        za, ra = a
        zb, rb = b
        if ra != rb:
            return -1 if ra else 1
        if ra:
            return -1 if za < zb else (1 if za > zb else 0)
        xa, xb = za.re(), zb.re()
        if xa != xb:
            return -1 if xa < xb else 1
        ya, yb = za.im(), zb.im()
        return -1 if ya < yb else (1 if ya > yb else 0)
    return [z for z, r in sorted(zip(roots, real), key=functools.cmp_to_key(cmp))]

def evaluate(s, R, rational_pi=False):
    """the value of the expression s (Python syntax) in the field R: a field
    element, or a bool for a relation. With rational_pi, sin, cos, tan at
    rational multiples of pi and exp at rational multiples of pi*I are
    computed as algebraic numbers (sin_pi, ...; for qqbar)."""
    tree = ast.parse(_expand_big(s), mode="eval").body
    memo = {}

    def num(node):
        # exact rational value of a constant subexpression, or None
        k = id(node)
        if k not in memo:
            memo[k] = _num(node)
        return memo[k]

    def _num(node):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ("fibonacci", "factorial"):
            n = num(node.args[0])
            if n is None or n.denominator != 1:
                return None
            n = int(n)
            if node.func.id == "factorial":
                return Fraction(math.factorial(n))
            a, b = 0, 1
            for _ in range(n):
                a, b = b, a + b
            return Fraction(a)
        if isinstance(node, ast.BinOp):
            a, b = num(node.left), num(node.right)
            if a is None or b is None:
                return None
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Div): return a / b if b != 0 else None
            if isinstance(node.op, ast.Pow) and b.denominator == 1 and abs(b) < 100000:
                return a ** int(b) if (a != 0 or b >= 0) else None
            return None
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            v = num(node.operand)
            return None if v is None else -v
        return _rational(node) if isinstance(node, ast.Constant) else None

    def frac(f):
        return R(f.numerator) / f.denominator if f.denominator != 1 else R(f.numerator)

    def ev(node, sub):
        if not sub:
            f = num(node)
            if f is not None:
                return frac(f)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return R(node.value)
        if isinstance(node, ast.BinOp):
            if isinstance(node.op, ast.Pow):
                e = None if sub else num(node.right)
                b = ev(node.left, sub)
                if e is not None:
                    if e.denominator == 1:
                        return b ** int(e)
                    if e == Fraction(1, 2):
                        return b.sqrt()
                    return b ** frac(e)
                return b ** ev(node.right, sub)
            a, b = ev(node.left, sub), ev(node.right, sub)
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(node.op, ast.Div): return a / b
            raise Unsupported(type(node.op).__name__)
        if isinstance(node, ast.UnaryOp):
            v = ev(node.operand, sub)
            if isinstance(node.op, ast.USub): return -v
            if isinstance(node.op, ast.UAdd): return v
            raise Unsupported("unary operator")
        if isinstance(node, ast.Name):
            n = node.id
            if n in sub: return sub[n]
            if n == "pi": return R.pi()
            if n == "I": return R.i()
            if n == "E": return R(1).exp()
            if n == "GoldenRatio": return (1 + R(5).sqrt()) / 2
            if n == "Catalan": return R.catalan()
            if n == "EulerGamma": return R.euler()
            raise Unsupported("name " + n)
        if isinstance(node, ast.Call):
            if isinstance(node.func, ast.Attribute) and node.func.attr == "subs":
                var, val = node.args
                return ev(node.func.value, dict(sub, **{var.id: ev(val, sub)}))
            fn = node.func.id
            args = node.args
            if fn in _RELATIONS:
                return _RELATIONS[fn](ev(args[0], sub), ev(args[1], sub))
            if rational_pi and fn in ("sin", "cos", "tan", "exp") and not sub:
                pm = _pi_multiple(args[0])
                if pm is not None and pm[1] == (fn == "exp"):
                    q = frac(pm[0])
                    return {"sin": R.sin_pi, "cos": R.cos_pi, "tan": R.tan_pi, "exp": R.exp_pi_i}[fn](q)
            if fn in _FUNCS:
                return getattr(ev(args[0], sub), _FUNCS[fn])()
            if fn == "expand":
                return ev(args[0], sub)
            if fn == "real_root":
                b = ev(args[0], sub); n = int(num(args[1]))
                if n % 2 == 1 and b < 0:
                    return -((-b) ** (R(1) / n))
                return b ** (R(1) / n)
            if fn == "cbrt":
                return ev(args[0], sub) ** (R(1) / 3)
            if fn == "LambertW":
                if len(args) == 1:
                    return R.lambertw(ev(args[0], sub))
                return R.lambertw(ev(args[0], sub), int(num(args[1])))
            if fn == "zeta":
                if len(args) == 1:
                    return R.zeta(ev(args[0], sub))
                return R.hurwitz_zeta(ev(args[0], sub), ev(args[1], sub))
            if fn == "polylog":
                return R.polylog(ev(args[0], sub), ev(args[1], sub))
            if fn == "polygamma":
                return R.polygamma(ev(args[0], sub), ev(args[1], sub))
            if fn == "RootOf":
                return rootof(args[0], int(num(args[1])))
            raise Unsupported("function " + fn)
        raise Unsupported(type(node).__name__)

    def rootof(ptree, k):
        from flint_ctypes import PolynomialRing_gr_poly
        c = _poly_coeffs(ptree)
        den = math.lcm(*[x.denominator for x in c])
        P = PolynomialRing_gr_poly(R)([int(x * den) for x in c])
        rts = list(P.roots()[0])
        # (the real roots come first: the nonreal ones are sorted only
        # when the index is beyond them)
        real = [z for z in rts if z == z.conj()]
        if k < len(real):
            return _exact_sort(real)[k]
        return _exact_sort(rts)[k]

    return ev(tree, {})

def check(case, R, rational_pi=False):
    """True if the field gets the case right, False if it gets it wrong
    (an exception if it cannot decide)"""
    cid, cat, expect, expr, approx = case[:5]
    x = evaluate(expr, R, rational_pi)
    if expect in ("true", "false"):
        if not isinstance(x, bool):
            raise TypeError("a relation was expected")
        return x == (expect == "true")
    if expect == "zero":
        return (x == 0) is True
    if expect == "nonzero":
        return (x == 0) is False
    P = _poly_coeffs(ast.parse(expect[8:], mode="eval").body)
    y = R(0)
    for c in reversed(P):
        y = y * x + (R(c.numerator) / c.denominator)
    return (y == 0) is True


# ---------------------------------------------------------------------------
# program cases: matrices and polynomial roots
# ---------------------------------------------------------------------------

def _poly_ring(R):
    from flint_ctypes import PolynomialRing_gr_poly
    return PolynomialRing_gr_poly(R)

def _roots(p):
    """the roots of p with multiplicities, in exact CRootOf order"""
    rts, mults = p.roots()
    rts, mults = list(rts), [int(m) for m in mults]
    order = _exact_sort(rts)
    return order, [mults[[i for i, r in enumerate(rts) if r is z][0]] for z in order]

def _select(roots, approx, R, radius="1/1000"):
    """the unique root within radius of approx (a decimal string), by exact
    comparisons |r - approx| < radius in the field"""
    a = R(Fraction(approx).numerator) / Fraction(approx).denominator
    rad = R(Fraction(radius).numerator) / Fraction(radius).denominator
    found = [r for r in roots if abs(r - a) < rad]
    if len(found) != 1:
        raise ValueError("%d roots near %s" % (len(found), approx))
    return found[0]

def prog_cayley_hamilton(R):
    from flint_ctypes import Mat
    A = Mat(R, 2, 2)([[5, R.pi()], [1, -1]]) ** 4
    return A.charpoly()(A) == Mat(R, 2, 2)([[0, 0], [0, 0]])

def prog_sage37927(R):
    from flint_ctypes import Mat
    I = R.i(); v1 = -I; v2 = -R(2).sqrt()
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
    return (M * v).is_zero() is True and X.ncols() == 1

def prog_matexp_log(R):
    from flint_ctypes import Mat
    M = Mat(R)([[1, -1], [-1, -1]])
    return M.exp().log() == M

def prog_roots_product(R):
    x = _poly_ring(R).gen()
    h = (x - R(3).sqrt()) ** 2 * (2 + R(2).sqrt() * x + x ** 2)
    rts, mults = _roots(h)
    p = _poly_ring(R)([1])
    for r, m in zip(rts, mults):
        p = p * (x - r) ** m
    return p == h and mults == [2, 1, 1] and rts[0] == R(3).sqrt()

def prog_x5(R):
    x = _poly_ring(R).gen()
    rts, _ = _roots(x ** 5 - x - 1)
    return (len(rts) == 5 and all(r ** 5 - r - 1 == 0 for r in rts) and sum(rts, R(0)) == 0
            and rts[0] > R(11673) / 10000 and rts[0] < R(11674) / 10000)

def prog_disc(R):
    P = _poly_ring(R)
    x = P.gen()
    def disc1(b, c, d): return b**2*c**2 - 4*b**3*d - 4*c**3 + 18*b*c*d - 27*d**2
    def disc2(s1, s2, s3): return ((s1 - s2)*(s1 - s3)*(s2 - s3))**2
    ok = True
    for p in [x*(x - 2)*(x - 4), x*(x - 2)*(x - 4) + 1,
              (x - R(2).sqrt())*(x - R(2)**(R(1)/3))*(x - R(3).sqrt())]:
        d, c, b = p[0], p[1], p[2]
        rts, _ = _roots(p)
        ok = ok and (disc1(b, c, d) == disc2(rts[0], rts[1], rts[2]))
    return ok

def _34gon(R):
    sqrt = lambda v: R(v).sqrt() if not hasattr(v, "sqrt") else v.sqrt()
    rt17 = sqrt(17); rt2 = sqrt(2)
    eps = sqrt(17 + rt17); epss = sqrt(17 - rt17)
    delta = rt17 - 1
    alpha = sqrt(34 + 6*rt17 + rt2*delta*epss - 8*rt2*eps)
    x = rt2*sqrt(15 + rt17 + rt2*(alpha + epss))/8
    y = rt2*sqrt(epss**2 - rt2*(alpha + epss))/8
    return x, y

def prog_34gon(R):
    x, y = _34gon(R)
    cx, cy = R(1), R(0)
    for i in range(34):
        cx, cy = x*cx - y*cy, x*cy + y*cx
    return cx == 1 and cy == 0

def prog_34gon_vs_poly(R):
    x, y = _34gon(R)
    t = _poly_ring(R).gen()
    p = 256*t**8 - 128*t**7 - 448*t**6 + 192*t**5 + 240*t**4 - 80*t**3 - 40*t**2 + 8*t + 1
    x2 = _select(list(p.roots()[0]), "0.98297", R)
    y2 = (1 - x2**2).sqrt()
    return x == x2 and y == y2

def prog_arprec(R):
    t = _poly_ring(R).gen()
    p = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1
    a = _select(list(p.roots()[0]), "1.17628", R)
    lhs = a**630 - 1
    num = (a**315 - 1) * (a**210 - 1) * (a**126 - 1)**2 * (a**90 - 1) * (a**3 - 1)**3 * (a**2 - 1)**5 * (a - 1)**3
    den = (a**35 - 1) * (a**15 - 1)**2 * (a**14 - 1)**2 * (a**5 - 1)**6 * a**68
    return lhs == num / den

def prog_roots_pi_sqrt2(R):
    t = _poly_ring(R).gen()
    pi, s2 = R.pi(), R(2).sqrt()
    rts, _ = _roots(t**2 - (pi + s2)*t + pi*s2)
    return len(rts) == 2 and rts[0] == s2 and rts[1] == pi

def prog_roots_pi_pm_sqrt2(R):
    t = _poly_ring(R).gen()
    pi = R.pi()
    rts, _ = _roots(t**2 - 2*pi*t + pi**2 - 2)
    return (len(rts) == 2 and rts[1] - rts[0] == 2*R(8).sqrt()/2
            and all((r - pi)**2 == 2 for r in rts) and rts[0] > 1)

def prog_roots_e_sqrt2_i(R):
    t = _poly_ring(R).gen()
    e, s2, i = R(1).exp(), R(2).sqrt(), R.i()
    rts, _ = _roots((t - e)*(t - s2)*(t - i))
    return len(rts) == 3 and rts[0] == s2 and rts[1] == e and rts[2] == i

def prog_quintic_power_sums(R):
    t = _poly_ring(R).gen()
    pi = R.pi()
    rts = list((t**5 - pi*t - 1).roots()[0])     # (no ordering needed)
    p = [sum((r**k for r in rts), R(0)) for k in range(1, 6)]
    return len(rts) == 5 and p[0] == 0 and p[1] == 0 and p[2] == 0 and p[3] == 4*pi and p[4] == 5

def prog_quintic_real_roots(R):
    # x^5 - pi x - 1 has three real roots, the smallest in (-3/2, -1)
    t = _poly_ring(R).gen()
    rts, _ = _roots(t**5 - R.pi()*t - 1)
    nreal = sum(r == r.conj() for r in rts)
    return nreal == 3 and rts[0] > -R(3)/2 and rts[0] < -1

def prog_cubic_sqrt_pi(R):
    # the roots of x^3 - sqrt(pi) x^2 - x + sqrt(pi) are sqrt(pi), -1, 1
    t = _poly_ring(R).gen()
    sp = R.pi().sqrt()
    rts, _ = _roots(t**3 - sp*t**2 - t + sp)
    return len(rts) == 3 and rts[0] == -1 and rts[1] == 1 and rts[2] == sp

_PS = "gr_tower profile"
# (id, category, algebraic, function, source, what it stresses, LaTeX)
PROGRAMS = [
    ("prog.cayley_hamilton", "matrix", False, prog_cayley_hamilton, _NB,
     "characteristic polynomial over Q(pi)", r"\chi_A(A) = 0,\; A = \begin{pmatrix} 5 & \pi \\ 1 & -1 \end{pmatrix}^{4}"),
    ("prog.sage37927", "matrix", True, prog_sage37927, "Sage trac #37927",
     "nullspace of a 14 x 10 matrix over Q(i, sqrt 2)", r"\dim \ker M = 1,\; M \in \mathbb{Q}(i, \sqrt{2})^{14 \times 10}"),
    ("prog.matexp_log", "matrix", False, prog_matexp_log, _CI + "23",
     "matrix exponential and logarithm", r"\log(\exp(M)) = M,\; M = \begin{pmatrix} 1 & -1 \\ -1 & -1 \end{pmatrix}"),
    ("prog.roots_product", "roots", True, prog_roots_product, _FT + " (qqbar roots)",
     "roots with multiplicities over Q(sqrt 2, sqrt 3)", r"\prod_k (x - r_k)^{m_k} = (x-\sqrt{3})^2 (x^2 + \sqrt{2} x + 2)"),
    ("prog.sage_x5", "roots", True, prog_x5, _SQ,
     "the roots of x^5 - x - 1", r"r^5 - r - 1 = 0,\; \textstyle\sum r = 0,\; 1.1673 < r_0 < 1.1674"),
    ("prog.sage_disc", "roots", True, prog_disc, _SQ,
     "cubic discriminants from coefficients and from roots", r"\operatorname{disc}(p) = \prod_{i<j} (r_i - r_j)^2"),
    ("prog.sage_34gon", "roots", True, prog_34gon, _SQ,
     "rotate 34 times by a nested radical angle", r"(x + i y)^{34} = 1"),
    ("prog.sage_34gon_vs_poly", "roots", True, prog_34gon_vs_poly, _SQ + " ('currently infinitely long')",
     "nested radical against a root of a degree-8 polynomial, selected by |r - 0.98297| < 1/1000",
     r"x = r,\; |r - 0.98297| < 10^{-3},\; 256 r^8 - 128 r^7 - \cdots + 1 = 0"),
    ("prog.sage_arprec", "roots", True, prog_arprec, _SQ + " (ARPREC)",
     "a degree-10 unit to the power 630, selected by |a - 1.17628| < 1/1000",
     r"\alpha^{630} - 1 = \frac{(\alpha^{315}-1)(\alpha^{210}-1)(\alpha^{126}-1)^2 \cdots}{(\alpha^{35}-1)(\alpha^{15}-1)^2 \cdots \alpha^{68}}"),
    ("prog.roots_pi_sqrt2", "roots", False, prog_roots_pi_sqrt2, _PS,
     "roots of a quadratic over Q(pi, sqrt 2): the discriminant is a square", r"x^2 - (\pi + \sqrt{2}) x + \pi \sqrt{2} = 0 \;\Rightarrow\; x \in \{\sqrt{2}, \pi\}"),
    ("prog.roots_pi_pm_sqrt2", "roots", False, prog_roots_pi_pm_sqrt2, _PS,
     "roots over Q(pi): pi +- sqrt 2", r"x^2 - 2 \pi x + \pi^2 - 2 = 0 \;\Rightarrow\; x = \pi \pm \sqrt{2}"),
    ("prog.roots_e_sqrt2_i", "roots", False, prog_roots_e_sqrt2_i, _PS,
     "roots of a cubic over Q(e, sqrt 2, i), real and complex", r"(x - e)(x - \sqrt{2})(x - i) \;\Rightarrow\; \{\sqrt{2}, e, i\}"),
    ("prog.cubic_sqrt_pi", "roots", False, prog_cubic_sqrt_pi, _PS,
     "a cubic over Q(sqrt(pi)) with rational roots", r"x^3 - \sqrt{\pi} x^2 - x + \sqrt{\pi} = 0 \;\Rightarrow\; x \in \{-1, 1, \sqrt{\pi}\}"),
    ("prog.quintic_power_sums", "roots", False, prog_quintic_power_sums, _PS,
     "power sums of the roots of an irreducible quintic over Q(pi)", r"r^5 - \pi r - 1 = 0:\; \textstyle\sum r^4 = 4\pi,\; \sum r^5 = 5"),
    ("prog.quintic_real_roots", "roots", False, prog_quintic_real_roots, _PS,
     "real roots of an irreducible quintic over Q(pi): conjugation and ordering of implicit algebraic numbers",
     r"\#\{r \in \mathbb{R} : r^5 - \pi r - 1 = 0\} = 3,\; -\tfrac{3}{2} < r_0 < -1"),
]


# ---------------------------------------------------------------------------
# Cardano and Ferrari: radical root formulas of cubics and quartics
# ---------------------------------------------------------------------------

_W3 = ["1", "(-1 + sqrt(3)*I)/2", "(-1 - sqrt(3)*I)/2"]

def cardano(p, q, k):
    """the k-th Cardano root (k = 0, 1, 2) of y^3 + p y + q (p != 0), as an
    expression: u w - p/(3 u w) with u the principal cube root of
    -q/2 + sqrt(q^2/4 + p^3/27) and w = 1, (-1 + sqrt(-3))/2, (-1 - sqrt(-3))/2"""
    return "(u*%s - (%s)/(3*u*%s)).subs(u, (-(%s)/2 + sqrt((%s)**2/4 + (%s)**3/27))**(1/3))" % (
        _W3[k], p, _W3[k], q, q, p)

def cardano_general(b, c, d, k):
    """the k-th root of y^3 + b y^2 + c y + d, through y = z - b/3"""
    p = "(%s) - (%s)**2/3" % (c, b)
    q = "2*(%s)**3/27 - (%s)*(%s)/3 + (%s)" % (b, b, c, d)
    return "(%s - (%s)/3)" % (cardano(p, q, k), b)

def ferrari(p, q, r, k, mk=0):
    """the k-th root (k = 0..3) of x^4 + p x^2 + q x + r (q != 0), with m
    the mk-th Cardano root of the resolvent m^3 + p m^2 + (p^2/4 - r) m - q^2/8:
    with s = sqrt(2m), the roots of x^2 -+ s x + p/2 + m +- q/(2s)"""
    m = cardano_general(p, "(%s)**2/4 - (%s)" % (p, r), "-(%s)**2/8" % q, mk)
    s1, s2 = [("", "+"), ("", "-"), ("-", "+"), ("-", "-")][k]
    inner = "-2*(%s) - 2*m %s 2*(%s)/s" % (p, "-" if s1 == "" else "+", q)
    return "((%ss %s sqrt(%s))/2).subs(s, sqrt(2*m)).subs(m, %s)" % (s1, s2, inner, m)

# (name, cubic y^3 + p y + q or quartic x^4 + p x^2 + q x + r, coefficients,
#  closed forms of the formula roots in order, or None); the order was
#  determined numerically (50 digits)
RADICAL_POLYS = [
    ('c_rat', 'cubic', ('-7', '6'), ['2', '-3', '1']),
    ('c_cos9', 'cubic', ('-3', '1'), ['2*cos(2*pi/9)', '2*cos(8*pi/9)', '2*cos(4*pi/9)']),
    ('c_cbrt', 'cubic', ('-6', '-6'), ['2**(1/3) + 4**(1/3)', '2**(1/3)*(-1 - sqrt(3)*I)/2 + 4**(1/3)*(-1 + sqrt(3)*I)/2', '2**(1/3)*(-1 + sqrt(3)*I)/2 + 4**(1/3)*(-1 - sqrt(3)*I)/2']),
    ('c_sqrt2', 'cubic', ('-sqrt(2)', '-1'), None),
    ('c_gauss', 'cubic', ('1 + I', 'sqrt(3)'), None),
    ('c_pi', 'cubic', ('-pi', '-1'), None),
    ('c_e', 'cubic', ('1', '-E'), None),
    ('c_double', 'cubic', ('-3*pi**2', '2*pi**3'), ['pi', '-2*pi', 'pi']),
    ('c_pi_sqrt2', 'cubic', ('-pi**2 - sqrt(2)*pi - 2', 'sqrt(2)*pi**2 + 2*pi'), ['pi', '-pi - sqrt(2)', 'sqrt(2)']),
    ('q_rat', 'quartic', ('-25', '60', '-36'), ['3', '2', '1', '-6']),
    ('q_generic', 'quartic', ('0', '1', '1'), None),
    ('q_alg', 'quartic', ('-sqrt(2)', '1', 'sqrt(3)'), None),
    ('q_sqrts', 'quartic', ('-6 - sqrt(6) - sqrt(3) - sqrt(2)', '2*sqrt(6) + 5 + 3*sqrt(3) + 4*sqrt(2)', '-3*sqrt(2) - 2*sqrt(3) - sqrt(6)'), ['sqrt(3)', 'sqrt(2)', '1', '-1 - sqrt(2) - sqrt(3)']),
    ('q_pi', 'quartic', ('0', '-pi', '-1'), None),
    ('q_pi_roots', 'quartic', ('-pi**2 - 3*pi - 7', '6 + 9*pi + 3*pi**2', '2*pi*(-pi - 3)'), ['pi', '2', '1', '-pi - 3']),
]

_RAD_NOTES = {
    'c_rat': 'casus irreducibilis with rational roots: complex cube roots must cancel',
    'c_cos9': 'casus irreducibilis: the roots are 2 cos(2 pi k/9)',
    'c_cbrt': 'one real root cbrt(2) + cbrt(4)',
    'c_sqrt2': 'coefficients in Q(sqrt 2)',
    'c_gauss': 'complex algebraic coefficients',
    'c_pi': 'three real roots over Q(pi): radicals of complex numbers over Q(pi)',
    'c_e': 'coefficients in Q(e)',
    'c_double': 'a double root: zero discriminant, cube root of -pi^3',
    'c_pi_sqrt2': 'casus irreducibilis over Q(pi, sqrt 2) with roots pi, sqrt 2, -pi - sqrt 2',
    'q_rat': 'rational roots through the resolvent cubic',
    'q_generic': 'Galois group S4',
    'q_alg': 'coefficients in Q(sqrt 2, sqrt 3)',
    'q_sqrts': 'roots 1, sqrt 2, sqrt 3, -1 - sqrt 2 - sqrt 3: denesting through the resolvent',
    'q_pi': 'coefficients in Q(pi)',
    'q_pi_roots': 'roots pi, 1, 2, -pi - 3: denesting over Q(pi)',
}

def _rad_poly(kind, c):
    if kind == "cubic":
        return "x**3 + (%s)*x + (%s)" % c
    return "x**4 + (%s)*x**2 + (%s)*x + (%s)" % c

def _rad_roots(kind, c):
    if kind == "cubic":
        return [cardano(c[0], c[1], k) for k in range(3)]
    return [ferrari(c[0], c[1], c[2], k) for k in range(4)]

def _subs_roots(expr, roots):
    for i, r in enumerate(roots):
        expr = "(%s).subs(r%d, %s)" % (expr, i, r)
    return expr

def radical_formula_cases():
    """expression cases: each Cardano/Ferrari root is a root of the
    polynomial; Vieta's formulas; the roots equal their closed forms"""
    src = "gr_tower profile (Cardano, Ferrari)"
    out = []
    for name, kind, c, known in RADICAL_POLYS:
        P, rs, n = _rad_poly(kind, c), _rad_roots(kind, c), (3 if kind == "cubic" else 4)
        note = _RAD_NOTES[name]
        what = "Cardano" if kind == "cubic" else "Ferrari"
        try:
            ptex = to_latex(P)
        except Exception:
            ptex = None
        for k in ((0, 1) if kind == "cubic" else (0, 2)):
            tex = None if ptex is None else r"p(r_%d) = 0,\; p = %s,\; r_%d \text{ by %s}" % (k, ptex, k, what)
            out.append(("rad.%s.poly%d" % (name, k), "radicals", "zero", "(%s).subs(x, %s)" % (P, rs[k]),
                        None, src, "%s root %d of the polynomial; %s" % (what, k, note), tex))
        def neg(t):
            return t[1:] if t.startswith("-") and t[1:].isdigit() else ("-" + t if t.isdigit() else "-(%s)" % t)
        lt = lambda t: to_latex(t) or "?"
        if kind == "cubic":
            e2 = ("r0*r1 + r0*r2 + r1*r2 - (%s)" % c[0], r"r_0 r_1 + r_0 r_2 + r_1 r_2 = %s" % lt(c[0]))
            en = ("r0*r1*r2 + (%s)" % c[1], r"r_0 r_1 r_2 = %s" % lt(neg(c[1])))
        else:
            e2 = ("r0*r1 + r0*r2 + r0*r3 + r1*r2 + r1*r3 + r2*r3 - (%s)" % c[0], r"e_2(r_0, \ldots, r_3) = %s" % lt(c[0]))
            en = ("r0*r1*r2*r3 - (%s)" % c[2], r"r_0 r_1 r_2 r_3 = %s" % lt(c[2]))
        for tag, (e, t) in (("vieta2", e2), ("vieta%d" % n, en)):
            tex = None if ptex is None else r"%s,\; p = %s \text{ (%s roots)}" % (t, ptex, what)
            out.append(("rad.%s.%s" % (name, tag), "radicals", "zero", _subs_roots(e, rs), None, src,
                        "Vieta's formula for the %s roots; %s" % (what, note), tex))
        if known:
            for k, kv in enumerate(known):
                try:
                    ktex = to_latex(kv)
                except Exception:
                    ktex = None
                tex = None if ptex is None or ktex is None else r"r_%d = %s,\; p = %s \text{ (%s)}" % (k, ktex, ptex, what)
                out.append(("rad.%s.root%d" % (name, k), "radicals", "zero", "%s - (%s)" % (rs[k], kv), None, src,
                            "%s root %d equals its closed form; %s" % (what, k, note), tex))
    return out

def _radical_program(kind, c):
    def f(R):
        rational_pi = isinstance(R, __import__("flint_ctypes").ComplexAlgebraicField_qqbar)
        cs = [evaluate(t, R, rational_pi) for t in c]
        t = _poly_ring(R).gen()
        p = t**3 + cs[0]*t + cs[1] if kind == "cubic" else t**4 + cs[0]*t**2 + cs[1]*t + cs[2]
        tower = list(p.roots()[0])
        formula = [evaluate(r, R, rational_pi) for r in _rad_roots(kind, c)]
        # each formula root equals exactly one root found by root finding,
        # with the multiplicities matching
        mults = [int(m) for m in p.roots()[1]]
        used = [0] * len(tower)
        for z in formula:
            hits = [i for i, w in enumerate(tower) if z == w]
            if len(hits) != 1:
                return False
            used[hits[0]] += 1
        return used == mults
    return f

def radical_programs():
    out = []
    for name, kind, c, known in RADICAL_POLYS:
        P = _rad_poly(kind, c)
        try:
            ptex = to_latex(P)
        except Exception:
            ptex = None
        what = "Cardano" if kind == "cubic" else "Ferrari"
        tex = None if ptex is None else r"\{r_k \text{ by %s}\} = \operatorname{roots}(p),\; p = %s" % (what, ptex)
        out.append(("prog.rad_%s" % name, "radicals", is_algebraic(" + ".join(c)), _radical_program(kind, c),
                    "gr_tower profile (Cardano, Ferrari)",
                    "the %s roots against the roots found by root finding (exact equality, multiplicities); %s"
                    % (what, _RAD_NOTES[name]), tex, _sympy_radical_program(kind, c)))
    return out


# ---------------------------------------------------------------------------
# matrix functions (gr_mat_exp, gr_mat_log, gr_mat_sqrt, gr_mat_pow_scalar:
# through the Jordan form, gr_mat_func_jordan / gr_mat_func_param_jordan)
# ---------------------------------------------------------------------------

# (name, entries as expressions, what it is); the exp(log) and log(exp)
# identities need: invertible (log), eigenvalues with |Im| < pi (log exp)
MATRICES = [
    ("jordan2", [["2", "1"], ["0", "2"]], "a Jordan block"),
    ("jordan3", [["3", "1", "0"], ["0", "3", "1"], ["0", "0", "3"]], "a 3 x 3 Jordan block"),
    ("rot", [["0", "-1"], ["1", "0"]], "eigenvalues +-i"),
    ("quad", [["1", "2"], ["3", "4"]], "eigenvalues (5 +- sqrt 33)/2, one negative"),
    ("cubic", [["0", "0", "1"], ["1", "0", "1"], ["0", "1", "0"]], "companion of x^3 - x - 1 (one real, two complex eigenvalues)"),
    ("spd4", [["2", "1", "0", "1"], ["1", "3", "1", "0"], ["0", "1", "4", "1"], ["1", "0", "1", "5"]],
     "symmetric positive definite 4 x 4, eigenvalues of degree 4"),
    ("tridiag5", [["4", "1", "0", "0", "0"], ["1", "4", "1", "0", "0"], ["0", "1", "4", "1", "0"],
                  ["0", "0", "1", "4", "1"], ["0", "0", "0", "1", "4"]], "eigenvalues 4 + 2 cos(k pi/6)"),
    ("block4", [["2", "1", "0", "0"], ["0", "2", "0", "0"], ["0", "0", "1", "-1"], ["0", "0", "1", "1"]],
     "a Jordan block and a rotation-scaling block"),
    ("pi_sqrt2", [["pi", "sqrt(2)"], ["1", "I"]], "entries pi, sqrt 2, i"),
    ("trans3", [["1", "pi", "0"], ["0", "1", "E"], ["0", "0", "2"]], "defective, transcendental entries"),
    ("gauss", [["1 + I", "2"], ["sqrt(3)", "-I"]], "complex algebraic entries"),
    ("exp_entries", [["exp(1)", "1"], ["log(2)", "exp(-1)"]], "entries e, 1/e, log 2"),
]

_MAT_IDENTITIES = {
    "explog": (r"\exp(\log M) = M", "exp(log M) = M", False),
    "logexp": (r"\log(\exp M) = M", "log(exp M) = M (eigenvalues with |Im| < pi)", False),
    "det": (r"\det \exp M = e^{\operatorname{tr} M}", "det exp(M) = exp(trace M)", False),
    "inv": (r"\exp(M) \exp(-M) = I", "exp(M) exp(-M) = I", False),
    "double": (r"\exp(2M) = \exp(M)^2", "exp(2M) = exp(M)^2", False),
    "commute": (r"\exp(M + M^2) = \exp(M) \exp(M^2)", "exp(M + M^2) = exp(M) exp(M^2) (commuting)", False),
    "sqrt": (r"\sqrt{M}^{\,2} = M", "the principal square root squared", True),
    "cbrt": (r"(M^{1/3})^3 = M", "the principal cube root cubed", True),
    "powsum": (r"M^{1/2} M^{1/3} = M^{5/6}", "rational powers add", True),
    "powpi": (r"M^{\pi} M^{-\pi} = I", "transcendental powers", False),
}

_MAT_TESTS = {
    "jordan2": ["explog", "logexp", "det", "inv", "double", "commute", "sqrt", "cbrt", "powsum", "powpi"],
    "jordan3": ["explog", "logexp", "det", "inv", "sqrt", "cbrt", "powsum"],
    "rot": ["explog", "logexp", "det", "inv", "double", "sqrt", "cbrt"],
    "quad": ["explog", "logexp", "det", "inv", "sqrt", "cbrt", "powsum"],
    "cubic": ["explog", "logexp", "det", "inv", "sqrt", "cbrt"],
    "spd4": ["explog", "logexp", "det", "inv", "sqrt"],
    "tridiag5": ["explog", "logexp", "det", "sqrt"],
    "block4": ["explog", "logexp", "det", "inv", "double", "sqrt", "cbrt", "powsum", "powpi"],
    "pi_sqrt2": ["explog", "logexp", "det", "inv", "sqrt"],
    "trans3": ["explog", "logexp", "det", "inv", "commute", "sqrt", "powpi"],
    "gauss": ["explog", "logexp", "det", "inv", "sqrt", "cbrt"],
    "exp_entries": ["explog", "det", "inv", "sqrt"],
}

def _mat_eval(rows, R):
    from flint_ctypes import Mat
    rational_pi = isinstance(R, __import__("flint_ctypes").ComplexAlgebraicField_qqbar)
    return Mat(R)([[evaluate(e, R, rational_pi) for e in row] for row in rows])

def _mat_identity(test, rows):
    def f(R):
        from flint_ctypes import Mat
        M = _mat_eval(rows, R)
        n = len(rows)
        I = Mat(R)([[1 if i == j else 0 for j in range(n)] for i in range(n)])
        if test == "explog":
            return M.log().exp() == M
        if test == "logexp":
            return M.exp().log() == M
        if test == "det":
            return M.exp().det() == M.trace().exp()
        if test == "inv":
            return M.exp() * (-M).exp() == I
        if test == "double":
            return (2 * M).exp() == M.exp() ** 2
        if test == "commute":
            return (M + M * M).exp() == M.exp() * (M * M).exp()
        if test == "sqrt":
            return M.sqrt() ** 2 == M
        if test == "cbrt":
            return (M ** (R(1) / 3)) ** 3 == M
        if test == "powsum":
            return M ** (R(1) / 2) * M ** (R(1) / 3) == M ** (R(5) / 6)
        if test == "powpi":
            return M ** R.pi() * M ** (-R.pi()) == I
        raise ValueError(test)
    return f

def _sympy_mat_zero(D):
    out = True
    for e in D:
        z = _sympy_zero(e)
        if z is False:
            return False
        if z is None:
            out = None
    return out

def _sympy_mat_identity(test, rows):
    def f():
        import sympy as S
        M = S.Matrix([[S.sympify(e) for e in row] for row in rows])
        n = M.shape[0]
        I = S.eye(n)
        if test == "explog":
            D = M.log().exp() - M
        elif test == "logexp":
            D = M.exp().log() - M
        elif test == "det":
            return _sympy_zero(M.exp().det() - S.exp(M.trace()))
        elif test == "inv":
            D = M.exp() * (-M).exp() - I
        elif test == "double":
            D = (2 * M).exp() - M.exp() ** 2
        elif test == "commute":
            D = (M + M * M).exp() - M.exp() * (M * M).exp()
        elif test == "sqrt":
            D = (M ** S.Rational(1, 2)) ** 2 - M
        elif test == "cbrt":
            D = (M ** S.Rational(1, 3)) ** 3 - M
        elif test == "powsum":
            D = M ** S.Rational(1, 2) * M ** S.Rational(1, 3) - M ** S.Rational(5, 6)
        elif test == "powpi":
            D = M ** S.pi * M ** (-S.pi) - I
        return _sympy_mat_zero(D)
    return f

def _rotation_program(theta):
    def f(R):
        from flint_ctypes import Mat
        t = evaluate(theta, R)
        J = Mat(R)([[0, -1], [1, 0]])
        c, s = t.cos(), t.sin()
        return (t * J).exp() == Mat(R)([[c, -s], [s, c]])
    def g():
        import sympy as S
        t = S.sympify(theta)
        D = (t * S.Matrix([[0, -1], [1, 0]])).exp() - S.Matrix([[S.cos(t), -S.sin(t)], [S.sin(t), S.cos(t)]])
        return _sympy_mat_zero(D)
    return f, g

def matrix_programs():
    src = "gr_tower profile (matrix functions)"
    out = []
    for name, rows, what in MATRICES:
        alg = all(is_algebraic(e) for row in rows for e in row)
        try:
            mtex = r"\begin{pmatrix} %s \end{pmatrix}" % r" \\ ".join(" & ".join(to_latex(e) or "?" for e in row) for row in rows)
        except Exception:
            mtex = "M"
        for test in _MAT_TESTS[name]:
            tex, note, algebraic_ok = _MAT_IDENTITIES[test]
            out.append(("prog.mat_%s_%s" % (name, test), "matrix", alg and algebraic_ok, _mat_identity(test, rows), src,
                        "%s; %s" % (note, what), r"%s,\; M = %s" % (tex, mtex), _sympy_mat_identity(test, rows)))
    for theta, tname in (("1", "1"), ("sqrt(2)", "sqrt2"), ("2*pi/7", "2pi7"), ("pi + 1/3", "pi13")):
        f, g = _rotation_program(theta)
        out.append(("prog.mat_rotation_%s" % tname, "matrix", False, f, src,
                    "the exponential of a rotation generator against cos and sin",
                    r"\exp\begin{pmatrix} 0 & -\theta \\ \theta & 0 \end{pmatrix} = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix},\; \theta = %s" % (to_latex(theta) or theta), g))
    return out


# ---------------------------------------------------------------------------
# SymPy versions of program cases (True/False for a right/wrong answer,
# None if SymPy cannot decide)
# ---------------------------------------------------------------------------

def _sympy_all(zs):
    out = True
    for z in zs:
        if z is False:
            return False
        if z is None:
            out = None
    return out

def _sympy_match(formula, roots):
    """each formula root equals exactly one root (by equals), multiset-wise"""
    roots = list(roots)
    unknown = False
    for z in formula:
        hit = None
        for j, w in enumerate(roots):
            e = _sympy_zero(z - w)
            if e is True:
                hit = j
                break
            if e is None:
                unknown = True
        if hit is None:
            return None if unknown else False
        roots.pop(hit)
    return True

def _sympy_roots_list(P, x):
    """the roots of P (with multiplicity) by SymPy: CRootOf for rational
    polynomials, otherwise roots() (radical formulas), None if incomplete"""
    import sympy as S
    poly = S.Poly(P, x)
    if all(c.is_Rational for c in poly.all_coeffs()):
        return poly.all_roots()
    rs = S.roots(poly)
    if sum(rs.values()) != poly.degree():
        return None
    return [r for r, m in rs.items() for _ in range(m)]

def sympy_prog_roots_pi_sqrt2():
    import sympy as S
    x = S.Symbol("x")
    rs = _sympy_roots_list(x**2 - (S.pi + S.sqrt(2))*x + S.pi*S.sqrt(2), x)
    return None if rs is None else _sympy_match([S.sqrt(2), S.pi], rs)

def sympy_prog_roots_pi_pm_sqrt2():
    import sympy as S
    x = S.Symbol("x")
    rs = _sympy_roots_list(x**2 - 2*S.pi*x + S.pi**2 - 2, x)
    return None if rs is None else _sympy_match([S.pi - S.sqrt(2), S.pi + S.sqrt(2)], rs)

def sympy_prog_roots_e_sqrt2_i():
    import sympy as S
    x = S.Symbol("x")
    rs = _sympy_roots_list(S.expand((x - S.E)*(x - S.sqrt(2))*(x - S.I)), x)
    return None if rs is None else _sympy_match([S.sqrt(2), S.E, S.I], rs)

def sympy_prog_cubic_sqrt_pi():
    import sympy as S
    x = S.Symbol("x")
    sp = S.sqrt(S.pi)
    rs = _sympy_roots_list(x**3 - sp*x**2 - x + sp, x)
    return None if rs is None else _sympy_match([-1, 1, sp], rs)

def sympy_prog_x5():
    import sympy as S
    x = S.Symbol("x")
    rs = S.Poly(x**5 - x - 1, x).all_roots()
    ok = [_sympy_zero(r**5 - r - 1) for r in rs] + [_sympy_zero(S.Add(*rs))]
    real = S.Gt(rs[0], S.Rational(11673, 10000)) and S.Lt(rs[0], S.Rational(11674, 10000))
    return _sympy_all(ok + [bool(real) if real in (S.true, S.false) else None])

def sympy_prog_roots_product():
    import sympy as S
    x = S.Symbol("x")
    h = S.expand((x - S.sqrt(3))**2 * (2 + S.sqrt(2)*x + x**2))
    rs = _sympy_roots_list(h, x)
    if rs is None:
        return None
    return _sympy_zero(S.expand(S.prod([x - r for r in rs]) - h)) if len(rs) == 4 else False

def sympy_prog_disc():
    import sympy as S
    x = S.Symbol("x")
    out = []
    for p in [x*(x - 2)*(x - 4), x*(x - 2)*(x - 4) + 1, (x - S.sqrt(2))*(x - S.cbrt(2))*(x - S.sqrt(3))]:
        p = S.Poly(S.expand(p), x)
        d, c, b = [p.coeff_monomial(x**k) for k in range(3)]
        rs = _sympy_roots_list(p.as_expr(), x)
        if rs is None or len(rs) != 3:
            return None
        d1 = b**2*c**2 - 4*b**3*d - 4*c**3 + 18*b*c*d - 27*d**2
        d2 = ((rs[0] - rs[1])*(rs[0] - rs[2])*(rs[1] - rs[2]))**2
        out.append(_sympy_zero(d1 - d2))
    return _sympy_all(out)

def sympy_prog_34gon():
    import sympy as S
    sq = S.sqrt
    rt17, rt2 = sq(17), sq(2)
    eps, epss = sq(17 + rt17), sq(17 - rt17)
    alpha = sq(34 + 6*rt17 + rt2*(rt17 - 1)*epss - 8*rt2*eps)
    x = rt2*sq(15 + rt17 + rt2*(alpha + epss))/8
    y = rt2*sq(epss**2 - rt2*(alpha + epss))/8
    z = S.expand((x + S.I*y)**34)
    return _sympy_zero(z - 1)

def sympy_prog_34gon_vs_poly():
    import sympy as S
    t = S.Symbol("t")
    sq = S.sqrt
    rt17, rt2 = sq(17), sq(2)
    eps, epss = sq(17 + rt17), sq(17 - rt17)
    alpha = sq(34 + 6*rt17 + rt2*(rt17 - 1)*epss - 8*rt2*eps)
    x = rt2*sq(15 + rt17 + rt2*(alpha + epss))/8
    p = 256*t**8 - 128*t**7 - 448*t**6 + 192*t**5 + 240*t**4 - 80*t**3 - 40*t**2 + 8*t + 1
    near = [r for r in S.Poly(p, t).all_roots() if S.Abs(r - S.Rational(98297, 100000)) < S.Rational(1, 1000)]
    if len(near) != 1:
        return False
    return _sympy_zero(x - near[0])

def sympy_prog_arprec():
    import sympy as S
    t = S.Symbol("t")
    p = t**10 + t**9 - t**7 - t**6 - t**5 - t**4 - t**3 + t + 1
    near = [r for r in S.Poly(p, t).all_roots() if S.Abs(r - S.Rational(117628, 100000)) < S.Rational(1, 1000)]
    if len(near) != 1:
        return False
    a = near[0]
    num = (a**315 - 1) * (a**210 - 1) * (a**126 - 1)**2 * (a**90 - 1) * (a**3 - 1)**3 * (a**2 - 1)**5 * (a - 1)**3
    den = (a**35 - 1) * (a**15 - 1)**2 * (a**14 - 1)**2 * (a**5 - 1)**6 * a**68
    return _sympy_zero(a**630 - 1 - num / den)

def sympy_prog_cayley_hamilton():
    import sympy as S
    A = S.Matrix([[5, S.pi], [1, -1]])**4
    lam = S.Symbol("lam")
    cp = A.charpoly(lam)
    C = S.zeros(2, 2)
    for k, c in enumerate(reversed(cp.all_coeffs())):
        C += c * A**k
    return _sympy_mat_zero(C)

def sympy_prog_matexp_log():
    import sympy as S
    M = S.Matrix([[1, -1], [-1, -1]])
    return _sympy_mat_zero(M.exp().log() - M)

def sympy_prog_sage37927():
    import sympy as S
    I = S.I; v1 = -I; v2 = -S.sqrt(2)
    M = S.Matrix([[0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
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
    X = M.nullspace(iszerofunc=lambda e: _sympy_zero(e))
    v = S.Matrix([-108, 0, 0, 1, 0, 12, 0, -60, 0, 64])
    return _sympy_all([_sympy_mat_zero(M * v), len(X) == 1])

def _sympy_radical_program(kind, c):
    def f():
        import sympy as S
        x = S.Symbol("x")
        P = S.sympify(_rad_poly(kind, c))
        formula = [S.sympify(r) for r in _rad_roots(kind, c)]
        rs = _sympy_roots_list(S.expand(P), x)
        return None if rs is None else _sympy_match(formula, rs)
    return f

SYMPY_PROGRAMS = {
    "prog.roots_pi_sqrt2": sympy_prog_roots_pi_sqrt2,
    "prog.roots_pi_pm_sqrt2": sympy_prog_roots_pi_pm_sqrt2,
    "prog.roots_e_sqrt2_i": sympy_prog_roots_e_sqrt2_i,
    "prog.cubic_sqrt_pi": sympy_prog_cubic_sqrt_pi,
    "prog.sage_x5": sympy_prog_x5,
    "prog.roots_product": sympy_prog_roots_product,
    "prog.sage_disc": sympy_prog_disc,
    "prog.sage_34gon": sympy_prog_34gon,
    "prog.sage_34gon_vs_poly": sympy_prog_34gon_vs_poly,
    "prog.sage_arprec": sympy_prog_arprec,
    "prog.cayley_hamilton": sympy_prog_cayley_hamilton,
    "prog.matexp_log": sympy_prog_matexp_log,
    "prog.sage37927": sympy_prog_sage37927,
}


# ---------------------------------------------------------------------------
# benchmarks from algebraic number theory, real algebraic geometry and
# translation surfaces (qqbar-checkable unless noted)
# ---------------------------------------------------------------------------

def _kronecker(D, a):
    """the Kronecker symbol (D/a) for a > 0"""
    res = 1
    while a % 2 == 0:
        a //= 2
        if D % 2 == 0:
            return 0
        if D % 8 in (3, 5):
            res = -res
    D %= a
    t, n = 1, a
    while D:
        while D % 2 == 0:
            D //= 2
            if n % 8 in (3, 5):
                t = -t
        D, n = n, D
        if D % 4 == 3 and n % 4 == 3:
            t = -t
        D %= n
    return res * t if n == 1 else 0

def _fundamental_unit(D):
    """the fundamental unit of the real quadratic field of discriminant D
    (as an expression), by search (small D only)"""
    d = D if D % 4 == 1 else D // 4
    y = 1
    while True:
        for s in ((-4, 4) if D % 4 == 1 else (-1, 1)):
            x2 = d * y * y + s
            x = math.isqrt(x2) if x2 >= 0 else -1
            if x >= 0 and x * x == x2:
                return ("(%d + %d*sqrt(%d))/2" % (x, y, d)) if D % 4 == 1 else ("(%d + %d*sqrt(%d))" % (x, y, d))
        y += 1

def _primitive_root(p):
    fs = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    return next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in fs))

# discriminant: class number (checked numerically to 50 digits)
_CLASS_NUMBERS = {5: 1, 8: 1, 12: 1, 13: 1, 17: 1, 21: 1, 24: 1, 28: 1, 29: 1, 40: 2, 60: 2, 65: 2,
                  136: 2, 229: 3, 316: 3, 401: 5}

def number_theory_cases():
    out = []
    # the analytic class number formula for real quadratic fields:
    # eps^h = prod_{0 < a < D/2, (a, D) = 1} sin(pi a/D)^(-chi(a))
    src = "Dirichlet's class number formula (real quadratic fields)"
    for D, h in sorted(_CLASS_NUMBERS.items()):
        eps = _fundamental_unit(D)
        num, den = [], []
        for a in range(1, (D + 1) // 2):
            if math.gcd(a, D) == 1:
                c = _kronecker(D, a)
                (den if c == 1 else num).append("sin(%d*pi/%d)" % (a, D))
        prod = "(%s)/(%s)" % ("*".join(num) or "1", "*".join(den) or "1")
        out.append(("nt.cnf%d" % D, "number theory", "zero", "%s - (%s)**%d" % (prod, eps, h), None, src,
                    "a quotient of %d sines in Q(zeta_%d) is the unit %s to the power h = %d" % (len(num) + len(den), 4 * D, eps, h)))
        if D in (5, 12, 13, 40, 229):
            terms = []
            for a in range(1, (D + 1) // 2):
                if math.gcd(a, D) == 1:
                    terms.append(("+" if _kronecker(D, a) == 1 else "-") + "log(sin(%d*pi/%d))" % (a, D))
            out.append(("nt.cnf%d_log" % D, "number theory", "zero", "%s + %d*log(%s)" % (" ".join(terms).lstrip("+"), h, eps), None,
                        src, "the logarithmic form: multiplicative relations among cyclotomic units (transcendental)"))
    # cubic Gauss and Jacobi sums: G(chi)^3 = p J(chi, chi)
    src = "Gauss/Jacobi sums of cubic characters (Ireland-Rosen, ch. 8)"
    for p in (7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 151, 307):
        g = _primitive_root(p)
        ind, x = {}, 1
        for k in range(p - 1):
            ind[x] = k
            x = x * g % p
        G = " + ".join("exp(2*pi*I*%d/%d)" % ((p * (ind[k] % 3) + 3 * k) % (3 * p), 3 * p) for k in range(1, p))
        c = [0, 0, 0]
        for k in range(2, p):
            c[(ind[k] + ind[(1 - k) % p]) % 3] += 1
        A, B = c[0] - c[2], c[1] - c[2]
        out.append(("nt.cubic_gauss%d" % p, "number theory", "zero",
                    "(%s)**3 - %d*(%d + %d*exp(2*pi*I/3))" % (G, p, A, B), None, src,
                    "a sum of %d roots of unity of order %d cubed is p times a Eisenstein integer" % (p - 1, 3 * p)))
    # Gauss's cubic period polynomial: x^3 + x^2 - (p-1)/3 x - (p(L+3) - 1)/27
    src = "Gauss's cubic period polynomial (4p = L^2 + 27 M^2, L = 1 mod 3)"
    for p in (7, 13, 31, 103, 307, 997, 9973):
        g = _primitive_root(p)
        L = next(L for L in range(-2 * math.isqrt(p) - 3, 2 * math.isqrt(p) + 3)
                 if L % 3 == 1 and (4 * p - L * L) > 0 and (4 * p - L * L) % 27 == 0 and math.isqrt((4 * p - L * L) // 27) ** 2 == (4 * p - L * L) // 27)
        H = sorted({pow(g, 3 * k, p) for k in range((p - 1) // 3)})
        eta = " + ".join("2*cos(2*pi*%d/%d)" % (a, p) for a in H if a < p / 2)
        poly = "x**3 + x**2 - %d*x - %d" % ((p - 1) // 3, (p * (L + 3) - 1) // 27)
        out.append(("nt.cubic_period%d" % p, "number theory", "minpoly:" + poly, eta, None, src,
                    "a period of %d cosines is a root of a cubic" % len([a for a in H if a < p / 2])))
    # quadratic Gauss sums sum_k exp(2 pi i k^2/n)
    src = "Gauss's evaluation of quadratic Gauss sums"
    for n in (100, 101, 102, 103, 1000, 1001):
        s = " + ".join("exp(2*pi*I*%d/%d)" % (k * k % n, n) for k in range(n))
        val = ["(1 + I)*sqrt(%d)", "sqrt(%d)", "0", "I*sqrt(%d)"][n % 4]
        val = val % n if "%d" in val else val
        out.append(("nt.quad_gauss%d" % n, "number theory", "zero", "%s - %s" % (s, val), None, src,
                    "%d roots of unity, n = %d mod 4" % (n, n % 4)))
    # norms of cyclotomic units: prod_{(k, n) = 1} (1 - zeta_n^k)
    src = "norms of 1 - zeta_n (Washington, Introduction to cyclotomic fields, Prop. 2.8)"
    for n, v in ((343, 7), (512, 2), (210, 1), (2310, 1)):
        prod = "*".join("(1 - exp(2*pi*I*%d/%d))" % (k, n) for k in range(1, n) if math.gcd(k, n) == 1)
        out.append(("nt.cyclo_norm%d" % n, "number theory", "zero", "%s - %d" % (prod, v), None, src,
                    "a product of phi(%d) = %d factors" % (n, len([k for k in range(1, n) if math.gcd(k, n) == 1]))))
    # Sage QQbar doctests
    src = _SQ
    out.append(("sage.re_quintic", "algebraic",
                "minpoly:x**10 + 3/16*x**6 + 11/32*x**5 - 1/64*x**2 + 1/128*x - 1/1024",
                "re(RootOf(x**5 - x - 1, 4))", None, src + " (long time)", "the real part of a nonreal quintic root: degree 10"))
    out.append(("sage.im_quintic", "algebraic",
                "minpoly:x**20 - 5/8*x**16 - 95/256*x**12 - 625/1024*x**10 - 5/512*x**8 - 1875/8192*x**6 + 25/4096*x**4 - 625/32768*x**2 + 2869/1048576",
                "im(RootOf(x**5 - x - 1, 4))", None, src + " (long time, 10 s in 2013)", "the imaginary part: degree 20"))
    w = ["1", "exp(2*pi*I/3)", "exp(4*pi*I/3)", "1"] + ["exp(2*pi*I*%d/5)" % k for k in range(1, 5)]
    P = {k: "(" + " + ".join("(%s)**%d" % (x, k) for x in w) + ")" for k in range(1, 5)}
    schur = ("((P1**3 + 3*P1*P2 + 2*P3)/6*(P1**2 + P2)/2 - (P1**4 + 6*P1**2*P2 + 3*P2**2 + 8*P1*P3 + 6*P4)/24*P1)"
             ".subs(P1, %s).subs(P2, %s).subs(P3, %s).subs(P4, %s)" % (P[1], P[2], P[3], P[4]))
    out.append(("sage.schur", "cyclotomic", "zero", schur, None, src + " (long time)",
                "the Schur polynomial s_(3,2) of 8 variables at roots of unity in Q(zeta_15) (Jacobi-Trudi)"))
    # real algebraic geometry: Mignotte's polynomials x^d - 2(101 x - 1)^2,
    # two real roots at distance sqrt(2) 101^(-d/2-1) (to 20 digits)
    src = "Mignotte polynomials (root separation; cf. Tsigaridas-Emiris, cs/0604066)"
    for d in (20, 50, 100):
        M = "x**%d - 2*(101*x - 1)**2" % d
        sep = "(RootOf(%s, 2) - RootOf(%s, 1))" % (M, M)
        ref = "sqrt(2)/101**%d" % (d // 2 + 1)
        out.append(("ra.mignotte%d_lower" % d, "real algebraic", "true", "Gt(%s, 99/100*%s)" % (sep, ref), None, src,
                    "the two close real roots of a degree-%d polynomial: separation %.1e" % (d, 2 ** 0.5 * 101.0 ** (-(d // 2 + 1)))))
        out.append(("ra.mignotte%d_upper" % d, "real algebraic", "true", "Lt(%s, 101/100*%s)" % (sep, ref), None, src, "upper bound"))
    # translation surfaces
    out.append(("fs.ngon17", "cyclotomic", "zero",
                "cos(2*pi/17) + (-a**14/2 + 15*a**12/2 - 45*a**10 + 275*a**8/2 - 225*a**6 + 189*a**4 - 70*a**2 + 15/2).subs(a, 2*sin(2*pi/17))",
                None, "sage-flatsurf polygons.regular_ngon(17) doctest", "a vertex of the regular 17-gon in the field Q(2 sin(2 pi/17)) of degree 16"))
    return out

def _veech_program(n):
    def f(R):
        from flint_ctypes import Mat
        rational_pi = isinstance(R, __import__("flint_ctypes").ComplexAlgebraicField_qqbar)
        lam = evaluate("2*cos(pi/%d)" % n, R, rational_pi)
        ST = Mat(R)([[0, -1], [1, 0]]) * Mat(R)([[1, lam], [0, 1]])
        return ST ** n == Mat(R)([[-1, 0], [0, -1]]) and ST ** (n // 2 if n % 2 == 0 else n) != Mat(R)([[1, 0], [0, 1]])
    def g():
        import sympy as S
        lam = 2 * S.cos(S.pi / n)
        ST = S.Matrix([[0, -1], [1, 0]]) * S.Matrix([[1, lam], [0, 1]])
        return _sympy_mat_zero(ST ** n + S.eye(2))
    return f, g

def _quartic_sqrt2_program():
    # sage #18242: a root of x^4 + x^3 + sqrt(2) x + 1 has minimal polynomial
    # x^8 + 2x^7 + x^6 + 2x^4 + 2x^3 - 2x^2 + 1
    def f(R):
        t = _poly_ring(R).gen()
        rts = list((t**4 + t**3 + R(2).sqrt()*t + 1).roots()[0])
        return len(rts) == 4 and all(r**8 + 2*r**7 + r**6 + 2*r**4 + 2*r**3 - 2*r**2 + 1 == 0 for r in rts)
    return f

def _degenerate_sort_program():
    # Sage QQbar doctest: p1((x - 1)^2) has 16 roots, all with real part 1
    def f(R):
        t = _poly_ring(R).gen()
        p1 = lambda u: u**8 + 74*u**7 + 2300*u**6 + 38928*u**5 + 388193*u**4 + 2295312*u**3 + 7613898*u**2 + 12066806*u + 5477001
        rts, _ = _roots(p1((t - 1)**2))
        return len(rts) == 16 and all(r.re() == 1 for r in rts) and all(rts[i].im() < rts[i + 1].im() for i in range(15))
    def g():
        import sympy as S
        x = S.Symbol("x")
        p1 = lambda u: u**8 + 74*u**7 + 2300*u**6 + 38928*u**5 + 388193*u**4 + 2295312*u**3 + 7613898*u**2 + 12066806*u + 5477001
        rs = S.Poly(S.expand(p1((x - 1)**2)), x).all_roots()
        return _sympy_all([_sympy_zero(S.re(r) - 1) for r in rs])
    return f, g

def benchmark_programs():
    out = []
    for n in (5, 7, 8, 9, 12, 17, 30):
        f, g = _veech_program(n)
        out.append(("prog.veech%d" % n, "number theory", True, f, "Veech (1989): the double regular n-gon",
                    "the Veech group relation (ST)^n = -I with T = [[1, 2 cos(pi/%d)], [0, 1]]" % n,
                    r"(ST)^{%d} = -I,\; S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix},\; T = \begin{pmatrix} 1 & 2\cos(\pi/%d) \\ 0 & 1 \end{pmatrix}" % (n, n), g))
    out.append(("prog.sage18242", "roots", True, _quartic_sqrt2_program(), "Sage issue #18242",
                "minimal polynomial of the roots of a quartic over Q(sqrt 2) (31 s in Sage before the fix)",
                r"r^4 + r^3 + \sqrt{2} r + 1 = 0 \Rightarrow r^8 + 2r^7 + r^6 + 2r^4 + 2r^3 - 2r^2 + 1 = 0"))
    f, g = _degenerate_sort_program()
    out.append(("prog.sage_degenerate_sort", "roots", True, f, _SQ + " (long time)",
                "16 roots with equal real parts: sorting needs exact equality of real parts",
                r"p_1((x-1)^2):\; \operatorname{Re} r_k = 1,\; \operatorname{Im} r_0 < \cdots < \operatorname{Im} r_{15}", g))
    return out


# ---------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------

def _case(t):
    if callable(t[3]):
        cid, cat, algebraic, func, source, note, tex = t[:7]
        return dict(id=cid, category=cat, kind="program", expect="true", expr=None, approx=None,
                    source=source, note=note, tex=tex, func=func, algebraic=algebraic,
                    sympy_func=(t[7] if len(t) > 7 else None))
    cid, cat, expect, expr, approx, source, note = t[:7]
    return dict(id=cid, category=cat, kind="expr", expect=expect, expr=expr, approx=approx,
                source=source, note=note, tex=(t[7] if len(t) > 7 else None), func=None, algebraic=None,
                sympy_func=None)

def all_cases(families=True, big=True):
    cs = list(CASES) + radical_formula_cases() + number_theory_cases() + (family_cases() if families else []) + (big_cases() if big else [])
    cs = [_case(c) for c in cs if big or "BIG:" not in c[3]]
    progs = [p if len(p) > 7 or p[0] not in SYMPY_PROGRAMS else tuple(p) + (SYMPY_PROGRAMS[p[0]],)
             for p in PROGRAMS + radical_programs() + matrix_programs() + benchmark_programs()]
    return cs + [_case(p) for p in progs]

ENGINES = ("tower", "ca", "qqbar", "sympy")
STATUSES = ("ok", "wrong", "undecided", "unable", "error", "timeout", "n/a")


# ---------------------------------------------------------------------------
# engines
# ---------------------------------------------------------------------------

def _field(name):
    import flint_ctypes as F
    return {"tower": F.ComplexField_tower, "ca": F.ComplexField_ca, "qqbar": F.ComplexAlgebraicField_qqbar}[name]()

def _sympy_zero(e):
    z = e.is_zero
    if z is None:
        z = e.equals(0)
    return z if isinstance(z, bool) else None

def sympy_check(case):
    """True/False for a right/wrong answer by SymPy, None if undecided"""
    import sympy as S
    loc = {n: S.Symbol(n) for n in ("x", "v", "a")}
    e = S.sympify(case["expr"], locals=loc)
    expect = case["expect"]
    if expect in ("true", "false"):
        want = (expect == "true")
        if e is S.true or e is S.false:
            return bool(e) == want
        if isinstance(e, (S.Eq, S.Ne)):
            z = _sympy_zero(e.lhs - e.rhs)
            if z is None:
                return None
            return (z if isinstance(e, S.Eq) else not z) == want
        return None
    if expect.startswith("minpoly:"):
        e = S.sympify(expect[8:], locals=loc).subs(loc["x"], e)
    z = _sympy_zero(e)
    if z is None:
        return None
    return z if expect != "nonzero" else not z

def _run_case(case, engine):
    """(status, seconds, detail); runs in the subprocess"""
    import flint_ctypes as F
    t = time.time()
    try:
        if engine == "sympy":
            r = case["sympy_func"]() if case["kind"] == "program" else sympy_check(case)
        else:
            R = _field(engine)
            if case["kind"] == "program":
                r = case["func"](R)
            else:
                r = check((case["id"], case["category"], case["expect"], case["expr"], case["approx"]),
                          R, rational_pi=(engine == "qqbar"))
        dt = time.time() - t
        if r is None:
            return ("undecided", dt, "")
        return ("ok" if r is True else "wrong", dt, "" if isinstance(r, bool) else repr(r)[:200])
    except F.Undecidable as ex:
        return ("undecided", time.time() - t, str(ex)[:300])
    except F.FlintUnableError as ex:
        return ("unable", time.time() - t, str(ex)[:300])
    except F.FlintDomainError as ex:
        # qqbar: a transcendental value the syntactic test did not catch
        return ("n/a" if engine == "qqbar" else "error", time.time() - t, str(ex)[:300])
    except NotImplementedError as ex:
        return ("unable", time.time() - t, ("%s: %s" % (type(ex).__name__, ex))[:300])
    except Exception as ex:
        return ("error", time.time() - t, ("%s: %s" % (type(ex).__name__, ex))[:300])

def _worker(case, engine, conn):
    sys.setrecursionlimit(200000)
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    out = []
    threading.stack_size(1 << 29)     # deep recursion in the evaluator (qqbar, big cases)
    th = threading.Thread(target=lambda: out.append(_run_case(case, engine)))
    th.start(); th.join()
    conn.send(out[0] if out else ("error", 0.0, "no result"))
    conn.close()

def applicable(case, engine):
    if engine == "sympy":
        return case["kind"] == "expr" or case.get("sympy_func") is not None
    if engine == "qqbar":
        return case["algebraic"] if case["kind"] == "program" else is_algebraic(case["expr"])
    return True

def _fmt(res):
    st, dt = res["status"], res.get("time")
    if st == "ok":
        return "%.4f" % dt
    if st in ("n/a", "timeout"):
        return st
    return "%s %.3f" % (st.upper() if st == "wrong" else st, dt or 0.0)

def run(cases, engines=("tower", "ca"), timeout=60.0, jobs=1, verbose=True):
    """run the cases on the engines; fills case["results"][engine]"""
    if "sympy" in engines:
        import sympy        # imported once, before forking
    ctx = mp.get_context("fork")
    for c in cases:
        c["results"] = {}
        # (the big expressions are generated here, outside the timings)
        if c["expr"] and "BIG:" in c["expr"]:
            c["expr"] = _expand_big(c["expr"])
    queue = []
    for i, c in enumerate(cases):
        for e in engines:
            if applicable(c, e):
                queue.append((i, e))
            else:
                c["results"][e] = dict(status="n/a", time=None, detail="")
    running, printed = [], 0
    if verbose:
        print("case | " + " | ".join(engines)); sys.stdout.flush()
    while queue or running:
        while queue and len(running) < max(1, jobs):
            i, e = queue.pop(0)
            a, b = ctx.Pipe(duplex=False)
            p = ctx.Process(target=_worker, args=(cases[i], e, b))
            p.start(); b.close()
            running.append((p, a, i, e, time.time()))
        time.sleep(0.002)
        still = []
        for (p, a, i, e, t0) in running:
            res = None
            if a.poll():
                try:
                    st, dt, det = a.recv()
                    res = dict(status=st, time=dt, detail=det)
                except EOFError:
                    res = dict(status="error", time=time.time() - t0, detail="crashed")
                p.join()
            elif not p.is_alive():
                p.join()
                res = dict(status="error", time=time.time() - t0, detail="crashed (exit code %s)" % p.exitcode)
            elif time.time() - t0 > timeout:
                p.kill(); p.join()
                res = dict(status="timeout", time=timeout, detail="")
            if res is None:
                still.append((p, a, i, e, t0))
            else:
                a.close()
                cases[i]["results"][e] = res
        running = still
        while verbose and printed < len(cases) and len(cases[printed]["results"]) == len(engines):
            c = cases[printed]
            print(c["id"] + " | " + " | ".join(_fmt(c["results"][e]) for e in engines))
            sys.stdout.flush()
            printed += 1
    return cases

def summary(cases, engines):
    """per engine: counts by status, total and median time of the ok runs"""
    out = {}
    for e in engines:
        rs = [c["results"][e] for c in cases if e in c["results"]]
        cnt = {s: sum(r["status"] == s for r in rs) for s in STATUSES}
        ts = sorted(r["time"] for r in rs if r["status"] == "ok")
        out[e] = dict(counts=cnt, total=sum(ts), median=(ts[len(ts) // 2] if ts else None),
                      applicable=len(rs) - cnt["n/a"])
    return out

def print_summary(cases, engines):
    S = summary(cases, engines)
    print()
    print("%-8s " % "" + " ".join("%9s" % s for s in STATUSES) + "   time(ok)")
    for e in engines:
        print("%-8s " % e + " ".join("%9d" % S[e]["counts"][s] for s in STATUSES) + "   %.2f s" % S[e]["total"])
    wrong = [(c["id"], e) for c in cases for e in engines if c["results"].get(e, {}).get("status") == "wrong"]
    if wrong:
        print("WRONG:", ", ".join("%s (%s)" % w for w in wrong))


# ---------------------------------------------------------------------------
# LaTeX (via fexpr) and the HTML report
# ---------------------------------------------------------------------------

_FEXPR_FUNCS = {"exp": "Exp", "log": "Log", "sin": "Sin", "cos": "Cos", "tan": "Tan", "atan": "Atan",
                "asin": "Asin", "acos": "Acos", "sinh": "Sinh", "cosh": "Cosh", "tanh": "Tanh",
                "asinh": "Asinh", "acosh": "Acosh", "atanh": "Atanh", "abs": "Abs", "Abs": "Abs",
                "arg": "Arg", "re": "Re", "im": "Im", "conjugate": "Conjugate", "floor": "Floor",
                "ceiling": "Ceil", "erf": "Erf", "erfc": "Erfc", "erfi": "Erfi", "gamma": "Gamma",
                "sqrt": "Sqrt", "sign": "Sign", "digamma": "DigammaFunction", "LambertW": "LambertW",
                "polylog": "PolyLog", "fibonacci": "Fibonacci", "factorial": "Factorial",
                "Eq": "Equal", "Ne": "NotEqual", "Lt": "Less", "Le": "LessEqual", "Gt": "Greater",
                "Ge": "GreaterEqual", "real_root": "RealRoot", "RootOf": "RootOf"}

def to_fexpr(s):
    """the expression s (Python syntax) as an fexpr"""
    from flint_ctypes import fexpr as F
    def c(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, int):
            return F(n.value)
        if isinstance(n, ast.Name):
            return F({"pi": "Pi", "I": "NumberI", "E": "NumberE", "GoldenRatio": "GoldenRatio",
                      "Catalan": "CatalanConstant", "EulerGamma": "ConstGamma"}.get(n.id, n.id))
        if isinstance(n, ast.UnaryOp):
            return F("Neg")(c(n.operand)) if isinstance(n.op, ast.USub) else c(n.operand)
        if isinstance(n, ast.BinOp):
            # a + (-b) and a + (-b)*c print as subtractions
            if isinstance(n.op, (ast.Add, ast.Sub)):
                r = n.right
                neg = None
                if isinstance(r, ast.UnaryOp) and isinstance(r.op, ast.USub):
                    neg = r.operand
                elif isinstance(r, ast.BinOp) and isinstance(r.op, (ast.Mult, ast.Div)) and \
                        isinstance(r.left, ast.UnaryOp) and isinstance(r.left.op, ast.USub):
                    neg = ast.BinOp(r.left.operand, r.op, r.right)
                if neg is not None:
                    return F("Sub" if isinstance(n.op, ast.Add) else "Add")(c(n.left), c(neg))
            op = {ast.Add: "Add", ast.Sub: "Sub", ast.Mult: "Mul", ast.Div: "Div", ast.Pow: "Pow"}[type(n.op)]
            return F(op)(c(n.left), c(n.right))
        if isinstance(n, ast.Call):
            if isinstance(n.func, ast.Attribute) and n.func.attr == "subs":
                return F("Where")(c(n.func.value), F("Def")(c(n.args[0]), c(n.args[1])))
            fn = n.func.id
            args = [c(a) for a in n.args]
            if fn == "expand":
                return args[0]
            if fn == "cbrt":
                return F("Pow")(args[0], F("Div")(1, 3))
            if fn == "zeta":
                return F("RiemannZeta")(*args) if len(args) == 1 else F("HurwitzZeta")(*args)
            if fn == "polygamma":
                return F("DigammaFunction")(args[1], args[0])
            return F(_FEXPR_FUNCS[fn])(*args)
        raise Unsupported(type(n).__name__)
    return c(ast.parse(s, mode="eval").body)

def to_latex(s, limit=700):
    """LaTeX for the expression s via fexpr, or None if it is too large"""
    if s is None or len(s) > limit:
        return None
    try:
        t = to_fexpr(s).latex()
    except Exception:
        return None
    for a, b in ((r"\operatorname{ConstGamma}", r"\gamma"), (r"\operatorname{RealRoot}", r"\operatorname{realroot}")):
        t = t.replace(a, b)
    return t

def _meta(engines, timeout, jobs):
    import platform, subprocess
    m = dict(date=time.strftime("%Y-%m-%d %H:%M"), engines=list(engines), timeout=timeout, jobs=jobs,
             python=platform.python_version(), machine=platform.machine())
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        m["commit"] = subprocess.run(["git", "-C", here, "log", "-1", "--format=%h %s"],
                                     capture_output=True, text=True).stdout.strip()
    except Exception:
        m["commit"] = ""
    if "sympy" in engines:
        import sympy
        m["sympy"] = sympy.__version__
    return m

def to_json(cases, meta):
    out = []
    for c in cases:
        d = {k: v for k, v in c.items() if k not in ("func", "sympy_func")}
        if d["expr"] and len(d["expr"]) > 20000:
            d["expr"] = d["expr"][:2000] + " ... (%d characters)" % len(d["expr"])
        out.append(d)
    return dict(meta=meta, cases=out)

_CSS = r"""
:root {
  --bg: #f7f6f2; --panel: #ffffff; --fg: #1c1d1f; --muted: #62656b; --faint: #9a9da3;
  --rule: #dddbd3; --rule2: #ecebe5; --accent: #0d5f73; --accent-bg: #e3eff2;
  --ok: #1d7a4b; --ok-bg: #e2f2e8; --wrong: #c22630; --wrong-bg: #fbe3e3;
  --undecided: #9a6400; --undecided-bg: #f8eed9; --unable: #b0502a; --unable-bg: #f7e6dd;
  --error: #9b2f7b; --error-bg: #f5e2ef; --timeout: #46608a; --timeout-bg: #e3e8f1;
  --na: #a7a9ae; --na-bg: transparent; --code-bg: #f1f0ea;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #131518; --panel: #1b1e22; --fg: #e6e6e3; --muted: #a0a3a8; --faint: #6c7076;
  --rule: #33373d; --rule2: #262a2f; --accent: #6cc3d6; --accent-bg: #16323a;
  --ok: #52c98a; --ok-bg: #173225; --wrong: #ff6f73; --wrong-bg: #3d1a1c;
  --undecided: #e5ad48; --undecided-bg: #352a14; --unable: #ee8a5c; --unable-bg: #3a2318;
  --error: #e27cc6; --error-bg: #361a2e; --timeout: #8ea8d6; --timeout-bg: #1f2737;
  --na: #5c6066; --code-bg: #22262b; color-scheme: dark } }
:root[data-theme="dark"] {
  --bg: #131518; --panel: #1b1e22; --fg: #e6e6e3; --muted: #a0a3a8; --faint: #6c7076;
  --rule: #33373d; --rule2: #262a2f; --accent: #6cc3d6; --accent-bg: #16323a;
  --ok: #52c98a; --ok-bg: #173225; --wrong: #ff6f73; --wrong-bg: #3d1a1c;
  --undecided: #e5ad48; --undecided-bg: #352a14; --unable: #ee8a5c; --unable-bg: #3a2318;
  --error: #e27cc6; --error-bg: #361a2e; --timeout: #8ea8d6; --timeout-bg: #1f2737;
  --na: #5c6066; --code-bg: #22262b; color-scheme: dark }
* { box-sizing: border-box }
body { background: var(--bg); color: var(--fg); margin: 0;
  font: 14px/1.5 "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif }
main { max-width: 1240px; margin: 0 auto; padding-inline: 20px; padding-block: 28px 64px;
  display: flex; flex-direction: column; gap: 36px }
h1 { font-size: 26px; line-height: 1.2; margin: 0; font-weight: 600; letter-spacing: -0.01em; text-wrap: balance }
h1 code { font: inherit; font-family: "IBM Plex Mono", ui-monospace, monospace; font-weight: 500 }
h2 { font-size: 13px; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted);
  font-weight: 600; margin: 0 0 12px }
.sub { color: var(--muted); margin-top: 6px; display: flex; flex-wrap: wrap; gap: 4px 16px }
.sub span b { color: var(--fg); font-weight: 500 }
code, .mono { font-family: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12.5px }
.scroll { overflow-x: auto; min-width: 0 }
table { border-collapse: collapse; font-variant-numeric: tabular-nums }
th { text-align: left; font-weight: 500; color: var(--muted); font-size: 12px; padding: 6px 10px;
  border-bottom: 1px solid var(--rule); white-space: nowrap }
td { padding: 6px 10px; border-bottom: 1px solid var(--rule2); vertical-align: top }
td.n, th.n { text-align: right }
.dot { display: inline-block; width: 8px; height: 8px; border-radius: 50%; margin-right: 6px;
  vertical-align: 1px; background: currentColor }
.s-ok { color: var(--ok) } .s-wrong { color: var(--wrong) } .s-undecided { color: var(--undecided) }
.s-unable { color: var(--unable) } .s-error { color: var(--error) } .s-timeout { color: var(--timeout) }
.s-na { color: var(--na) }
.zero { color: var(--faint) }
.bar { display: flex; height: 8px; width: 160px; border-radius: 2px; overflow: hidden; background: var(--rule2) }
.bar i { display: block; height: 100% }
.bar .b-ok { background: var(--ok) } .bar .b-wrong { background: var(--wrong) }
.bar .b-undecided { background: var(--undecided) } .bar .b-unable { background: var(--unable) }
.bar .b-error { background: var(--error) } .bar .b-timeout { background: var(--timeout) }
.mini { display: flex; align-items: center; gap: 8px; white-space: nowrap }
.mini .bar { width: 64px; height: 6px }
.alert { border: 1px solid var(--wrong); background: var(--wrong-bg); border-radius: 6px; padding: 12px 16px }
.alert h2 { color: var(--wrong); margin-bottom: 6px }
.alert ul { margin: 0; padding-left: 18px }
.good { border-left: 3px solid var(--ok); padding: 4px 12px; color: var(--muted) }
.legend { display: flex; flex-wrap: wrap; gap: 4px 16px; color: var(--muted); font-size: 12.5px; margin-top: 10px }
.grid2 { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 520px), 1fr)); gap: 36px }
.grid2 > * { min-width: 0 }
.controls { display: flex; flex-wrap: wrap; gap: 8px 12px; align-items: center; margin-bottom: 12px }
.controls input, .controls select { font: inherit; color: var(--fg); background: var(--panel);
  border: 1px solid var(--rule); border-radius: 5px; padding: 5px 8px }
.controls input { width: min(100%, 260px) }
.controls .count { color: var(--muted); margin-left: auto }
#cases td.id { white-space: nowrap }
#cases td.id .cid { font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 12.5px; font-weight: 500 }
#cases td.id .cat { display: block; color: var(--muted); font-size: 12px }
#cases td.f { max-width: 420px; min-width: 200px }
#cases td.f .fx { overflow-x: auto; max-width: 420px; padding-bottom: 2px }
#cases td.f .big { color: var(--faint); font-style: italic }
#cases td.r { white-space: nowrap; font-size: 12.5px }
#cases td.r.wrong { background: var(--wrong-bg) }
#cases tr.row { cursor: pointer }
#cases tr.row:hover td { background: var(--rule2) }
#cases tr.row:hover td.r.wrong { background: var(--wrong-bg) }
#cases tr.detail td { background: var(--panel); padding: 12px 16px 16px }
.det { display: grid; gap: 6px 16px; grid-template-columns: max-content minmax(0, 1fr); font-size: 13px }
.det dt { color: var(--muted) }
.det dd { margin: 0; min-width: 0 }
.det pre { margin: 0; white-space: pre-wrap; word-break: break-all; background: var(--code-bg);
  padding: 8px 10px; border-radius: 4px; max-height: 220px; overflow: auto }
.expect { font-size: 12px; color: var(--muted); white-space: nowrap }
math { font-size: 1.08em }
footer { color: var(--faint); font-size: 12.5px }
"""

_JS = r"""
(function () {
  function render() {
    if (!window.katex) return;
    document.querySelectorAll('[data-tex]').forEach(function (el) {
      try { katex.render(el.getAttribute('data-tex'), el, {output: 'mathml', throwOnError: false, displayMode: false}); }
      catch (e) {}
    });
  }
  render();
  var q = document.getElementById('q'), cat = document.getElementById('cat'), st = document.getElementById('st');
  var rows = Array.prototype.slice.call(document.querySelectorAll('#cases tr.row'));
  var ref = document.getElementById('report').getAttribute('data-ref');
  function keep(r) {
    var s = (q.value || '').toLowerCase();
    if (s && r.getAttribute('data-search').indexOf(s) < 0) return false;
    if (cat.value && r.getAttribute('data-cat') !== cat.value) return false;
    var res = JSON.parse(r.getAttribute('data-res')), v = st.value, ks = Object.keys(res);
    var app = ks.filter(function (k) { return res[k] !== 'n/a'; });
    var oks = app.filter(function (k) { return res[k] === 'ok'; });
    if (v === 'wrong') return ks.some(function (k) { return res[k] === 'wrong'; });
    if (v === 'disagree') return oks.length > 0 && oks.length < app.length;
    if (v === 'refnot') return res[ref] !== 'ok';
    if (v === 'refonly') return res[ref] === 'ok' && oks.length < app.length;
    if (v === 'allok') return oks.length === app.length;
    if (v === 'noneok') return oks.length === 0;
    return true;
  }
  function apply() {
    var n = 0;
    rows.forEach(function (r) {
      var k = keep(r); r.hidden = !k; if (k) n++;
      var d = r.nextElementSibling; if (!k && d) d.hidden = true;
    });
    document.getElementById('count').textContent = n + ' of ' + rows.length + ' cases';
  }
  [q, cat, st].forEach(function (el) { el.addEventListener('input', apply); });
  rows.forEach(function (r) {
    r.addEventListener('click', function () { var d = r.nextElementSibling; d.hidden = !d.hidden; });
  });
  apply();
})();
"""

def _bar(counts, total, cls="bar"):
    if not total:
        return '<span class="%s"></span>' % cls
    parts = "".join('<i class="b-%s" style="width:%.3f%%"></i>' % (s, 100.0 * counts.get(s, 0) / total)
                    for s in ("ok", "wrong", "undecided", "unable", "error", "timeout") if counts.get(s))
    return '<span class="%s">%s</span>' % (cls, parts)

def _sec(t):
    if t is None:
        return "&ndash;"
    if t < 0.001:
        return "%.0f&thinsp;&micro;s" % (t * 1e6)
    if t < 1:
        return "%.1f&thinsp;ms" % (t * 1e3)
    return "%.2f&thinsp;s" % t

def _sclass(s):
    return "s-" + ("na" if s == "n/a" else s)

def html_report(data, fragment=False):
    """the HTML report (a string) for the results in data (as written by to_json);
    with fragment, without the document skeleton (for embedding)"""
    E = html.escape
    meta, cases = data["meta"], data["cases"]
    engines = meta["engines"]
    ref = engines[0]
    S = summary(cases, engines)
    out = []
    w = out.append
    if not fragment:
        w('<!doctype html><html lang="en"><head><meta charset="utf-8">'
          '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">')
    w('<title>gr_tower profile</title>')
    w('<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>')
    w('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap">')
    w('<style>%s</style>' % _CSS)
    w(('' if fragment else '</head><body>') + '<main id="report" data-ref="%s">' % E(ref))

    # header
    w('<header><h1>Exact numbers profile: <code>%s</code></h1><div class="sub">' % E(" · ".join(engines)))
    w('<span><b>%d</b> cases</span>' % len(cases))
    w('<span>time limit <b>%g&thinsp;s</b></span>' % meta["timeout"])
    if meta.get("jobs", 1) > 1:
        w('<span><b>%d</b> parallel jobs</span>' % meta["jobs"])
    w('<span>%s</span>' % E(meta.get("date", "")))
    if meta.get("commit"):
        w('<span>FLINT <span class="mono">%s</span></span>' % E(meta["commit"][:60]))
    if meta.get("sympy"):
        w('<span>SymPy %s</span>' % E(meta["sympy"]))
    w('</div></header>')

    # wrong answers
    wrong = [(c, e) for c in cases for e in engines if c["results"].get(e, {}).get("status") == "wrong"]
    if wrong:
        w('<section class="alert"><h2>Wrong answers (%d)</h2><ul>' % len(wrong))
        for c, e in wrong:
            w('<li><b>%s</b> on <span class="mono">%s</span> (expected %s)%s</li>'
              % (E(e), E(c["id"]), E(c["expect"].split(":")[0]), (": " + E(c["note"])) if c["note"] else ""))
        w('</ul></section>')
    else:
        w('<div class="good">No engine gave a wrong answer.</div>')

    # engines
    w('<section><h2>Engines</h2><div class="scroll"><table><thead><tr><th>engine</th>')
    for s in STATUSES:
        w('<th class="n"><span class="dot %s"></span>%s</th>' % (_sclass(s), s))
    w('<th>of applicable</th><th class="n">time on ok</th><th class="n">median</th></tr></thead><tbody>')
    for e in engines:
        cnt = S[e]["counts"]
        w('<tr><td><b>%s</b></td>' % E(e))
        for s in STATUSES:
            v = cnt[s]
            w('<td class="n%s">%d</td>' % ("" if v else " zero", v))
        w('<td>%s</td><td class="n">%s</td><td class="n">%s</td></tr>'
          % (_bar(cnt, S[e]["applicable"]), _sec(S[e]["total"]), _sec(S[e]["median"])))
    w('</tbody></table></div></section>')

    # head to head and categories
    w('<div class="grid2">')
    w('<section><h2>Statuses</h2><div class="legend" style="margin:0;flex-direction:column;gap:4px">'
      '<span><span class="dot s-ok"></span>ok: the right answer</span>'
      '<span><span class="dot s-wrong"></span>wrong: a wrong answer (a near miss decided zero, a zero decided nonzero, a false relation)</span>'
      '<span><span class="dot s-undecided"></span>undecided: the zero test or comparison cannot decide</span>'
      '<span><span class="dot s-unable"></span>unable: an operation cannot be computed</span>'
      '<span><span class="dot s-error"></span>error: any other failure</span>'
      '<span><span class="dot s-timeout"></span>timeout: over the time limit</span>'
      '<span><span class="dot s-na"></span>n/a: not applicable (qqbar: transcendental; SymPy: programs without a SymPy version)</span></div></section>')
    if len(engines) > 1:
        w('<section><h2>Against %s</h2><div class="scroll"><table><thead><tr><th>engine</th>'
          '<th class="n">both ok</th><th class="n">only %s</th><th class="n">only other</th>'
          '<th class="n">time ratio<br>geometric mean</th><th class="n">median</th></tr></thead><tbody>' % (E(ref), E(ref)))
        for e in engines[1:]:
            both, only_r, only_e, ratios = 0, 0, 0, []
            for c in cases:
                a, b = c["results"].get(ref, {}), c["results"].get(e, {})
                if b.get("status") == "n/a":
                    continue
                ra, rb = a.get("status") == "ok", b.get("status") == "ok"
                if ra and rb:
                    both += 1
                    ratios.append(max(b["time"], 1e-4) / max(a["time"], 1e-4))
                elif ra:
                    only_r += 1
                elif rb:
                    only_e += 1
            ratios.sort()
            gm = math.exp(sum(map(math.log, ratios)) / len(ratios)) if ratios else None
            md = ratios[len(ratios) // 2] if ratios else None
            f = lambda x: "&ndash;" if x is None else ("%.3g&times;" % x)
            w('<tr><td><b>%s</b></td><td class="n">%d</td><td class="n">%d</td><td class="n">%d</td>'
              '<td class="n">%s</td><td class="n">%s</td></tr>' % (E(e), both, only_r, only_e, f(gm), f(md)))
        w('</tbody></table></div><div class="legend"><span>ratio = time of the engine / time of %s, '
          'on the cases both get right (times below 0.1&thinsp;ms counted as 0.1&thinsp;ms)</span></div></section>' % E(ref))
    w('</div>')
    cats = []
    for c in cases:
        if c["category"] not in cats:
            cats.append(c["category"])
    w('<section><h2>Categories</h2><div class="scroll"><table><thead><tr><th>category</th><th class="n">cases</th>')
    for e in engines:
        w('<th>%s</th>' % E(e))
    w('</tr></thead><tbody>')
    for k in sorted(cats):
        cs = [c for c in cases if c["category"] == k]
        w('<tr><td>%s</td><td class="n">%d</td>' % (E(k), len(cs)))
        for e in engines:
            cnt = {s: sum(c["results"].get(e, {}).get("status") == s for c in cs) for s in STATUSES}
            app = len(cs) - cnt["n/a"]
            if app == 0:
                w('<td class="zero">n/a</td>')
            else:
                w('<td><span class="mini">%s<span>%d/%d</span></span></td>' % (_bar(cnt, app, "bar"), cnt["ok"], app))
        w('</tr>')
    w('</tbody></table></div></section>')

    # cases
    w('<section><h2>Cases</h2><div class="controls">'
      '<input id="q" type="search" placeholder="Search id, source, note" aria-label="search">'
      '<select id="cat" aria-label="category"><option value="">all categories</option>')
    for k in sorted(cats):
        w('<option>%s</option>' % E(k))
    w('</select><select id="st" aria-label="status filter"><option value="">all results</option>'
      '<option value="wrong">a wrong answer</option><option value="disagree">engines disagree</option>'
      '<option value="refnot">%s not ok</option><option value="refonly">%s ok, another not</option>'
      '<option value="allok">all ok</option><option value="noneok">none ok</option></select>'
      '<span class="count" id="count"></span></div>' % (E(ref), E(ref)))
    w('<div class="scroll"><table id="cases"><thead><tr><th>case</th><th>formula</th><th>expect</th>')
    for e in engines:
        w('<th>%s</th>' % E(e))
    w('</tr></thead><tbody>')
    for c in cases:
        res = {e: c["results"].get(e, {}).get("status", "n/a") for e in engines}
        search = " ".join([c["id"], c["category"], c.get("source") or "", c.get("note") or ""]).lower()
        w('<tr class="row" data-cat="%s" data-search="%s" data-res="%s">'
          % (E(c["category"]), E(search), E(json.dumps(res))))
        w('<td class="id"><span class="cid">%s</span><span class="cat">%s</span></td>' % (E(c["id"]), E(c["category"])))
        tex = c.get("tex") or to_latex(c["expr"])
        if tex:
            w('<td class="f"><div class="fx"><span data-tex="%s">%s</span></div></td>'
              % (E(tex), E(c["expr"] if c["expr"] else tex)))
        elif c["expr"]:
            w('<td class="f"><span class="big">%d-character expression</span></td>' % len(c["expr"]))
        else:
            w('<td class="f"></td>')
        ex = c["expect"]
        w('<td class="expect">%s</td>' % E("minpoly" if ex.startswith("minpoly:") else ex))
        for e in engines:
            r = c["results"].get(e, {})
            s = r.get("status", "n/a")
            t = r.get("time")
            txt = "n/a" if s == "n/a" else (s if s != "ok" else "") + (" " if s != "ok" else "") + \
                  ("" if s in ("timeout",) or t is None else _sec(t))
            w('<td class="r%s"><span class="%s"><span class="dot"></span></span>%s</td>'
              % (" wrong" if s == "wrong" else "", _sclass(s), txt.strip()))
        w('</tr><tr class="detail" hidden><td colspan="%d"><dl class="det">' % (3 + len(engines)))
        if c["expr"]:
            w('<dt>expression</dt><dd><pre>%s</pre></dd>' % E(c["expr"]))
        if ex.startswith("minpoly:"):
            w('<dt>minimal polynomial</dt><dd class="mono">%s</dd>' % E(ex[8:]))
        if c.get("approx"):
            w('<dt>value</dt><dd class="mono">%s</dd>' % E(c["approx"]))
        if c.get("note"):
            w('<dt>stresses</dt><dd>%s</dd>' % E(c["note"]))
        if c.get("source"):
            w('<dt>source</dt><dd>%s</dd>' % E(c["source"]))
        for e in engines:
            d = c["results"].get(e, {}).get("detail")
            if d:
                w('<dt>%s</dt><dd class="mono">%s</dd>' % (E(e), E(d)))
        w('</dl></td></tr>')
    w('</tbody></table></div></section>')
    w('<footer>Generated by <span class="mono">src/python/gr_tower_profile.py</span>. '
      'Click a row for the expression, source and messages.</footer></main>')
    w('<script src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.js"></script>')
    w('<script>%s</script>' % _JS + ('' if fragment else '</body></html>'))
    return "".join(out)


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip(),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pattern", nargs="?", default="", help="run the cases whose id contains this string")
    ap.add_argument("--engines", default="tower,ca", help="comma-separated, from tower,ca,qqbar,sympy")
    ap.add_argument("--category", default=None, help="run one category")
    ap.add_argument("--timeout", type=float, default=60.0, help="seconds per case and engine")
    ap.add_argument("--jobs", type=int, default=1, help="subprocesses at a time")
    ap.add_argument("--no-families", action="store_true", help="leave out the scalable families")
    ap.add_argument("--no-big", action="store_true", help="leave out the big expressions")
    ap.add_argument("--json", default=None, help="write the results to this file")
    ap.add_argument("--html", default=None, help="write an HTML report to this file")
    ap.add_argument("--from-json", default=None, help="read results from this file instead of running")
    ap.add_argument("--fragment", action="store_true", help="HTML without the document skeleton")
    ap.add_argument("--list", action="store_true", help="list the cases and exit")
    args = ap.parse_args(argv)

    if args.from_json:
        data = json.load(open(args.from_json))
        if args.html:
            open(args.html, "w").write(html_report(data, args.fragment))
        return
    engines = [e for e in args.engines.split(",") if e]
    for e in engines:
        if e not in ENGINES:
            ap.error("unknown engine %s" % e)
    cases = all_cases(families=not args.no_families, big=not args.no_big)
    cases = [c for c in cases if args.pattern in c["id"] and (args.category is None or c["category"] == args.category)]
    if args.list:
        for c in cases:
            print("%-28s %-18s %-8s %s" % (c["id"], c["category"], c["expect"].split(":")[0],
                                           (c["expr"] or c["tex"])[:100]))
        return
    run(cases, engines, args.timeout, args.jobs)
    print_summary(cases, engines)
    data = to_json(cases, _meta(engines, args.timeout, args.jobs))
    if args.json:
        json.dump(data, open(args.json, "w"), indent=1)
    if args.html:
        open(args.html, "w").write(html_report(data, args.fragment))

if __name__ == "__main__":
    main()
