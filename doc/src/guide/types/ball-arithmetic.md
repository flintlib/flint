---
title: Ball arithmetic
summary: Rigorous real and complex arithmetic with balls: the inclusion principle, printing, overestimation, predicates, precision and accuracy, the guess-and-verify strategy, and the conventions specific to the ball types.
order: 2
---

FLINT computes rigorously with real and complex numbers in *ball
arithmetic* (midpoint-radius interval arithmetic). A real number is
represented by a ball $[m \pm r]$, the set of real numbers $t$ with
$|t - m| \le r$, where the midpoint $m$ is an arbitrary-precision binary
floating-point number (`arf_t`) and the radius $r$ a low-precision one
(`mag_t`, with a 30-bit mantissa). Every operation returns a ball that
contains all possible exact results for inputs in the input balls, so a
computation produces a rigorous enclosure of the exact result with no
error analysis by the user: the radius of the output proves the accuracy
of the midpoint. Since only the midpoint has full precision, a ball
operation at high precision costs about the same as the corresponding
floating-point operation, half the cost of inf-sup interval arithmetic.
Van der Hoeven's report [@Hoe2009] is a good introduction to ball
arithmetic, [@Joh2017] describes the design of FLINT's implementation
(developed as the separate library Arb until FLINT 3.0), and Rump's
survey [@Rum2010] covers rigorous numerical methods in general.

The modules [arb](arb) (real balls) and [acb](acb) (complex balls, as
pairs of real balls) provide arithmetic and elementary functions; on
them are built polynomials and power series ([arb_poly](arb_poly),
[acb_poly](acb_poly)), matrices ([arb_mat](arb_mat),
[acb_mat](acb_mat)), special functions ([arb_hypgeom](arb_hypgeom),
[acb_hypgeom](acb_hypgeom), [acb_elliptic](acb_elliptic),
[acb_modular](acb_modular), [acb_theta](acb_theta),
[acb_dirichlet](acb_dirichlet)), root finding and integration
([arb_calc](arb_calc), [acb_calc](acb_calc)), and generic rings
([](gr_domains#real-complex)). The `dfloat` and `decimal` modules
provide further ball types with less functionality
([](guide/types/choosing#real-complex)). [](arb/getting-started) is a
practical introduction to the programming interface; the algorithms
are described in [](algorithms/ball-arithmetic).

## Ball semantics {#semantics}

A ball implementation of a function $f : A \to B$ is a function $F$ that
maps sets $X \subseteq A$ to sets $F(X) \subseteq B$ such that

$$
x \in X \implies f(x) \in F(X).
$$

This *inclusion principle* is the contract of every function in FLINT's
ball arithmetic; the documentation simply writes $f(x)$ for the ball
implementation. Functions that take a precision argument, such as
`arb_add`, round the midpoint of the output to `prec` bits and add the
rounding error to the radius; functions that are always exact, such as
`arb_neg`, have no precision argument.

## A first program {#first}

An `arb_t` is initialised with `arb_init` (which sets it to 0) and freed
with `arb_clear`; output arguments come first, then the inputs, then
the precision in bits.

```c example name=ball_first
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y;
    arb_init(x);
    arb_init(y);

    arb_set_ui(x, 3);          /* x = 3, exactly */
    arb_const_pi(y, 128);      /* y = pi, to 128 bits */
    arb_sub(y, y, x, 53);      /* y = y - x, to 53 bits */

    arb_printn(y, 20, 0);  flint_printf("\n");
    arb_printd(y, 20);     flint_printf("\n");
    arb_print(y);          flint_printf("\n");

    arb_clear(x);
    arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
[0.1415926535897932 +/- 5.48e-17]
0.14159265358979322702 +/- 2.7756e-17
(1275353243890785 * 2^-53) +/- (536870913 * 2^-84)
```

`arb_printn` prints a decimal ball that is guaranteed to contain the
binary ball, showing only digits of the midpoint that are correct up to
one unit in the last place; this is the function to use for output.
`arb_printd` prints decimal approximations of the midpoint and radius,
which are not guaranteed to enclose the ball, and `arb_print` the exact
binary representation. `arb_get_str` and `arb_set_str` convert to and
from strings ([](arb/printing)). `flint_cleanup` frees cached data
([](guide/conventions/memory#cleanup)).

A common mistake is to create a ball from a `double` constant such as
`2.3`, which is in fact the binary number
2.29999999999999982236431605997495353221893310546875: `arb_set_d`
converts it exactly, giving a ball that does not contain $23/10$. To
enclose $23/10$, use one of

```c
arb_set_str(x, "2.3", prec);

arb_set_ui(x, 23);
arb_div_ui(x, x, 10, prec);

fmpq_set_si(q, 23, 10);      /* q is an fmpq_t */
arb_set_fmpq(x, q, prec);
```

## Overestimation {#overestimation}

The computed enclosure is usually larger than the smallest possible
one, because of rounding errors and cancellation in the individual
steps, and because each operation treats its inputs as independent (the
*dependency problem*), so that correlated errors add up instead of
cancelling:

```c example name=ball_dependency
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y, z;
    arb_init(x); arb_init(y); arb_init(z);

    arb_set_str(x, "[3 +/- 0.1]", 53);

    arb_add_ui(y, x, 1, 53);
    arb_sub(z, y, x, 53);           /* (x + 1) - x */
    flint_printf("(x + 1) - x = "); arb_printd(z, 5); flint_printf("\n");

    arb_sub(z, x, x, 53);           /* x - x */
    flint_printf("x - x       = "); arb_printd(z, 5); flint_printf("\n");

    arb_set_str(x, "[0 +/- 1]", 53);
    arb_mul(z, x, x, 53);           /* x^2 for x in [-1, 1] */
    flint_printf("x * x       = "); arb_printd(z, 5); flint_printf("\n");

    arb_clear(x); arb_clear(y); arb_clear(z);
    flint_cleanup();
    return 0;
}
```
```output
(x + 1) - x = 1.0000 +/- 0.20000
x - x       = 0 +/- 0.20000
x * x       = 0 +/- 1.0000
```

In unfavourable cases the radius grows exponentially with the number of
steps. When the inputs are exact or can be computed to any accuracy, an
overestimation by a factor $2^n$ can usually be compensated by $n$
extra bits of precision, but an iteration that converges in
floating-point arithmetic may fail to converge in ball arithmetic.
Techniques that make floating-point code numerically stable also give
tighter enclosures; ball arithmetic is not a substitute for numerical
analysis.

## Predicates {#predicates}

A ball predicate returns 1 when the property is certainly true, and 0
when it is false *or unknown*. To show that a property is false, test
the negated property. The following code is therefore *wrong* in
general:

```c
if (arb_is_positive(x))
{
    ...  /* x > 0 */
}
else
{
    ...  /* wrong: x <= 0 is not known here */
}
```

and the correct version has three branches:

```c
if (arb_is_positive(x))
{
    ...  /* x > 0 */
}
else if (arb_is_nonpositive(x))
{
    ...  /* x <= 0 */
}
else
{
    ...  /* the sign of x is not known */
}
```

For balls, $x \le y$ means that it holds for all points of the balls
(`arb_le`). Some predicates, such as `arb_overlaps` and `arb_contains`,
are predicates on balls as sets; the negation of `arb_contains_zero` is
`arb_is_nonzero`, while the pointwise test "$x = 0$" is `arb_is_zero`.
[](arb/predicates) lists all predicates and their treatment of NaN and
infinities.

## Precision and accuracy {#precision}

The precision argument sets the *working precision*, not the accuracy
of the result, which depends on the function, the input and the
conditioning of the problem. For a ball $[m \pm r]$ approximating $x$,
the absolute error is $r$ and the relative error $r / \max(0, |m| - r)$;
`arb_rel_accuracy_bits` estimates $-\log_2$ of the relative error
cheaply, and `arb_get_mag` and `arb_get_rad_arb` give the radius.

When a fixed expression is evaluated, the error of the result is
generically bounded by $2^{O(1) - \mathit{prec}}$, where the $O(1)$ term
depends on the expression and on the implementation. For numerically
well-behaved expressions, adding 10 or 20 *guard bits* to the target
accuracy suffices. Otherwise, increase the precision in exponential
steps until the result is accurate enough (Ziv's strategy); a partial
result also estimates the $O(1)$ term: if 120 bits of precision give 80
bits of accuracy, a little over 140 bits will probably give 100 bits.
The next program evaluates $(1 + 1/n)^n - e$ for $n = 10^{15}$, which
loses about 51 bits to cancellation:

```c example name=ball_precision_loop
#include <flint/flint.h>
#include <flint/arb.h>

/* f = (1 + 1/n)^n - e for n = 10^15, which suffers from cancellation */
void f(arb_t res, slong prec)
{
    arb_t t, e;
    arb_init(t); arb_init(e);
    arb_set_ui(t, 1000000000000000);   /* n = 10^15 */
    arb_inv(t, t, prec);
    arb_log1p(t, t, prec);             /* log(1 + 1/n) */
    arb_mul_ui(t, t, 1000000000000000, prec);
    arb_exp(t, t, prec);
    arb_const_e(e, prec);
    arb_sub(res, t, e, prec);
    arb_clear(t); arb_clear(e);
}

int main(void)
{
    arb_t y;
    slong prec, goal = 53;

    arb_init(y);
    for (prec = goal + 10; ; prec *= 2)
    {
        f(y, prec);
        flint_printf("prec = %4wd: accuracy %3wd bits, f = ", prec,
                     arb_rel_accuracy_bits(y));
        arb_printn(y, 16, 0); flint_printf("\n");
        if (arb_rel_accuracy_bits(y) >= goal)
            break;
    }
    arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
prec =   63: accuracy   8 bits, f = [-1.36e-15 +/- 2.24e-18]
prec =  126: accuracy  71 bits, f = [-1.359140914229521e-15 +/- 3.72e-31]
```

(Computing $1 + 1/n$ first and then its logarithm, instead of using
`arb_log1p`, would lose another 50 bits.)

Simple operations, such as arithmetic and most elementary functions,
compute $f(m)$ to `prec` bits of relative accuracy, up to a few units in
the last place, and bound the propagated error separately
([](algorithms/ball-arithmetic/error-propagation)); their enclosures are
close to the best possible. More complex functions, such as many higher
transcendental functions, are evaluated by sequences of ball operations
that may depend on `prec` (for example, a series truncated where the
remainder is below $2^{-\mathit{prec}}$). Their enclosures can be far
from optimal, and convergence as $\mathit{prec} \to \infty$ is not
guaranteed in general, although it holds in practice. The inclusion
principle is the contract; the tightness of the enclosures is a matter
of quality of implementation.

## Guess and verify {#guess-and-verify}

The precision loop above is an instance of the strategy typical of
ball arithmetic: instead of determining in advance the precision needed
for a correct result, guess one, use ball arithmetic to verify that the
result is accurate enough, and repeat with a higher precision
otherwise. A floating-point computation with unknown rounding errors is
like a Monte Carlo algorithm, which may return a wrong answer; ball
arithmetic turns it into a Las Vegas algorithm, which never returns a
wrong answer but may fail to terminate.

The next program sums the Taylor series of $\sin x$ naively. Ball
arithmetic handles the rounding errors of the finite sum; the tail is
bounded by hand: from the first term $t$ with $|t| \le
2^{-\mathit{prec}}$ on, the terms alternate and decrease, so the tail is
at most $|t|$, which is added to the radius with `arb_add_error_mag`.
For $x = 2016.1$ the largest term is about $10^{873}$ while the sum is
less than 1, so about 2900 bits are lost to cancellation. The program
doubles the precision until the sign of $\sin x$ is known:

```c example name=ball_sin_naive
#include <flint/flint.h>
#include <flint/arb.h>

/* sin(x) by its Taylor series. We stop at the first term t with
   |t| <= 2^-prec; from there on the terms alternate and decrease,
   so the tail is bounded by |t|, which is added to the radius. */
void arb_sin_naive(arb_t res, const arb_t x, slong prec)
{
    arb_t s, t, x2;
    mag_t bound;
    slong k;

    arb_init(s); arb_init(t); arb_init(x2); mag_init(bound);

    arb_mul(x2, x, x, prec);
    arb_set(t, x);                     /* t = x^(2k+1) / (2k+1)! */

    for (k = 0; ; k++)
    {
        arb_get_mag(bound, t);
        if (mag_cmp_2exp_si(bound, -prec) <= 0)
        {
            arb_add_error_mag(s, bound);
            break;
        }

        if (k % 2 == 0)
            arb_add(s, s, t, prec);
        else
            arb_sub(s, s, t, prec);

        arb_mul(t, t, x2, prec);
        arb_div_ui(t, t, (2 * k + 2) * (2 * k + 3), prec);
    }

    arb_swap(res, s);
    arb_clear(s); arb_clear(t); arb_clear(x2); mag_clear(bound);
}

int main(void)
{
    arb_t x, y;
    slong prec;

    arb_init(x); arb_init(y);

    for (prec = 64; ; prec *= 2)
    {
        arb_set_str(x, "2016.1", prec);
        arb_sin_naive(y, x, prec);
        flint_printf("Using %5wd bits, sin(x) = ", prec);
        arb_printn(y, 10, 0); flint_printf("\n");
        if (!arb_contains_zero(y))      /* the sign is known: stop */
            break;
    }

    arb_set_str(x, "2016.1", 64);
    arb_sin(y, x, 64);
    flint_printf("arb_sin, 64 bits:  ");
    arb_printn(y, 20, 0); flint_printf("\n");

    arb_clear(x); arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
Using    64 bits, sin(x) = [+/- 4.49e+859]
Using   128 bits, sin(x) = [+/- 2.44e+840]
Using   256 bits, sin(x) = [+/- 7.15e+801]
Using   512 bits, sin(x) = [+/- 6.18e+724]
Using  1024 bits, sin(x) = [+/- 4.61e+570]
Using  2048 bits, sin(x) = [+/- 2.57e+262]
Using  4096 bits, sin(x) = [-0.7190842207 +/- 1.20e-11]
arb_sin, 64 bits:  [-0.7190842207119598 +/- 5.36e-17]
```

At 4096 bits the sign is known (the radius printed there is the error of
rounding to 10 digits). The built-in `arb_sin` reduces the argument
modulo $\pi/4$ first, so 64 bits suffice. Only infinite processes such
as series need a hand-made bound; [](tutorials/numerics/series) treats
tail bounds systematically.

The loop does not terminate if the quantity is exactly zero, such as
$\sin \pi$:

```c example name=ball_sin_pi
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y;
    slong prec;

    arb_init(x); arb_init(y);

    /* sin(pi) = 0: no precision suffices to determine its sign */
    for (prec = 64; prec <= 4096; prec *= 4)
    {
        arb_const_pi(x, prec);
        arb_sin(y, x, prec);
        flint_printf("%5wd bits: sin(pi)   = ", prec);
        arb_printn(y, 10, 0); flint_printf("\n");
    }

    /* sin(pi x) at the exact input x = 1 is exact */
    arb_one(x);
    arb_sin_pi(y, x, 64);
    flint_printf("  64 bits: sin(pi*1) = "); arb_printn(y, 10, 0); flint_printf("\n");

    arb_clear(x); arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
   64 bits: sin(pi)   = [+/- 1.59e-19]
  256 bits: sin(pi)   = [+/- 2.83e-77]
 1024 bits: sin(pi)   = [+/- 2.03e-308]
 4096 bits: sin(pi)   = [+/- 2.74e-1233]
  64 bits: sin(pi*1) = 0
```

The sign of a nonzero number can always be found by computing it
accurately enough, but an exact zero cannot be recognised numerically
unless the whole computation is exact, as with `arb_sin_pi` at an
integer. A loop that doubles the precision therefore needs a stopping
criterion suited to the application. Deciding such questions exactly is
the purpose of the [exact real numbers](guide/types/choosing#exact).

## The polynomial time guarantee {#polynomial-time}

FLINT provides a soft guarantee that the time used to evaluate a ball
function depends polynomially on `prec` and on the bit size of the input,
uniformly in the numerical value of the input, so that the user can
bound the resources of a computation by its precision. For example, an
accurate $\sin x$ for $x = 2^{2^n}$ needs a multiple of $2\pi$ to
$2^n$ bits, and $\exp x$ has an exponent of $2^n$ bits. Internal work
parameters (working precision, number of terms) are therefore capped by
polynomial, usually linear, functions of `prec`, and a crude bound is
returned when the cap is exceeded:

```c example name=ball_huge_input
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y;
    arb_init(x); arb_init(y);

    /* x = 2^(2^40): an exact number with a 41-bit exponent */
    arb_one(x);
    arb_mul_2exp_si(x, x, WORD(1) << 40);

    arb_sin(y, x, 64);
    flint_printf("sin(x)  = "); arb_printn(y, 10, 0); flint_printf("\n");
    arb_exp(y, x, 64);
    flint_printf("exp(x)  = "); arb_printn(y, 10, 0); flint_printf("\n");
    arb_neg(x, x);
    arb_exp(y, x, 64);
    flint_printf("exp(-x) = "); arb_printn(y, 10, 0); flint_printf("\n");

    arb_clear(x); arb_clear(y);
    flint_cleanup();
    return 0;
}
```
```output
sin(x)  = [+/- 1.00]
exp(x)  = [+/- inf]
exp(-x) = [+/- 1.69e-102435199438739363750012109250103232700]
```

The last bound is $2^{-2^{128}}$. Such bounds are often enough: in
$\log x + \sin(x) \exp(-x)$ for large $x$, the second term is
negligible, and $[\pm 1] \cdot [\pm 2^{-N}]$ suffices to evaluate the
expression accurately. The guarantee is soft: inputs that are most
likely chosen deliberately are not cut off; for example, the cost of
the Riemann zeta function $\zeta(\sigma + it)$ grows linearly with
$|t|$.

## Conventions {#conventions}

The ball modules follow the conventions of all modules
([](guide/conventions)); the following points are specific to them.

### Precision limits {#limits}

Precisions are `slong` values and are not checked against overflow
([](guide/correctness/limits)). There is no hard upper limit, but
$2^{36} \approx 10^{11}$ bits on 64-bit systems and $2^{24}$ bits on
32-bit systems can be considered safe (GMP integers, used internally,
have at most about $2^{37}$ bits). Loops that double the precision need
an upper bound, since the result may never become accurate. The minimum
precision is 2 bits.

Exponents are not limited: midpoints and radii have arbitrary-precision
exponents, so $2^{2^{100}}$ is an ordinary ball. Conversions whose exact
result has a size that depends on the exponent are the caller's
responsibility: converting $2^{2^{100}}$ to an integer with
`arf_get_fmpz` can return a wrong value or crash.

### Aliasing and algebraic simplification {#aliasing}

Inputs and outputs may be aliased unless documented otherwise
([](guide/conventions/types-and-naming#aliasing)). The rule concerns
whole objects of the same type: in `arb_poly_evaluate`, computing
$y = f(x)$, the output `y` must not be a coefficient of `f`, and an
`arb_srcptr` array of length 1 is not an `arb_t` for aliasing purposes.

Identical input pointers give the implementation *permission* to
simplify algebraically: `acb_mul(z, x, x, prec)` encloses $\{t^2 : t \in
x\}$, not the generally larger $\{t u : t, u \in x\}$. To enclose
$f(x, y)$ for two independent values in $[-1, 1]$, pass two distinct
balls, not the same ball twice. Simplification is allowed but not
guaranteed: `arb_add(z, x, x, prec)`, `arb_sub(z, x, x, prec)` and
currently `arb_mul(z, x, x, prec)` treat their inputs as independent;
use `arb_mul_2exp_si(z, x, 1)` to double a ball and `arb_zero` to get
zero.

### Cached constants {#caches}

Constants such as $\pi$ and $\log 2$ are read from tables at low
precision and computed and cached per thread at higher precision, as
are Bernoulli numbers and the logarithms of small primes used by the
elementary functions. This does not affect results, which are always
rounded to the requested precision. The caches are freed by
`flint_cleanup` ([](guide/conventions/memory#cleanup),
[](guide/conventions/threading)).

The use of hardware floating-point arithmetic in the ball modules is
described in [](guide/correctness/guarantees#floating-point).
