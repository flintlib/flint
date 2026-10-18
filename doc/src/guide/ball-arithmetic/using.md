---
title: Using ball arithmetic
summary: The inclusion principle, a first program, binary and decimal conversion, overestimation, predicates, and a worked example that turns a series into a rigorous function with the guess-and-verify strategy.
order: 2
---

This page introduces computing with real balls (`arb_t`) through
examples. The same principles apply to complex balls (`acb_t`) and to
polynomials and matrices of balls. The full interface is documented on
the module page [arb](arb); [](arb/getting-started) is a shorter
practical introduction.

## Ball semantics {#semantics}

Let $f : A \to B$ be a function. A ball implementation of $f$ is a
function $F$ that maps sets $X \subseteq A$ to sets $F(X) \subseteq B$
such that

$$
x \in X \implies f(x) \in F(X).
$$

In other words, $F(X)$ is an *enclosure* of $\{f(x) : x \in X\}$. This is
the *inclusion principle*, the fundamental contract of every function in
FLINT's ball arithmetic. The documentation simply writes $f(x)$ for the
ball implementation, with the understanding that the input is a set of
points and the output an enclosure.

The sets used are balls $[m \pm r] = [m - r, m + r]$ whose midpoint $m$
and radius $r$ are binary floating-point numbers $u \cdot 2^v$, $u, v \in
\mathbb{Z}$, together with the special values $\pm\infty$ and NaN. Any
ball whose midpoint is a real number can be rounded to a ball whose
midpoint has at most `prec` bits, by rounding the midpoint and adding the
rounding error to the radius. Functions that take a precision argument,
such as `arb_add`, round their output in this way; functions that are
always exact, such as `arb_neg`, have no precision argument.

## A first program {#first}

The interface follows the conventions of GMP and the rest of FLINT. An
`arb_t` must be initialised with `arb_init` (which sets it to 0) and freed
with `arb_clear`; it is passed by reference; and the output arguments of a
function come first, then the inputs, then the precision in bits.

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

`flint_cleanup` frees the data that FLINT has cached in the calling
thread (at 128 bits, $\pi$ is read from a table, but higher precisions
cache it); it is not required, but gives clean output in memory
checkers. See [](guide/conventions/memory#cleanup).

## Binary and decimal {#decimal}

Balls are stored in binary, and conversion to decimal generally needs
rounding. The three printing functions used above show the same ball in
different ways:

- `arb_printn` prints a *decimal ball* that is guaranteed to contain the
  binary ball. By default it shows only digits of the midpoint that are
  correct up to one unit in the last place, and the printed radius
  includes the error of the decimal conversion. This is the function to
  use for output.
- `arb_printd` prints decimal approximations of the midpoint and the
  radius. The printed interval is not guaranteed to contain the ball.
- `arb_print` prints the exact internal representation, with binary
  exponents.

`arb_get_str` returns the string printed by `arb_printn`, and
`arb_set_str` reads such strings, and plain decimal numbers, back
rigorously. The guide [](arb/printing) gives the details.

A common mistake is to create a ball from a `double` constant such as
`2.3`, which is in fact the binary number
2.29999999999999982236431605997495353221893310546875. `arb_set_d` converts
it exactly, giving a ball that does not contain $23/10$. To enclose the
rational number $23/10$, use one of

```c
arb_set_str(x, "2.3", prec);

arb_set_ui(x, 23);
arb_div_ui(x, x, 10, prec);

fmpq_set_si(q, 23, 10);      /* q is an fmpq_t */
arb_set_fmpq(x, q, prec);
```

## Quality of enclosures {#overestimation}

The main problem of ball and interval arithmetic is *overestimation*:
the computed enclosure is usually larger than the smallest possible
enclosure of the result. It comes from rounding errors and from
cancellation in the individual steps, and techniques that make
floating-point code numerically stable also make ball arithmetic produce
tighter enclosures.

A particular source of overestimation is the *dependency problem*. Each
operation treats its inputs as independent numbers, so correlated errors
add up instead of cancelling. In floating-point arithmetic, $(x + 1) - x$
gives 1 up to a small rounding error even if $x$ itself has a large error,
because that error cancels. In ball arithmetic, with $x = [3 \pm 0.1]$,
the result is $[1 \pm 0.2]$: the radius has doubled.

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

Even $x - x$ is not 0, and $x \cdot x$ for $x \in [-1, 1]$ contains
negative numbers. In unfavourable cases the radius grows exponentially
with the number of steps.

If all inputs of a computation are exact or are constants that can be
computed to any accuracy (such as $\pi$), an overestimation by a factor
$2^n$ can usually be compensated by working with $n$ extra bits of
precision, at a cost polynomial in $n$. In some situations, however,
overestimation makes an algorithm exponentially slow or prevents it from
converging: an iterative root-finding method that converges in
floating-point arithmetic may fail to converge in ball arithmetic.
Ball arithmetic is not a substitute for numerical analysis, and some
experimentation may be needed to find out whether, and how, it can be
used for a given problem.

## Predicates {#predicates}

A ball implementation of a predicate $f : \mathbb{R} \to
\{\mathrm{true}, \mathrm{false}\}$ needs a third value meaning that the
answer could be either. In FLINT, predicates return the `int` 1 when the
answer is certainly true, and 0 when it is false *or unknown*. To find
out that a predicate is certainly false, one must test whether the
negated predicate is certainly true. The following code is therefore
*wrong* in general:

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

In the same way, $x \le y$ for balls $X$ and $Y$ means that $x \le y$
for all $x \in X$ and $y \in Y$ (this is `arb_le`). Some predicates, such
as `arb_overlaps` and `arb_contains`, are predicates on balls as sets
rather than ball versions of pointwise predicates. Others are
complementary: `arb_contains_zero` tests whether the ball contains the
point 0, so its negation is `arb_is_nonzero`, while the pointwise
predicate "$x = 0$" is `arb_is_zero`:

```c
if (arb_is_zero(x))
{
    ...  /* x = 0 */
}
else if (arb_is_nonzero(x))     /* equivalently: !arb_contains_zero(x) */
{
    ...  /* x != 0 */
}
else
{
    ...  /* unknown */
}
```

The guide [](arb/predicates) lists all predicates and their exact
meaning, including the treatment of NaN and infinities.

## A worked example: the sine function {#example}

We implement $\sin x$ naively by its Taylor series $\sum_{k=0}^{\infty}
(-1)^k x^{2k+1}/(2k+1)!$. The infinite series must be split into a finite
sum, which ball arithmetic evaluates rigorously, and a tail, which has to
be bounded by hand. We stop at the first term $t$ with $|t| \le
2^{-\mathit{prec}} < 1$. From that point on the terms alternate and
decrease in magnitude, so the tail is bounded by $|t|$, and we add this
bound to the radius of the sum with `arb_add_error_mag`. All other errors
are handled automatically, so the output is a ball guaranteed to contain
$\sin x$.

The largest term of the series for $\sin x$ has magnitude about $x^x/x!$,
or $10^{873}$ for $x = 2016.1$, and the final sum is less than 1, so about
873 digits (2900 bits) are lost to cancellation. The program does not know
this in advance. It *guesses* a precision, uses ball arithmetic to
*verify* whether the result is accurate enough to determine the sign of
$\sin x$, and doubles the precision until it is:

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

At 4096 bits the sign is known; the radius $1.2 \cdot 10^{-11}$ printed
there is the error of rounding the midpoint to 10 decimal digits, not the
error of the computation. The built-in `arb_sin` reduces the argument
modulo $\pi/4$ before summing a series, so it needs only 64 bits for the
same input. When a computation can be expressed with built-in functions,
all error bounds are automatic; only infinite processes such as series
need a hand-made bound. As an exercise, the reader can improve the naive
function by subtracting a suitable multiple of $2\pi$ from $x$ first
(using `arb_const_pi`, `arb_div` and `arf_get_fmpz`), after which 64 bits
suffice for $\sin(2016.1)$ and with a little more work also for
$\sin(e^{2016.1})$.

## Guess and verify {#guess-and-verify}

The loop above follows the *guess-and-verify* strategy that is typical of
ball arithmetic: instead of determining a priori the precision needed for
a correct result, guess a precision, use ball arithmetic to verify that
the result is accurate enough, and repeat with a higher precision (or
report failure) otherwise. If rounding errors are thought of as random, a
floating-point computation is like a Monte Carlo algorithm, which may
return a wrong answer; ball arithmetic turns it into a Las Vegas
algorithm, which never returns a wrong answer but may fail to terminate.

The loop does not terminate if the sign to be determined is that of an
expression that is exactly zero, such as $\sin \pi$:

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

The sign of a nonzero real number can always be found by computing it
accurately enough, but the sign of an expression that is exactly zero
cannot be decided numerically, unless the whole computation happens to be
exact, as with `arb_sin_pi` at an exact integer. A program must choose a
stopping criterion that suits its application: stopping when $|\sin x| <
10^{-10000}$ is clear would give useful information about the input $x =
\pi$, although it would not prove that $\sin \pi = 0$. Deciding such
questions exactly is the purpose of [exact real numbers](guide/exact).
