---
title: Computing π and other constants
summary: Precision versus accuracy, reading and printing balls, the built-in constants, and an AGM iteration for π with an explicit truncation error bound; the example programs pi and pi_agm.
order: 1
---

The simplest rigorous computation is a single call to a function that
returns a mathematical constant. This tutorial uses it to introduce the
basic notions of ball arithmetic: the working precision, the accuracy of
a result, and how balls are printed. It then computes $\pi$ with an
iteration of our own, where the error that ball arithmetic does not see
(the truncation error of the iteration) has to be added by hand.

## Precision and accuracy {#precision}

Every function in [](arb) takes a *working precision* `prec` in bits.
It controls the precision of the floating-point operations performed
internally, and most functions aim at a result with a relative error of
about $2^{-\text{prec}}$. The result is a ball `arb_t` with a midpoint and
a radius, and the true value is guaranteed to lie in the ball. The
*accuracy* of the result is measured on the result itself:
`arb_rel_accuracy_bits` returns about $-\log_2$ of the relative radius.
It is usually a little smaller than the precision, because of rounding
errors, and it can be much smaller when there is cancellation or when
the input balls are wide.

For $d$ decimal digits, a precision of $d \log_2 10 \approx 3.32 d$ bits
plus a few guard bits suffices for the constants below. The following
program computes several constants to 40 digits and $\pi$ to 10000
digits.

```c example name=numerics_constants_builtin
#include <flint/arb.h>

int main(void)
{
    arb_t x;
    slong digits = 40;
    slong prec = digits * 3.3219280948873623 + 5;   /* log2(10) bits per digit */

    arb_init(x);

    arb_const_pi(x, prec);
    flint_printf("prec = %wd bits, accuracy = %wd bits\n",
        prec, arb_rel_accuracy_bits(x));
    flint_printf("pi    = "); arb_printn(x, digits, 0); flint_printf("\n");
    flint_printf("pi    = "); arb_printn(x, 15, 0); flint_printf("\n");
    flint_printf("pi    ~ "); arb_printn(x, 15, ARB_STR_NO_RADIUS); flint_printf("\n");

    arb_const_e(x, prec);
    flint_printf("e     = "); arb_printn(x, digits, 0); flint_printf("\n");
    arb_const_log2(x, prec);
    flint_printf("log 2 = "); arb_printn(x, digits, 0); flint_printf("\n");
    arb_const_euler(x, prec);
    flint_printf("gamma = "); arb_printn(x, digits, 0); flint_printf("\n");
    arb_const_catalan(x, prec);
    flint_printf("G     = "); arb_printn(x, digits, 0); flint_printf("\n");
    arb_zeta_ui(x, 3, prec);
    flint_printf("z(3)  = "); arb_printn(x, digits, 0); flint_printf("\n");

    /* 10000 digits, printed condensed */
    arb_const_pi(x, 10000 * 3.3219280948873623 + 5);
    arb_printn(x, 10000, ARB_STR_CONDENSE * 10);
    flint_printf("\n");

    arb_clear(x);
    flint_cleanup();
    return 0;
}
```

```output
prec = 137 bits, accuracy = 136 bits
pi    = [3.141592653589793238462643383279502884197 +/- 1.81e-40]
pi    = [3.14159265358979 +/- 3.24e-15]
pi    ~ 3.14159265358979
e     = [2.718281828459045235360287471352662497757 +/- 2.60e-40]
log 2 = [0.6931471805599453094172321214581765680755 +/- 3.58e-42]
gamma = [0.5772156649015328606065120900824024310422 +/- 4.14e-41]
G     = [0.9159655941772190150546035149323841107741 +/- 5.21e-41]
z(3)  = [1.202056903159594285399738161511449990765 +/- 2.58e-41]
[3.1415926535{...9979 digits...}5525637568 +/- 1.72e-10000]
```

## Reading the output {#printing}

`arb_printn(x, n, flags)` prints a decimal ball that contains the binary
ball `x`, with at most $n$ digits of the midpoint. The conversion is
exact: the printed radius is the true radius plus the error of rounding
the midpoint to decimal, rounded up. The number of digits is chosen so
that the printed midpoint is correct to within one unit in the last
place. Thus `[3.14159265358979 +/- 3.24e-15]` says that
$3.14159265358979 - 3.24 \cdot 10^{-15} \le \pi \le 3.14159265358979 +
3.24 \cdot 10^{-15}$; the radius is large compared with the binary radius
of about $10^{-41}$ only because the midpoint was rounded to 15 digits.

With `ARB_STR_NO_RADIUS`, only the midpoint is printed, still rounded so
that it is correct to within one unit in the last place; this is the
right format for tables of digits. `ARB_STR_CONDENSE * m` replaces long
runs of digits by the first and last $m$ digits and a count, which is
how the example programs print millions of digits. The same strings are
available from `arb_get_str`, and `arb_printd` prints the midpoint and
radius as floating-point numbers without any guarantee on the digits,
which is useful for debugging.

When the radius is large compared with the midpoint, there may be no
correct digit at all, and the ball is printed as `[+/- 2.72]`, meaning
only that $|x| \le 2.72$. This is the typical symptom of insufficient
precision or of a formula that loses accuracy, and later tutorials show
how to deal with it.

## How the constants are computed {#algorithms}

The constants are computed by fast convergent series evaluated with
binary splitting [@HP1998] (see
[](algorithms/integers/binary-splitting)), which makes the cost
quasi-linear in the precision:

- $\pi$ (`arb_const_pi`): the Chudnovsky series [@CC1988] (below
  about 4600 bits, a precomputed table is rounded instead);
- $e$ (`arb_const_e`): $\sum_k 1/k!$;
- $\log 2$ (`arb_const_log2`), Catalan's constant `arb_const_catalan`
  and $\zeta(3)$ (`arb_const_apery`, used by `arb_zeta_ui` for $n = 3$):
  short hypergeometric series with a fast rate of convergence;
- Euler's constant $\gamma$ (`arb_const_euler`): the Brent–McMillan
  algorithm [@BM1980] with the error bound of [@BJ2013].

The details are in [](algorithms/special-functions/constants). Each
constant is cached: the first call at a given precision computes it, and
later calls at the same or lower precision only round the cached value.
The cache is local to the thread and is freed by `flint_cleanup`, which
is why the programs in these tutorials call it before returning (this keeps
memory checkers quiet; it is not required for correctness).

## Our own π: the AGM iteration {#agm}

The arithmetic-geometric mean gives an iteration for $\pi$ that converges
quadratically: the number of correct digits roughly doubles in each step
[@Sal1976], [@Bre1976agm]. Start with $a_0 = 1$, $b_0 = 1/\sqrt{2}$,
$s_0 = 2 b_0 - \tfrac12$ and set
$$
a_k = \frac{a_{k-1} + b_{k-1}}{2}, \qquad
b_k = \sqrt{a_{k-1} b_{k-1}}, \qquad
s_k = s_{k-1} - 2^k (a_k - b_k)^2, \qquad
\pi_k = \frac{(a_k + b_k)^2}{s_k}.
$$
Then $\pi_k > \pi$ and $\pi_k - \pi \le 2^{k+8} e^{-\pi 2^{k+1}}$ (this is
the bound used by `examples/pi_agm.c`, which follows Salamin's error
analysis).

If we evaluate the recurrence in ball arithmetic, the ball for $\pi_k$
contains the exact value of $\pi_k$, not $\pi$: the rounding errors are
tracked, but the fact that $\pi_k$ is only an approximation of $\pi$ is
invisible to the arithmetic. This *truncation error* has to be added to
the radius explicitly, here with `arb_add_error_2exp_si`. Forgetting this
step is the most common mistake in rigorous numerics; the result then
looks precise but is wrong. The program below takes enough steps for a
truncation error below $2^{-1000}$, prints each $\pi_k$ to 40 digits (by
midpoint only, to show the convergence), adds the error bound, and
checks the result against `arb_const_pi`.

```c example name=numerics_constants_agm
#include <flint/arb.h>

/* Upper bound for log2 of the truncation error after n AGM steps:
   2^(n+8) exp(-pi 2^(n+1)) <= 2^e with e = n + 8 - (pi/log 2) 2^(n+1).
   We use 4.53236 < pi/log(2) = 4.5323601... and add 1 for safety
   against rounding in the double arithmetic. */
static slong agm_error_exponent(slong n)
{
    return (n + 8) - 4.53236 * (double) (UWORD(1) << (n + 1)) + 1;
}

int main(void)
{
    arb_t a, b, s, t, pi, ref;
    slong k, n, prec = 1000;

    arb_init(a); arb_init(b); arb_init(s); arb_init(t);
    arb_init(pi); arb_init(ref);

    /* the number of steps such that the truncation error is < 2^-prec */
    for (n = 1; agm_error_exponent(n) > -prec; n++) ;

    arb_one(a);                         /* a = 1 */
    arb_sqrt_ui(b, 2, prec);
    arb_mul_2exp_si(b, b, -1);          /* b = sqrt(2)/2 */
    arb_mul_2exp_si(s, b, 1);
    arb_set_d(t, 0.5);
    arb_sub(s, s, t, prec);             /* s = 2b - 1/2 */

    for (k = 1; k <= n; k++)
    {
        arb_add(t, a, b, prec);
        arb_mul_2exp_si(t, t, -1);      /* a' = (a + b)/2 */
        arb_mul(b, a, b, prec);
        arb_sqrt(b, b, prec);           /* b' = sqrt(a b) */
        arb_swap(a, t);
        arb_sub(t, a, b, prec);
        arb_sqr(t, t, prec);
        arb_mul_2exp_si(t, t, k);
        arb_sub(s, s, t, prec);         /* s' = s - 2^k (a' - b')^2 */

        arb_add(pi, a, b, prec);
        arb_sqr(pi, pi, prec);
        arb_div(pi, pi, s, prec);       /* (a + b)^2 / s */
        flint_printf("step %wd: truncation error < 2^%wd, ", k, agm_error_exponent(k));
        arb_printn(pi, 40, ARB_STR_NO_RADIUS);
        flint_printf("\n");
    }

    /* account for the truncation error */
    arb_add_error_2exp_si(pi, agm_error_exponent(n));

    arb_const_pi(ref, prec);
    flint_printf("overlaps arb_const_pi: %d\n", arb_overlaps(pi, ref));
    flint_printf("accuracy: %wd bits\n", arb_rel_accuracy_bits(pi));
    arb_printn(pi, 300, ARB_STR_CONDENSE * 10);
    flint_printf("\n");

    arb_clear(a); arb_clear(b); arb_clear(s); arb_clear(t);
    arb_clear(pi); arb_clear(ref);
    flint_cleanup();
    return 0;
}
```

```output
step 1: truncation error < 2^-8, 3.141680293297653293918070424560009382796
step 2: truncation error < 2^-25, 3.141592653895446496002914758818043486109
step 3: truncation error < 2^-60, 3.141592653589793238466360602706631321758
step 4: truncation error < 2^-132, 3.141592653589793238462643383279502884197
step 5: truncation error < 2^-276, 3.141592653589793238462643383279502884197
step 6: truncation error < 2^-565, 3.141592653589793238462643383279502884197
step 7: truncation error < 2^-1144, 3.141592653589793238462643383279502884197
overlaps arb_const_pi: 1
accuracy: 994 bits
[3.1415926535{...278 digits...}2602491413 +/- 3.64e-299]
```

Some remarks on the program:

- The error bound is computed in `double` arithmetic. This is safe here
  because it is rounded in the safe direction by construction (a lower
  bound for $\pi/\log 2$ and an extra 1 in the exponent), and because
  the exponent is an integer far from the limits of `double`. In general,
  bounds that enter a proof should be computed with [](mag) or in ball
  arithmetic.
- `arb_overlaps` is the right test for agreement of two enclosures of
  the same number: two correct balls for $\pi$ must overlap, so a result
  of 0 proves a bug. A result of 1 does not prove that both are correct,
  but agreement between two independent methods is strong evidence.
- The 994 accurate bits out of 1000 reflect the accumulated rounding
  errors of 7 steps. `examples/pi_agm.c` uses a working precision of
  $p + 5 + 2 \log_2 p$ for a target of $p$ bits to compensate, and rounds
  the final division to $p$ bits.
- `arb_const_pi` does not use the AGM: although both methods have
  quasi-linear complexity, the Chudnovsky series is faster in practice.

## The example programs {#programs}

`examples/pi.c` computes $\pi$, or another constant, to $n$ digits by
calling the corresponding function at a precision of $n \log_2 10 + 5$
bits, and prints the result condensed with `ARB_STR_CONDENSE`:

```
> build/examples/pi 1000000
precision = 3321933 bits... cpu/wall(s): 0.208 0.279
virt/peak/rss/peak: 35.9 MB / 35.9 MB / 18.7 MB / 18.7 MB
[3.14159265358979323846{...999959 digits...}42209010610577945815 +/- 1.38e-1000000]
```

The options are `-constant name` with `name` one of `pi`, `e`, `log2`,
`euler`, `catalan`, `khinchin`, `zeta3` and `reciprocal_fibonacci`;
`-threads n` to use $n$ threads (binary splitting is parallelised); and
`-condense m` to print the first and last $m$ digits (the default is
20; `-condense 0` prints all digits):

```
> build/examples/pi 100000 -constant euler
precision = 332197 bits... cpu/wall(s): 0.309 0.555
virt/peak/rss/peak: 26.7 MB / 27.1 MB / 9.06 MB / 9.06 MB
[0.57721566490153286060{...99960 digits...}56340316331811093897 +/- 3.80e-100001]
```

`examples/pi_agm.c` implements the AGM iteration above at a precision
with guard bits and prints the truncation error bound of each step:

```
> build/examples/pi_agm 100000
precision = 332197 bits...
Iteration 1 / 16: error bound 0.00391
Iteration 2 / 16: error bound 2.98e-8
Iteration 3 / 16: error bound 8.67e-19
...
Iteration 15 / 16: error bound 4.82e-89409
Iteration 16 / 16: error bound 1.38e-178824
cpu/wall(s): 0.038 0.04
virt/peak/rss/peak: 26.0 MB / 26.0 MB / 8.59 MB / 8.59 MB
[3.14159265358979323846{...99959 digits...}76742080565549362465 +/- 3.91e-100000]
```

The last step is taken because the bound of step 15, about
$2^{-297000}$, is still larger than $2^{-p}$ for the working precision
$p \approx 332000$; the final truncation error is far below the rounding
error, so the radius of the result is dominated by rounding.
