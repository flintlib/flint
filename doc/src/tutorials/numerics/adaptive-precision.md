---
title: Choosing the precision: the logistic map and Hilbert matrices
summary: Restarting with doubled precision until the result is accurate, how fast accuracy is lost in an iteration and why the formula matters, and the determinant of ill-conditioned Hilbert matrices; the example programs logistic and hilbert_matrix.
order: 2
---

The working precision of a ball computation is chosen in advance, but
the accuracy of the result is only known afterwards. When the accuracy
is insufficient, the standard remedy is to repeat the whole computation
at a higher precision, typically doubling it each time. Because the cost
of arithmetic grows at least linearly with the precision, the failed
attempts cost at most about as much as the successful one, so this
strategy is within a small constant factor of knowing the right
precision in advance. This tutorial applies it to two problems where the
needed precision is hard to predict: a chaotic iteration and the
determinant of an ill-conditioned matrix.

## The logistic map {#logistic}

The logistic map $x_{k+1} = r x_k (1 - x_k)$ is chaotic for most $r$
between about $3.57$ and $4$: two orbits that start close together
separate exponentially fast. In floating-point arithmetic, the computed
$x_k$ has no correct digits after a few dozen steps, although nothing in
the output shows it. In ball arithmetic, the radius grows instead, and
`arb_rel_accuracy_bits` shows how many bits are left.

The program below first iterates at a fixed precision of 128 bits and
prints the accuracy every 10 steps, for two ways of writing the same
map. It then computes $x_{100}$ for $x_0 = 1/2$, $r = 15/4$ to 10 digits
(36 bits) by restarting with doubled precision as soon as the accuracy
drops below the goal; checking during the iteration avoids running a
hopeless attempt to the end.

```c example name=numerics_logistic
#include <flint/arb.h>

/* one step x -> r x (1 - x) of the logistic map */
void logistic_step(arb_t x, const arb_t r, int centered, slong prec)
{
    arb_t t;
    arb_init(t);
    if (!centered)
    {
        arb_sub_ui(t, x, 1, prec);
        arb_neg(t, t);
        arb_mul(x, x, t, prec);         /* x (1 - x) */
    }
    else
    {
        arb_set_d(t, 0.5);
        arb_sub(t, x, t, prec);
        arb_sqr(t, t, prec);
        arb_set_d(x, 0.25);
        arb_sub(x, x, t, prec);         /* 1/4 - (x - 1/2)^2 */
    }
    arb_mul(x, x, r, prec);
    arb_clear(t);
}

int main(void)
{
    arb_t x, y, r;
    slong i, n = 100, prec, goal = 36;     /* 36 bits: about 10 digits */

    arb_init(x); arb_init(y); arb_init(r);
    arb_set_d(r, 3.75);

    /* how the accuracy decays at a fixed precision */
    arb_set_d(x, 0.5);
    arb_set_d(y, 0.5);
    flint_printf("step   accurate bits at prec = 128\n");
    for (i = 1; i <= 50; i++)
    {
        logistic_step(x, r, 0, 128);
        logistic_step(y, r, 1, 128);
        if (i % 10 == 0)
            flint_printf("%4wd   %4wd   %4wd\n", i,
                arb_rel_accuracy_bits(x), arb_rel_accuracy_bits(y));
    }

    /* restart with doubled precision until x_n is accurate enough */
    for (prec = 64; ; prec *= 2)
    {
        arb_set_d(x, 0.5);
        for (i = 0; i < n; i++)
        {
            logistic_step(x, r, 0, prec);
            if (arb_rel_accuracy_bits(x) < goal)
                break;
        }

        if (i == n)
        {
            flint_printf("prec = %wd: success\n", prec);
            break;
        }
        flint_printf("prec = %wd: lost accuracy at step %wd\n", prec, i);
    }

    flint_printf("x_%wd = ", n);
    arb_printn(x, 10, 0);
    flint_printf("\n");

    arb_clear(x); arb_clear(y); arb_clear(r);
    flint_cleanup();
    return 0;
}
```

```output
step   accurate bits at prec = 128
  10    118    121
  20     99    115
  30     79    109
  40     59    100
  50     40     94
prec = 64: lost accuracy at step 18
prec = 128: lost accuracy at step 53
prec = 256: success
x_100 = [0.8882939923 +/- 1.60e-11]
```

### Why the formula matters {#dependency}

With the formula $r x (1 - x)$ (second column), about 2 bits are lost per
step; with the equivalent $r(\tfrac14 - (x - \tfrac12)^2)$ (third column),
only about 0.6. The reason is the *dependency problem* of interval and
ball arithmetic: each operation bounds its error for all values in its
input balls independently, without knowing that $x$ and $1 - x$ are
correlated. If $x = [m \pm \varepsilon]$, then $x (1 - x)$ is computed
as $[m \pm \varepsilon] \cdot [1 - m \pm \varepsilon]$, whose radius is about $(|m| + |1 - m|) \varepsilon =
\varepsilon$ for $0 < m < 1$, so each step multiplies the radius by about
$r = 3.75$, or 1.9 bits. The true sensitivity of one step is the
derivative $|r (1 - 2m)|$, which is often much smaller. In the centred
form, `arb_sqr` treats its input as a single variable and returns the
radius $2 |m - \tfrac12| \varepsilon + \varepsilon^2$ for the square of
$x - \tfrac12$, so the radius grows by the factor $|r(1 - 2m)|$ per step
(up to second-order terms and rounding), which is the rate at which
nearby orbits actually separate.

The general lesson: when a ball is much wider than expected, look for a
variable that occurs several times in a formula, and rewrite the formula
so that it occurs once, or so that the dependent part is handled by a
single operation (`arb_sqr` for a square; `arb_expm1` instead of
`arb_exp` followed by a subtraction; library functions that evaluate a
whole expression, such as the special functions, rather than a formula
assembled from elementary operations). The precision can compensate for
any constant overestimation factor per step, but at the cost of more
bits: here, about three times as many.

### Inexact input {#inexact-input}

Here $x_0$ and $r$ are exact binary numbers, so the computed ball
contains the exact $x_{100}$. If the parameters were decimal numbers such
as $r = 3.7$, which are not exactly representable,
`arb_set_str(r, "3.7", prec)` would give a ball of radius about $2^{-\text{prec}}$
containing $3.7$, and the result would contain $x_{100}$ for all $r$ in
that ball, in particular for $r = 3.7$. Each restart must then recreate
the input at the new precision, as `examples/logistic.c` does, otherwise
the input error limits the accuracy no matter how high the working
precision.

## The determinant of the Hilbert matrix {#hilbert}

The Hilbert matrix $H_n = (1/(i + j - 1))_{i,j=1}^{n}$ is the classical
example of an ill-conditioned matrix: its condition number grows like
$e^{3.5 n}$, and its determinant
$$
\det H_n = \frac{c_n^4}{c_{2n}}, \qquad c_n = \prod_{i=1}^{n-1} i!,
$$
decreases roughly like $4^{-n^2}$. To compute $\det H_n$ in ball
arithmetic we create $H_n$ at precision $p$ with `arb_mat_hilbert` (the
entries are balls containing $1/(i+j-1)$, so the computed determinant
contains $\det H_n$) and call `arb_mat_det`. Gaussian elimination loses
about $\log_2$ of the condition number in bits, so the result is
useless until the precision exceeds that. The program doubles the
precision until the determinant has 30 accurate bits, and then checks the
result against the exact value computed with [](fmpq).

```c example name=numerics_hilbert_det
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>

/* exact determinant of the n x n Hilbert matrix: c(n)^4 / c(2n),
   where c(n) = 1! 2! ... (n-1)! */
void hilbert_det_exact(fmpq_t res, slong n)
{
    fmpz_t c, f;
    slong i;
    fmpz_init(c); fmpz_init(f);
    fmpz_one(c);
    fmpz_one(f);
    fmpq_one(res);
    for (i = 1; i < 2 * n; i++)
    {
        fmpz_mul_ui(f, f, i);                   /* f = i! */
        if (i < n)
            fmpz_mul(c, c, f);                  /* c = c(n) at the end */
        fmpz_mul(fmpq_denref(res), fmpq_denref(res), f);
    }
    fmpz_pow_ui(fmpq_numref(res), c, 4);
    fmpq_canonicalise(res);
    fmpz_clear(c); fmpz_clear(f);
}

int main(void)
{
    arb_mat_t A;
    arb_t det;
    fmpq_t exact;
    slong n = 40, prec;

    arb_mat_init(A, n, n);
    arb_init(det);
    fmpq_init(exact);

    for (prec = 32; ; prec *= 2)
    {
        arb_mat_hilbert(A, prec);           /* entries 1/(i+j+1), i, j >= 0 */
        arb_mat_det(det, A, prec);
        flint_printf("prec = %3wd: ", prec);
        arb_printn(det, 10, 0);
        flint_printf("\n");
        if (arb_rel_accuracy_bits(det) >= 30)
            break;
    }

    hilbert_det_exact(exact, n);
    flint_printf("contains the exact value: %d\n", arb_contains_fmpq(det, exact));

    arb_mat_clear(A);
    arb_clear(det);
    fmpq_clear(exact);
    flint_cleanup();
    return 0;
}
```

```output
prec =  32: [+/- 1.53e-128]
prec =  64: [+/- 1.19e-256]
prec = 128: [+/- 3.42e-492]
prec = 256: [1.097111417e-932 +/- 3.27e-944]
contains the exact value: 1
```

The failed attempts are not useless: `[+/- 1.53e-128]` proves that
$|\det H_{40}| \le 1.53 \cdot 10^{-128}$. A ball that contains zero
cannot determine the sign, though, and for many applications (is the
matrix singular? is a quadratic form definite?) the first precision at
which the ball excludes zero is what matters; `examples/hilbert_matrix.c`
stops there. For large matrices at moderate precision, `arb_mat_det`
uses a preconditioned method rather than Gaussian elimination in ball
arithmetic; see [](arb_mat) for the available methods.

## The example programs {#programs}

`examples/logistic.c` computes the $n$-th iterate of the logistic map,
`logistic n [x0] [r] [digits]` with defaults $x_0 = 0.5$, $r = 3.75$ and
10 digits, with the restart strategy above (and the formula
$r x (1 - x)$). The inputs are parsed with `arb_set_str` at each
precision:

```
> build/examples/logistic 1234 0.1 3.99 30
Trying prec=64 bits...ran out of accuracy at step 0
Trying prec=128 bits...ran out of accuracy at step 10
Trying prec=256 bits...ran out of accuracy at step 76
Trying prec=512 bits...ran out of accuracy at step 205
Trying prec=1024 bits...ran out of accuracy at step 461
Trying prec=2048 bits...ran out of accuracy at step 974
Trying prec=4096 bits...success!
cpu/wall(s): 0.006 0.007
x_1234 = [0.256445391958651410579677945635 +/- 3.92e-31]
```

At 64 bits the goal of 30 digits (103 bits) is not reached even at step
0. For `logistic 10000` the program needs 32768 bits.

`examples/hilbert_matrix.c` computes $\det H_n$, doubling the precision
from 20 bits until the ball does not contain zero:

```
> build/examples/hilbert_matrix 200
prec=20: [+/- 1.32e-335]
prec=40: [+/- 1.63e-545]
prec=80: [+/- 1.30e-933]
prec=160: [+/- 3.62e-1926]
prec=320: [+/- 1.81e-4129]
prec=640: [+/- 3.84e-8838]
prec=1280: [2.955454297e-23924 +/- 8.29e-23935]
success!
```

With `-eig`, it computes instead the smallest eigenvalue of $H_n$: it
computes approximate eigenvalues and eigenvectors with
`acb_mat_approx_eig_qr`, certifies them with `acb_mat_eig_simple` (which
proves that the eigenvalues are simple and encloses each one), and since
$H_n$ is symmetric sets the imaginary parts to zero and sorts:

```
> build/examples/hilbert_matrix -eig 10
prec=20: nan
prec=40: nan
prec=80: nan
prec=160: [1.093153819e-13 +/- 3.80e-23]
success!
```

Here `nan` means that the eigenvalues could not be isolated at that
precision; the output is an indeterminate ball rather than a wrong
number.
