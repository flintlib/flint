---
title: The Riemann zeta function, its zeros and L-functions
summary: Certified zeros of the Riemann zeta function on the critical line, the Hardy Z-function, Laurent and Taylor coefficients, and Dirichlet L-functions with acb_dirichlet; the example programs zeta_zeros, lvalue, lcentral and keiper_li.
order: 9
---

The Riemann zeta function is the most studied special function in
computational number theory, and [](acb_dirichlet) contains the tools to
compute with it rigorously: values and derivatives of $\zeta(s)$ and of
Dirichlet $L$-functions $L(s, \chi)$ anywhere in the complex plane, and
the zeros on the critical line, which are isolated with a proof that none
is missed. The algorithms are described in
[](algorithms/special-functions/zeta).

## Zeros on the critical line {#zeros}

The nontrivial zeros of $\zeta(s)$ lie in the critical strip
$0 < \operatorname{Re}(s) < 1$, and the Riemann hypothesis asserts that
they lie on the line $\operatorname{Re}(s) = 1/2$. On this line,
$$
Z(t) = e^{i \theta(t)} \zeta(\tfrac12 + it)
$$
is real for real $t$, where $\theta(t)$ is the Riemann–Siegel theta
function, so the zeros $\tfrac12 + i t$ on the line correspond to the
real zeros of the Hardy $Z$-function, which can be found from sign
changes. Sign changes only show that there are *at least* so many zeros;
to know that all zeros in the strip up to a given height are on the line
and have been found, one compares with the number of zeros $N(T)$ with
$0 < \operatorname{Im}(s) \le T$ in the whole strip, which is known
through the argument principle. In practice this is done with Gram
points (where $\theta(t)$ is a multiple of $\pi$, so that sign changes
are expected between consecutive ones) and Turing's method [@Tur1953],
[@Leh1970], [@Tru2014], which certifies $N(T)$ from the signs of $Z$ at a
few Gram points beyond $T$, without evaluating $\zeta$ off the line.

`acb_dirichlet_hardy_z_zeros(res, n, len, prec)` uses this to compute
the zeros number $n, n + 1, \ldots, n + \text{len} - 1$ (counted from 1 by
increasing height; the zeros are simple, which is part of what the
computation verifies), as balls `arb_t` of about `prec` bits of
accuracy containing the imaginary parts. It isolates them as just
described, evaluating $Z(t)$ with the Riemann–Siegel formula, and refines
each one by sign-checked secant and bisection steps. At large heights
($10^{11}$ and beyond) it may instead use Platt's method [@Pla2017],
which isolates thousands of consecutive zeros at once with a fast
multipoint evaluation of $Z$, when a cost model predicts that this is
faster.

## Values, derivatives and L-functions {#values}

The program below computes three zeros near the bottom of the critical
strip and three near the millionth zero, evaluates $Z$ and $Z'$ at the
first zero, computes the first Taylor coefficients of
$\zeta(s) - 1/(s - 1)$ at $s = 1$, and evaluates the Dirichlet
$L$-function of the nontrivial character modulo 4 at $s = 2$, which is
Catalan's constant. The two last results are checked against the
constants of [](arb).

```c example name=numerics_zeta
#include <flint/fmpz.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_dirichlet.h>

int main(void)
{
    arb_ptr zeros;
    acb_ptr z;
    acb_t s;
    arb_t t;
    fmpz_t n;
    dirichlet_group_t G;
    dirichlet_char_t chi;
    slong i, prec = 128;

    zeros = _arb_vec_init(3);
    z = _acb_vec_init(3);
    acb_init(s);
    arb_init(t);
    fmpz_init(n);

    /* zeros number 1, 2, 3 and 10^6, 10^6 + 1, 10^6 + 2 */
    fmpz_one(n);
    acb_dirichlet_hardy_z_zeros(zeros, n, 3, prec);
    for (i = 0; i < 3; i++)
    {
        flint_printf("rho_%wd = 1/2 + i ", i + 1);
        arb_printn(zeros + i, 30, 0);
        flint_printf("\n");
    }
    fmpz_set_ui(n, 1000000);
    acb_dirichlet_hardy_z_zeros(zeros, n, 3, prec);
    for (i = 0; i < 3; i++)
    {
        flint_printf("rho_%wd = 1/2 + i ", 1000000 + i);
        arb_printn(zeros + i, 30, 0);
        flint_printf("\n");
    }

    /* Z(t) and Z'(t) at the first zero */
    fmpz_one(n);
    acb_dirichlet_hardy_z_zeros(zeros, n, 1, prec);
    acb_set_arb(s, zeros + 0);
    acb_dirichlet_hardy_z(z, s, NULL, NULL, 2, prec);
    flint_printf("Z(t_1)  = "); acb_printn(z + 0, 10, 0); flint_printf("\n");
    flint_printf("Z'(t_1) = "); acb_printn(z + 1, 30, 0); flint_printf("\n");

    /* Taylor coefficients of zeta(s) - 1/(s-1) at s = 1 */
    acb_one(s);
    acb_dirichlet_l_jet(z, s, NULL, NULL, 1, 3, prec);
    for (i = 0; i < 3; i++)
    {
        flint_printf("[(s-1)^%wd] = ", i);
        acb_printn(z + i, 30, 0);
        flint_printf("\n");
    }
    arb_const_euler(t, prec);
    flint_printf("equals Euler's constant: %d\n", arb_overlaps(acb_realref(z), t));

    /* L(2, chi) for the nontrivial character mod 4 is Catalan's constant */
    dirichlet_group_init(G, 4);
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, 3);
    acb_set_ui(s, 2);
    acb_dirichlet_l(z, s, G, chi, prec);
    flint_printf("L(2, chi_4) = "); acb_printn(z, 30, 0); flint_printf("\n");
    arb_const_catalan(t, prec);
    flint_printf("equals Catalan's constant: %d\n", arb_overlaps(acb_realref(z), t));
    dirichlet_char_clear(chi);
    dirichlet_group_clear(G);

    _arb_vec_clear(zeros, 3);
    _acb_vec_clear(z, 3);
    acb_clear(s);
    arb_clear(t);
    fmpz_clear(n);
    flint_cleanup();
    return 0;
}
```

```output
rho_1 = 1/2 + i [14.1347251417346937904572519836 +/- 3.76e-29]
rho_2 = 1/2 + i [21.0220396387715549926284795939 +/- 3.10e-30]
rho_3 = 1/2 + i [25.0108575801456887632137909926 +/- 3.72e-29]
rho_1000000 = 1/2 + i [600269.677012444955521233914270 +/- 4.91e-25]
rho_1000001 = 1/2 + i [600270.301090711698660098258485 +/- 1.76e-25]
rho_1000002 = 1/2 + i [600270.747870594366133557973222 +/- 1.38e-25]
Z(t_1)  = [+/- 8.06e-37]
Z'(t_1) = [0.793160433356506116013897565274 +/- 3.53e-31]
[(s-1)^0] = [0.577215664901532860606512090082 +/- 4.03e-31]
[(s-1)^1] = [0.0728158454836767248605863758749 +/- 1.32e-33]
[(s-1)^2] = [-0.00484518159643615924226519301761 +/- 3.74e-33]
equals Euler's constant: 1
L(2, chi_4) = [0.915965594177219015054603514932 +/- 3.85e-31]
equals Catalan's constant: 1
```

Some points worth noting:

- The zero numbered $10^6$ is computed directly, without finding the
  zeros below it: the isolation only needs Gram points near the target
  and the count $N(T)$ from Turing's method.
- $Z(t_1)$ evaluated on the ball for the first zero is a ball around
  zero, as it must be; its tiny radius reflects the accuracy of the ball
  for $t_1$.
- `acb_dirichlet_hardy_z` and `acb_dirichlet_l_jet` compute `len` Taylor
  coefficients, $f(s), f'(s), f''(s)/2!, \ldots$ With the `deflate` flag,
  `acb_dirichlet_l_jet` removes the pole at $s = 1$, so the coefficients
  above are those of the Laurent expansion
  $\zeta(s) = \frac{1}{s - 1} + \sum_{k \ge 0} \frac{(-1)^k \gamma_k}{k!}
  (s - 1)^k$, with the Stieltjes constants $\gamma_0 = \gamma$,
  $\gamma_1 = -0.0728\ldots$, $\gamma_2 = -0.0096\ldots$ (see also
  `acb_dirichlet_stieltjes`, which computes $\gamma_k$ directly, also for
  large $k$ [@JB2018]).
- Passing `NULL` for the group and the character selects the Riemann
  zeta function. A Dirichlet character is given by a group
  `dirichlet_group_t` for the modulus $q$ and a character
  `dirichlet_char_t`; `dirichlet_char_log(chi, G, m)` selects the
  character with Conrey label $(q, m)$, here $\chi_4(3) = -1$. See
  [](dirichlet) for characters and [](acb_dirichlet) for the
  $L$-functions.

## The example programs {#programs}

`examples/zeta_zeros.c` prints consecutive zeros of $Z(t)$, with their
index:

```
zeta_zeros [-n n] [-count n] [-prec n] [-digits n] [-threads n] [-platt] [-noplatt] [-v] [-verbose] [-h] [-help]
```

It starts at the zero number `-n` (default 1) and prints `-count` zeros
(without `-count` it runs until interrupted). The precision is given in
bits with `-prec` or in digits with `-digits`; by default it is
$64 + \lceil \log_2 n \rceil$ bits, enough to separate zeros at that
height comfortably. By default the program lets
`acb_dirichlet_hardy_z_zeros` choose the method; `-noplatt` forces the
Riemann–Siegel method and `-platt` the large height method (valid from
the 10000th zero, with the precision then used as a working precision),
which is best for many consecutive zeros at a large height. `-threads`
uses several threads, and `-v` prints the parameters chosen.

```
> build/examples/zeta_zeros -n 1 -count 5 -digits 30
1	14.1347251417346937904572519836
2	21.0220396387715549926284795939
3	25.0108575801456887632137909926
4	30.4248761258595132103118975306
5	32.9350615877391896906623689641
cpu/wall(s): 0.008 0.008
virt/peak/rss/peak: 24.8 MB / 24.8 MB / 9.30 MB / 9.30 MB
> build/examples/zeta_zeros -n 1048449114 -count 2 -digits 25
1048449114	388858886.0022851217767971
1048449115	388858886.0023936897027167
cpu/wall(s): 0.114 0.133
virt/peak/rss/peak: 26.0 MB / 26.9 MB / 11.4 MB / 11.9 MB
```

The second pair is a close pair of zeros (a Lehmer pair) with a distance
of about $10^{-4}$, more than a thousand times smaller than the average
spacing $2\pi / \log(t / 2\pi) \approx 0.35$ at that height; telling them apart requires evaluating $Z$ near its tiny
local extremum between them, which the program does without being told.
The printed digits are correct (the digits are printed with
`ARB_STR_NO_RADIUS`, so the last digit may be off by one unit).

`examples/lvalue.c` evaluates $L(s, \chi)$, or with `-z` the
$Z$-function, and optionally its Taylor series, at a point $s = a + bi$:

```
lvalue [-character q n] [-re a] [-im b] [-prec p] [-z] [-deflate] [-len l]
```

The character is given by its Conrey label $(q, n)$ (default $(1, 0)$,
the Riemann zeta function), the point by `-re` and `-im` (default
$s = 1/2$), and the precision in bits (default 53). With `-len l` it
prints $l$ Taylor coefficients, and with `-deflate` it removes the pole
at $s = 1$ as above:

```
> build/examples/lvalue -character 4 3 -re 2 -prec 128
L(s) = [0.91596559417721901505460351493238411077 +/- 7.86e-39]
> build/examples/lvalue -re 1 -deflate -len 4
L(s) = [0.577215664901532861 +/- 7.54e-19]
L'(s) = [0.072815845483676725 +/- 5.70e-19]
[x^2] L(s+x) = [-0.004845181596436159 +/- 8.03e-19]
[x^3] L(s+x) = [-0.000342305736717224 +/- 6.97e-19]
> build/examples/lvalue -z -re 14.134725 -len 3
Z(s) = [-1.1241835e-7 +/- 3.84e-15]
Z'(s) = [0.7931604148843 +/- 5.24e-14]
[x^2] Z(s+x) = [0.0651645864916 +/- 6.45e-14]
```

(each followed by timing and memory lines).

`examples/lcentral.c` computes the central values $L(1/2, \chi)$ for all
primitive characters $\chi$ modulo each $q$ in a range,
`lcentral [--quiet] [--check] [--prec bits] qmin qmax` (default 100
bits). It is much faster than calling `acb_dirichlet_l` for each
character: all $L(1/2, \chi)$ modulo $q$ are obtained from the Hurwitz
zeta values $\zeta(1/2, k/q)$ by a discrete Fourier transform
(`acb_dirichlet_l_vec_hurwitz`), and the Hurwitz zeta values for all
moduli in the range are computed from a common precomputation
(`acb_dirichlet_hurwitz_precomp_init_num`). With `--check`, the program
verifies that no computed central value contains zero (it is
conjectured that $L(1/2, \chi) \ne 0$ for all primitive $\chi$); the check can only fail when the precision is too low, in which case
the program reports `Value could be zero!` and aborts.

```
> build/examples/lcentral 1 5
3,2: [0.48086755769682862618122006324 +/- 7.42e-30] + [+/- 2.45e-30]*I
4,3: [0.66769145718960917665869092930 +/- 1.62e-30] + [+/- 9.34e-31]*I
5,2: [0.76374788011728687822451215264 +/- 2.90e-30] + [0.21696476751886069363858659310 +/- 3.82e-30]*I
5,4: [0.23175094750401575588338366176 +/- 2.12e-30] + [+/- 1.44e-30]*I
5,3: [0.76374788011728687822451215264 +/- 2.90e-30] + [-0.21696476751886069363858659310 +/- 3.82e-30]*I
```

For real characters the value is real, but the computation through the
complex Fourier transform gives a ball for the imaginary part that only
contains zero.

`examples/keiper_li.c` computes the Keiper–Li coefficients
$\lambda_0, \ldots, \lambda_n$: with
$\xi(s) = \tfrac12 s (s - 1) \pi^{-s/2} \Gamma(s/2) \zeta(s)$, they are
defined by $\lambda_n = \frac{1}{(n-1)!} \frac{d^n}{ds^n}
\left[s^{n-1} \log \xi(s)\right]_{s=1}$ for $n \ge 1$, and Li's criterion
states that the Riemann hypothesis holds if and only if $\lambda_n > 0$
for all $n \ge 1$ [@Kei1992], [@Li1997]. The program computes the power
series of $\log \xi(s)$ at $s = 0$ from those of $\zeta(s)$
(`_arb_poly_zeta_series`), $\log \Gamma(1 + s/2)$
(`_arb_poly_lgamma_series`) and $\log(1 - s)$, and obtains the
$\lambda_n$ with a binomial transform (`_arb_poly_binomial_transform`).
Since there is heavy cancellation, the default working precision is
about $1.1 n + 50$ bits. It was used for the record computations of
[@Joh2013]:

```
keiper_li n [-prec prec] [-threads num_threads] [-out out_file]
```

`-out` writes all coefficients, as exact midpoints and radii, to a file.
For $n = 1000$ the program prints the first and last coefficients with
their error bounds, from $\lambda_1 = 0.0230957\ldots$ to
$\lambda_{1000} = 2.32605\ldots$.

The zeros can also be found with the general root isolation of
[](arb_calc), as in `real_roots 0` (see
[](tutorials/numerics/real-roots)), which is slower but shows the method
in its simplest form, and counted with the argument principle and
numerical integration, as in `integrals -i 22` (see
[](tutorials/numerics/integrals)).
