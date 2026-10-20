---
title: acb_modular.h
header: acb_modular.h
description: The modular group, Jacobi theta functions, the Dedekind eta function, modular forms and Hilbert class polynomials.
category: Special functions
order: 4
---

This module evaluates modular forms and the Jacobi theta functions of a
variable $\tau$ in the upper half-plane
$\mathbb{H} = \{\tau \in \mathbb{C} : \operatorname{Im}(\tau) > 0\}$ with rigorous
error bounds: the theta functions $\theta_1, \ldots, \theta_4$ and their
derivatives, the Dedekind eta function, Klein's $j$-invariant, the
modular lambda function, the discriminant $\Delta$ and the Eisenstein
series. It also provides arithmetic in the modular group
$\mathrm{PSL}_2(\Z)$, reduction to the fundamental domain, and Hilbert
class polynomials.

Every evaluation first moves $\tau$ to (or close to) the fundamental
domain with a modular transformation, so that $\operatorname{Im}(\tau) \ge
\sqrt3/2$ and $|q| \le e^{-\pi\sqrt3/2} \approx 0.066$ with
$q = e^{\pi i \tau}$, and then sums a $q$-series with an addition
sequence or rectangular splitting [@EHJ2016]. At precision $p$ this needs
$O(\sqrt{p})$ terms, and the cost does not depend on how close $\tau$
was to the real line. The algorithms are described in
[](algorithms/special-functions/elliptic-modular) and
[](algorithms/special-functions/theta-functions); see also [@Joh2019].
Theta functions in any genus are in [`acb_theta`](acb_theta), elliptic
functions and integrals in [`acb_elliptic`](acb_elliptic).

::: summary algorithms/special-functions/elliptic-modular :::

## Conventions {#conventions}

- $\tau$ is a complex ball, normally in $\mathbb{H}$. If the ball for
  $\tau$ is not certainly in the upper half-plane (it touches the real
  line), functions defined only on $\mathbb{H}$ return balls with
  infinite radius or indeterminate values. The same holds for the
  summation functions that take $q$ and require $|q| < 1$.
- The variable $q$ is $e^{2\pi i \tau}$ for modular forms
  (`acb_modular_eta_sum`) and $e^{\pi i \tau}$ for theta functions
  (`acb_modular_theta_sum`, `acb_modular_theta_const_sum`); each function
  says which.
- The theta functions are those of [](#theta) and [@NIST2012, §20.2],
  with the variable $z$ multiplied by $\pi$: here
  $\theta_3(z, \tau) = \sum_{n \in \Z} e^{\pi i (n^2 \tau + 2nz)}$, so
  the period in $z$ is 1, and derivatives are taken with respect to this
  $z$.
- Outputs may alias inputs unless stated otherwise.

## The modular group {#psl2z}

::: type
psl2z_struct
psl2z_t

An element of the modular group $\mathrm{PSL}_2(\Z)$, stored as an
integer matrix $\begin{pmatrix} a & b \\ c & d \end{pmatrix}$ with
$ad - bc = 1$, in the canonical form $c \ge 0$, and $d > 0$ if $c = 0$.
The members `a`, `b`, `c`, `d` are of type `fmpz`. A `psl2z_t` is an
array of length 1, so that it is passed by reference.
:::

::: function
void psl2z_init(psl2z_t g)
void psl2z_clear(psl2z_t g)

Initialise `g` to the identity, and free it.
:::

::: function
void psl2z_swap(psl2z_t f, psl2z_t g)
void psl2z_set(psl2z_t f, const psl2z_t g)
void psl2z_one(psl2z_t g)

Swap `f` and `g`, set `f` to a copy of `g`, and set `g` to the identity.
:::

::: function
int psl2z_is_one(const psl2z_t g)
int psl2z_equal(const psl2z_t f, const psl2z_t g)
int psl2z_is_correct(const psl2z_t g)

Return whether `g` is the identity, whether `f` and `g` are equal, and
whether `g` is a valid element in canonical form ($ad - bc = 1$,
$c \ge 0$, and $d > 0$ if $c = 0$).
:::

::: function
void psl2z_print(const psl2z_t g)
void psl2z_fprint(FILE * file, const psl2z_t g)

Print `g` as `[a b; c d]` to standard output or to `file`.
:::

::: function
void psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g)
void psl2z_inv(psl2z_t h, const psl2z_t g)

Set `h` to the product $fg$, respectively the inverse $g^{-1}$, in
canonical form. Aliasing is allowed.
:::

::: function
void psl2z_randtest(psl2z_t g, flint_rand_t state, slong bits)

Sets `g` to a random element with entries of at most about `bits` bits
(at least 1): random $a$ and $d$ are divided by their GCD, completed with
Bézout coefficients, and the signs are canonicalised.
:::

## Modular transformations and the fundamental domain {#fundamental-domain}

The standard fundamental domain of $\mathrm{PSL}_2(\Z)$ acting on
$\mathbb{H}$ is $\{\tau : |\operatorname{Re}(\tau)| \le 1/2, |\tau| \ge 1\}$.
The functions below find a transformation that maps a point to the
domain, up to a tolerance $\varepsilon$; an approximate reduction is
enough for evaluation, since points slightly outside the domain still
have $|q|$ bounded well below 1. All evaluation functions of this module
use $1 - \varepsilon = 63/64$.

::: function
void acb_modular_transform(acb_t w, const psl2z_t g, const acb_t z, slong prec)

Sets `w` to $gz = (az + b)/(cz + d)$.
:::

::: function
void acb_modular_fundamental_domain_approx_d(psl2z_t g, double x, double y, double one_minus_eps)
void acb_modular_fundamental_domain_approx_arf(psl2z_t g, const arf_t x, const arf_t y, const arf_t one_minus_eps, slong prec)

Set `g` to a transformation that should map $x + yi$ (with $y > 0$) to
the fundamental domain or just outside it, by the classical reduction
algorithm in `double` or `arf` arithmetic at precision `prec`. The output
is always a valid element of $\mathrm{PSL}_2(\Z)$, but it is not
guaranteed to reduce the point; the caller must check the result.

::: details
The iteration alternately applies $z \gets z + b$ (with $b$ the nearest
integer to $-\operatorname{Re}(z)$) and $z \gets -1/z$, and stops when
$|x| \le 1/2$ and $x^2 + y^2 \ge 1 - \varepsilon$, where
$1 - \varepsilon$ is `one_minus_eps`. It also stops if too many steps
have been taken without convergence, or if the numbers become too large
or too small for the working precision.
:::
:::

::: function
void acb_modular_fundamental_domain_approx(acb_t w, psl2z_t g, const acb_t z, const arf_t one_minus_eps, slong prec)

Sets `g` to a transformation that should map the midpoint of `z` to the
fundamental domain or just outside it, and `w` to $gz$ evaluated in ball
arithmetic. The output `g` is always valid, but the caller must check
that `w` is close enough to the domain.

::: details
If $\operatorname{Im}(z) \ge 1$, only the translation is applied. Otherwise
`acb_modular_fundamental_domain_approx_d` is tried first (when
$2^{-40} < \operatorname{Im}(z)$ and $|\operatorname{Re}(z)| < 2^{40}$),
and if its result is not good enough,
`acb_modular_fundamental_domain_approx_arf` is used at precision `prec`.
If $z$ is not finite, not in the upper half-plane, or has real part
larger than $2^{\mathit{prec}}$, `g` is set to the identity.
:::
:::

::: function
int acb_modular_is_in_fundamental_domain(const acb_t z, const arf_t tol, slong prec)

Returns nonzero if it is certain that $|\operatorname{Re}(z)| \le 1/2 + \varepsilon$
and $|z| \ge 1 - \varepsilon$, where $\varepsilon$ is `tol`, and zero if
this is false or cannot be decided.
:::

## Jacobi theta functions {#theta}

The four Jacobi theta functions are

$$
\begin{aligned}
\theta_1(z,\tau) &= -i \sum_{n=-\infty}^{\infty} (-1)^n e^{\pi i [(n + 1/2)^2 \tau + (2n + 1) z]}
                 = 2 q_{1/4} \sum_{n=0}^{\infty} (-1)^n q^{n(n+1)} \sin((2n+1) \pi z), \\
\theta_2(z,\tau) &= \sum_{n=-\infty}^{\infty} e^{\pi i [(n + 1/2)^2 \tau + (2n + 1) z]}
                 = 2 q_{1/4} \sum_{n=0}^{\infty} q^{n(n+1)} \cos((2n+1) \pi z), \\
\theta_3(z,\tau) &= \sum_{n=-\infty}^{\infty} e^{\pi i [n^2 \tau + 2n z]}
                 = 1 + 2 \sum_{n=1}^{\infty} q^{n^2} \cos(2n \pi z), \\
\theta_4(z,\tau) &= \sum_{n=-\infty}^{\infty} (-1)^n e^{\pi i [n^2 \tau + 2n z]}
                 = 1 + 2 \sum_{n=1}^{\infty} (-1)^n q^{n^2} \cos(2n \pi z),
\end{aligned}
$$

where $q = e^{\pi i \tau}$ and $q_{1/4} = e^{\pi i \tau / 4}$. Many
authors write $q^{1/4}$ for $q_{1/4}$, but the principal fourth root
$(q)^{1/4}$ differs from $q_{1/4}$ in general, so this module writes
$q^k$ only for integer $k$. In the notation of [`acb_theta`](acb_theta),
$(\theta_{0,0}, \theta_{0,1}, \theta_{1,0}, \theta_{1,1}) = (\theta_3, \theta_4, \theta_2, -\theta_1)$.
[](acb_modular/theta) explains how the functions below are
evaluated and gives an example.

::: function
void acb_modular_theta(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, slong prec)
void acb_modular_theta_notransform(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, slong prec)

Set `theta1`, ..., `theta4` to $\theta_1(z, \tau), \ldots, \theta_4(z, \tau)$.
`acb_modular_theta` moves $\tau$ to the fundamental domain, reduces $z$
modulo the transformed lattice when its imaginary part is large, and
calls `acb_modular_theta_sum`. The `notransform` version sums the series
at $\tau$ directly, which is only efficient when $\operatorname{Im}(\tau)$
is not small.
:::

::: function
void acb_modular_theta_jet(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau, slong len, slong prec)
void acb_modular_theta_jet_notransform(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau, slong len, slong prec)

Set each output vector to the first `len` Taylor coefficients of
$\theta_i(z + x, \tau)$ at $x = 0$, so that entry $k$ is
$\theta_i^{(k)}(z, \tau)/k!$, with derivatives with respect to $z$. The
`notransform` version does not move $\tau$ to the fundamental domain or
reduce $z$.
:::

::: function
void _acb_modular_theta_series(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec)
void acb_modular_theta_series(acb_poly_t theta1, acb_poly_t theta2, acb_poly_t theta3, acb_poly_t theta4, const acb_poly_t z, const acb_t tau, slong len, slong prec)

Set the outputs to $\theta_i(z(x), \tau)$ truncated to length `len`, for
a power series $z(x)$, by composing the theta jets with the nonconstant
part of $z$. Any of the outputs may be `NULL`, in which case that
function is not stored.
:::

## Theta transformations and summation {#theta-sum}

These are the building blocks of the theta function evaluation. They are
documented for users who want to sum the series directly, for example to
avoid recomputing $q$.

::: function
void acb_modular_theta_transform(int * R, int * S, int * C, const psl2z_t g)

Computes the data of the transformation formula expressing theta
functions at $\tau$ through theta functions at $\tau' = g\tau$: integers
$R_i$, $S_i$ (arrays of length 4) and $C \in \{0, 1\}$ such that for
$i = 0, 1, 2, 3$

$$
\theta_{1+i}(z,\tau) = e^{\pi i R_i / 4} \cdot A \cdot B \cdot \theta_{1+S_i}(z',\tau'),
$$

where $z' = z$ and $A = B = 1$ if $C = 0$, and if $C = 1$

$$
z' = \frac{-z}{c \tau + d}, \quad
A = \sqrt{\frac{i}{c \tau + d}}, \quad
B = \exp\left(-\pi i c \frac{z^2}{c \tau + d}\right).
$$

::: details
$A$ is well defined with the principal square root since
$A^2 = i/(c\tau + d)$ lies in the right half-plane.

If $c = 0$, then $C = 0$: $\theta_i(z, \tau) = e^{-\pi i b / 4} \theta_i(z, \tau+b)$
for $i = 1, 2$, while $\theta_3$ and $\theta_4$ are unchanged when $b$ is
even and swap places when $b$ is odd.

For $c > 0$, $C = 1$ and the general transformations are given by
Rademacher [@Rad1973] (whose formulas contain typos and have the
transformed and untransformed variables exchanged). With
$\theta_{m,n}$ defined for $m, n \in \Z$ by

$$
\theta_{0,0} = \theta_3, \quad \theta_{0,1} = \theta_4, \quad
\theta_{1,0} = \theta_2, \quad \theta_{1,1} = i \theta_1, \quad
\theta_{m+2,n} = (-1)^n \theta_{m,n}, \quad
\theta_{m,n+2} = \theta_{m,n},
$$

we have

$$
\begin{aligned}
\theta_1(z,\tau) &= \varepsilon_1 A B \, \theta_1(z', \tau'), \\
\theta_2(z,\tau) &= \varepsilon_2 A B \, \theta_{1-c,1+a}(z', \tau'), \\
\theta_3(z,\tau) &= \varepsilon_3 A B \, \theta_{1+d-c,1-b+a}(z', \tau'), \\
\theta_4(z,\tau) &= \varepsilon_4 A B \, \theta_{1+d,1-b}(z', \tau'),
\end{aligned}
$$

where the $\varepsilon_i$ are 8th roots of unity. With
$\varepsilon(a,b,c,d) = \exp(\pi i R(a,b,c,d) / 12)$ the 24th root of
unity of the eta transformation (`acb_modular_epsilon_arg`),

$$
\begin{aligned}
\varepsilon_1(a,b,c,d) &= \exp(\pi i [R(-d,b,c,-a) + 1] / 4), \\
\varepsilon_2(a,b,c,d) &= \exp(\pi i [-R(a,b,c,d) + (5+(2-c)a)] / 4), \\
\varepsilon_3(a,b,c,d) &= \exp(\pi i [-R(a,b,c,d) + (4+(c-d-2)(b-a))] / 4), \\
\varepsilon_4(a,b,c,d) &= \exp(\pi i [-R(a,b,c,d) + (3-(2+d)b)] / 4).
\end{aligned}
$$
:::
:::

::: function
void acb_modular_theta_sum(acb_ptr theta1, acb_ptr theta2, acb_ptr theta3, acb_ptr theta4, const acb_t w, int w_is_unit, const acb_t q, slong len, slong prec)

Sets the outputs to the first `len` $\ge 1$ Taylor coefficients of
$\theta_1(z+x,\tau)/q_{1/4}$, $\theta_2(z+x,\tau)/q_{1/4}$,
$\theta_3(z+x,\tau)$ and $\theta_4(z+x,\tau)$, given $w = e^{\pi i z}$ and
$q = e^{\pi i \tau}$, by summing truncated theta series with a rigorous
bound for the tail. A nonzero `w_is_unit` asserts that $|w| = 1$ (real
$z$). Aliasing between inputs and outputs is not allowed.

::: details
The factor $q_{1/4}$ is left out of $\theta_1$ and $\theta_2$ because it
often cancels in quotients; the caller computes $q_{1/4} = e^{\pi i\tau/4}$
and $q = q_{1/4}^4$ and multiplies it back. The function is meant for
$|q| \ll 1$, after a modular transformation; it accepts any $q$ but
returns useless enclosures when the series does not converge quickly.

The four series are summed together. For $k = 0, 1, 2, \ldots$ the
exponents of $q$ are $\lfloor (k+2)^2/4 \rfloor = 1, 2, 4, 6, 9, \ldots$
(squares and pronic numbers interleaved), alternately for
$(\theta_3, \theta_4)$ and $(\theta_1, \theta_2)$, and the exponents of
$w$ are $\pm (k+2)$:

$$
\begin{array}{llll}
& \theta_1, \theta_2 & q^0 & (w^1 \pm w^{-1}) \\
k = 0  & \theta_3, \theta_4 & q^1 & (w^2 \pm w^{-2}) \\
k = 1  & \theta_1, \theta_2 & q^2 & (w^3 \pm w^{-3}) \\
k = 2  & \theta_3, \theta_4 & q^4 & (w^4 \pm w^{-4}) \\
k = 3  & \theta_1, \theta_2 & q^6 & (w^5 \pm w^{-5}) \\
k = 4  & \theta_3, \theta_4 & q^9 & (w^6 \pm w^{-6}) \\
\end{array}
$$

The powers of $q$ are computed with the addition sequence of
`acb_modular_addseq_theta`. Writing the even cosine terms as
$w^{2n} + w^{-2n}$, the odd cosine terms as $w (w^{2n} + w^{-2n-2})$ and
the sine terms as $w (w^{2n} - w^{-2n-2})$, only even powers of $w$ and
$w^{-1}$ are needed. The truncation bound is given in
[](algorithms/special-functions/theta-functions#genus-one).
:::
:::

::: function
void acb_modular_addseq_theta(slong * exponents, slong * aindex, slong * bindex, slong num)

Constructs an addition sequence for the first `num` squares and pronic
numbers interleaved (excluding zero), $1, 2, 4, 6, 9, 12, 16, 20, 25, 30, \ldots$:
entry $k$ of `exponents` is the $k$-th number, and it equals the sum of
the entries at `aindex[k]` and `bindex[k]`, or twice the first plus the
second.
:::

::: function
void acb_modular_theta_const_sum(acb_t theta2, acb_t theta3, acb_t theta4, const acb_t q, slong prec)
void acb_modular_theta_const_sum_basecase(acb_t theta2, acb_t theta3, acb_t theta4, const acb_t q, slong N, slong prec)
void acb_modular_theta_const_sum_rs(acb_t theta2, acb_t theta3, acb_t theta4, const acb_t q, slong N, slong prec)

Compute the theta constants $\theta_2(0, \tau)/q_{1/4}$,
$\theta_3(0, \tau)$ and $\theta_4(0, \tau)$ from $q = e^{\pi i \tau}$
without a modular transformation. The `basecase` and `rs` versions sum
the truncated series $\sum q^{k(k+1)}$, $\sum q^{k^2}$ and
$\sum (-1)^k q^{k^2}$ over the $k \in \Z$ with exponent less than $N$,
with a short addition sequence and with rectangular splitting
respectively [@EHJ2016], without a tail bound. `acb_modular_theta_const_sum` chooses $N$ from $|q|$ and `prec`,
uses the basecase for $N < 1800$ and adds the truncation error.
:::

## The Dedekind eta function {#eta}

::: function
void acb_modular_eta(acb_t z, const acb_t tau, slong prec)

Sets `z` to the Dedekind eta function
$\eta(\tau) = e^{\pi i \tau/12} \prod_{n \ge 1} (1 - e^{2\pi i n\tau})$,
by moving $\tau$ to the fundamental domain, summing the pentagonal number
series with `acb_modular_eta_sum`, and applying the transformation
formula.
:::

::: function
int acb_modular_epsilon_arg(const psl2z_t g)

Given $g = (a, b; c, d)$, returns an integer $R$ such that
$\varepsilon(a,b,c,d) = e^{\pi i R / 12}$ is the 24th root of unity in
the transformation formula of the eta function,

$$
\eta\left(\frac{a\tau+b}{c\tau+d}\right) = \varepsilon(a,b,c,d) \sqrt{c\tau+d} \, \eta(\tau).
$$
:::

::: function
void acb_modular_eta_sum(acb_t eta, const acb_t q, slong prec)

Sets `eta` to $\prod_{n \ge 1} (1 - q^n) = \sum_{k \in \Z} (-1)^k q^{k(3k-1)/2}$,
the eta function without the factor $e^{\pi i \tau / 12}$, given
$q = e^{2\pi i \tau}$, by summing the pentagonal number series with a
rigorous tail bound. Like `acb_modular_theta_sum`, it is meant for
$|q| \ll 1$. The series is summed with a short addition sequence, or with
rectangular splitting from 400 terms [@EHJ2016].
:::

::: function
void acb_modular_addseq_eta(slong * exponents, slong * aindex, slong * bindex, slong num)

Constructs an addition sequence for the first `num` generalised
pentagonal numbers (excluding zero), $1, 2, 5, 7, 12, 15, 22, 26, 35, 40, \ldots$,
in the format of `acb_modular_addseq_theta`.
:::

## Modular forms and functions {#forms}

The functions below move $\tau$ to the fundamental domain, compute the
theta constants $\theta_2, \theta_3, \theta_4$ at $z = 0$ with
`acb_modular_theta_const_sum`, and apply the transformation law of the
form. When $\operatorname{Re}(\tau)$ is a multiple of $1/2$ and the value
is real, the imaginary part of the output is set to exactly zero.

::: function
void acb_modular_j(acb_t z, const acb_t tau, slong prec)

Sets `z` to Klein's $j$-invariant, normalised so that $j(i) = 1728$,
computed on the fundamental domain as
$j(\tau) = 32 (\theta_2^8+\theta_3^8+\theta_4^8)^3 / (\theta_2 \theta_3 \theta_4)^8$
with $\theta_i = \theta_i(0, \tau)$. The output is indeterminate unless
$\operatorname{Im}(\tau)$ is certainly positive.
:::

::: function
void acb_modular_lambda(acb_t r, const acb_t tau, slong prec)

Sets `r` to the modular lambda function
$\lambda(\tau) = \theta_2^4(0,\tau) / \theta_3^4(0,\tau)$, which is
invariant under the subgroup $\Gamma(2)$ of transformations
$(a, b; c, d)$ with $a, d$ odd and $b, c$ even.
:::

::: function
void acb_modular_delta(acb_t z, const acb_t tau, slong prec)

Sets `z` to the modular discriminant $\Delta(\tau) = \eta(\tau)^{24}$,
computed as $(\theta_2 \theta_3 \theta_4 / 2)^8$, which transforms as
$\Delta(g\tau) = (c\tau+d)^{12} \Delta(\tau)$. The discriminant is
sometimes defined with an extra factor $(2\pi)^{12}$, which is omitted
here.
:::

::: function
void acb_modular_eisenstein(acb_ptr r, const acb_t tau, slong len, slong prec)

Sets the `len` entries of `r` to the Eisenstein series
$G_4(\tau), G_6(\tau), G_8(\tau), \ldots$, where

$$
G_{2k}(\tau) = \sum_{(m, n) \ne (0, 0)} \frac{1}{(m+n\tau )^{2k}},
\qquad
G_{2k}\left(\frac{a\tau+b}{c\tau+d}\right) = (c\tau+d)^{2k} G_{2k}(\tau).
$$

$G_4$ and $G_6$ are computed from theta constants on the fundamental
domain, the others by the convolution recurrence of
[@Apostol1997, Theorem 1.13].

::: details
With $a = \theta_2^4$, $b = \theta_3^4$, $c = \theta_4^4$ (theta
constants at the reduced $\tau$),

$$
G_4 = \frac{\pi^4}{90} (a^2 + b^2 + c^2), \qquad
G_6 = \frac{\pi^6}{945} (b^3 + c^3 - 3a^2(b + c)) ,
$$

and with $c_k = (2k-1) G_{2k}$, the recurrence is

$$
(2k+1)(k-3) c_k = 3 \sum_{j=2}^{k-2} c_j c_{k-j}, \qquad k \ge 4.
$$
:::
:::

## Hilbert class polynomials {#class-polynomials}

::: function
void acb_modular_hilbert_class_poly(fmpz_poly_t res, slong D)

Sets `res` to the Hilbert class polynomial $H_D$ of the negative
discriminant $D$, the minimal polynomial over $\mathbb{Q}$ of $j$ at a CM point
of discriminant $D$. If $D \ge 0$ or $D \equiv 2, 3 \pmod 4$, `res` is
set to the zero polynomial. See [](acb_modular/class-polynomials).
:::

## Elliptic integrals and functions {#elliptic}

These wrappers are kept for backwards compatibility; use the functions of
[`acb_elliptic`](acb_elliptic) instead.

::: function
void acb_modular_elliptic_k(acb_t k, const acb_t m, slong prec)
void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, slong len, slong prec)
void acb_modular_elliptic_e(acb_t res, const acb_t m, slong prec)
void acb_modular_elliptic_p(acb_t r, const acb_t z, const acb_t tau, slong prec)
void acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, slong len, slong prec)

Call `acb_elliptic_k`, `acb_elliptic_k_jet`, `acb_elliptic_e`,
`acb_elliptic_p` and `acb_elliptic_p_jet` respectively.
:::

## Internal functions {#internal}

::: function
void acb_modular_fill_addseq(slong * tab, slong len)

Builds an addition sequence for a set of integers that is assumed to be
fairly dense. On input, entry $i$ of `tab` is $-1$ if $i$ is to be part
of the sequence and 0 otherwise; on output, a marked entry $i$ holds a
value $j$ such that $j$ and $i - j$ are also in the sequence (doubling is
preferred). Entries 0 and 1 are ignored; 1 is always assumed to be in the
sequence.
:::

::: function
void _acb_modular_mul(acb_t z, acb_t tmp1, acb_t tmp2, const acb_t x, const acb_t y, slong wprec, slong prec)

Sets `z` to $xy$ at precision `wprec`, first rounding the inputs to
`wprec` bits (using the temporaries `tmp1`, `tmp2`) when `prec` exceeds
1024. Used for the powers of $q$ in the series summation, where terms of
small magnitude need less precision.
:::
