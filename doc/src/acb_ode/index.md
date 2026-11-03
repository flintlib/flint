---
title: acb_ode.h
header: acb_ode.h
description: Linear differential equations with polynomial coefficients.
category: Real and complex numbers
---

::: warning
This module is experimental; the API is subject to change.
:::

This module provides functions for solving linear homogeneous differential
equations

$$
a_r(x) y^{(r)}(x) + \dots + a_1(x) y'(x) + a_0(x) y(x) = 0
$$

with polynomial coefficients $a_i(x) \in \mathbb C[x]$.
When working locally at 0, it is often convenient to rewrite the equation as

$$
\tilde a_r(x) (\theta^r y)(x) + \dots + \tilde a_1(x) (\theta y)(x) + \tilde a_0(x) y(x) = 0
$$

where $\theta = x \frac{\mathrm d}{\mathrm d x}$ is the *Euler derivation*.
The integer $r$ is called the *order* of the equation.

Some high-level functions take as input differential operators represented as
`gr_ore_poly_t` objects and may accept both of the above forms, as well as
various coefficient types.
Lower-level functions typically expect pointers to vectors of `acb_poly_struct`
objects instead. Unless stated otherwise, these correspond to coefficient
vectors with respect to the Euler derivation.

This module currently supports evaluating solutions given by initial conditions
at an *ordinary* or a *regular singular* point.
Initial conditions at a point are coordinates in a given *local basis* of the
solution space.
At an ordinary point, the equation admits an $r$-dimensional vector space of
analytic solutions, all of valuation less than $r$.
At a regular singular point, it admits an $r$-dimensional space of solutions
of the form

$$
\sum_\lambda x^\lambda \sum_{k=0}^{K-1} y_{\lambda,k}(x) \frac{\log(x)^k}{k!}
$$

where the $y_{\lambda,k}$ are analytic and $\lambda$ ranges over the local
*exponents* of the differential equation, which are the roots of the local
*indicial polynomial*.
A maximal set of exponents that differ by integers is called a *group*
(this has nothing to do with the algebraic structure of the same name).
We call the exponent of minimal real part in a group the *leader* of the group,
and the difference of an exponent with its group leader its *shift*.
Solutions that can be written with a single common $\lambda$ are also said to
belong to the same group.

There is currently no support for inhomogeneous equations.

## Evaluation of fundamental matrices {#evaluation-of-fundamental-matrices}

::: type
acb_ode_basis_t

Constants used to identify families of local bases of solutions.
To specify a basis, one needs a constant of this type, an expansion
point, and possibly some additional information.

Currently the following two families are supported.
In both cases, the basis is the concatenation of blocks corresponding to
groups of local exponents.
Within each block, the columns are ordered by increasing shift relative to
the group leader (that is, by increasing real part of the exponent), and,
for a given shift, by decreasing degree in $\log(x)$ of the coefficient
of the leading power of $x$.
The order of the blocks is specified separately.

Let $E$ denote the set of pairs $(\alpha, k)$ where $\alpha$ is a local
exponent of multiplicity $m$ and $k$ ranges from $0$ to $m - 1$.
At a regular singular point, the cardinality of $E$ is equal to the order of
the differential equation.

- `ACB_ODE_BASIS_ECHELON`: The basis is made of a solution of leading term
  $x^{\alpha} \log(x)^k/k!$ for each $(\alpha, k) \in E$.
  The series expansion of the basis element attached to
  $(\alpha, k)$ involves no term in $x^{\beta} \log(x)^\ell$
  with $(\beta, \ell) \in E \setminus \{ (\alpha, k) \}$.

  This is the default when no basis choice can be specified.

- `ACB_ODE_BASIS_CASCADE`: The basis contains a solution of leading term
  $x^{\alpha} \log(x)^{m-1}/(m-1)!$
  for each exponent $\alpha$, with $m$ the multiplicity of $\alpha$.
  The series expansions of these basis elements involve no term in
  $x^{\beta} \log(x)^\ell$
  with $(\beta, \ell) \in E \setminus \{ (\alpha, m-1) \}$.
  The remaining elements of the basis are obtained by repeatedly changing
  each $\log(x)^k/k!$ to $\log(x)^{k-1}/(k-1)!$ and dropping the terms not
  involving $\log(x)$ in previously constructed elements.

More basis families may be added in the future, including some not sharing the
block structure of the previous two.
:::

::: function
int acb_ode_fundamental_matrix(acb_mat_t mat, const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop, const acb_ode_exponents_t expos, acb_srcptr lcrt, const acb_t pt, acb_ode_basis_t basis, ulong flags, slong prec)

Set `mat` to the value at `pt` of the fundamental matrix in the
basis `basis` of the differential operator `dop`.

This function is currently limited to operators `dop` of algebra type
`ORE_ALGEBRA_DERIVATIVE` or
`ORE_ALGEBRA_EULER_DERIVATIVE`.

The origin may be an ordinary point or a regular singular point of `dop`.

The evaluation point `pt` must lie within an open disk centered at the
origin and containing no singular point of `dop` except, possibly, for the
origin itself.

The `expos` parameter corresponds to the exponent structure of `dop` at the
origin in the format returned by `acb_ode_exponents`.
In general, it cannot be automatically determined from `dop` and must be
provided by the caller (see `acb_ode_exponents_t`).
When the exponent structure can be computed by `acb_ode_exponents`
(typically, over an exact field of constants with a complete zero test),
though, `expos` may be set to `NULL`.
Note however that, in this case, the order of the exponent groups is
unspecified.

Similarly, the user should provide in `lcrt` enclosures of the nonzero roots
of the leading coefficient of `dop`, with multiple roots repeated. However,
when the ground field supports computing complex ball roots using
`gr_poly_roots_other`, the `lcrt` argument may be set to `NULL`.

Each of the first $r$ columns of the output matrix `mat`, where $r$ is the
order of `dop`, is set to a vector of initial coefficients of the Taylor
expansion at `pt` of the corresponding element of the basis of type `basis`
of the solution space, with solution groups ordered as in `expos`.
(Extraneous columns are left untouched. Note that when `expos` is omitted
on input, the order of the columns is not fully specified.)
The row dimension of `mat` determines the number of Taylor coefficients to
be evaluated.
The current implementation is not designed to evaluate a number of
Taylor coefficients exceeding the order of `dop` efficiently.

The following `flags` are supported:

- `ACB_ODE_ORDINARY`: Promises that the origin is an ordinary point.

The return value is a `gr`-style status code. It may happen, however, that
the return value is `GR_SUCCESS` but the output matrix is indeterminate.

::: details
The example below illustrates the basic usage of this function with an exact
operator as input. See also the documentation of `acb_ode_exponents_t` for
examples of how to construct a suitable exponents object in cases where it
cannot automatically be derived from the operator.
```c example name=acb_ode_fundamental_matrix
#include "flint/acb_types.h"
#include "flint/acb.h"
#include "flint/acb_mat.h"
#include "flint/acb_ode.h"
#include "flint/gr.h"
#include "flint/gr_ore_poly.h"

int main(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ptr dop;

    int status = GR_SUCCESS;

    slong prec = 50;

    /* Differential operator */

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ctx_init_gr_ore_poly(Dop, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
    status |= gr_ctx_set_gen_name(Pol, "z");
    status |= gr_ctx_set_gen_name(Dop, "Tz");

    GR_TMP_INIT(dop, Dop);

    status |= gr_ore_poly_set_str(dop,  "(z^2 + 1)*Tz^2 + (z^2 - 1)*Tz", Dop);

    flint_printf("dop = %{gr}\n", dop, Dop);

    /* Evaluation point */

    acb_t pt;
    acb_init(pt);
    acb_const_pi(pt, prec + 4);
    arb_one(acb_imagref(pt));
    acb_inv(pt, pt, prec + 4);

    /* Output matrix */

    acb_mat_t mat;
    slong dop_order = gr_ore_poly_length(dop, Dop) - 1;
    acb_mat_init(mat, dop_order, dop_order);

    /* Numerical solution */

    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt, 0, 0, prec);

    acb_mat_printd(mat, prec * 0.3);

    acb_clear(pt);
    acb_mat_clear(mat);
    GR_TMP_CLEAR(dop, Dop);

    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(ZZ);

    flint_cleanup_master();

    return status;
}
```
```output
dop = (-1 + z^2)*Tz + (1 + z^2)*Tz^2
[(1.00000000000000 + 0j)  +/-  (0, 0j), (0.283455752470505 - 0.0850599850774542j)  +/-  (5.16e-16, 1.37e-16j)]
[(0 + 0j)  +/-  (0, 0j), (0.927899890214499 + 0.0459003554793228j)  +/-  (1.46e-15, 4.89e-16j)]
```
:::
:::

::: function
void _acb_ode_fundamental_matrix_vec(acb_mat_struct * mat, const acb_poly_struct * dop, slong dop_len, const acb_ode_exponents_t expos, acb_srcptr lcrt, acb_ode_bound_t bound, acb_srcptr pts, slong npts, acb_ode_basis_t basis, acb_ode_sum_worker_t sum_worker, void * worker_data, ulong flags, slong prec)
int acb_ode_fundamental_matrix_vec(acb_mat_struct * mat, const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop, const acb_ode_exponents_t expos, acb_srcptr lcrt, acb_srcptr pts, slong npts, acb_ode_basis_t basis, ulong flags, slong prec)

Like `acb_ode_fundamental_matrix` but evaluate the fundamental matrix at
multiple points.
The output parameter `mat` should point to an array of `npts` matrices,
not necessarily with the same number of rows.

The underscore method assumes that `dop` has order at least one and does
not allow the `expos` and `lcrt` parameters to be omitted.
Additionally, the caller must provide an object `bound` pre-filled using
`acb_ode_bound_precompute` and a summation callback `sum_worker`
with an arbitrary user data pointer `worker_data`.
The user data is accessible to `sum_worker` in the `data` field of its
`sum` argument.
:::

## Series solutions {#series-solutions}

::: function
int acb_ode_solution_series(gr_poly_vec_t series, const gr_ore_poly_t dop, gr_ore_poly_ctx_t Dop, const acb_ode_group_t group, const gr_vec_t ini, slong n, ulong flags, gr_ctx_t CC)

Write to `series` a vector of polynomials such that
$x^\lambda \sum_k \mathit{series}[k] \log(x)^k/k!$
is equal to the series expansion, truncated to order `n`, of the solution
of `dop` belonging to the group `group` of leader $\lambda$ and whose
coordinates in the block corresponding to that group of the
`ACB_ODE_BASIS_ECHELON` basis are given by `ini`.
Here `ini` should be a vector of elements of `CC` and `series` should be a
vector of polynomials over `CC`.
The supported values for `CC` are currently limited to complex ball fields.
The `flags` parameter may be set to `ACB_ODE_ORDINARY` to indicate that the
origin is known to be an ordinary point.
:::

## Solution groups and exponents {#solution-groups-and-exponents}

::: type
acb_ode_shift_struct
acb_ode_shift_t
acb_ode_group_struct
acb_ode_group_t
acb_ode_exponents_struct
acb_ode_exponents_t

Types used to represent the structure of a local basis of solutions, that
is, the roots of the indicial polynomial and the integer differences between
these roots.

An `acb_ode_exponents_struct` contains an array `groups` of length `ngroups`
of `acb_ode_group_struct` objects.
Each `acb_ode_group_struct` contains an exponent `leader` of type
`acb_t` and an array of length `nshifts` of `acb_ode_shift_struct`
objects.
An `acb_ode_shift_struct` is a pair of an integer shift `n` and a
multiplicity `mult`.
The shifts are sorted in increasing order.

::: details
```c example name=acb_ode_exponents
#include "flint/acb_types.h"
#include "flint/acb.h"
#include "flint/acb_mat.h"
#include "flint/acb_ode.h"
#include "flint/gr.h"
#include "flint/gr_ore_poly.h"

int
fundamental_matrix(const char * dop_str,
                   const acb_ode_exponents_struct * expos,
                   double pt_d)
{
    gr_ctx_t CC, Pol, Dop;
    gr_ptr dop;

    slong prec = 30;

    int status = GR_SUCCESS;

    gr_ctx_init_complex_acb(CC, prec);
    gr_ctx_init_gr_poly(Pol, CC);
    gr_ctx_init_gr_ore_poly(Dop, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);

    GR_TMP_INIT(dop, Dop);

    status |= gr_ctx_set_gen_name(Pol, "z");
    status |= gr_ctx_set_gen_name(Dop, "Tz");
    status |= gr_ore_poly_set_str(dop, dop_str, Dop);

    flint_printf("%{gr}\n", dop, Dop);

    slong dop_order = gr_ore_poly_length(dop, Dop) - 1;

    acb_mat_t mat;
    acb_mat_init(mat, dop_order, dop_order);

    acb_t pt;
    acb_init(pt);
    acb_set_d(pt, pt_d);

    GR_MUST_SUCCEED(acb_ode_fundamental_matrix(mat, dop, Dop, expos, NULL, pt, 0, 0, prec));

    flint_printf("@(x=%{acb}) =\n%{acb_mat}\n\n", pt, mat);

    acb_mat_clear(mat);
    GR_TMP_CLEAR(dop, Dop);
    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);
    gr_ctx_clear(CC);
    acb_clear(pt);

    return status;
}

void
apery(void)
{
    flint_printf("Apery:\n");
    acb_ode_shift_struct shift[1] = {{ .n = 0, .mult = 3 }};
    acb_ode_group_struct grp[1] = {{ .nshifts = 1, .shifts = shift }};
    acb_init(grp->leader);
    acb_zero(grp->leader);
    acb_ode_exponents_struct expos[1] = {{ .ngroups = 1, .groups = grp }};

    fundamental_matrix(
            "(z^2 - 34*z + 1)*Tz^3 + (3*z^2 - 51*z)*Tz^2 + (3*z^2 - 27*z)*Tz + z^2 - 5*z",
            expos,
            0.015625);

    acb_clear(grp->leader);
}

void
multiple_shifts(void)
{
    flint_printf("Multiple shifts:\n");
    const char * dop = "Tz^6 - 6*Tz^5 + 12*Tz^4 - 10*Tz^3 + 3*Tz^2 + z^2";

    acb_ode_shift_struct shift[3] = {
        { .n = 0, .mult = 2 },
        { .n = 1, .mult = 3 },
        { .n = 3, .mult = 1 },
    };
    acb_ode_group_struct grp[1] = {{ .nshifts = 3, .shifts = shift }};
    acb_init(grp->leader);
    acb_zero(grp->leader);
    acb_ode_exponents_struct expos[1] = {{ .ngroups = 1, .groups = grp }};

    fundamental_matrix(dop, expos, 2.);

    acb_clear(grp->leader);
}

void
whittaker(void)
{
    flint_printf("Whittaker:\n");
    slong prec = 30;

    acb_t kappa, mu, half;
    acb_init(kappa);
    acb_init(mu);
    acb_init(half);

    const char * dop = "4*Tz^2 - 4*Tz - z^2 + 8*z - 11";
    acb_set_si(kappa, 2);
    acb_set_si(mu, 3);
    acb_sqrt(mu, mu, prec);
    acb_set_d(half, .5);

    acb_ode_shift_struct shift[1] = { { .n = 0, .mult = 1 } };
    acb_ode_group_struct grp[2] = {
        { .nshifts = 1, .shifts = shift },
        { .nshifts = 1, .shifts = shift },
    };
    acb_init(grp[0].leader);
    acb_sub(grp[0].leader, half, mu, prec);
    acb_init(grp[1].leader);
    acb_add(grp[1].leader, half, mu, prec);
    acb_ode_exponents_struct expos[1] = {{ .ngroups = 2, .groups = grp }};

    fundamental_matrix(dop, expos, 2.);

    acb_clear(grp[0].leader);
    acb_clear(grp[1].leader);
    acb_clear(kappa);
    acb_clear(mu);
    acb_clear(half);
}

int
main(void)
{
    apery();
    multiple_shifts();
    whittaker();

    flint_cleanup_master();

    return 0;
}
```
```output
Apery:
(-5.000000000*z + z^2) + (-27.00000000*z + 3.000000000*z^2)*Tz + (-51.00000000*z + 3.000000000*z^2)*Tz^2 + (1 - 34.00000000*z + z^2)*Tz^3
@(x=0.0156250) =
[[[8.47877 +/- 4.00e-7] + [+/- 4.87e-18] * i, [-4.32824 +/- 2.79e-6] + [+/- 1.82e-18] * i, [1.10474 +/- 7.86e-7] + [+/- 3.53e-19] * i], [[-294.557 +/- 2.94e-4] + [+/- 1.98e-14] * i, [56.9120 +/- 2.98e-5] + [+/- 7.35e-15] * i, [9.32690 +/- 1.03e-6] + [+/- 1.42e-15] * i], [[9538.55 +/- 2.53e-3] + [+/- 3.94e-11] * i, [-1947.59 +/- 2.18e-3] + [+/- 1.47e-11] * i, [266.409 +/- 1.12e-4] + [+/- 2.82e-12] * i]]

Multiple shifts:
(z^2) + 3.000000000*Tz^2 - 10.00000000*Tz^3 + 12.00000000*Tz^4 - 6.000000000*Tz^5 + Tz^6
@(x=2) =
[[[-1.56921 +/- 4.41e-6] + [+/- 2.73e-20] * i, [1.99074 +/- 3.49e-6] + [+/- 1.62e-20] * i, [0.313250 +/- 1.66e-7] + [+/- 1.30e-22] * i, [1.52686 +/- 1.35e-6] + [+/- 2.48e-18] * i, [1.92285 +/- 7.79e-7] + [+/- 1.89e-18] * i, [7.99000 +/- 9.49e-7] + [+/- 6.28e-17] * i], [[-1.22252 +/- 4.74e-6] + [+/- 1.72e-19] * i, [0.981490 +/- 2.90e-7] + [+/- 9.68e-20] * i, [0.594985 +/- 2.24e-7] + [+/- 9.03e-22] * i, [1.98603 +/- 3.79e-6] + [+/- 1.36e-17] * i, [0.828657 +/- 2.06e-7] + [+/- 1.04e-17] * i, [11.9750 +/- 3.32e-6] + [+/- 3.46e-16] * i], [[-0.268112 +/- 8.15e-8] + [+/- 5.01e-19] * i, [0.236121 +/- 3.99e-7] + [+/- 2.67e-19] * i, [0.208879 +/- 4.84e-8] + [+/- 2.95e-21] * i, [0.444069 +/- 1.45e-7] + [+/- 3.40e-17] * i, [-0.127383 +/- 4.17e-7] + [+/- 2.60e-17] * i, [5.97500 +/- 4.97e-6] + [+/- 8.64e-16] * i], [[0.100534 +/- 4.28e-7] + [+/- 8.97e-19] * i, [-0.00462277 +/- 1.63e-9] + [+/- 4.44e-19] * i, [-0.0717298 +/- 4.27e-8] + [+/- 5.95e-21] * i, [-0.000419269 +/- 5.43e-10] + [+/- 5.10e-17] * i, [-0.0351200 +/- 4.70e-8] + [+/- 3.90e-17] * i, [0.987504 +/- 1.34e-7] + [+/- 1.30e-15] * i], [[-0.0192704 +/- 2.57e-9] + [+/- 1.11e-18] * i, [-0.000576132 +/- 1.26e-10] + [+/- 4.99e-19] * i, [0.00227053 +/- 3.62e-9] + [+/- 8.34e-21] * i, [0.00976466 +/- 4.08e-9] + [+/- 5.10e-17] * i, [-0.00172191 +/- 4.47e-9] + [+/- 3.90e-17] * i, [-0.00312293 +/- 3.33e-9] + [+/- 1.30e-15] * i], [[0.00723142 +/- 1.07e-9] + [+/- 9.94e-19] * i, [5.14228e-7 +/- 3.32e-13] + [+/- 4.00e-19] * i, [0.000366028 +/- 6.34e-10] + [+/- 8.64e-21] * i, [-0.00324687 +/- 4.53e-9] + [+/- 3.57e-17] * i, [0.000179376 +/- 1.61e-10] + [+/- 2.73e-17] * i, [-0.000311880 +/- 2.05e-11] + [+/- 9.07e-16] * i]]

Whittaker:
(-11.00000000 + 8.000000000*z - z^2) - 4.000000000*Tz + 4.000000000*Tz^2
@(x=2) =
[[[1.99998 +/- 9.67e-7] + [+/- 2.31e-14] * i, [1.96034 +/- 4.73e-6] + [+/- 5.47e-16] * i], [[-0.382362 +/- 5.52e-7] + [+/- 2.49e-13] * i, [1.35729 +/- 4.85e-6] + [+/- 5.61e-15] * i]]
```
:::
:::

::: function
void acb_ode_exponents_init(acb_ode_exponents_t expos)
void acb_ode_exponents_clear(acb_ode_exponents_t expos)

Initialize, respectively clear an `acb_ode_exponents_t` object.
:::

::: function
int acb_ode_exponents(acb_ode_exponents_t expos, const gr_ore_poly_t dop, gr_ctx_t Dop, gr_ctx_t CC)

Set `expos` to the exponent structure of the operator `dop` at the origin,
computed over `CC`. The only currently supported `CC` are complex ball fields.
:::

::: function
void acb_ode_exponents_ordinary(acb_ode_exponents_t expos, slong order)

Set `expos` to the exponent structure of an operator of order `order` at
an ordinary point.
:::

::: function
void acb_ode_exponents_randtest(acb_ode_exponents_t expos, flint_rand_t state, slong len, slong disp, slong prec, slong mag_bits)

Set `expos` to a random exponent structure with `len` exponents in total
(counting multiplicities) and a maximum dispersion of `disp`.
:::

::: function
slong acb_ode_exponents_length(const acb_ode_exponents_t expos)

Return the total number of solutions across all groups in `expos`, counted
with multiplicities.
:::

::: function
void acb_ode_exponents_println(const acb_ode_exponents_t expos)

Print an exponent structure.
:::

::: function
void acb_ode_group_init(acb_ode_group_t group, slong len)

Initialize `group` with space for up to `len` shifts.
:::

::: function
void acb_ode_group_clear(acb_ode_group_t group)

Clear `group`.
:::

::: function
void acb_ode_group_set(acb_ode_group_t dest, const acb_ode_group_t src)

Set `dest` (which must have space for `src->nshifts` shifts) to a copy of
`src`.
:::

::: function
slong acb_ode_group_length(const acb_ode_group_t group)

Return the number of solutions of `group`, i.e., the sum of
multiplicities of its shifts.
:::

::: function
slong acb_ode_group_multiplicity (const acb_ode_group_t group, slong n)

Return the multiplicity of `n` as a shift in group `group`.
:::

::: function
slong acb_ode_group_nlogs (const acb_ode_group_t group, slong n)

Return a bound on the length of the term of index `n` in any solution of
the group `group` when that term is viewed as a polynomial in $\log(x)$.
:::

::: function
void acb_ode_indicial_polynomial_from_exponents(acb_poly_t ind, const acb_ode_exponents_t expos, slong prec)

Set `ind` to the indicial polynomial of an equation with exponents `expos`.
:::

## Partial sums of logarithmic series solutions {#partial-sums-of-logarithmic-series-solutions}

This section describes low-level functions for computing series expansions
and/or values of partial sums or derivatives of partial sums of solutions.
These functions can handle several solutions from the same group (and possibly,
in the future, from conjugate groups) in parallel, and can evaluate them at
zero, one, or more points.

The general usage pattern is as follows:

1. initialize a new context object of type `acb_ode_sum_t`;
2. configure it by setting a differential operator
   (`acb_ode_sum_set_diffop`), solution group (`acb_ode_sum_set_group` or
   similar) and some vectors of initial values (`acb_ode_sum_set_ini_*`), in
   this order;
3. optionally set some evaluation points and other options;
4. call `acb_ode_sum_precompute`;
5. call one or more summation functions such as `acb_ode_sum_divconquer` to
   compute partial sums of solutions;
6. extract the value and/or coefficients of the solution.

Steps 3 to 6 can be repeated with the same `acb_ode_sum_t` object, but changing
the parameters set at step 2 without first clearing the context object is not
supported in general.

::: type
acb_ode_sum_struct
acb_ode_sum_t

An `acb_ode_sum_t` object `sum` stores, among other things:

- the terms in $x^{\lambda+i}$ where $\lambda$ is the group leader in the
  solutions being computed, for $i$ between `sum->n0` and `sum->n - 1`,
- some additional coefficients, associated with indices `i >= sum->n` and
  collectively referred to as the *residual buffer*.

The following options can be set by writing to the corresponding fields
after initializing the object:

- `ulong flags`:

  - `ACB_ODE_WANT_SERIES`: Keep all series coefficients.

  - `ACB_ODE_APPROX`: Skip the computation of rigorous error bounds.

- `slong wp`: Working precision for the computation of series coefficients and
  other intermediate results subject to overestimation due to dependencies.

- `slong sum_wp`: Working precision for the accumulation of partial sums.
:::

::: type
acb_ode_sol_struct
acb_ode_sol_t

Represents an individual logarithmic series solution in the process of being
computed.
Mainly makes sense as a component of an `acb_ode_sum_t` object or
given such an object as context.
:::

::: function
void acb_ode_sum_init(acb_ode_sum_t sum, slong dop_len, slong npts, slong nsols, slong nder)

Initialize a context object suitable for computing the derivatives of
order 0 to `nder - 1` of `nsols` solutions of an equation of order
`dop_len - 1` at `npts` points.
When only the coefficients of the series are of interest, `npts` can be set
to zero.
:::

::: function
void acb_ode_sum_clear(acb_ode_sum_t sum)

Clear the context object.
:::

::: function
void acb_ode_sum_zero(acb_ode_sum_t sum)

Reset the summation. Defining data such as the differential operator is
kept, but all partial sums are set to zero.
:::

::: function
void acb_ode_sum_fit_length(acb_ode_sum_t sum, slong len)

Ensure that `sum` has space for `len` terms, counting both partial sum
terms and entries of the residual buffer.
:::

::: function
void _acb_ode_sum_discard_head(acb_ode_sum_struct * sum, slong n0)

Discard the coefficients of index less than `n0` of the truncated series
stored in `sum`, increasing `sum->n0` accordingly.
:::

::: function
void acb_ode_sum_set_diffop(acb_ode_sum_t sum, const acb_poly_struct * dop, slong dop_len, const mag_t cvrad)

Set the differential operator whose solutions are to be computed.
The operator must be nonzero.
The `cvrad` parameter is a lower bound on the radius of convergence of the
solutions and will limit the magnitude of the evaluation points.
:::

::: function
void acb_ode_sum_set_group(acb_ode_sum_t sum, const acb_ode_group_t grp)

Set the solution group.
:::

::: function
void acb_ode_sum_set_ini_echelon(acb_ode_sum_t sum)
void acb_ode_sum_set_ini_highest(acb_ode_sum_t sum)

Set initial values defining the solutions to be computed.
The solution group must have been set before.
The `echelon` version sets initial values corresponding to the elements of
the echelonized basis of solutions (as defined above) belonging to the
`sum` object's group.
The `highest` version sets initial values for a single one of these basis
elements per shift, namely the one with the highest power of $\log(x)$ in
the leading term.
:::

::: function
void acb_ode_sum_set_points(acb_ode_sum_t sum, acb_srcptr pts, slong npts)

Set the evaluation points.
:::

::: function
slong acb_ode_sum_max_nlogs(acb_ode_sum_struct * sum)
slong acb_ode_sum_max_nlogs_xn(acb_ode_sum_struct * sum, slong n)

Return the maximum length of any of the truncated solutions stored in `sum`
viewed as a polynomial in $\log(x)$, respectively of the terms in $x^n$ of any
of those series.
:::

::: function
void acb_ode_sum_precompute(acb_ode_sum_t sum)

Perform various precomputations.
:::

::: type
acb_ode_sum_worker_t

Pointer to a worker function
`void worker(acb_ode_sum_t sum, slong nterms, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)`.
:::

::: function
void acb_ode_sum_divconquer(acb_ode_sum_t sum, slong nterms, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

Sum the series expansions of the solutions contained in `sum`, updating the
`sum` object with their coefficients and/or partial sums.
The `nterms` parameter can be used to limit the number of terms (-1 means
unlimited).
Unless `nterms` is set to a nonnegative value, the caller must supply two
objects `bound` and `gbound` precomputed using `acb_ode_bound_precompute` and
`acb_ode_group_bound_precompute`.

This function restarts the summation from $n=0$ (possibly reusing some internal
data) even if `sum` already contains a partial sum from a previous run.
:::

All `_acb_ode_sum_*` functions below assume that
`acb_ode_sum_precompute` has been called before.
The numeric precision of the computations they perform is controlled by two
parameters, `sum->wp` and `sum->sum_wp` (see `acb_ode_sum_t`), which
generally must be set by the caller but may be updated by the function.

The `_acb_ode_sum_forward_*` functions expect an object *sum* in which a
certain number of initial entries of the residual buffer are equal to the
corresponding coefficients of the image by the differential operator of the
solutions truncated at order `sum->n` (the current partial sums).
They extend the partial sums by a possibly different number of new terms,
destroying the corresponding entries of the residual buffer (if any).
They do not otherwise update the residual buffer to match the new partial sum.
Additionally, these functions may leave *sum* in a state that requires
`_acb_ode_sum_fix` to be called before passing *sum* to a function other
than `_acb_ode_sum_*`.

::: function
void _acb_ode_sum_forward_1(acb_ode_sum_struct * sum)

Add one term to `sum`, consuming the entry of index `sum->n` of the
residual buffer.
:::

::: function
int _acb_ode_sum_forward_divconquer(acb_ode_sum_struct * sum, slong high, slong block_len, slong stride, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

Extend the sum until the partial sum of order `high` is reached, some error
tolerance based on `prec` is met, or the working precision and/or tail bound
quality are determined to be insufficient, whichever happens first.

The first min(`block_len`, `high` - `sum->n`) entries of the residual buffer
must be set to the corresponding residual coefficients and are consumed.
Unless the `ACB_ODE_WANT_SERIES` flag is set, series coefficients
before the last `block_len` computed are discarded.

The `bound` and `gbound` parameters may both be `NULL`, in which case the
function stops exactly at term `high` and returns zero.
Alternatively, they may be set to objects precomputed using respectively
`acb_ode_bound_precompute` and `acb_ode_group_bound_precompute`.
In this case, the function returns nonzero if it detects either that the sum
has converged (within a tolerance based on `prec`) or that the current
working precision is insufficient to significantly improve the current
estimate.

The value of `block_len` must be at least the maximum degree of the
coefficients of the differential operator written in terms of the Euler
derivation. The algorithm computes strides of `stride` terms (checking for
convergence and making some internal adjustments after each stride) by
chaining blocks of `block_len` terms that are computed using a
divide-and-conquer method.
:::

::: function
void _acb_ode_sum_update_residuals(acb_ode_sum_struct * sum, slong low, slong mid, slong high)

Add to the entries of index `mid` to `high - 1` of the residual buffer
the coefficients of $x^i$ (`mid <= i < high`) in the image by the
differential operator of the terms in $x^i$ (`low <= i < mid`) of the
current partial sum.
:::

::: function
void _acb_ode_sum_fix(acb_ode_sum_struct * sum)

Restore some invariants of `sum` broken by the `sum_forward_*` functions.
After a sequence of calls to those functions, `_acb_ode_sum_fix` must be
called before accessing the computed sum.
:::

::: function
double _acb_ode_sum_update_tail_bounds(acb_ode_sum_struct * sum, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

Update `sum` with rigorous upper bounds on the magnitude of the tails of
the solutions.

Return a non-rigorous estimate of the sum over all solutions being computed
of the gain of accuracy in bits with respect to the bounds previously stored
in `sum`.
It can happen that some of the tail bounds worsen when taking into account
more terms of the series while the bounds previously computed using fewer
terms are known to remain valid.
In this case, the corresponding series contribute negatively to the return
value of the function but the existing tighter bounds are not overwritten.

The entries of index `sum->n` to `sum->n + d - 1` where `d`
is the maximum degree in `x` of the coefficients of the differential
operator must be preloaded with the corresponding residual coefficients.
The caller must supply two objects `bound` and `gbound` precomputed using
`acb_ode_bound_precompute` and `acb_ode_group_bound_precompute`.
:::

::: function
int _acb_ode_sum_done(acb_ode_sum_struct * sum, slong stride, acb_ode_bound_t bound, acb_ode_group_bound_t gbound, slong prec)

Return nonzero when either the sum has converged (within a tolerance based
on `prec`) or the current working precision is insufficient to significantly
improve the current estimate.
Update some heuristic convergence estimates stored in `sum`, assuming that
`stride` terms have been added since the last update.
When returning nonzero, also update the rigorous tail bounds if they have
improved (unless the `ACB_ODE_APPROX` flag is set in `sum`).

This function assumes the same preconditions as
`_acb_ode_sum_update_tail_bounds`.
Additionally, it requires `sum->n >= sum->n0 + d` where `d` is the maximum
degree of the coefficients of the differential operator written in terms of
the Euler derivation.
(In other words, it is illegal to call this function after the coefficient
of index `n - d` has been discarded.)

If `sum->nder` is zero, this function always returns zero.
:::

::: function
void acb_ode_sol_init(acb_ode_sol_t sol, slong nshifts, slong nlogs, slong npts, slong nder)

Initialize an `acb_ode_sol_t` object suitable for a solution
belonging to a solution group with `nshifts` shifts, involving $\log(x)$ to
powers less than `nlogs`, the derivatives of order up to `nder - 1` of
which are to be evaluated at `npts` points.
For rigorous error bounds, `nlogs` should be greater than or equal to the
length as a polynomial in $\log(x)$ of the mathematical solution even if
only a partial sum of lower degree in $\log(x)$ is to be computed.
:::

::: function
void acb_ode_sol_clear(acb_ode_sol_t sol)

Clear `sol`.
:::

::: function
void acb_ode_sol_fit_length(acb_ode_sol_t sol, slong len)

Ensure that the `series` field of `sol` has space for `len` coefficients.
:::

::: function
void acb_ode_sol_zero(acb_ode_sol_t sol)

Zero the series coefficients and partial sums (without touching the
initial conditions), reset bounds and termination flag.
:::

::: function
acb_ptr acb_ode_sol_sum_ptr(const acb_ode_sol_t sol, slong p, slong k, slong i)

Return a pointer to the coefficient of index `i` in the Taylor expansion
at the evaluation point of index `p` of the partial sum of the coefficient
of $\log(x)^k/k!$.
:::

::: function
void acb_ode_sol_set_ini(acb_ode_sol_t sol, acb_srcptr ini, const acb_ode_group_t group)

Set the initial values stored in `sol` as required for computing the solution of
coordinates `ini` in the echelonized basis of solutions from `group`.
:::

::: function
void _acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo, acb_srcptr sums, slong stride, mag_srcptr tb, const acb_t pt, slong nlogs, slong ord, slong nfrobshifts, slong prec)
void acb_ode_sol_jet(acb_poly_struct * val, const acb_t expo, const acb_ode_sol_t sol, slong p, const acb_t pt, slong ord, slong nfrobshifts, slong prec)

Compute the jets of order `ord` at the evaluation point `pt` of
the Frobenius shifts of order 0 to `nfrobshifts` - 1 of the solution `sol`,
based on the partial sums currently stored in the `sol` object.
The integer `p` is the index of the partial sum at `pt` in `sol`.
The Frobenius shift of order `j` of $y(x)$ is defined as the series obtained
by replacing $\log(x)^k/k!$ for $k \geq j$ by $\log(x)^{k-j}/(k-j)!$ in the
series expansion and dropping the terms with $k < j$.

The underscore version assumes that `ord` > 0, that the length with respect
to $\log(x)$ of the complete series solution (as opposed to the current
truncation) is at most `nlogs`, and that `sums + k * stride` points to the
jet of the coefficient of $\log(x)^k/k!$ for all `k < nlogs`.
:::

::: function
slong acb_ode_sol_estimate_terms(mag_t est, const acb_ode_sol_t sol, slong off, slong len, const mag_t radpow)

Set `est` to a heuristic estimate of the magnitude of the coefficients at
positions `off` to `off + len - 1` in the coefficient buffer of `sol`,
multiplied by `radpow`.
Return an estimate of the relative accuracy in bits of these coefficients.
:::

::: function
slong acb_ode_sol_estimate_sums(mag_t mag, mag_t rad, const acb_ode_sol_t sol)

Set `mag` to a heuristic estimate of the magnitude of the partial sums
encoded by `sol` and set `rad` to an estimate of their radii.
Return an estimate of the minimum relative accuracy in bits of the sums.
:::

## Functions for working with differential operators {#functions-for-working-with-differential-operators}

::: function
void acb_ode_local_diffop_ctx_init(gr_ctx_t Pol, gr_ore_poly_ctx_t Euler, gr_ctx_t Dop, gr_ctx_t Cst)

Given a ring `Dop` of differential operators over a univariate polynomial ring,
initialize `Pol` to a univariate polynomial ring over `Cst` and `Euler` to a
ring of differential operators in the Euler derivation over `Pol`, with
generator names derived from those of `Dop`.
:::

::: function
int acb_ode_local_diffop(gr_ore_poly_t ldop, ulong * flags, const gr_ore_poly_t dop, gr_ctx_t Dop, gr_ore_poly_ctx_t Euler)

Set `ldop` to an operator in `Euler` equal to `dop` up to multiplication from
the left by a power of `x`, removing common power-of-$x$ factors from the
coefficients of the result. Update the flag set `* flags` by setting the flags
`ACB_ODE_NONZERO`, `ACB_ODE_REGULAR`, and `ACB_ODE_ORDINARY` when the operator
is determined, respectively, to be nonzero, to have a regular (in the sense of
ordinary or regular singular) point at the origin, and to have an ordinary point
at the origin.
:::

::: function
int acb_ode_lc_roots(acb_ptr lcrt, const gr_ore_poly_t dop, gr_ctx_t Dop, int include_zero, gr_ctx_t CC)

Write to `lcrt` enclosures of the roots of the leading coefficient of `dop`
(belonging to `Dop`), with multiple roots repeated. (For an operator in the
standard or the Euler derivation and whose coefficients have no nontrivial
common factor, this is the same as its singular points.) Roots at the origin are
only listed if the `include_zero` flag is set. The return value is a `gr`-style
status code.
:::

::: function
void _acb_ode_apply_diffop_basecase_precomp(acb_poly_struct * g, slong goff, acb_srcptr weights, slong weights_nlogs, const acb_poly_struct * f, slong foff, slong flen, slong nlogs, slong start, slong len, slong prec)
void acb_ode_apply_diffop_basecase(acb_poly_struct * g, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, const acb_poly_struct * f, slong nlogs, slong start, slong len, slong prec)
void _acb_ode_apply_diffop_polmul(acb_poly_struct * g, slong goff, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, const acb_poly_struct * f, slong foff, slong flen, slong nlogs, slong start, slong len, slong prec)
void acb_ode_apply_diffop_polmul(acb_poly_struct * g, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, const acb_poly_struct * f, slong nlogs, slong start, slong len, slong prec)

Compute a slice of the coefficients of the image of a series by a
differential operator.

More precisely, the operator `dop` is applied to the logarithmic series

$$
f = x^{\mathit{expo} + \mathit{offset}} \sum_{p,k} f[k,p] x^p \log(x)^k/k!,
$$

resulting in a logarithmic series

$$
x^{\mathit{expo} + \mathit{offset}}
    \sum_{m',k'} g[k',m'] x^{m'} \log(x)^{k'}/k'!
$$

whose terms of index $m' = start + p'$ for $0 \leq p' < \mathit{len}$ are
written into `g`.

The non-underscore versions represent `f` and `g` as arrays of polynomials
corresponding to the coefficients of $\log(x)^k/k!$.

The underscore versions read the coefficients of $f[k]$ starting at offset
`foff` in the coefficient array of the element of index $k$ of the array `f`
and *add* the coefficients of $g[k]$ to the entries starting at offset
`goff` in the coefficient array of the element of index $k$ of `g`.

The `expo` parameter can be `NULL`; this is understood as zero.

The `polmul` version uses fast polynomial multiplication. The `basecase`
version uses a quadratic-time (w.r.t. `len`) algorithm. The `precomp`
version is similar to the `basecase` version but takes as input, instead of
the differential operator, an array of weights that can be precomputed using
`_acb_ode_apply_diffop_basecase_weights`.
:::

::: function
void _acb_ode_apply_diffop_basecase_weights(acb_ptr weights, const acb_poly_struct * dop, slong dop_len, acb_srcptr expo, slong offset, slong flen, slong nlogs, slong start, slong len, slong prec)

Set `weights` to an array of `len` × `nlogs` × `flen` weights that can be
used to apply the differential operator `dop` to multiple series using
`_acb_ode_apply_diffop_basecase_precomp`.
See that function for the meaning of the remaining parameters.
The `expo` parameter can be `NULL`; this is understood as zero.
:::

## Other summation utilities {#other-summation-utilities}

::: function
void _acb_ode_poly_negdivrevhigh(acb_ptr res, acb_srcptr a, acb_srcptr cst, acb_srcptr b, slong len, slong prec)

Solve the equation

$$
a(1/X) \cdot \mathit{res}(X) + \mathit{cst} \cdot b(X) = (\text{terms of negative degree}),
$$

where `a`, `res`, `b` are polynomials of length `len`. The `cst` pointer may
be null; this is understood as `cst = 1`.
:::

::: function
void acb_ode_poly_taylor_shift_aps_trunc(acb_poly_t g, const acb_poly_t f, acb_srcptr a, slong n, slong len, slong prec)

Shift the variable by `a + n`, truncate the result to length `len`.
:::

## Random generation {#random-generation}

The following function generates differential operators with a high probability
of being regular singular with a nontrivial exponent structure.

::: function
void acb_ode_randtest_acb(acb_poly_struct * dop, acb_struct * lcroots, acb_ode_exponents_t expos, flint_rand_t state, slong disp, slong lcdeg, slong clen, slong len, slong prec)

Write to `dop` a random differential operator, to `lcroots` the roots of
its leading coefficient, and to `expos` its exponent structure.
The generated operator has length `len`, a leading coefficient of degree
`lcdeg`, remaining coefficients of length at most `clen`, and exponent dispersion at most `disp`.
:::

## Error bounds {#error-bounds}

The functions in this section operate simultaneously on one or more logarithmic
series solutions of the differential equation $L(y) = 0$ belonging to the same
*group*, that is, of the form

$$
y(x) = x^\lambda \sum_{k=0}^{K-1} y_k(x) \frac{\log(x)^k}{k!}
     = x^\lambda \sum_{k=0}^{K-1} \sum_{j=0}^{\infty} c_{k,j} x^j \frac{\log(x)^k}{k!}.
$$

with the same $\lambda$.
Given one such solution, the main goal is to bound the tails

$$
y_{k,n:}(a) = \sum_{j=n}^{\infty} c_{k,j} a^j, \qquad 0 \leq k < K,
$$

of the power series $y_k$ truncated to some order $n$ and evaluated at
$a \in \mathbb C$.

The tail bound is computed from the *normalized residual* $q(x)$ defined as
follows. Let $Q_0$ be the monic indicial polynomial of $L$ at the origin.
The equation

$$
(Q_0(x\cdot\mathrm d/\mathrm dx))(q) = L(y_{n:})
$$

has a unique solution of the form

$$
q(x) = x^{\lambda+n} \sum_{k} q_{k}(x) \frac{\log(x)^k}{k!}
     = x^{\lambda+n} \sum_{k,j} q_{k,j} x^j \frac{\log(x)^k}{k!}
$$

where the sums are finite and $q_{k,j} = 0$ whenever $\lambda + n + j$ is a root
of multiplicity strictly more than $k$ of $Q_0$.

All functions in this section assume that the leading coefficient with respect
to the Euler derivation of the differential equation of interest does not
vanish at the origin. (For an equation with coprime coefficients, this is
equivalent to the origin being an ordinary point or a regular singular point.)

::: type
acb_ode_bound_struct
acb_ode_bound_t

Stores precomputed bounds independent of the particular solution and the
truncation order.
:::

::: type
acb_ode_group_bound_struct
acb_ode_group_bound_t

Stores additional precomputed bounds common to solutions belonging to the
same group.
:::

::: function
void acb_ode_bound_init(acb_ode_bound_t bound)
void acb_ode_bound_clear(acb_ode_bound_t bound)
:::

::: function
void acb_ode_group_bound_init(acb_ode_group_bound_t gbound)
void acb_ode_group_bound_clear(acb_ode_group_bound_t gbound)
:::

::: function
void acb_ode_bound_precompute(acb_ode_bound_t bound, const acb_poly_struct * dop, slong dop_len, acb_srcptr lcrt, slong pol_part_len, slong prec)

Fill `bound` with data used when computing bounds on the tails of series
solutions of a given equation.

The input consists of the operator `dop` of length `dop_len`,
a vector `lcrt` containing enclosures of the roots of the leading
coefficient (with multiple roots repeated),
a minimum number `pol_part_len` of initial terms of the local series
expansion of `dop` to be considered individually,
and a working precision `prec`.

As a general rule, increasing `pol_part_len` leads to sharper bounds but
makes both the precomputation and subsequent evaluation of the bound slower.
Values of `prec` beyond 30–50 bits are rarely useful except for equations
with extremely large coefficients and/or clusters of singular points.

When this function is called repeatedly with the same destination object,
`dop` must correspond to the same operator (possibly represented at
different precision) for all calls.
In addition, changing `pol_part_len` invalidates any
`acb_ode_group_bound_t` objects derived from `bound`.
:::

::: function
void acb_ode_group_bound_precompute(acb_ode_group_bound_t gbound, const acb_poly_struct * dop, slong dop_len, const acb_ode_exponents_t expos, slong grp, const acb_ode_bound_t bound, slong prec)

Fill `gbound` with data used for computing bounds on solutions from a given
group.

The input consists of the operator `dop` of length `dop_len`,
its local exponent structure `expos`,
the index `grp` of the group in the exponent structure,
a group-independent bound data object computed with
`acb_ode_bound_precompute`,
and a numeric working precision `prec`.
:::

::: function
void acb_ode_tail_bound_jet(arb_poly_t res, const acb_ode_group_t group, const acb_ode_bound_t bound, const acb_ode_group_bound_t gbound, slong n, slong nlogs, const arb_poly_t nres_maj, const arb_t rad, slong ord, slong prec)
void acb_ode_tail_bound_jet_precomp(arb_poly_t res, const acb_ode_bound_t bound, slong n, const arb_poly_t itg_pol, const arb_poly_t itg_num, const arb_poly_t nres_maj, const arb_t rad, slong ord, slong prec)

Compute upper bounds on the first `ord` coefficients of the Taylor
expansion with respect to $z$ of $y_{k,n:}(a + z)$ where $|a|$ is bounded
by `rad`.
The truncation order $n$ is assumed to be larger than any integer $s$ such that
$Q_0(\lambda + s) = 0$.
The returned bounds are valid uniformly for all `k`, for all evaluation
points $a$ in the disk of radius `rad` around zero, and for all
solutions $y$ that belong to the given `group` and satisfy the following
conditions:

- The polynomial `nres_maj` is a common majorant of the
  coefficients $q_k(x)$ of the normalized residual, i.e., for all $j$
  and $k$, the absolute value $|q_{k,j}|$ is bounded by the coefficient
  of $x^j$ in `nres_maj`.
- The degree in $\log(x)$ of the truncation at order $x^n$ of $y(x)$ is
  strictly smaller than `nlogs`. In other words, $c_{k,j} = 0$ when $j < n$
  and $k \geq nlogs$.

The caller must supply global and group-specific bound objects `bound`,
`gbound` computed as described above.

The `precomp` variant is suitable for repeated calls with the same `group`,
the same value of `nlogs`, and values of `n` bounded from below.
It takes two polynomials `itg_pol` and `itg_num` obtained by calling
`acb_ode_bound_precompute_integrand` with $n_0 \leq n$.
Typically, values of $n_0$ closer to `n` yield sharper bounds, but the effect
becomes less pronounced as `n` grows.
:::

::: function
void acb_ode_bound_precompute_integrand(arb_poly_t itg_pol, arb_poly_t itg_num, const acb_ode_group_t group, const acb_ode_bound_t bound, const acb_ode_group_bound_t gbound, slong n0, slong nlogs, slong prec)

Compute polynomials `itg_pol` and `itg_num` that can be passed to
`acb_ode_tail_bound_jet_precomp` to obtain bounds on the tails of
solutions from the group `group` at truncation orders >= `n0`.

The remaining input arguments have the same meaning as the corresponding
arguments of `acb_ode_tail_bound_jet`.
:::

## Bounds on rational sequences {#bounds-on-rational-sequences}

The functions in this section compute upper bounds of the expression

$$
F_T(n) = n \cdot \sum_{t=0}^{T-1} \left| [X^t] \frac{p(n+X)}{X^{-m(n)} q(n+X)} \right|
$$

where $q$ is a monic polynomial with $q(0) = 0$, $m(n)$ denotes the multiplicity
of $n$ as a root of $q$, and $p$ is a polynomial of strictly smaller degree
than $q$.
All functions work with a vector of numerators $p$ and a single denominator $q$.

The function and parameter names reflect that, in the intended application,
$q$ is the monic indicial polynomial of the differential equation being solved,
recentered at one of its roots.

Some functions accept a working precision *prec* for intermediate computations
even though their outputs are represented in fixed precision.

::: type
acb_ode_ind_lbound_struct
acb_ode_ind_lbound_t

A precomputed lower bound on $|q(k)|/k^r$ where $q$ is a monic polynomial
of degree $r$ with $q(0)=0$.
The data structure contains an integer `length` and a vector `r` of
length `length` of structures each containing a root `root` of $q$,
its multiplicity `mult`, an integer `n_min` such that `|1-root/k|` is
nondecreasing for `k >= n_min`, and a lower bound `global_lbound` of
`|1-root/k|^mult` valid for all `k >= 1` except the integer roots of $q$.
:::

::: function
void acb_ode_ind_lbound_init(acb_ode_ind_lbound_t ind_lbound)
void acb_ode_ind_lbound_clear(acb_ode_ind_lbound_t ind_lbound)

Initialize, respectively clear `ind_lbound`.
:::

::: function
void acb_ode_ind_lbound_precompute(acb_ode_ind_lbound_t ind_lbound, const acb_ode_exponents_t expos, slong grp, slong prec)

Fill `ind_lbound` with data corresponding to $q(X) = Q(\lambda + X)$ where
$Q$ is the monic polynomial whose roots are given by `expos` and $\lambda$
is the leader of the group of index `grp`.
:::

::: function
void acb_ode_ind_lbound_eval(mag_t res, const acb_ode_ind_lbound_t ind_lbound, slong n, slong prec)

Set `res` to a lower bound of $|q(k)|/k^r$ for all $k \geq n$ that are not
roots of $q$, with $q$ as above.
:::

In the following descriptions $\tau(n) = \sum_{k=0}^{n} m(k)$.

::: type
acb_ode_stairs_struct
acb_ode_stairs_t
:::

::: function
void acb_ode_stairs_init(acb_ode_stairs_t stairs)
void acb_ode_stairs_clear(acb_ode_stairs_t stairs)

Initialize, respectively clear `stairs`.
:::

::: function
void acb_ode_stairs_precompute(acb_ode_stairs_t stairs, const acb_poly_struct * num, slong len, const acb_poly_t ind, const acb_ode_group_t group, const acb_ode_ind_lbound_t ind_lbound, slong prec)

Fill `stairs` with data expected by `acb_ode_bound_rat_vec`.
The parameters have the same meaning as in that function.
The caller must supply an `ind_lbound` structure precomputed as described
above.
:::

::: function
void acb_ode_bound_rat_vec(mag_ptr res, const acb_poly_struct * num, slong len, const acb_poly_t ind, const acb_ode_group_t group, const acb_ode_ind_lbound_t ind_lbound, const acb_ode_stairs_t stairs, slong n0, slong ord, slong prec)

For each $p$ in the vector `num` of length `len`, compute a bound on
$F_{\rho(n)}(n)$ valid for all $n \geq n_0$, where $\rho(n) = \mathit{ord}$
as long as there is no integer root of $q$ between $n_0$ and $n$
(inclusive), and $\rho(n) = \tau(n)$ otherwise.
The parameter `ind` specifies the common denominator $q$.
The integer `n0` is assumed to be nonnegative.

The caller must supply the integer roots of $q$ in the parameter `group`
(the `leader` field is ignored), as well as two data structures `ind_lbound`
and `stairs` precomputed as described above.
:::

::: function
void acb_ode_bound_rat_ordinary_vec(mag_ptr res, const acb_poly_struct * num, slong len, const acb_poly_t ind, const acb_ode_ind_lbound_t ind_lbound, slong n0, slong ord, slong prec)

Similar but assume $m(n_0) = 0$ and compute a bound on $F_{T}(n)$ valid
for all $n \geq n_0$ with $m(n) = 0$.
The value of $T$ is given by the parameter `ord`.
:::

::: function
void acb_ode_bound_rat_ref_vec(mag_ptr res, const acb_poly_struct * num, slong len, const acb_poly_t ind, slong n, slong mult, slong ord, slong prec)

Similar but assume that $m(n)$ is equal to `mult` and compute a bound on
$F_{T}(n)$ at the given (nonnegative) `n`.
:::

## Miscellaneous bounds and estimates {#miscellaneous-bounds-and-estimates}

::: function
void _acb_ode_solution_growth(mag_t order, mag_t base, const acb_poly_struct * dop, slong dop_len)

Estimate parameters `order` and `base` such that the coefficient sequences
of series solutions of `dop` are
$\mathrm O(\mathit{base}^n/n!^{1/\mathit{order}})$.
Assumes that the leading coefficient of `dop` is a constant.
The outputs are currently not rigorous bounds.
:::

::: function
slong acb_ode_choose_prec(slong * rec_prec, const acb_poly_struct * dop, slong dop_len, mag_srcptr rad, mag_srcptr cvrad, slong tgt_prec)

Set `* rec_prec` to the default working precision for computing the coefficients
of a series solution of `dop` that converges on a disk of radius `cvrad` when
that solution is to be evaluated on points of magnitude `rad` to a target
precision `tgt_prec`, and return the default working precision at which to
accumulate the sum.
:::

::: type
acb_ode_cvest_struct
acb_ode_cvest_t

Stores running estimates of the rate of decay of the coefficients of a series,
the rate at which their computation is losing precision, and other related
quantities.
:::

::: function
void acb_ode_cvest_init(acb_ode_cvest_t cvest)
void acb_ode_cvest_clear(acb_ode_cvest_t cvest)

Initialize, respectively clear `cvest`.
:::

::: function
void acb_ode_cvest_update(acb_ode_cvest_t cvest, const acb_ode_cvest_t old, mag_t est, slong accuracy, slong stride, slong prec, slong work_prec)

Compute new estimates based on those in `old`, assuming that `stride` additional
terms have been computed, that the last few terms of the sum have magnitude
about `est` and are accurate to about `accuracy` bits, that the current working
precision for the computation of the coefficients is `work_prec` and that the
target precision is `prec`.
:::
