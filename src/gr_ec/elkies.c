/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Elkies' step: for a prime l at which E has a rational l-isogeny, the
    kernel polynomial h of that isogeny, of degree (l - 1)/2 instead of the
    (l^2 - 1)/2 of the division polynomial.

    l is such a prime -- an Elkies prime -- exactly when Phi_l(j(E), Y) has
    a root j~ in F_q, which is the j-invariant of the isogenous curve. From
    j~ and the partial derivatives of Phi_l at (j, j~), the equations of
    Ramanujan for the derivatives of the Eisenstein series give, in the
    normalisation a = -3 E4, b = -2 E6 (theta = q d/dq):

        theta j   = -j E6 / E4
        theta j~  = -Phi_X theta j / Phi_Y
        r         = E6~/E4~ = -theta j~ / (l j~)
        E4~       = j~ r^2 / (j~ - 1728),   E6~ = r E4~
        A~ = -3 l^4 E4~,   B~ = -2 l^6 E6~

    for the curve of the normalised isogeny; differentiating
    Phi(j(tau), j(l tau)) = 0 a second time fixes E2 - l E2~, and the sum
    of the x-coordinates of the (l - 1)/2 kernel points up to sign is

        s1 = -(l/2) (E2 - l E2~).

    Both were derived here from those equations, and checked: the first by
    j(A~, B~) = j~, the second, including the constant -l/2, against curves
    whose kernel is rational so that its x-coordinates are roots of psi_l.

    The remaining power sums come from Velu's formula,

        wp~(z) = wp(z) + sum over Q in C \ {0} of (wp(z + Q) - wp(Q)),

    whose coefficient of z^(2k) reads

        c~_k - c_k = (2 / (2k)!) sum over Q in C* up to sign of wp^(2k)(Q),

    with c_k the Laurent coefficients of wp and wp^(2k) = D_k(wp) for the
    polynomials D_0 = x, D_{k+1} = 4 f D_k'' + 2 f' D_k', f = x^3 + a x + b
    (wp'^2 = 4 f(wp) and wp'' = 2 f'(wp)).
    D_k has degree k + 1 and leading coefficient (2k + 1)!, so equation k
    gives s_{k+1} from the lower power sums, and Newton's identities turn
    s_1, ..., s_d into h.

    The characteristic has to exceed 2 l + 3 or so for the divisions; point
    counting uses l far below the characteristic.
*/

#include "fmpz.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_vec.h"
#include "fmpz_vec.h"
#include "gr_ec.h"
#include "impl.h"

/* Laurent coefficients c_1, ..., c_n of wp for y^2 = x^3 + a x + b */
static int
_wp_coeffs(gr_ptr c, slong n, gr_srcptr a, gr_srcptr b, gr_ctx_t R)
{
    slong sz = R->sizeof_elem, k, m;
    gr_ptr t;
    int status = GR_SUCCESS;

#define C(i) GR_ENTRY(c, (i) - 1, sz)

    GR_TMP_INIT(t, R);

    /* c_1 = g2/20 = -a/5, c_2 = g3/28 = -b/7 */
    if (n >= 1)
        status |= gr_div_si(C(1), a, -5, R);
    if (n >= 2)
        status |= gr_div_si(C(2), b, -7, R);

    /* c_k = 3 / ((k - 2)(2k + 3)) sum_{m=1}^{k-2} c_m c_{k-1-m} */
    for (k = 3; k <= n && status == GR_SUCCESS; k++)
    {
        status |= gr_zero(C(k), R);

        for (m = 1; m <= k - 2; m++)
        {
            status |= gr_mul(t, C(m), C(k - 1 - m), R);
            status |= gr_add(C(k), C(k), t, R);
        }

        status |= gr_mul_ui(C(k), C(k), 3, R);
        status |= gr_div_ui(C(k), C(k), (ulong) ((k - 2) * (2 * k + 3)), R);
    }

#undef C

    GR_TMP_CLEAR(t, R);

    return status;
}

/*
    The value and partial derivatives up to order two of
    Phi(X, Y) = sum_k X^k phi_k(Y), k < n, at (x, y).
*/
static int
_bivariate_partials(gr_ptr P0, gr_ptr PX, gr_ptr PY, gr_ptr PXX, gr_ptr PXY,
        gr_ptr PYY, const gr_poly_struct * phi, slong n, gr_srcptr x,
        gr_srcptr y, gr_ctx_t R)
{
    slong k;
    gr_ptr v0, v1, v2, xk, xk1, xk2, w;
    gr_poly_t g;
    int status = GR_SUCCESS;

    GR_TMP_INIT4(v0, v1, v2, w, R);
    GR_TMP_INIT3(xk, xk1, xk2, R);
    gr_poly_init(g, R);

    status |= gr_zero(P0, R);
    status |= gr_zero(PX, R);
    status |= gr_zero(PY, R);
    status |= gr_zero(PXX, R);
    status |= gr_zero(PXY, R);
    status |= gr_zero(PYY, R);

    /* xk = x^k, xk1 = k x^(k-1), xk2 = k (k - 1) x^(k-2) */
    status |= gr_one(xk, R);
    status |= gr_zero(xk1, R);
    status |= gr_zero(xk2, R);

    for (k = 0; k < n && status == GR_SUCCESS; k++)
    {
        status |= gr_poly_evaluate(v0, phi + k, y, R);
        status |= gr_poly_derivative(g, phi + k, R);
        status |= gr_poly_evaluate(v1, g, y, R);
        status |= gr_poly_derivative(g, g, R);
        status |= gr_poly_evaluate(v2, g, y, R);

        status |= gr_mul(w, xk, v0, R);
        status |= gr_add(P0, P0, w, R);
        status |= gr_mul(w, xk, v1, R);
        status |= gr_add(PY, PY, w, R);
        status |= gr_mul(w, xk, v2, R);
        status |= gr_add(PYY, PYY, w, R);
        status |= gr_mul(w, xk1, v0, R);
        status |= gr_add(PX, PX, w, R);
        status |= gr_mul(w, xk1, v1, R);
        status |= gr_add(PXY, PXY, w, R);
        status |= gr_mul(w, xk2, v0, R);
        status |= gr_add(PXX, PXX, w, R);

        /* step to k + 1: (x^k)'' = x (x^(k-1))'' + 2 (x^(k-1))' etc. */
        status |= gr_mul(xk2, xk2, x, R);
        status |= gr_mul_two(w, xk1, R);
        status |= gr_add(xk2, xk2, w, R);
        status |= gr_mul(xk1, xk1, x, R);
        status |= gr_add(xk1, xk1, xk, R);
        status |= gr_mul(xk, xk, x, R);
    }

    gr_poly_clear(g, R);
    GR_TMP_CLEAR3(xk, xk1, xk2, R);
    GR_TMP_CLEAR4(v0, v1, v2, w, R);

    return status;
}

/*
    The normalised isogenous curve and the kernel's s1, from a root jt of
    Phi_l(j, Y). GR_UNABLE when the formulas degenerate -- j or j~ in
    {0, 1728}, or a vanishing partial derivative, which happens at the
    singular points of the modular curve.
*/
static int
_elkies_isogeny(gr_ptr At, gr_ptr Bt, gr_ptr s1, gr_srcptr a, gr_srcptr b,
        gr_srcptr j, gr_srcptr jt, const gr_poly_struct * phi, ulong l,
        gr_ctx_t R)
{
    gr_ptr PX, PY, PXX, PXY, PYY, E4, E6, tj, tjt, E4t, E6t, r, u, w, v0, rhs, K;
    int status = GR_SUCCESS;

    GR_TMP_INIT5(PX, PY, PXX, PXY, PYY, R);
    GR_TMP_INIT5(E4, E6, tj, tjt, E4t, R);
    GR_TMP_INIT5(E6t, r, u, w, v0, R);
    GR_TMP_INIT2(rhs, K, R);

    /* the special j-invariants make E4 or E6 vanish on one side */
    if (gr_is_zero(a, R) != T_FALSE || gr_is_zero(b, R) != T_FALSE
            || gr_is_zero(jt, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    status |= gr_sub_ui(u, jt, 1728, R);

    if (status != GR_SUCCESS || gr_is_zero(u, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* the partials of Phi = sum X^k phi_k(Y) at (X, Y) = (j, j~) */
    status |= _bivariate_partials(v0, PX, PY, PXX, PXY, PYY, phi, l + 2,
            j, jt, R);

    if (status != GR_SUCCESS || gr_is_zero(PX, R) != T_FALSE
            || gr_is_zero(PY, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* E4 = -a/3, E6 = -b/2, theta j = -j E6/E4 */
    status |= gr_div_si(E4, a, -3, R);
    status |= gr_div_si(E6, b, -2, R);
    status |= gr_div(tj, E6, E4, R);
    status |= gr_mul(tj, tj, j, R);
    status |= gr_neg(tj, tj, R);

    /* theta j~ = -PX theta j / PY */
    status |= gr_mul(tjt, PX, tj, R);
    status |= gr_div(tjt, tjt, PY, R);
    status |= gr_neg(tjt, tjt, R);

    /* r = -theta j~ / (l j~), E4~ = j~ r^2 / (j~ - 1728), E6~ = r E4~ */
    status |= gr_mul_ui(u, jt, l, R);
    status |= gr_div(r, tjt, u, R);
    status |= gr_neg(r, r, R);

    status |= gr_sqr(E4t, r, R);
    status |= gr_mul(E4t, E4t, jt, R);
    status |= gr_sub_ui(u, jt, 1728, R);
    status |= gr_div(E4t, E4t, u, R);
    status |= gr_mul(E6t, r, E4t, R);

    if (status != GR_SUCCESS || gr_is_zero(E4t, R) != T_FALSE
            || gr_is_zero(E6t, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* A~ = -3 l^4 E4~, B~ = -2 l^6 E6~ */
    status |= gr_set_ui(u, l, R);
    status |= gr_sqr(w, u, R);
    status |= gr_sqr(w, w, R);                  /* l^4 */
    status |= gr_mul(At, E4t, w, R);
    status |= gr_mul_si(At, At, -3, R);
    status |= gr_mul(w, w, u, R);
    status |= gr_mul(w, w, u, R);               /* l^6 */
    status |= gr_mul(Bt, E6t, w, R);
    status |= gr_mul_si(Bt, Bt, -2, R);

    /*
        (K/6)(E2 - l E2~) = -[ PXX tj^2 + 2 PXY tj tjt + PYY tjt^2
                               + PX j (2/3 E6^2/E4^2 + E4/2)
                               + PY l^2 j~ (2/3 E6~^2/E4~^2 + E4~/2) ]
        with K = PY l j~ E6~/E4~
    */
    status |= gr_sqr(rhs, tj, R);
    status |= gr_mul(rhs, rhs, PXX, R);
    status |= gr_mul(u, tj, tjt, R);
    status |= gr_mul(u, u, PXY, R);
    status |= gr_mul_two(u, u, R);
    status |= gr_add(rhs, rhs, u, R);
    status |= gr_sqr(u, tjt, R);
    status |= gr_mul(u, u, PYY, R);
    status |= gr_add(rhs, rhs, u, R);

    status |= gr_div(u, E6, E4, R);
    status |= gr_sqr(u, u, R);
    status |= gr_mul_ui(u, u, 2, R);
    status |= gr_div_ui(u, u, 3, R);
    status |= gr_div_ui(w, E4, 2, R);
    status |= gr_add(u, u, w, R);
    status |= gr_mul(u, u, j, R);
    status |= gr_mul(u, u, PX, R);
    status |= gr_add(rhs, rhs, u, R);

    status |= gr_div(u, E6t, E4t, R);
    status |= gr_sqr(u, u, R);
    status |= gr_mul_ui(u, u, 2, R);
    status |= gr_div_ui(u, u, 3, R);
    status |= gr_div_ui(w, E4t, 2, R);
    status |= gr_add(u, u, w, R);
    status |= gr_mul(u, u, jt, R);
    status |= gr_mul(u, u, PY, R);
    status |= gr_mul_ui(u, u, l * l, R);
    status |= gr_add(rhs, rhs, u, R);

    status |= gr_div(K, E6t, E4t, R);
    status |= gr_mul(K, K, jt, R);
    status |= gr_mul(K, K, PY, R);
    status |= gr_mul_ui(K, K, l, R);

    /* E2 - l E2~ = -6 rhs / K, and s1 = -(l/2)(E2 - l E2~) = 3 l rhs / K */
    status |= gr_mul_ui(s1, rhs, 3 * l, R);
    status |= gr_div(s1, s1, K, R);

cleanup:
    GR_TMP_CLEAR2(rhs, K, R);
    GR_TMP_CLEAR5(E6t, r, u, w, v0, R);
    GR_TMP_CLEAR5(E4, E6, tj, tjt, E4t, R);
    GR_TMP_CLEAR5(PX, PY, PXX, PXY, PYY, R);

    return status;
}

/*
    The same from a root g of Mueller's canonical polynomial Phi^c_l(X, j)
    (see modular_poly.c), with f = l^s (eta(l tau)/eta(tau))^(2s) and
    theta f / f = (s/12) D where D = l E2~ - E2:

        theta j = -j E6/E4,   theta f = -Phi_J theta j / Phi_X,
        D = 12 theta f / (s g),                      s1 = l D / 2,

    differentiating Phi^c(f, j) = 0 once more, where the E2 of theta^2 f
    and of theta^2 j cancel,

        l^2 E4~ = (s + 1) D^2 + E4
                  + 144 / (s g Phi_X) (U - Phi_J theta j V),
        U = Phi_XX theta f^2 + 2 Phi_XJ theta f theta j + Phi_JJ theta j^2,
        V = 2 E6/(3 E4) + E4^2/(2 E6),

    then Delta~ = g^(12/s) Delta / l^12 and j~ = E4~^3 / Delta~, and from
    Phi^c(l^s/f, j~) = 0, which holds because l^s/f(tau) = f(-1/(l tau)),

        E6~ = -Phi_X(f*, j~) f* s D E4~ / (12 l j~ Phi_J(f*, j~)),  f* = l^s/g.

    Derived here like the classical ones, and checked the same way. E6~ is
    then also held to E6~^2 = E4~^3 - 1728 Delta~, which costs nothing and
    catches any degenerate root the conditions above miss.
*/
static int
_elkies_isogeny_canonical(gr_ptr At, gr_ptr Bt, gr_ptr s1, gr_srcptr a,
        gr_srcptr b, gr_srcptr j, gr_srcptr g, const gr_poly_struct * phi,
        ulong l, ulong s, gr_ctx_t R)
{
    gr_ptr P0, PX, PY, PXX, PXY, PYY, E4, E6, Dl, tj, tf, D, U, V, E4t, E6t,
        Dlt, jt, fs, u, w;
    int status = GR_SUCCESS;

    GR_TMP_INIT5(P0, PX, PY, PXX, PXY, R);
    GR_TMP_INIT5(PYY, E4, E6, Dl, tj, R);
    GR_TMP_INIT5(tf, D, U, V, E4t, R);
    GR_TMP_INIT5(E6t, Dlt, jt, fs, u, R);
    GR_TMP_INIT(w, R);

    if (gr_is_zero(a, R) != T_FALSE || gr_is_zero(b, R) != T_FALSE
            || gr_is_zero(g, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    status |= _bivariate_partials(P0, PX, PY, PXX, PXY, PYY, phi, l + 2,
            g, j, R);

    if (status != GR_SUCCESS || gr_is_zero(PX, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* E4 = -a/3, E6 = -b/2, Delta = (E4^3 - E6^2)/1728 */
    status |= gr_div_si(E4, a, -3, R);
    status |= gr_div_si(E6, b, -2, R);
    status |= gr_pow_ui(Dl, E4, 3, R);
    status |= gr_sqr(u, E6, R);
    status |= gr_sub(Dl, Dl, u, R);
    status |= gr_div_ui(Dl, Dl, 1728, R);

    /* theta j, theta f, D */
    status |= gr_div(tj, E6, E4, R);
    status |= gr_mul(tj, tj, j, R);
    status |= gr_neg(tj, tj, R);

    status |= gr_mul(tf, PY, tj, R);
    status |= gr_div(tf, tf, PX, R);
    status |= gr_neg(tf, tf, R);

    status |= gr_mul_ui(u, g, s, R);
    status |= gr_div(D, tf, u, R);
    status |= gr_mul_ui(D, D, 12, R);

    /* U and V */
    status |= gr_sqr(U, tf, R);
    status |= gr_mul(U, U, PXX, R);
    status |= gr_mul(u, tf, tj, R);
    status |= gr_mul(u, u, PXY, R);
    status |= gr_mul_two(u, u, R);
    status |= gr_add(U, U, u, R);
    status |= gr_sqr(u, tj, R);
    status |= gr_mul(u, u, PYY, R);
    status |= gr_add(U, U, u, R);

    status |= gr_div(V, E6, E4, R);
    status |= gr_mul_ui(V, V, 2, R);
    status |= gr_div_ui(V, V, 3, R);
    status |= gr_sqr(u, E4, R);
    status |= gr_div(u, u, E6, R);
    status |= gr_div_ui(u, u, 2, R);
    status |= gr_add(V, V, u, R);

    /* l^2 E4~ */
    status |= gr_mul(u, PY, tj, R);
    status |= gr_mul(u, u, V, R);
    status |= gr_sub(U, U, u, R);
    status |= gr_mul(u, g, PX, R);
    status |= gr_mul_ui(u, u, s, R);
    status |= gr_div(U, U, u, R);
    status |= gr_mul_ui(U, U, 144, R);

    status |= gr_sqr(E4t, D, R);
    status |= gr_mul_ui(E4t, E4t, s + 1, R);
    status |= gr_add(E4t, E4t, E4, R);
    status |= gr_add(E4t, E4t, U, R);
    status |= gr_div_ui(E4t, E4t, l * l, R);

    /* Delta~ = g^(12/s) Delta / l^12, j~ = E4~^3 / Delta~ */
    status |= gr_pow_ui(Dlt, g, 12 / s, R);
    status |= gr_mul(Dlt, Dlt, Dl, R);
    status |= gr_set_ui(u, l, R);
    status |= gr_pow_ui(u, u, 12, R);
    status |= gr_div(Dlt, Dlt, u, R);

    if (status != GR_SUCCESS || gr_is_zero(Dlt, R) != T_FALSE
            || gr_is_zero(E4t, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    status |= gr_pow_ui(jt, E4t, 3, R);
    status |= gr_div(jt, jt, Dlt, R);

    /* the partials of Phi^c at (f*, j~) */
    status |= gr_set_ui(fs, l, R);
    status |= gr_pow_ui(fs, fs, s, R);
    status |= gr_div(fs, fs, g, R);

    status |= _bivariate_partials(P0, PX, PY, PXX, PXY, PYY, phi, l + 2,
            fs, jt, R);

    if (status != GR_SUCCESS || gr_is_zero(PY, R) != T_FALSE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* E6~ = -PX* f* s D E4~ / (12 l j~ PY*) */
    status |= gr_mul(E6t, PX, fs, R);
    status |= gr_mul(E6t, E6t, D, R);
    status |= gr_mul(E6t, E6t, E4t, R);
    status |= gr_mul_ui(E6t, E6t, s, R);
    status |= gr_mul(u, jt, PY, R);
    status |= gr_mul_ui(u, u, 12 * l, R);
    status |= gr_div(E6t, E6t, u, R);
    status |= gr_neg(E6t, E6t, R);

    /* E6~^2 = E4~^3 - 1728 Delta~ */
    status |= gr_sqr(u, E6t, R);
    status |= gr_pow_ui(w, E4t, 3, R);
    status |= gr_sub(w, w, u, R);
    status |= gr_mul_ui(u, Dlt, 1728, R);

    if (status != GR_SUCCESS || gr_is_zero(E6t, R) != T_FALSE
            || gr_equal(w, u, R) != T_TRUE)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* A~ = -3 l^4 E4~, B~ = -2 l^6 E6~, s1 = l D / 2 */
    status |= gr_set_ui(u, l, R);
    status |= gr_pow_ui(w, u, 4, R);
    status |= gr_mul(At, E4t, w, R);
    status |= gr_mul_si(At, At, -3, R);
    status |= gr_pow_ui(w, u, 6, R);
    status |= gr_mul(Bt, E6t, w, R);
    status |= gr_mul_si(Bt, Bt, -2, R);

    status |= gr_mul_ui(s1, D, l, R);
    status |= gr_div_ui(s1, s1, 2, R);

cleanup:
    GR_TMP_CLEAR(w, R);
    GR_TMP_CLEAR5(E6t, Dlt, jt, fs, u, R);
    GR_TMP_CLEAR5(tf, D, U, V, E4t, R);
    GR_TMP_CLEAR5(PYY, E4, E6, Dl, tj, R);
    GR_TMP_CLEAR5(P0, PX, PY, PXX, PXY, R);

    return status;
}

/*
    The kernel polynomial of degree d = (l - 1)/2 of the normalised isogeny
    from y^2 = x^3 + a x + b to y^2 = x^3 + At x + Bt, given the sum s1 of
    its roots.
*/
static int
_kernel_polynomial(gr_poly_t h, gr_srcptr a, gr_srcptr b, gr_srcptr At,
        gr_srcptr Bt, gr_srcptr s1, ulong l, gr_ctx_t R)
{
    slong sz = R->sizeof_elem, d = (l - 1) / 2, k, m;
    gr_ptr c, ct, s, sig, t, fact;
    gr_poly_t D, Dp, Dpp, f, fp, tmp;
    int status = GR_SUCCESS;

    c = gr_heap_init_vec(FLINT_MAX(d, 1), R);
    ct = gr_heap_init_vec(FLINT_MAX(d, 1), R);
    s = gr_heap_init_vec(d + 1, R);
    sig = gr_heap_init_vec(d + 1, R);
    GR_TMP_INIT2(t, fact, R);

    gr_poly_init(D, R);
    gr_poly_init(Dp, R);
    gr_poly_init(Dpp, R);
    gr_poly_init(f, R);
    gr_poly_init(fp, R);
    gr_poly_init(tmp, R);

#define SS(i) GR_ENTRY(s, (i), sz)
#define SG(i) GR_ENTRY(sig, (i), sz)
#define CC(i) GR_ENTRY(c, (i) - 1, sz)
#define CT(i) GR_ENTRY(ct, (i) - 1, sz)

    status |= gr_set_ui(SS(0), d, R);

    if (d >= 1)
        status |= gr_set(SS(1), s1, R);

    if (d >= 2)
    {
        status |= _wp_coeffs(c, d - 1, a, b, R);
        status |= _wp_coeffs(ct, d - 1, At, Bt, R);
    }

    /* f = x^3 + a x + b, then 2 f' = 6 x^2 + 2 a = wp'', and D_0 = x */
    status |= gr_poly_zero(f, R);
    status |= gr_poly_set_coeff_si(f, 3, 1, R);
    status |= gr_poly_set_coeff_scalar(f, 1, a, R);
    status |= gr_poly_set_coeff_scalar(f, 0, b, R);
    status |= gr_poly_derivative(fp, f, R);
    status |= gr_poly_mul_si(fp, fp, 2, R);

    status |= gr_poly_zero(D, R);
    status |= gr_poly_set_coeff_si(D, 1, 1, R);

    /* (2k)! as it goes */
    status |= gr_one(fact, R);

    for (k = 1; k <= d - 1 && status == GR_SUCCESS; k++)
    {
        /* D_k = 4 f D_{k-1}'' + 2 f' D_{k-1}' */
        status |= gr_poly_derivative(Dp, D, R);
        status |= gr_poly_derivative(Dpp, Dp, R);
        status |= gr_poly_mul(D, f, Dpp, R);
        status |= gr_poly_mul_si(D, D, 4, R);
        status |= gr_poly_mul(tmp, fp, Dp, R);
        status |= gr_poly_add(D, D, tmp, R);

        status |= gr_mul_ui(fact, fact, (2 * k - 1) * (2 * k), R);

        /* sum_{m<=k} [D_k]_m s_m = (c~_k - c_k) (2k)!/2 - [D_k]_{k+1} s_{k+1} */
        status |= gr_sub(t, CT(k), CC(k), R);
        status |= gr_mul(t, t, fact, R);
        status |= gr_div_ui(t, t, 2, R);

        for (m = 0; m <= k; m++)
        {
            gr_ptr dm;
            GR_TMP_INIT(dm, R);
            status |= gr_poly_get_coeff_scalar(dm, D, m, R);
            status |= gr_mul(dm, dm, SS(m), R);
            status |= gr_sub(t, t, dm, R);
            GR_TMP_CLEAR(dm, R);
        }

        /* divide by the leading coefficient (2k + 1)! */
        {
            gr_ptr lc;
            GR_TMP_INIT(lc, R);
            status |= gr_poly_get_coeff_scalar(lc, D, k + 1, R);
            status |= gr_div(SS(k + 1), t, lc, R);
            GR_TMP_CLEAR(lc, R);
        }
    }

    /* Newton: k sigma_k = sum_{i=1}^k (-1)^(i-1) sigma_{k-i} s_i */
    status |= gr_one(SG(0), R);

    for (k = 1; k <= d && status == GR_SUCCESS; k++)
    {
        status |= gr_zero(SG(k), R);

        for (m = 1; m <= k; m++)
        {
            status |= gr_mul(t, SG(k - m), SS(m), R);

            if (m % 2)
                status |= gr_add(SG(k), SG(k), t, R);
            else
                status |= gr_sub(SG(k), SG(k), t, R);
        }

        status |= gr_div_ui(SG(k), SG(k), k, R);
    }

    /* h = sum_k (-1)^k sigma_k x^(d-k) */
    status |= gr_poly_zero(h, R);

    for (k = 0; k <= d && status == GR_SUCCESS; k++)
    {
        status |= gr_set(t, SG(k), R);

        if (k % 2)
            status |= gr_neg(t, t, R);

        status |= gr_poly_set_coeff_scalar(h, d - k, t, R);
    }

#undef SS
#undef SG
#undef CC
#undef CT

    gr_poly_clear(tmp, R);
    gr_poly_clear(fp, R);
    gr_poly_clear(f, R);
    gr_poly_clear(Dpp, R);
    gr_poly_clear(Dp, R);
    gr_poly_clear(D, R);

    GR_TMP_CLEAR2(t, fact, R);
    gr_heap_clear_vec(sig, d + 1, R);
    gr_heap_clear_vec(s, d + 1, R);
    gr_heap_clear_vec(ct, FLINT_MAX(d, 1), R);
    gr_heap_clear_vec(c, FLINT_MAX(d, 1), R);

    return status;
}

/*
    Phi(X, j) as a polynomial in X, for Phi = sum_k X^k phi_k(Y), k < n
    (the classical polynomial is symmetric, so this is also Phi_l(j, Y)).
*/
static int
_modular_poly_at_j(gr_poly_t g, const gr_poly_struct * phi, slong n,
        gr_srcptr j, gr_ctx_t R)
{
    slong k;
    int status = GR_SUCCESS;

    gr_poly_fit_length(g, n, R);

    for (k = 0; k < n; k++)
        status |= gr_poly_evaluate(GR_ENTRY(g->coeffs, k, R->sizeof_elem),
                phi + k, j, R);

    _gr_poly_set_length(g, n, R);
    _gr_poly_normalise(g, R);

    return status;
}

/* the roots of g in F_q: those of gcd(X^q - X, g); also X^q mod g */
static int
_rational_roots(gr_vec_t roots, gr_poly_t Xq, const gr_poly_t g,
        const gr_poly_preinv_t P, const fmpz_t q, gr_ctx_t R)
{
    gr_poly_t t, split;
    fmpz_vec_t mult;
    int status = GR_SUCCESS;

    gr_poly_init(t, R);
    gr_poly_init(split, R);
    fmpz_vec_init(mult, 0);

    status |= gr_poly_preinv_powmod_x_fmpz(Xq, q, P, R);
    status |= gr_poly_set_coeff_si(split, 1, 1, R);
    status |= gr_poly_sub(t, Xq, split, R);
    status |= gr_poly_gcd(split, t, g, R);

    if (status == GR_SUCCESS && gr_poly_length(split, R) >= 2)
        status |= gr_poly_roots(roots, mult, split, 0, R);
    else
        gr_vec_set_length(roots, 0, R);

    fmpz_vec_clear(mult);
    gr_poly_clear(split, R);
    gr_poly_clear(t, R);

    return status;
}

/*
    At an Atkin prime, Frobenius acts on the l + 1 isogenies -- the roots
    of Phi^c_l(X, j) -- without a fixed point, which forces its matrix to
    have conjugate eigenvalues outside F_l; their ratio has some order r,
    every orbit then has exactly r elements, and Phi^c_l(X, j) splits into
    irreducible factors all of degree r. r pins t^2/q modulo l down to a
    handful of values (atkin.c).

    That reasoning needs the roots to be distinct, so r is only reported
    when g is squarefree and its distinct-degree factorisation is a single
    degree dividing l + 1; *r = 0 otherwise.
*/
static int
_atkin_order(ulong * r, const gr_poly_t g, const gr_poly_t Xq,
        const gr_poly_preinv_t P, ulong l, const fmpz_t q, gr_ctx_t R)
{
    gr_poly_t dg, c;
    gr_poly_vec_t fac;
    fmpz_vec_t degs;
    int status = GR_SUCCESS;

    *r = 0;

    gr_poly_init(dg, R);
    gr_poly_init(c, R);
    gr_poly_vec_init(fac, 0, R);
    fmpz_vec_init(degs, 0);

    status |= gr_poly_derivative(dg, g, R);
    status |= gr_poly_gcd(c, g, dg, R);

    if (status == GR_SUCCESS && gr_poly_length(c, R) == 1)
    {
        status |= _gr_poly_factor_distinct_deg_with_frob(fac, degs, g, P, Xq, q, R);

        if (status == GR_SUCCESS && fac->length == 1
                && gr_poly_length(fac->entries, R) == (slong) l + 2
                && fmpz_cmp_ui(degs->entries, 1) > 0
                && fmpz_abs_fits_ui(degs->entries)
                && (l + 1) % fmpz_get_ui(degs->entries) == 0)
            *r = fmpz_get_ui(degs->entries);
    }

    fmpz_vec_clear(degs);
    gr_poly_vec_clear(fac, R);
    gr_poly_clear(c, R);
    gr_poly_clear(dg, R);

    return status;
}

/*
    If l is an Elkies prime for E, sets h to the kernel polynomial of a
    rational l-isogeny and *elkies to 1; otherwise *elkies to 0.

    l is Elkies exactly when the canonical polynomial Phi^c_l(X, j) has a
    root in F_q (each of its roots is the value of f at one of the l + 1
    isogenies, as each root of Phi_l(j, Y) is the j-invariant of one);
    taking the gcd with X^q - X first keeps the root finding to a
    polynomial of degree 1, 2 or l + 1, and decides Elkies against Atkin on
    its own. When the canonical formulas degenerate for every root, the
    classical polynomial is tried as well.

    Returns GR_UNABLE when l is Elkies but every formula degenerates, and
    when E is not a short model; the caller then falls back to Schoof's
    step for this l.
*/
int
_gr_ec_elkies_kernel(gr_poly_t h, int * elkies, ulong * atkin_r, ulong l,
        const fmpz_t q, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_poly_struct * phi;
    gr_poly_t g, Xq;
    gr_poly_preinv_t P;
    gr_vec_t roots;
    gr_ptr j, At, Bt, s1;
    ulong s;
    slong k, i;
    int status = GR_SUCCESS, done = 0, pass;

    *elkies = 0;

    if (atkin_r != NULL)
        *atkin_r = 0;

    if (gr_ec_ctx_model(ctx) != GR_EC_SHORT_WEIERSTRASS || l < 3 || l % 2 == 0)
        return GR_UNABLE;

    phi = flint_malloc((l + 2) * sizeof(gr_poly_struct));

    for (k = 0; k < (slong) l + 2; k++)
        gr_poly_init(phi + k, R);

    gr_poly_init(g, R);
    gr_poly_init(Xq, R);
    gr_poly_preinv_init(P, R);
    gr_vec_init(roots, 0, R);
    GR_TMP_INIT4(j, At, Bt, s1, R);

    status |= gr_ec_ctx_j_invariant(j, ctx);

    /* pass 0: canonical polynomial; pass 1: classical */
    for (pass = 0; pass < 2 && !done && status == GR_SUCCESS; pass++)
    {
        if (pass == 0)
            status |= _gr_ec_modular_polynomial_canonical(phi, &s, l, R);
        else
            status |= _gr_ec_modular_polynomial(phi, l, R);

        status |= _modular_poly_at_j(g, phi, l + 2, j, R);

        if (status != GR_SUCCESS || gr_poly_length(g, R) < 2)
        {
            status = GR_UNABLE;
            break;
        }

        status |= gr_poly_preinv_set(P, g, R);
        status |= _rational_roots(roots, Xq, g, P, q, R);

        if (status != GR_SUCCESS)
            break;

        if (roots->length == 0)
        {
            /* Atkin; the canonical pass is the only one that gets here */
            if (atkin_r != NULL)
                status |= _atkin_order(atkin_r, g, Xq, P, l, q, R);
            break;
        }

        *elkies = 1;

        /* any root will do; the first one whose formulas do not degenerate */
        for (i = 0; i < roots->length && !done; i++)
        {
            gr_srcptr r = GR_ENTRY(roots->entries, i, R->sizeof_elem);
            int st;

            if (pass == 0)
                st = _elkies_isogeny_canonical(At, Bt, s1, GR_EC_A4(ctx),
                        GR_EC_A6(ctx), j, r, phi, l, s, R);
            else
                st = _elkies_isogeny(At, Bt, s1, GR_EC_A4(ctx),
                        GR_EC_A6(ctx), j, r, phi, l, R);

            if (st == GR_SUCCESS && _kernel_polynomial(h, GR_EC_A4(ctx),
                        GR_EC_A6(ctx), At, Bt, s1, l, R) == GR_SUCCESS)
                done = 1;
        }
    }

    if (status == GR_SUCCESS && *elkies && !done)
        status = GR_UNABLE;

    GR_TMP_CLEAR4(j, At, Bt, s1, R);
    gr_vec_clear(roots, R);
    gr_poly_preinv_clear(P, R);
    gr_poly_clear(Xq, R);
    gr_poly_clear(g, R);

    for (k = 0; k < (slong) l + 2; k++)
        gr_poly_clear(phi + k, R);

    flint_free(phi);

    return status;
}

/*
    Frobenius acts on the kernel of a rational isogeny as multiplication by
    an eigenvalue lambda: phi(P) = lambda P for the generic point P of the
    kernel. Found by the same discrete logarithm Schoof uses, but modulo a
    polynomial of degree (l - 1)/2 instead of (l^2 - 1)/2.
*/
static int
_elkies_step(ulong * lambda, ulong l, const fmpz_t q,
        gr_ec_tors_ctx_struct * C, gr_ec_ctx_t FLINT_UNUSED(ctx))
{
    gr_ec_tors_struct P, phiP;
    int status = GR_SUCCESS, found = 0;

    _gr_ec_tors_init(&P, C);
    _gr_ec_tors_init(&phiP, C);

    status |= _gr_ec_tors_frobenius(&P, &phiP, q, C);

    if (status == GR_SUCCESS)
        status |= _gr_ec_tors_discrete_log(lambda, &found, &phiP, &P, l, C);

    /* Frobenius is invertible, so lambda = 0 would be an error */
    if (status == GR_SUCCESS && (!found || *lambda == 0))
        status = GR_UNABLE;

    _gr_ec_tors_clear(&phiP, C);
    _gr_ec_tors_clear(&P, C);

    return status;
}

/*
    t mod l at an Elkies prime, from the eigenvalue: lambda is a root of
    X^2 - t X + q modulo l, so t = lambda + q / lambda. Sets *elkies to 0
    and does nothing else when l is an Atkin prime.
*/
int
_gr_ec_elkies_trace(ulong * tl, int * elkies, ulong * atkin_r, ulong l,
        const fmpz_t q, gr_ec_ctx_t ctx)
{
    gr_poly_t h;
    ulong lambda, ql;
    int status;

    gr_poly_init(h, GR_EC_ELEM_CTX(ctx));

    status = _gr_ec_elkies_kernel(h, elkies, atkin_r, l, q, ctx);

    if (status == GR_SUCCESS && *elkies)
    {
        status = _gr_ec_tors_solve(&lambda, h, l, q, _elkies_step, ctx);

        if (status == GR_SUCCESS)
        {
            ql = fmpz_fdiv_ui(q, l);
            *tl = n_addmod(lambda, n_mulmod2(ql, n_invmod(lambda, l), l), l);
        }
    }

    gr_poly_clear(h, GR_EC_ELEM_CTX(ctx));

    return status;
}
