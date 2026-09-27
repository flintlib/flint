/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    The classical modular polynomial Phi_l(X, Y), reduced into the base ring.

    Phi_l is the polynomial whose roots in X, for Y = j(tau), are the l + 1
    values j(l tau) and j((tau + k)/l), k = 0, ..., l - 1. Over Z its
    coefficients are enormous -- thousands of bits already for l around 50
    -- but point counting only ever wants it modulo p, so it is computed
    here directly in the base ring, from q-expansions, and the large
    integers never appear.

    Write q = e^(2 pi i tau), J = q j = E4^3 / prod (1 - q^n)^24, and
    split Phi_l(X, j(q)) = G(X) (X - j(q^l)) where G has the l roots
    j(zeta^k q^(1/l)). The power sums of those roots are integral:

        P_m = sum_k j(zeta^k q^(1/l))^m = l * (coefficients of j^m at
              exponents divisible by l, with q^(l t) read as q^t),

    so Newton's identities give the coefficients e_k of G, which are power
    series in q for k < l and have a single pole q^(-1) for k = l (a product
    of fewer than l roots has pole order below one, and pole orders here
    are integers). The coefficient of X^(l+1-k) in Phi_l is then

        a_k = (-1)^k (e_k + e_{k-1} j(q^l)),

    a Laurent series with pole order at most l + 1, and it equals a
    polynomial of degree at most l + 1 in j. Since j^d = q^(-d) + ..., the
    exponents -(l+1), ..., 0 of a_k determine that polynomial by peeling
    off the leading term repeatedly. j(q^l) only contributes q^(-l) and 744
    to those exponents, so e_k is needed up to q^l, which needs the powers
    of J up to q^(l^2 + l).

    Everything is exact: the only divisions are by 1, ..., l in Newton's
    identities, which is why the characteristic has to exceed l + 1.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ec.h"
#include "impl.h"

/* the coefficient of q^i in f, zero for i < 0 */
static int
_coeff(gr_ptr res, const gr_poly_t f, slong i, gr_ctx_t R)
{
    if (i < 0)
        return gr_zero(res, R);

    return gr_poly_get_coeff_scalar(res, f, i, R);
}

/*
    J = q j(q) = E4^3 / prod (1 - q^n)^24 to length len, with
    E4 = 1 + 240 sum sigma_3(n) q^n.
*/
int
_gr_ec_j_qexp(gr_poly_t J, slong len, gr_ctx_t R)
{
    gr_poly_t E4, D;
    fmpz_poly_t Dz;
    ulong * sigma3;
    slong d, m;
    int status = GR_SUCCESS;

    sigma3 = flint_calloc(len, sizeof(ulong));

    for (d = 1; d < len; d++)
        for (m = d; m < len; m += d)
            sigma3[m] += (ulong) d * (ulong) d * (ulong) d;

    gr_poly_init(E4, R);
    gr_poly_init(D, R);
    fmpz_poly_init(Dz);

    gr_poly_fit_length(E4, len, R);
    status |= gr_one(E4->coeffs, R);

    for (m = 1; m < len; m++)
    {
        gr_ptr c = GR_ENTRY(E4->coeffs, m, R->sizeof_elem);
        status |= gr_set_ui(c, sigma3[m], R);
        status |= gr_mul_ui(c, c, 240, R);
    }

    _gr_poly_set_length(E4, len, R);
    _gr_poly_normalise(E4, R);

    fmpz_poly_eta_qexp(Dz, 24, len);
    status |= gr_poly_set_fmpz_poly(D, Dz, R);

    status |= gr_poly_inv_series(D, D, len, R);
    status |= gr_poly_mullow(J, E4, E4, len, R);
    status |= gr_poly_mullow(J, J, E4, len, R);
    status |= gr_poly_mullow(J, J, D, len, R);

    fmpz_poly_clear(Dz);
    gr_poly_clear(D, R);
    gr_poly_clear(E4, R);
    flint_free(sigma3);

    return status;
}

int
_gr_ec_modular_polynomial(gr_poly_struct * phi, ulong l, gr_ctx_t R)
{
    slong sz = R->sizeof_elem;
    slong L, T, k, i, d, m, t, W;
    gr_poly_struct * Jpow;
    gr_poly_struct * P;
    gr_poly_struct * e;
    gr_poly_t Ehat, S, tmp, Jl1;
    gr_ptr a, c, u, v;
    fmpz_t ch;
    int status = GR_SUCCESS;

    if (l < 2)
        return GR_DOMAIN;

    /* Newton divides by 1, ..., l */
    fmpz_init(ch);

    if (gr_ctx_fq_prime(ch, R) != GR_SUCCESS || fmpz_cmp_ui(ch, l + 1) <= 0)
    {
        fmpz_clear(ch);
        return GR_UNABLE;
    }

    fmpz_clear(ch);

    T = l;                          /* e_k is needed up to q^T */
    L = l * l + l + 1;              /* J^m up to q^(l^2 + l) */
    W = l + 2;                      /* exponents -(l+1), ..., 0 */

    Jpow = flint_malloc((l + 1) * sizeof(gr_poly_struct));
    P = flint_malloc(l * sizeof(gr_poly_struct));
    e = flint_malloc(l * sizeof(gr_poly_struct));

    for (m = 0; m <= (slong) l; m++)
        gr_poly_init(Jpow + m, R);
    for (m = 0; m < (slong) l; m++)
    {
        gr_poly_init(P + m, R);
        gr_poly_init(e + m, R);
    }

    gr_poly_init(Ehat, R);
    gr_poly_init(S, R);
    gr_poly_init(tmp, R);
    gr_poly_init(Jl1, R);

    a = gr_heap_init_vec((l + 2) * W, R);
    GR_TMP_INIT3(c, u, v, R);

    /* J^m for m = 0, ..., l */
    status |= gr_poly_one(Jpow + 0, R);
    status |= _gr_ec_j_qexp(Jpow + 1, L, R);

    for (m = 2; m <= (slong) l && status == GR_SUCCESS; m++)
        status |= gr_poly_mullow(Jpow + m, Jpow + m - 1, Jpow + 1, L, R);

    /* J^(l+1), only as far as q^(l+1) */
    status |= gr_poly_mullow(Jl1, Jpow + l, Jpow + 1, l + 2, R);

    if (status != GR_SUCCESS)
        goto cleanup;

    /* power sums P_m[t] = l J^m[l t + m], m = 1, ..., l - 1, t = 0, ..., T */
    for (m = 1; m < (slong) l; m++)
    {
        gr_poly_fit_length(P + m, T + 1, R);

        for (t = 0; t <= T; t++)
        {
            gr_ptr dst = GR_ENTRY(P[m].coeffs, t, sz);
            status |= _coeff(dst, Jpow + m, l * t + m, R);
            status |= gr_mul_ui(dst, dst, l, R);
        }

        _gr_poly_set_length(P + m, T + 1, R);
        _gr_poly_normalise(P + m, R);
    }

    /* e_0, ..., e_{l-1}, as power series to q^T */
    status |= gr_poly_one(e + 0, R);

    for (k = 1; k < (slong) l && status == GR_SUCCESS; k++)
    {
        status |= gr_poly_zero(e + k, R);

        for (i = 1; i <= k; i++)
        {
            status |= gr_poly_mullow(tmp, e + k - i, P + i, T + 1, R);

            if (i % 2)
                status |= gr_poly_add(e + k, e + k, tmp, R);
            else
                status |= gr_poly_sub(e + k, e + k, tmp, R);
        }

        status |= gr_set_ui(c, k, R);
        status |= gr_poly_div_scalar(e + k, e + k, c, R);
    }

    /*
        e_l = q^(-1) Ehat. The i = l term of Newton's identity is e_0 P_l,
        and P_l = q^(-1) Ph with Ph[t] = l J^l[l t].
    */
    status |= gr_poly_zero(S, R);

    for (i = 1; i < (slong) l; i++)
    {
        status |= gr_poly_mullow(tmp, e + l - i, P + i, T + 1, R);

        if (i % 2)
            status |= gr_poly_add(S, S, tmp, R);
        else
            status |= gr_poly_sub(S, S, tmp, R);
    }

    status |= gr_poly_shift_left(Ehat, S, 1, R);

    gr_poly_fit_length(tmp, T + 2, R);
    for (t = 0; t <= T + 1; t++)
    {
        gr_ptr dst = GR_ENTRY(tmp->coeffs, t, sz);
        status |= _coeff(dst, Jpow + l, l * t, R);
        status |= gr_mul_ui(dst, dst, l, R);
    }
    _gr_poly_set_length(tmp, T + 2, R);
    _gr_poly_normalise(tmp, R);

    if (l % 2)                          /* (-1)^(l-1) = +1 for odd l */
        status |= gr_poly_add(Ehat, Ehat, tmp, R);
    else
        status |= gr_poly_sub(Ehat, Ehat, tmp, R);

    status |= gr_set_ui(c, l, R);
    status |= gr_poly_div_scalar(Ehat, Ehat, c, R);

    if (status != GR_SUCCESS)
        goto cleanup;

    /*
        a_k at exponents i = -(l+1), ..., 0, stored in row k at column
        i + l + 1. With E_k(i) the coefficient of q^i in e_k,

            a_k(i) = (-1)^k (E_k(i) + E_{k-1}(i + l) + 744 E_{k-1}(i)).
    */
#define A(kk, ii) GR_ENTRY(a, (kk) * W + (ii) + (slong) l + 1, sz)
#define EK(res, kk, ii) \
    (((kk) < (slong) l) ? _coeff((res), e + (kk), (ii), R) \
                        : _coeff((res), Ehat, (ii) + 1, R))

    status |= gr_one(A(0, 0), R);

    for (k = 1; k <= (slong) l + 1 && status == GR_SUCCESS; k++)
    {
        for (i = -((slong) l + 1); i <= 0; i++)
        {
            gr_ptr dst = A(k, i);

            if (k <= (slong) l)
                status |= EK(dst, k, i);
            else
                status |= gr_zero(dst, R);

            status |= EK(u, k - 1, i + (slong) l);
            status |= gr_add(dst, dst, u, R);

            status |= EK(u, k - 1, i);
            status |= gr_mul_ui(u, u, 744, R);
            status |= gr_add(dst, dst, u, R);

            if (k % 2)
                status |= gr_neg(dst, dst, R);
        }
    }

    if (status != GR_SUCCESS)
        goto cleanup;

    /*
        Each a_k as a polynomial in j: the coefficient of q^(-d) is the
        coefficient of Y^d once the higher powers have been taken off, and
        j^d has coefficient J^d[i + d] at q^i.
    */
    for (k = 0; k <= (slong) l + 1 && status == GR_SUCCESS; k++)
    {
        gr_poly_struct * out = phi + (l + 1 - k);

        status |= gr_poly_zero(out, R);
        gr_poly_fit_length(out, l + 2, R);

        for (d = (slong) l + 1; d >= 1; d--)
        {
            const gr_poly_struct * Jd = (d <= (slong) l) ? Jpow + d : Jl1;

            status |= gr_set(c, A(k, -d), R);
            status |= gr_set(GR_ENTRY(out->coeffs, d, sz), c, R);

            for (i = -d; i <= 0; i++)
            {
                status |= _coeff(v, Jd, i + d, R);
                status |= gr_mul(v, v, c, R);
                status |= gr_sub(A(k, i), A(k, i), v, R);
            }
        }

        status |= gr_set(out->coeffs, A(k, 0), R);

        _gr_poly_set_length(out, l + 2, R);
        _gr_poly_normalise(out, R);
    }

#undef A
#undef EK

cleanup:
    GR_TMP_CLEAR3(c, u, v, R);
    gr_heap_clear_vec(a, (l + 2) * W, R);

    gr_poly_clear(Jl1, R);
    gr_poly_clear(tmp, R);
    gr_poly_clear(S, R);
    gr_poly_clear(Ehat, R);

    for (m = 0; m < (slong) l; m++)
    {
        gr_poly_clear(P + m, R);
        gr_poly_clear(e + m, R);
    }
    for (m = 0; m <= (slong) l; m++)
        gr_poly_clear(Jpow + m, R);

    flint_free(e);
    flint_free(P);
    flint_free(Jpow);

    return status;
}

/*
    Mueller's canonical modular polynomial Phi^c_l(X, J), reduced into the
    base ring: the minimal polynomial over Q(j) of

        f(tau) = l^s (eta(l tau) / eta(tau))^(2s),   s = 12 / gcd(12, l - 1),

    which is an invariant of Gamma_0(l). It has degree l + 1 in X like the
    classical one, but only degree v = s (l - 1)/12 in J (between (l - 1)/12
    and (l - 1)/2), which makes it much cheaper to compute and evaluate.

    With q = t^l, f = l^s q^v prod (1 - q^(l n))^(2s) / (1 - q^n)^(2s) has
    a zero of order v, and its l other conjugates are H(zeta^k t) for

        H(t) = t^(-v) B(t),   B(t) = prod (1 - t^n)^(2s) / (1 - t^(l n))^(2s),

    each with a pole of order v/l in q. As for Phi_l, their power sums

        P_m = l * (coefficient of t^(l u + v m) in B^m, at q^u)

    and Newton's identities give their elementary symmetric functions e_k,
    and the coefficient of X^(l+1-k) is a_k = (-1)^k (e_k + e_{k-1} f).
    Here e_{k-1} f only has positive exponents for k <= l, so a_k is read
    off e_k at the exponents -v, ..., 0, and a_{l+1} is the constant
    (-1)^(l+1) l^s e_l[q^(-v)].

    As for Phi_l, the characteristic has to exceed l + 1. phi receives the
    l + 2 coefficients of X^0, ..., X^(l+1), polynomials in J, and *s_out
    the exponent s.
*/
int
_gr_ec_modular_polynomial_canonical(gr_poly_struct * phi, ulong * s_out,
        ulong l, gr_ctx_t R)
{
    slong sz = R->sizeof_elem;
    slong L, v, s, k, i, d, m, u, Ld;
    gr_poly_struct * Bpow;
    gr_poly_struct * P;
    gr_poly_struct * e;
    gr_poly_struct * Jpow;
    gr_poly_t num, den, tmp;
    fmpz_poly_t z;
    gr_ptr a, c, w;
    fmpz_t ch;
    int status = GR_SUCCESS;

    if (l < 3 || !n_is_prime(l))
        return GR_DOMAIN;

    fmpz_init(ch);

    if (gr_ctx_fq_prime(ch, R) != GR_SUCCESS || fmpz_cmp_ui(ch, l + 1) <= 0)
    {
        fmpz_clear(ch);
        return GR_UNABLE;
    }

    fmpz_clear(ch);

    s = 12 / n_gcd(12, l - 1);
    v = s * (l - 1) / 12;
    L = l * v + l;                  /* B^m below t^(l v + l) */
    *s_out = s;

    Bpow = flint_malloc((l + 1) * sizeof(gr_poly_struct));
    P = flint_malloc((l + 1) * sizeof(gr_poly_struct));
    e = flint_malloc((l + 1) * sizeof(gr_poly_struct));
    Jpow = flint_malloc((v + 1) * sizeof(gr_poly_struct));

    for (m = 0; m <= (slong) l; m++)
    {
        gr_poly_init(Bpow + m, R);
        gr_poly_init(P + m, R);
        gr_poly_init(e + m, R);
    }
    for (d = 0; d <= v; d++)
        gr_poly_init(Jpow + d, R);

    gr_poly_init(num, R);
    gr_poly_init(den, R);
    gr_poly_init(tmp, R);
    fmpz_poly_init(z);

    a = gr_heap_init_vec(v + 1, R);
    GR_TMP_INIT2(c, w, R);

    /* B = prod (1 - t^n)^(2s) / prod (1 - t^(l n))^(2s) to length L */
    fmpz_poly_eta_qexp(z, 2 * s, L);
    status |= gr_poly_set_fmpz_poly(num, z, R);

    Ld = (L + l - 1) / l;
    fmpz_poly_eta_qexp(z, 2 * s, Ld);
    gr_poly_fit_length(den, (Ld - 1) * l + 1, R);
    status |= _gr_vec_zero(den->coeffs, (Ld - 1) * l + 1, R);
    for (i = 0; i < Ld && i < z->length; i++)
        status |= gr_set_fmpz(GR_ENTRY(den->coeffs, i * l, sz), z->coeffs + i, R);
    _gr_poly_set_length(den, (Ld - 1) * l + 1, R);
    _gr_poly_normalise(den, R);

    status |= gr_poly_inv_series(den, den, L, R);
    status |= gr_poly_mullow(Bpow + 1, num, den, L, R);
    status |= gr_poly_one(Bpow + 0, R);

    for (m = 2; m <= (slong) l && status == GR_SUCCESS; m++)
        status |= gr_poly_mullow(Bpow + m, Bpow + m - 1, Bpow + 1, L, R);

    if (status != GR_SUCCESS)
        goto cleanup;

    /*
        P_i has a pole, so e_k at q^0 needs e_{k-i} beyond q^0. By the pole
        orders, e_k is only needed on the exponents lo(k), ..., hi(k), with

            lo(k) = -floor(k v / l),   hi(k) = ceil((l - k) v / l),

        and computing it there needs e_{k-i} and P_i only on theirs, the
        same windows (lo(a) + lo(b) >= lo(a + b), and the tail of either
        factor above its window only reaches exponents above hi(k)). Each
        window holds about v + 1 coefficients, so the l^2/2 products of
        Newton's identities are of that length.
    */
#define LO(k) (-(slong) (((k) * v) / l))
#define HI(k) ((slong) ((((slong) l - (k)) * v + l - 1) / l))
#define WLEN(k) (HI(k) - LO(k) + 1)

    for (m = 1; m <= (slong) l; m++)
    {
        gr_poly_fit_length(P + m, WLEN(m), R);

        for (u = LO(m); u <= HI(m); u++)
        {
            gr_ptr dst = GR_ENTRY(P[m].coeffs, u - LO(m), sz);
            slong idx = l * u + v * m;

            if (idx < 0 || idx >= L)
                status |= gr_zero(dst, R);
            else
            {
                status |= gr_poly_get_coeff_scalar(dst, Bpow + m, idx, R);
                status |= gr_mul_ui(dst, dst, l, R);
            }
        }

        _gr_poly_set_length(P + m, WLEN(m), R);
    }

    /* e_0 = 1, window [0, v] */
    gr_poly_fit_length(e + 0, WLEN(0), R);
    status |= _gr_vec_zero(e[0].coeffs, WLEN(0), R);
    status |= gr_one(e[0].coeffs, R);
    _gr_poly_set_length(e + 0, WLEN(0), R);

    gr_poly_fit_length(tmp, v + 2, R);

    for (k = 1; k <= (slong) l && status == GR_SUCCESS; k++)
    {
        gr_poly_fit_length(e + k, WLEN(k), R);
        status |= _gr_vec_zero(e[k].coeffs, WLEN(k), R);
        _gr_poly_set_length(e + k, WLEN(k), R);

        for (i = 1; i <= k; i++)
        {
            /* e_{k-i} P_i, whose index 0 is the exponent lo(k-i) + lo(i) */
            slong base = LO(k - i) + LO(i), shift = base - LO(k);
            slong n = HI(k) - base + 1;
            slong la = FLINT_MIN(WLEN(k - i), n), lb = FLINT_MIN(WLEN(i), n);
            gr_srcptr pa = e[k - i].coeffs, pb = P[i].coeffs;

            if (n <= 0)
                continue;

            n = FLINT_MIN(n, la + lb - 1);

            if (la < lb)
            {
                FLINT_SWAP(slong, la, lb);
                FLINT_SWAP(gr_srcptr, pa, pb);
            }

            gr_poly_fit_length(tmp, n, R);
            status |= _gr_poly_mullow(tmp->coeffs, pa, la, pb, lb, n, R);

            if (i % 2)
                status |= _gr_vec_add(GR_ENTRY(e[k].coeffs, shift, sz),
                        GR_ENTRY(e[k].coeffs, shift, sz), tmp->coeffs, n, R);
            else
                status |= _gr_vec_sub(GR_ENTRY(e[k].coeffs, shift, sz),
                        GR_ENTRY(e[k].coeffs, shift, sz), tmp->coeffs, n, R);
        }

        status |= gr_set_ui(c, k, R);
        status |= gr_inv(c, c, R);
        status |= _gr_vec_mul_scalar(e[k].coeffs, e[k].coeffs, WLEN(k), c, R);
    }

    /* J = q j to length v + 1, and its powers, for the peeling */
    status |= gr_poly_one(Jpow + 0, R);
    if (v >= 1)
        status |= _gr_ec_j_qexp(Jpow + 1, v + 1, R);
    for (d = 2; d <= v && status == GR_SUCCESS; d++)
        status |= gr_poly_mullow(Jpow + d, Jpow + d - 1, Jpow + 1, v + 1, R);

    if (status != GR_SUCCESS)
        goto cleanup;

    /*
        a_k = (-1)^k e_k at the exponents -v, ..., 0, as a polynomial of
        degree <= v in J: the coefficient of q^(-d) is that of J^d once the
        higher powers are taken off, and J^d = q^(-d) (q j)^d.
    */
    for (k = 0; k <= (slong) l; k++)
    {
        gr_poly_struct * out = phi + (l + 1 - k);

        /* a[i + v] = coefficient of q^i, i = -v, ..., 0 */
        for (i = 0; i <= v; i++)
        {
            slong x = i - v - LO(k);      /* the exponent i - v in e_k */

            if (x < 0)
                status |= gr_zero(GR_ENTRY(a, i, sz), R);
            else
                status |= gr_poly_get_coeff_scalar(GR_ENTRY(a, i, sz), e + k, x, R);

            if (k % 2)
                status |= gr_neg(GR_ENTRY(a, i, sz), GR_ENTRY(a, i, sz), R);
        }

        status |= gr_poly_zero(out, R);
        gr_poly_fit_length(out, v + 1, R);

        for (d = v; d >= 1; d--)
        {
            status |= gr_set(c, GR_ENTRY(a, v - d, sz), R);
            status |= gr_set(GR_ENTRY(out->coeffs, d, sz), c, R);

            /* subtract c q^(-d) J^d from the exponents -d, ..., 0 */
            for (i = -d; i <= 0; i++)
            {
                status |= gr_poly_get_coeff_scalar(w, Jpow + d, i + d, R);
                status |= gr_mul(w, w, c, R);
                status |= gr_sub(GR_ENTRY(a, i + v, sz), GR_ENTRY(a, i + v, sz), w, R);
            }
        }

        status |= gr_set(out->coeffs, GR_ENTRY(a, v, sz), R);
        _gr_poly_set_length(out, v + 1, R);
        _gr_poly_normalise(out, R);
    }

    /* a_{l+1} = (-1)^(l+1) l^s e_l[q^(-v)], a constant; l is odd, and
       the window of e_l starts at lo(l) = -v */
    status |= gr_poly_get_coeff_scalar(c, e + l, 0, R);
    status |= gr_set_ui(w, l, R);
    status |= gr_pow_ui(w, w, s, R);
    status |= gr_mul(c, c, w, R);
    status |= gr_poly_set_scalar(phi + 0, c, R);

#undef LO
#undef HI
#undef WLEN

cleanup:
    GR_TMP_CLEAR2(c, w, R);
    gr_heap_clear_vec(a, v + 1, R);

    fmpz_poly_clear(z);
    gr_poly_clear(tmp, R);
    gr_poly_clear(den, R);
    gr_poly_clear(num, R);

    for (d = 0; d <= v; d++)
        gr_poly_clear(Jpow + d, R);
    for (m = 0; m <= (slong) l; m++)
    {
        gr_poly_clear(Bpow + m, R);
        gr_poly_clear(P + m, R);
        gr_poly_clear(e + m, R);
    }

    flint_free(Jpow);
    flint_free(e);
    flint_free(P);
    flint_free(Bpow);

    return status;
}
