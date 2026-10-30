/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Jets (Taylor expansions with respect to the parameter a) of the
    incomplete gamma function series used for rational a.

    Throughout, a = a0 + eps where a0 = p/q is the given rational number and
    eps is the formal jet variable; output vectors contain the coefficients
    of eps^0, ..., eps^(len-1). Error bounds are computed as majorant
    series: each coefficient of the truncation error is bounded separately,
    and the returned error is a bound valid simultaneously for all
    coefficients.
*/

#include <math.h>
#include "fmpq.h"
#include "arb.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"

/* Polynomial multiplication tuned for the shapes occurring in the binary
   splitting trees. The default algorithm selection in _arb_poly_mullow is
   based on prec, but the coefficients at the lower levels of the tree are
   much smaller than prec, where the classical algorithm is best; at the
   upper levels the coefficients are large and the block algorithm wins
   even for short series (it multiplies the integer polynomials with
   fmpz_poly, which uses classical multiplication with fft_small
   transforms for large coefficients: each coefficient is transformed
   once instead of once per coefficient product). */
static void
arb_poly_mullow_bsplit(arb_poly_t res, const arb_poly_t A, const arb_poly_t B, slong len, slong prec)
{
    slong lenA, lenB, m, n, bits, cutoff;

    lenA = FLINT_MIN(A->length, len);
    lenB = FLINT_MIN(B->length, len);

    if (lenA == 0 || lenB == 0)
    {
        arb_poly_zero(res);
        return;
    }

    n = FLINT_MIN(lenA + lenB - 1, len);
    m = FLINT_MIN(lenA, lenB);

    if (m == 1)
    {
        arb_poly_mullow(res, A, B, len, prec);
        return;
    }

    bits = FLINT_MAX(_arb_vec_bits(A->coeffs, lenA), _arb_vec_bits(B->coeffs, lenB));
    bits = FLINT_MIN(bits, prec);

    /* measured crossover between the classical and block algorithms */
    if (m <= 3)
        cutoff = 40000;
    else if (m == 4)
        cutoff = 24000;
    else if (m <= 7)
        cutoff = 12000;
    else
        cutoff = 40000 / m;

    if (res == A || res == B)
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        if (bits < cutoff)
            _arb_poly_mullow_classical(t->coeffs, A->coeffs, lenA, B->coeffs, lenB, n, prec);
        else
            _arb_poly_mullow_block(t->coeffs, A->coeffs, lenA, B->coeffs, lenB, n, prec);
        _arb_poly_set_length(t, n);
        _arb_poly_normalise(t);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(res, n);
        if (bits < cutoff)
            _arb_poly_mullow_classical(res->coeffs, A->coeffs, lenA, B->coeffs, lenB, n, prec);
        else
            _arb_poly_mullow_block(res->coeffs, A->coeffs, lenA, B->coeffs, lenB, n, prec);
        _arb_poly_set_length(res, n);
        _arb_poly_normalise(res);
    }
}

#define MAG_PREC 30

/* Helpers ****************************************************************/

/* upper bound for |x| */
static void
_mag_set_fmpq(mag_t res, const fmpq_t x)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz(res, fmpq_numref(x));
    mag_set_fmpz_lower(t, fmpq_denref(x));
    mag_div(res, res, t);
    mag_clear(t);
}

/* lower bound for |x| */
static void
_mag_set_fmpq_lower(mag_t res, const fmpq_t x)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz_lower(res, fmpq_numref(x));
    mag_set_fmpz(t, fmpq_denref(x));
    mag_div_lower(res, res, t);
    mag_clear(t);
}

/* upper bound for z^e e^(-z) where z > 0 is given by a ball */
static void
_mag_pow_fmpq_expinv(mag_t res, const arb_t z, const fmpq_t e, int with_exp)
{
    arb_t t, u;
    arb_init(t);
    arb_init(u);
    arb_set_round(t, z, MAG_PREC);
    if (!arb_is_positive(t))
    {
        mag_inf(res);
    }
    else
    {
        arb_pow_fmpq(u, t, e, MAG_PREC);
        if (with_exp)
        {
            arb_neg(t, t);
            arb_exp(t, t, MAG_PREC);
            arb_mul(u, u, t, MAG_PREC);
        }
        arb_get_mag(res, u);
    }
    arb_clear(t);
    arb_clear(u);
}

/* upper bound for |log(t)| for t in the ball z (requires z > 0) */
static void
_mag_abs_log_arb(mag_t res, const arb_t z)
{
    mag_t lo, hi, u;
    mag_init(lo);
    mag_init(hi);
    mag_init(u);
    arb_get_mag_lower(lo, z);
    arb_get_mag(hi, z);
    if (mag_is_zero(lo))
    {
        mag_inf(res);
    }
    else
    {
        mag_log(u, hi);         /* log(max(1, hi)) */
        mag_neg_log(res, lo);   /* -log(min(1, lo)) */
        mag_max(res, res, u);
    }
    mag_clear(lo);
    mag_clear(hi);
    mag_clear(u);
}

/* V <- V * (c + eps) */
static void
_mag_vec_mul_linear(mag_ptr V, const mag_t c, slong len)
{
    slong j;
    for (j = len - 1; j >= 1; j--)
    {
        mag_mul(V + j, V + j, c);
        mag_add(V + j, V + j, V + j - 1);
    }
    mag_mul(V, V, c);
}

/* V <- V / (c - eps), c is a lower bound for a positive number */
static void
_mag_vec_div_linear(mag_ptr V, const mag_t c, slong len)
{
    slong j;
    mag_div(V, V, c);
    for (j = 1; j < len; j++)
    {
        mag_add(V + j, V + j, V + j - 1);
        mag_div(V + j, V + j, c);
    }
}

/* res <- V * exp(L eps), res and V not aliased */
static void
_mag_vec_mul_exp(mag_ptr res, mag_srcptr V, const mag_t L, slong len)
{
    slong i, j;
    mag_t c, t;

    mag_init(c);
    mag_init(t);

    for (j = 0; j < len; j++)
        mag_zero(res + j);

    mag_one(c);  /* L^i / i! */
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            mag_mul(c, c, L);
            mag_div_ui(c, c, i);
        }

        for (j = i; j < len; j++)
        {
            mag_mul(t, V + j - i, c);
            mag_add(res + j, res + j, t);
        }
    }

    mag_clear(c);
    mag_clear(t);
}

static void
_mag_vec_max(mag_t res, mag_srcptr V, slong len)
{
    slong j;
    mag_zero(res);
    for (j = 0; j < len; j++)
        mag_max(res, res, V + j);
}

/* res = sum_{0 <= j < len} L^j / j! */
static void
_mag_exp_jet_sum(mag_t res, const mag_t L, slong len)
{
    slong j;
    mag_t c;
    mag_init(c);
    mag_one(c);
    mag_one(res);
    for (j = 1; j < len; j++)
    {
        mag_mul(c, c, L);
        mag_div_ui(c, c, j);
        mag_add(res, res, c);
    }
    mag_clear(c);
}

/* res = z^(a + eps) = z^a exp(eps log(z)); logz = log(z) may be NULL */
static void
_arb_pow_fmpq_jet(arb_ptr res, const arb_t z, const arb_t logz, const fmpq_t a, slong len, slong prec)
{
    slong j;

    arb_pow_fmpq(res, z, a, prec);

    if (len > 1)
    {
        arb_t t;
        arb_init(t);
        if (logz == NULL)
            arb_log(t, z, prec);
        else
            arb_set_round(t, logz, prec);
        for (j = 1; j < len; j++)
        {
            arb_mul(res + j, res + j - 1, t, prec);
            arb_div_ui(res + j, res + j, j, prec);
        }
        arb_clear(t);
    }
}

static void
_arb_poly_get_vec(arb_ptr res, const arb_poly_t f, slong len)
{
    slong j;
    for (j = 0; j < len; j++)
    {
        if (j < f->length)
            arb_set(res + j, f->coeffs + j);
        else
            arb_zero(res + j);
    }
}

/* f = c0 + c1 * eps (fmpz coefficients) */
static void
_arb_poly_set_linear_fmpz(arb_poly_t f, const fmpz_t c0, const fmpz_t c1)
{
    arb_poly_fit_length(f, 2);
    arb_set_fmpz(f->coeffs, c0);
    arb_set_fmpz(f->coeffs + 1, c1);
    _arb_poly_set_length(f, 2);
    _arb_poly_normalise(f);
}

static void
_arb_vec_add_error_mag(arb_ptr res, slong len, const mag_t err)
{
    slong j;
    for (j = 0; j < len; j++)
        arb_add_error_mag(res + j, err);
}

/* Asymptotic series at infinity *****************************************/

/*
    Gamma(a,z) - z^(a-1) e^(-z) sum_{k<N} (a-1)...(a-k) / z^k
        = (a-1)(a-2)...(a-N) Gamma(a-N, z).

    For z >= 1 and N >= a0 - 1, the coefficients of Gamma(a0-N+eps, z) are
    majorized by z^(a0-N-1) e^(-z) exp(log(z) eps) / (1 - eps/z), and the
    polynomial prefactor by prod (|a0-i| + eps).
*/
slong
_arb_hypgeom_gamma_upper_fmpq_inf_choose_N_jet(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol, slong len)
{
    slong N, ab;
    fmpz_t az;
    fmpq_t c, a1;
    mag_t E, zlow, zinv, L, cm, t;
    mag_ptr V, W, W2;

    if (len <= 1)
        return _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(err, a, z, abs_tol);

    fmpz_init(az);
    fmpq_init(c);
    fmpq_init(a1);
    mag_init(E);
    mag_init(zlow);
    mag_init(zinv);
    mag_init(L);
    mag_init(cm);
    mag_init(t);
    V = _mag_vec_init(len);
    W = _mag_vec_init(len);
    W2 = _mag_vec_init(len);

    fmpz_cdiv_q(az, fmpq_numref(a), fmpq_denref(a));
    arb_get_mag_lower(zlow, z);

    if (!fmpz_fits_si(az) || mag_cmp_2exp_si(zlow, 0) < 0)
    {
        mag_inf(err);
        N = -1;
        goto cleanup;
    }

    ab = fmpz_get_si(az);

    /* prefactor z^(a-1) e^(-z) */
    fmpq_sub_ui(a1, a, 1);
    _mag_pow_fmpq_expinv(E, z, a1, 1);

    mag_inv(zinv, zlow);
    _mag_abs_log_arb(L, z);

    if (mag_is_inf(E) || mag_is_inf(L))
    {
        mag_inf(err);
        N = -1;
        goto cleanup;
    }

    mag_one(V);

    for (N = 1; ; N++)
    {
        /* c = |a - N| */
        fmpq_sub_ui(c, a, N);
        _mag_set_fmpq(cm, c);

        /* stop if terms are increasing */
        if (mag_cmp(cm, zlow) > 0)
        {
            mag_inf(err);
            N = -1;
            break;
        }

        _mag_vec_mul_linear(V, cm, len);
        mag_mul(E, E, zinv);

        if (N >= ab - 1)
        {
            /* cheap test using the constant term */
            mag_mul(t, E, V);
            if (mag_cmp(t, abs_tol) >= 0)
                continue;

            _mag_vec_mul_exp(W, V, L, len);
            /* multiply by 1/(1 - eps/z) */
            mag_set(W2, W);
            for (slong j = 1; j < len; j++)
            {
                mag_mul(t, W2 + j - 1, zinv);
                mag_add(W2 + j, W + j, t);
            }

            _mag_vec_max(t, W2, len);
            mag_mul(t, t, E);

            if (mag_cmp(t, abs_tol) < 0)
            {
                mag_swap(err, t);
                break;
            }
        }
    }

cleanup:
    fmpz_clear(az);
    fmpq_clear(c);
    fmpq_clear(a1);
    mag_clear(E);
    mag_clear(zlow);
    mag_clear(zinv);
    mag_clear(L);
    mag_clear(cm);
    mag_clear(t);
    _mag_vec_clear(V, len);
    _mag_vec_clear(W, len);
    _mag_vec_clear(W2, len);

    return N;
}

static void
upper_bsplit_jet(arb_poly_t M, arb_poly_t S, arb_t Q, const fmpz_t ap, const fmpz_t aq, const arb_t z, slong na, slong nb, int cont, slong len, slong prec)
{
    if (nb - na == 0)
    {
        arb_poly_zero(M);
        arb_poly_zero(S);
        arb_one(Q);
    }
    else if (nb - na == 1)
    {
        fmpz_t t, u;
        fmpz_init(t);
        fmpz_init(u);

        /* M = -(a - (na + 1)) aq */
        fmpz_mul_ui(t, aq, na + 1);
        fmpz_sub(t, t, ap);
        fmpz_neg(u, aq);
        _arb_poly_set_linear_fmpz(M, t, u);
        arb_poly_truncate(M, len);

        arb_mul_fmpz(Q, z, aq, prec);
        arb_neg(Q, Q);
        arb_poly_set_arb(S, Q);

        fmpz_clear(t);
        fmpz_clear(u);
    }
    else
    {
        slong m;
        arb_poly_t M2, S2;
        arb_t Q2;

        m = na + (nb - na) / 2;

        arb_poly_init(M2);
        arb_poly_init(S2);
        arb_init(Q2);

        upper_bsplit_jet(M, S, Q, ap, aq, z, na, m, 1, len, prec);
        upper_bsplit_jet(M2, S2, Q2, ap, aq, z, m, nb, cont, len, prec);

        arb_poly_scalar_mul(S, S, Q2, prec);
        arb_poly_mullow_bsplit(S2, M, S2, len, prec);
        arb_poly_add(S, S, S2, prec);

        if (cont)
            arb_poly_mullow_bsplit(M, M, M2, len, prec);

        arb_mul(Q, Q, Q2, prec);

        arb_poly_clear(M2);
        arb_poly_clear(S2);
        arb_clear(Q2);
    }
}

void
_arb_hypgeom_gamma_upper_fmpq_inf_bsplit_jet(arb_ptr res, const fmpq_t a, const arb_t z, const arb_t logz, slong N, slong len, slong prec)
{
    arb_poly_t M, S;
    arb_t Q;
    arb_ptr T, U;
    fmpq_t a1;

    if (len == 1)
    {
        _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(res, a, z, N, prec);
        return;
    }

    arb_poly_init(M);
    arb_poly_init(S);
    arb_init(Q);
    T = _arb_vec_init(len);
    U = _arb_vec_init(len);
    fmpq_init(a1);

    N = FLINT_MAX(N, 0);

    upper_bsplit_jet(M, S, Q, fmpq_numref(a), fmpq_denref(a), z, 0, N, 0, len, prec);
    _arb_poly_get_vec(T, S, len);
    _arb_vec_scalar_div(T, T, len, Q, prec);

    /* z^(a-1+eps) e^(-z) */
    fmpq_sub_ui(a1, a, 1);
    _arb_pow_fmpq_jet(U, z, logz, a1, len, prec);
    arb_neg(Q, z);
    arb_exp(Q, Q, prec);
    _arb_vec_scalar_mul(U, U, len, Q, prec);

    _arb_poly_mullow(res, T, len, U, len, len, prec);

    arb_poly_clear(M);
    arb_poly_clear(S);
    arb_clear(Q);
    _arb_vec_clear(T, len);
    _arb_vec_clear(U, len);
    fmpq_clear(a1);
}

/* Lower series at zero ***************************************************/

/*
    gamma(a,z) = z^a e^(-z) sum_{k>=0} z^k / (a (a+1) ... (a+k)).

    The coefficients of 1/(a0+i+eps) are majorized by those of
    1/(|a0+i| - eps), and those of z^eps by exp(|log(z)| eps). Once
    a0 + N + 1 > z, the tail from index N is majorized by the term of index N
    times (c - eps)/(c - z - eps) = 1 + z/(c - z - eps), c = a0 + N + 1.
*/
slong
_arb_hypgeom_gamma_lower_fmpq_0_choose_N_jet(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol, slong len)
{
    slong N;
    fmpq_t c;
    mag_t E, zup, L, cm, d, t;
    mag_ptr V, W, W2;

    if (len <= 1)
        return _arb_hypgeom_gamma_lower_fmpq_0_choose_N(err, a, z, abs_tol);

    if (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0)
    {
        mag_inf(err);
        return -1;
    }

    fmpq_init(c);
    mag_init(E);
    mag_init(zup);
    mag_init(L);
    mag_init(cm);
    mag_init(d);
    mag_init(t);
    V = _mag_vec_init(len);
    W = _mag_vec_init(len);
    W2 = _mag_vec_init(len);

    arb_get_mag(zup, z);

    /* prefactor z^a e^(-z) */
    _mag_pow_fmpq_expinv(E, z, a, 1);
    _mag_abs_log_arb(L, z);

    if (mag_is_inf(E) || mag_is_inf(L))
    {
        mag_inf(err);
        N = -1;
        goto cleanup;
    }

    /* V = 1/(a+eps) */
    mag_one(V);
    _mag_set_fmpq_lower(cm, a);
    _mag_vec_div_linear(V, cm, len);

    for (N = 1; ; N++)
    {
        /* term of index N */
        fmpq_add_ui(c, a, N);
        _mag_set_fmpq_lower(cm, c);
        _mag_vec_div_linear(V, cm, len);
        mag_mul(E, E, zup);

        /* require a + N + 1 > z */
        fmpq_add_ui(c, a, N + 1);
        if (fmpz_sgn(fmpq_numref(c)) <= 0)
            continue;
        _mag_set_fmpq_lower(cm, c);
        mag_sub_lower(d, cm, zup);
        if (mag_is_zero(d))
            continue;

        /* cheap test using the constant term */
        mag_div(t, cm, d);
        mag_mul(t, t, V);
        mag_mul(t, t, E);
        if (mag_cmp(t, abs_tol) >= 0)
            continue;

        /* W2 = V * (1 + z/(d - eps)) */
        for (slong j = 0; j < len; j++)
            mag_set(W + j, V + j);
        _mag_vec_div_linear(W, d, len);
        for (slong j = 0; j < len; j++)
        {
            mag_mul(W + j, W + j, zup);
            mag_add(W + j, W + j, V + j);
        }

        _mag_vec_mul_exp(W2, W, L, len);
        _mag_vec_max(t, W2, len);
        mag_mul(t, t, E);

        if (mag_cmp(t, abs_tol) < 0)
        {
            mag_swap(err, t);
            break;
        }
    }

cleanup:
    fmpq_clear(c);
    mag_clear(E);
    mag_clear(zup);
    mag_clear(L);
    mag_clear(cm);
    mag_clear(d);
    mag_clear(t);
    _mag_vec_clear(V, len);
    _mag_vec_clear(W, len);
    _mag_vec_clear(W2, len);

    return N;
}

static void
lower_bsplit_jet(arb_t M, arb_poly_t S, arb_poly_t Q, const fmpz_t ap, const fmpz_t aq, const arb_t z, slong na, slong nb, int cont, slong len, slong prec)
{
    if (nb - na == 0)
    {
        arb_zero(M);
        arb_poly_zero(S);
        arb_poly_one(Q);
    }
    else if (nb - na == 1)
    {
        fmpz_t t;
        fmpz_init_set(t, ap);
        fmpz_addmul_ui(t, aq, na + 1);
        _arb_poly_set_linear_fmpz(S, t, aq);
        arb_poly_truncate(S, len);
        arb_poly_set(Q, S);
        arb_mul_fmpz(M, z, aq, prec);
        fmpz_clear(t);
    }
    else
    {
        slong m;
        arb_t M2;
        arb_poly_t S2, Q2;

        m = na + (nb - na) / 2;

        arb_init(M2);
        arb_poly_init(S2);
        arb_poly_init(Q2);

        lower_bsplit_jet(M, S, Q, ap, aq, z, na, m, 1, len, prec);
        lower_bsplit_jet(M2, S2, Q2, ap, aq, z, m, nb, cont, len, prec);

        arb_poly_mullow_bsplit(S, S, Q2, len, prec);
        arb_poly_scalar_mul(S2, S2, M, prec);
        arb_poly_add(S, S, S2, prec);

        if (cont)
            arb_mul(M, M, M2, prec);

        arb_poly_mullow_bsplit(Q, Q, Q2, len, prec);

        arb_clear(M2);
        arb_poly_clear(S2);
        arb_poly_clear(Q2);
    }
}

void
_arb_hypgeom_gamma_lower_fmpq_0_bsplit_jet(arb_ptr res, const fmpq_t a, const arb_t z, const arb_t logz, slong N, slong len, slong prec)
{
    arb_t M;
    arb_poly_t S, Q;
    arb_ptr T, U;
    slong j;

    if (len == 1)
    {
        _arb_hypgeom_gamma_lower_fmpq_0_bsplit(res, a, z, N, prec);
        return;
    }

    arb_init(M);
    arb_poly_init(S);
    arb_poly_init(Q);
    T = _arb_vec_init(len);
    U = _arb_vec_init(len);

    N = FLINT_MAX(N, 0);

    lower_bsplit_jet(M, S, Q, fmpq_numref(a), fmpq_denref(a), z, 0, N, 0, len, prec);

    if (S->length == 0)
    {
        _arb_vec_zero(res, len);
    }
    else
    {
        _arb_poly_div_series(T, S->coeffs, S->length, Q->coeffs, Q->length, len, prec);

        /* 1/(a+eps) */
        arb_set_fmpz(U, fmpq_denref(a));
        arb_div_fmpz(U, U, fmpq_numref(a), prec);
        for (j = 1; j < len; j++)
        {
            arb_mul(U + j, U + j - 1, U, prec);
            arb_neg(U + j, U + j);
        }
        _arb_poly_mullow(res, T, len, U, len, len, prec);

        /* z^(a+eps) e^(-z) */
        _arb_pow_fmpq_jet(U, z, logz, a, len, prec);
        arb_neg(M, z);
        arb_exp(M, M, prec);
        _arb_vec_scalar_mul(U, U, len, M, prec);
        _arb_poly_mullow(T, res, len, U, len, len, prec);
        _arb_vec_set(res, T, len);
    }

    arb_clear(M);
    arb_poly_clear(S);
    arb_poly_clear(Q);
    _arb_vec_clear(T, len);
    _arb_vec_clear(U, len);
}

/* Gamma function jet *****************************************************/

/*
    Gamma(a+eps) for rational a which is not a nonpositive integer, via
    gamma(a+eps,Z) + Gamma(a+eps,Z) with an integer Z, as in
    arb_hypgeom_gamma_fmpq: the lower series at Z needs about
    (1-alpha)/W((1-alpha)/(alpha e)) Z/alpha terms and the asymptotic series
    about (alpha-1)/W_{-1}((alpha-1)/(alpha e)) Z/alpha terms, at reduced
    precision since Gamma(a,Z) ~ Z^(a-1) e^(-Z) is small, where
    Z = alpha log(2) p. The asymptotic series can reach the tolerance
    2^-p |Gamma(a)| only if Z is at least about log(2) p / 2.

    If the asymptotic series does not converge to the tolerance (small
    precision, or large a), Z is increased; for large enough Z, the
    upper incomplete gamma function is negligible: its coefficients are
    int_Z^inf t^(a-1) log(t)^j/j! e^(-t) dt <= int_Z^inf t^a e^(-t) dt
    = Gamma(a+1, Z) <= 2 Z^a e^(-Z) for Z >= max(1, 2a).
*/

#ifndef GAMMA_JET_ALPHA
#define GAMMA_JET_ALPHA 0.55
#endif

void
_arb_hypgeom_gamma_fmpq_jet(arb_ptr res, const fmpq_t a, slong len, slong prec)
{
    double ad, lg, target, umag;
    slong Z, Zmin, N1, N2, wp, wp2;
    arb_t z, logz;
    arb_ptr U;
    mag_t tol, err, err2, T;
    fmpz_t ac;

    if (len <= 0)
        return;

    if (fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0)
    {
        _arb_vec_indeterminate(res, len);
        return;
    }

    if (len == 1)
    {
        arb_gamma_fmpq(res, a, prec);
        return;
    }

    ad = fmpq_get_d(a);

    if (!(fabs(ad) < 1e15))
    {
        _arb_vec_indeterminate(res, len);
        return;
    }

    arb_init(z);
    arb_init(logz);
    mag_init(tol);
    mag_init(err);
    mag_init(err2);
    mag_init(T);
    fmpz_init(ac);
    U = _arb_vec_init(len);

    /* log2 |Gamma(a)| (approximate) */
    lg = lgamma(ad) * 1.4426950408889634;
    if (!(lg > -1e15 && lg < 1e15))
        lg = 0;

    /* absolute tolerance for each of the two parts */
    target = floor(lg) - prec - 6;
    mag_one(tol);
    mag_mul_2exp_si(tol, tol, (slong) target);

    fmpz_cdiv_q(ac, fmpq_numref(a), fmpq_denref(a));

    Zmin = FLINT_MAX(2, 2 * (slong) ceil(fabs(ad)) + 1);
    Z = FLINT_MAX(Zmin, (slong) (GAMMA_JET_ALPHA * 0.693147 * (-target)) + 1);

    /* find Z for which the asymptotic series converges, or for which
       the upper incomplete gamma function is negligible */
    N2 = -1;
    while (1)
    {
        /* tail bound 2 Z^ceil(a) e^(-Z) */
        mag_set_ui_lower(T, Z);
        mag_expinv(T, T);
        mag_mul_2exp_si(T, T, 1);
        mag_set_ui(err, Z);
        mag_pow_fmpz(err, err, ac);
        mag_mul(T, T, err);

        if (mag_cmp(T, tol) <= 0)
            break;

        arb_set_ui(z, Z);
        N2 = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N_jet(err2, a, z, tol, len);

        if (N2 >= 0)
            break;

        Z += Z / 16 + 1;
    }

    arb_set_ui(z, Z);
    wp = prec + 10 + 2 * FLINT_BIT_COUNT(len);
    arb_log_ui(logz, Z, wp);

    N1 = _arb_hypgeom_gamma_lower_fmpq_0_choose_N_jet(err, a, z, tol, len);

    if (N1 < 0)
    {
        _arb_vec_indeterminate(res, len);
    }
    else
    {
        _arb_hypgeom_gamma_lower_fmpq_0_bsplit_jet(res, a, z, logz, N1, len, wp);

        if (N2 >= 0)
        {
            /* log2 of a bound for the coefficients of Gamma(a+eps, Z):
               Z^(a-1) e^(-Z) times exp(log(Z) eps) / (1 - eps/Z) */
            umag = (ad - 1) * log(Z) * 1.4426950408889634 - Z * 1.4426950408889634
                + len * 1.4426950408889634 * log(1 + log(Z)) + 2;
            wp2 = (slong) (umag - target) + 10 + 2 * FLINT_BIT_COUNT(len);
            wp2 = FLINT_MIN(wp2, wp);
            wp2 = FLINT_MAX(wp2, 30);

            _arb_hypgeom_gamma_upper_fmpq_inf_bsplit_jet(U, a, z, logz, N2, len, wp2);
            _arb_vec_add(res, res, U, len, wp);
            mag_add(err, err, err2);
        }
        else
        {
            mag_add(err, err, T);
        }

        _arb_vec_add_error_mag(res, len, err);
        _arb_vec_set_round(res, res, len, prec);
    }

    arb_clear(z);
    arb_clear(logz);
    mag_clear(tol);
    mag_clear(err);
    mag_clear(err2);
    mag_clear(T);
    fmpz_clear(ac);
    _arb_vec_clear(U, len);
}

/* Singular series ********************************************************/

/*
    Gamma(-n+eps, z) = [Gamma(-n+eps) - z^(-n+eps) (-z)^n / (n! eps)]
                       - z^(-n+eps) sum_{k != n} (-z)^k / (k! (k-n+eps))

    where the bracket equals

        (-1)^n / n! [Gamma(1+eps) / prod_{i=1}^n (1 - eps/i) - z^eps] / eps.

    The coefficients of 1/(k-n+eps) for k != n are bounded by 1 in absolute
    value, so the tail of the sum is majorized termwise by z^k/k! and the
    factor z^eps by exp(|log(z)| eps).
*/
slong
_arb_hypgeom_gamma_upper_singular_si_choose_N_jet(mag_t err, slong n, const arb_t z, const mag_t abs_tol, slong len)
{
    slong k;
    mag_t t, u, zm, F;

    mag_init(t);
    mag_init(u);
    mag_init(zm);
    mag_init(F);

    arb_get_mag(zm, z);

    /* F = sum_{j<len} L^j/j! bounds the coefficients of z^eps */
    _mag_abs_log_arb(u, z);
    _mag_exp_jet_sum(F, u, len);

    arb_get_mag_lower(t, z);
    mag_inv(t, t);
    mag_pow_ui(t, t, n);
    mag_mul(t, t, F);

    for (k = 1; ; k++)
    {
        mag_mul(t, t, zm);
        mag_div_ui(t, t, k);

        if (mag_cmp(t, abs_tol) < 0)
        {
            /* tail from index k: t * sum_m (z/(k+1))^m */
            mag_div_ui(u, zm, k + 1);
            mag_geom_series(u, u, 0);
            mag_mul(u, t, u);

            if (mag_cmp(u, abs_tol) < 0)
            {
                mag_swap(err, u);
                break;
            }
        }
    }

    mag_clear(t);
    mag_clear(u);
    mag_clear(zm);
    mag_clear(F);

    return k;
}

/*
    Binary splitting for sum_{k != n} h_k / (k - n + eps), h_k = (-z)^k / k!.
    The sum over [a,b) is T / (Q D) where Q = prod (k+1) and D = prod (k-n+eps)
    (omitting k = n), and P = (-z)^(b-a); only T and D are power series.
*/
static void
singular_bsplit_jet(arb_t P, arb_t Q, arb_poly_t D, arb_poly_t T, slong n, const arb_t z, slong na, slong nb, int cont, slong len, slong prec)
{
    if (nb - na == 0)
    {
        arb_one(P);
        arb_one(Q);
        arb_poly_one(D);
        arb_poly_zero(T);
    }
    else if (nb - na == 1)
    {
        slong k = na;

        arb_neg(P, z);
        arb_set_si(Q, k + 1);

        if (k == n)
        {
            arb_poly_one(D);
            arb_poly_zero(T);
        }
        else
        {
            arb_poly_fit_length(D, 2);
            arb_set_si(D->coeffs, k - n);
            arb_one(D->coeffs + 1);
            _arb_poly_set_length(D, 2);
            arb_poly_truncate(D, len);

            arb_poly_set_si(T, k + 1);
        }
    }
    else
    {
        slong m;
        arb_t P2, Q2;
        arb_poly_t D2, T2;

        m = na + (nb - na) / 2;

        arb_init(P2);
        arb_init(Q2);
        arb_poly_init(D2);
        arb_poly_init(T2);

        singular_bsplit_jet(P, Q, D, T, n, z, na, m, 1, len, prec);
        singular_bsplit_jet(P2, Q2, D2, T2, n, z, m, nb, cont, len, prec);

        /* T = Q2 T1 D2 + P1 T2 D1 */
        arb_poly_mullow_bsplit(T, T, D2, len, prec);
        arb_poly_scalar_mul(T, T, Q2, prec);
        arb_poly_mullow_bsplit(T2, T2, D, len, prec);
        arb_poly_scalar_mul(T2, T2, P, prec);
        arb_poly_add(T, T, T2, prec);

        arb_poly_mullow_bsplit(D, D, D2, len, prec);
        arb_mul(Q, Q, Q2, prec);

        if (cont)
            arb_mul(P, P, P2, prec);

        arb_clear(P2);
        arb_clear(Q2);
        arb_poly_clear(D2);
        arb_poly_clear(T2);
    }
}

void
_arb_hypgeom_gamma_upper_singular_si_bsplit_jet(arb_ptr res, slong n, const arb_t z, const arb_t logz, slong N, arb_srcptr gamma1, slong len, slong prec)
{
    arb_poly_t D, T;
    arb_t P, Q;
    arb_ptr G, U, V;
    arb_t t;
    fmpq_t c;
    slong i, j;

    if (len == 1)
    {
        _arb_hypgeom_gamma_upper_singular_si_bsplit(res, n, z, N, prec);
        return;
    }

    arb_poly_init(D);
    arb_poly_init(T);
    arb_init(P);
    arb_init(Q);
    G = _arb_vec_init(len + 1);
    U = _arb_vec_init(len + 1);
    V = _arb_vec_init(len + 1);
    arb_init(t);
    fmpq_init(c);

    N = FLINT_MAX(N, 0);

    /* G = Gamma(1+eps) / prod_{i=1}^n (1 - eps/i) - z^eps */
    if (gamma1 == NULL)
    {
        fmpq_one(c);
        _arb_hypgeom_gamma_fmpq_jet(G, c, len + 1, prec);
    }
    else
    {
        _arb_vec_set(G, gamma1, len + 1);
    }

    /* G <- G * i/(i - eps): G_j <- G_j + G_{j-1} / i, ascending in j */
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= len; j++)
        {
            arb_div_ui(t, G + j - 1, i, prec);
            arb_add(G + j, G + j, t, prec);
        }
    }

    fmpq_zero(c);
    _arb_pow_fmpq_jet(U, z, logz, c, len + 1, prec);
    _arb_vec_sub(G, G, U, len + 1, prec);

    /* bracket / eps, times (-1)^n / n! */
    arb_fac_ui(t, n, prec);
    if (n & 1)
        arb_neg(t, t);
    _arb_vec_scalar_div(V, G + 1, len, t, prec);

    /* the sum */
    singular_bsplit_jet(P, Q, D, T, n, z, 0, N, 0, len, prec);

    if (T->length != 0)
    {
        _arb_poly_div_series(G, T->coeffs, T->length, D->coeffs, D->length, len, prec);
        _arb_vec_scalar_div(G, G, len, Q, prec);

        /* times z^(-n+eps) */
        fmpq_set_si(c, -n, 1);
        _arb_pow_fmpq_jet(U, z, logz, c, len, prec);
        _arb_poly_mullow(res, G, len, U, len, len, prec);
        _arb_vec_sub(res, V, res, len, prec);
    }
    else
    {
        _arb_vec_set(res, V, len);
    }

    arb_poly_clear(D);
    arb_poly_clear(T);
    arb_clear(P);
    arb_clear(Q);
    _arb_vec_clear(G, len + 1);
    _arb_vec_clear(U, len + 1);
    _arb_vec_clear(V, len + 1);
    arb_clear(t);
    fmpq_clear(c);
}
