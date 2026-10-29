/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Jet version of _arb_gamma_upper_fmpq_step_bsplit: one bit-burst step
    for Gamma(a+eps, z) where a is rational and eps is the jet variable.

    Gamma(a, z0 + x) = Gamma(a, z0) - int_{z0}^{z0+x} f(t) dt with
    f(t) = t^(a-1) e^(-t). The Taylor coefficients c_k of f at z0 satisfy a
    three-term recurrence whose coefficients are linear in a; we run the
    same 3x3 matrix binary splitting as in the scalar case, with matrix
    entries that are power series in eps.

    Error bound: by Cauchy's estimate on the circle |t - z0| = R < z0, the
    coefficient of eps^j in c_k is bounded by M Lambda^j / j! / R^k where
    M bounds |t^(a-1) e^(-t)| and Lambda bounds |log(t)| on the circle.
    The scalar truncation bound is therefore multiplied by
    max_{j<len} Lambda^j / j!, and the step length R is chosen with this
    factor taken into account.
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

/* bound for the coefficients of t^(a-1+eps) e^(-t) on |t - z| <= 2^Rexp */
static void
taylor_M_jet(mag_t M, const arb_t a, const arb_t z, slong Rexp, slong len)
{
    arb_t t, u;
    mag_t lo, hi, L, Lt, R;

    arb_init(t);
    arb_init(u);
    mag_init(lo);
    mag_init(hi);
    mag_init(L);
    mag_init(Lt);
    mag_init(R);

    arb_sub_ui(u, a, 1, 53);
    arb_sgn(t, u);
    arb_mul_2exp_si(t, t, Rexp);
    arb_add(t, z, t, 53);
    arb_pow(t, t, u, 53);

    arb_one(u);
    arb_mul_2exp_si(u, u, Rexp);
    arb_sub(u, u, z, 53);
    arb_exp(u, u, 53);

    arb_mul(t, t, u, 53);

    arb_get_mag(M, t);

    if (len > 1 && mag_is_finite(M))
    {
        slong j;

        /* |log(t)| <= max(|log(z-R)|, |log(z+R)|) + |arg(t)|,
           |arg(t)| <= asin(R/z) <= (pi/2) R/z */
        mag_one(R);
        mag_mul_2exp_si(R, R, Rexp);

        arb_get_mag_lower(lo, z);
        arb_get_mag(hi, z);

        mag_div(Lt, R, lo);
        mag_mul_ui(Lt, Lt, 7);
        mag_mul_2exp_si(Lt, Lt, -2);   /* 7/4 > pi/2 */

        mag_sub_lower(lo, lo, R);
        mag_add(hi, hi, R);

        if (mag_is_zero(lo))
        {
            mag_inf(M);
        }
        else
        {
            mag_log(L, hi);
            mag_neg_log(lo, lo);
            mag_max(L, L, lo);
            mag_add(L, L, Lt);

            /* M *= max_{j<len} L^j/j! */
            mag_one(lo);
            mag_one(hi);
            for (j = 1; j < len; j++)
            {
                mag_mul(lo, lo, L);
                mag_div_ui(lo, lo, j);
                mag_max(hi, hi, lo);
            }

            mag_mul(M, M, hi);
        }
    }

    arb_clear(t);
    arb_clear(u);
    mag_clear(lo);
    mag_clear(hi);
    mag_clear(L);
    mag_clear(Lt);
    mag_clear(R);
}

/* choose N such that M * C^N / (1 - C) <= tol */
static slong
mag_geom_choose_N(const mag_t M, const mag_t C, const mag_t tol)
{
    mag_t t, u;
    slong N;

    if (mag_is_finite(M) && mag_is_zero(C))
        return 1;

    if (!mag_is_finite(M))
        return WORD_MAX;

    /* N = log(M / ((1 - C) tol)) / log(1/C) */
    mag_init(t);
    mag_init(u);

    mag_one(t);
    mag_sub_lower(t, t, C);
    mag_mul_lower(t, t, tol);
    mag_div(t, M, t);
    mag_log(t, t);

    mag_inv_lower(u, C);
    mag_log_lower(u, u);
    mag_div(t, t, u);

    if (mag_cmp_2exp_si(t, 60) > 0)
        N = WORD_MAX;
    else
    {
        N = ceil(mag_get_d(t));
        N = FLINT_MAX(N, 1);
    }

    mag_clear(t);
    mag_clear(u);

    return N;
}

static void
taylor_bound_jet(mag_t err, const arb_t a, const arb_t z, const mag_t x, slong Rexp, slong N, slong len)
{
    mag_t C, M;

    mag_init(C);
    mag_init(M);

    /* C = x / R */
    mag_mul_2exp_si(C, x, -Rexp);

    /* M R C^n / (1 - C) / N */
    mag_geom_series(err, C, N);

    if (!mag_is_inf(err))
    {
        taylor_M_jet(M, a, z, Rexp, len);
        mag_mul(err, err, M);

        mag_mul_2exp_si(err, err, Rexp);
        mag_div_ui(err, err, N);
    }

    mag_clear(C);
    mag_clear(M);
}

static slong
taylor_N_jet(const arb_t a, const arb_t z, const mag_t x, slong Rexp, const mag_t abs_tol, slong len)
{
    mag_t C, M;
    slong N;

    mag_init(C);
    mag_init(M);

    /* C = x / R */
    mag_mul_2exp_si(C, x, -Rexp);

    if (mag_cmp_2exp_si(C, 0) < 0)
    {
        taylor_M_jet(M, a, z, Rexp, len);
        mag_mul_2exp_si(M, M, Rexp);

        N = mag_geom_choose_N(M, C, abs_tol);
    }
    else
    {
        N = WORD_MAX;
    }

    mag_clear(C);
    mag_clear(M);

    return N;
}

static int
arb_hypgeom_gamma_upper_taylor_choose_jet(slong * res_N, mag_t err, const arb_t a, const arb_t z, const mag_t x, const mag_t abs_tol, slong len)
{
    slong N, New;
    mag_t zlow;
    slong Rexp;

    mag_init(zlow);
    arb_get_mag_lower(zlow, z);

    Rexp = 0;
    while (mag_cmp_2exp_si(zlow, Rexp + 1) < 0)
        Rexp--;

    N = taylor_N_jet(a, z, x, Rexp, abs_tol, len);

    while (N > 1 && mag_cmp_2exp_si(x, Rexp - 1) < 0)
    {
        New = taylor_N_jet(a, z, x, Rexp - 1, abs_tol, len);

        if (New <= N)
        {
            Rexp = Rexp - 1;
            N = New;
        }
        else
        {
            break;
        }
    }

    if (Rexp == 0)
    {
        while (N > 1 && mag_cmp_2exp_si(zlow, Rexp + 1) > 0)
        {
            New = taylor_N_jet(a, z, x, Rexp + 1, abs_tol, len);

            if (New <= N)
            {
                Rexp = Rexp + 1;
                N = New;
            }
            else
            {
                break;
            }
        }
    }

    mag_clear(zlow);

    if (N == WORD_MAX)
    {
        mag_inf(err);
        return 0;
    }

    *res_N = N;
    taylor_bound_jet(err, a, z, x, Rexp, N, len);

    return mag_cmp(err, abs_tol) <= 0;
}

/*
    The matrix for one step has the block structure

        [A B 0]
        [C D 0]
        [E F G]

    where G = Q (the scalar denominator). Only D depends on a.
*/
static void
gamma_upper_taylor_bsplit_jet(arb_poly_t A, arb_poly_t B, arb_poly_t C, arb_poly_t D,
    arb_poly_t E, arb_poly_t F, arb_t Q,
    const fmpz_t ap, const fmpz_t aq, const arb_t z0, const arb_t x,
    slong a, slong b, int cont, slong len, slong prec)
{
    if (b - a == 0)
    {
        arb_poly_one(A);
        arb_poly_zero(B);
        arb_poly_zero(C);
        arb_poly_one(D);
        arb_poly_zero(E);
        arb_poly_zero(F);
        arb_one(Q);
    }
    else if (b - a == 1)
    {
        slong n;
        fmpz_t t;
        arb_t u;

        n = a;
        fmpz_init(t);
        arb_init(u);

        /* Q = -z0*(n+1)*(n+2)*aq */
        fmpz_mul2_uiui(t, aq, n + 1, n + 2);
        arb_mul_fmpz(Q, z0, t, prec);
        arb_neg(Q, Q);

        arb_poly_zero(A);

        /* x Q */
        arb_mul(u, Q, x, prec);
        arb_poly_set_arb(B, u);

        /* aq n x */
        fmpz_mul_ui(t, aq, n);
        arb_mul_fmpz(u, x, t, prec);
        arb_poly_set_arb(C, u);

        /* x*(-(ap + aq*eps) + aq*(n + z0 + 1))*(n + 1) */
        arb_poly_fit_length(D, 2);
        arb_add_ui(D->coeffs, z0, n + 1, prec);
        arb_mul_fmpz(D->coeffs, D->coeffs, aq, prec);
        arb_sub_fmpz(D->coeffs, D->coeffs, ap, prec);
        arb_mul_ui(D->coeffs, D->coeffs, n + 1, prec);
        arb_mul(D->coeffs, D->coeffs, x, prec);
        fmpz_mul_ui(t, aq, n + 1);
        arb_mul_fmpz(D->coeffs + 1, x, t, prec);
        arb_neg(D->coeffs + 1, D->coeffs + 1);
        _arb_poly_set_length(D, 2);
        _arb_poly_normalise(D);
        arb_poly_truncate(D, len);

        arb_poly_set_arb(E, Q);
        arb_poly_zero(F);

        fmpz_clear(t);
        arb_clear(u);
    }
    else
    {
        arb_poly_t A2, B2, C2, D2, E2, F2, T, U;
        arb_t Q2;
        slong m;

        arb_poly_init(A2);
        arb_poly_init(B2);
        arb_poly_init(C2);
        arb_poly_init(D2);
        arb_poly_init(E2);
        arb_poly_init(F2);
        arb_poly_init(T);
        arb_poly_init(U);
        arb_init(Q2);

        m = a + (b - a) / 2;

        gamma_upper_taylor_bsplit_jet(A, B, C, D, E, F, Q, ap, aq, z0, x, a, m, 1, len, prec);
        gamma_upper_taylor_bsplit_jet(A2, B2, C2, D2, E2, F2, Q2, ap, aq, z0, x, m, b, cont, len, prec);

        /* E = E2 A + F2 C + Q2 E */
        arb_poly_scalar_mul(E, E, Q2, prec);
        arb_poly_mullow_bsplit(T, E2, A, len, prec);
        arb_poly_add(E, E, T, prec);
        arb_poly_mullow_bsplit(T, F2, C, len, prec);
        arb_poly_add(E, E, T, prec);

        /* F = E2 B + F2 D + Q2 F */
        arb_poly_scalar_mul(F, F, Q2, prec);
        arb_poly_mullow_bsplit(T, E2, B, len, prec);
        arb_poly_add(F, F, T, prec);
        arb_poly_mullow_bsplit(T, F2, D, len, prec);
        arb_poly_add(F, F, T, prec);

        if (cont)
        {
            /* [A B; C D] = [A2 B2; C2 D2] [A B; C D] */
            arb_poly_mullow_bsplit(T, A2, A, len, prec);
            arb_poly_mullow_bsplit(U, B2, C, len, prec);
            arb_poly_add(T, T, U, prec);           /* new A */

            arb_poly_mullow_bsplit(U, C2, A, len, prec);
            arb_poly_mullow_bsplit(A, D2, C, len, prec);
            arb_poly_add(C, U, A, prec);           /* new C */

            arb_poly_mullow_bsplit(U, A2, B, len, prec);
            arb_poly_mullow_bsplit(A, B2, D, len, prec);
            arb_poly_add(U, U, A, prec);           /* new B */

            arb_poly_mullow_bsplit(A, C2, B, len, prec);
            arb_poly_mullow_bsplit(B, D2, D, len, prec);
            arb_poly_add(D, A, B, prec);           /* new D */

            arb_poly_swap(A, T);
            arb_poly_swap(B, U);
        }

        arb_mul(Q, Q, Q2, prec);

        arb_poly_clear(A2);
        arb_poly_clear(B2);
        arb_poly_clear(C2);
        arb_poly_clear(D2);
        arb_poly_clear(E2);
        arb_poly_clear(F2);
        arb_poly_clear(T);
        arb_poly_clear(U);
        arb_clear(Q2);
    }
}

void
_arb_gamma_upper_fmpq_step_bsplit_jet(arb_ptr Gz1, const fmpq_t a, const arb_t z0, const arb_t logz0, const arb_t z1, arb_srcptr Gz0, const arb_t expmz0, const mag_t abs_tol, slong len, slong prec)
{
    arb_t x, Q, a_real, t;
    arb_poly_t A, B, C, D, E, F;
    arb_ptr T, U;
    mag_t xmag, err;
    slong N, j;
    fmpq_t a1;

    if (len == 1)
    {
        _arb_gamma_upper_fmpq_step_bsplit(Gz1, a, z0, z1, Gz0, expmz0, abs_tol, prec);
        return;
    }

    if (arb_is_zero(z0))
    {
        T = _arb_vec_init(len);
        mag_init(err);
        N = _arb_hypgeom_gamma_lower_fmpq_0_choose_N_jet(err, a, z1, abs_tol, len);
        if (N < 0)
        {
            _arb_vec_indeterminate(Gz1, len);
        }
        else
        {
            _arb_hypgeom_gamma_lower_fmpq_0_bsplit_jet(T, a, z1, NULL, N, len, prec);
            for (j = 0; j < len; j++)
                arb_add_error_mag(T + j, err);
            _arb_vec_sub(Gz1, Gz0, T, len, prec);
        }
        _arb_vec_clear(T, len);
        mag_clear(err);
        return;
    }

    mag_init(xmag);
    mag_init(err);
    arb_init(x);
    arb_init(Q);
    arb_init(a_real);
    arb_init(t);
    fmpq_init(a1);
    arb_poly_init(A);
    arb_poly_init(B);
    arb_poly_init(C);
    arb_poly_init(D);
    arb_poly_init(E);
    arb_poly_init(F);
    T = _arb_vec_init(len);
    U = _arb_vec_init(len);

    arb_sub(x, z1, z0, prec);
    arb_get_mag(xmag, x);

    arb_set_fmpq(a_real, a, 53);

    if (!arb_hypgeom_gamma_upper_taylor_choose_jet(&N, err, a_real, z0, xmag, abs_tol, len))
    {
        _arb_vec_indeterminate(Gz1, len);
    }
    else
    {
        gamma_upper_taylor_bsplit_jet(A, B, C, D, E, F, Q, fmpq_numref(a), fmpq_denref(a), z0, x, 0, N, 0, len, prec);

        /* E * Gz0 */
        _arb_poly_get_vec(T, E, len);
        _arb_poly_mullow(U, T, len, Gz0, len, len, prec);

        /* - F * z0^(a-1+eps) e^(-z0) */
        fmpq_sub_ui(a1, a, 1);
        _arb_pow_fmpq_jet(T, z0, logz0, a1, len, prec);
        _arb_vec_scalar_mul(T, T, len, expmz0, prec);

        if (F->length != 0)
        {
            _arb_poly_mullow(Gz1, F->coeffs, FLINT_MIN(F->length, len), T, len, len, prec);
            _arb_vec_sub(U, U, Gz1, len, prec);
        }

        _arb_vec_scalar_div(Gz1, U, len, Q, prec);

        for (j = 0; j < len; j++)
            arb_add_error_mag(Gz1 + j, err);
    }

    mag_clear(xmag);
    mag_clear(err);
    arb_clear(x);
    arb_clear(Q);
    arb_clear(a_real);
    arb_clear(t);
    fmpq_clear(a1);
    arb_poly_clear(A);
    arb_poly_clear(B);
    arb_poly_clear(C);
    arb_poly_clear(D);
    arb_poly_clear(E);
    arb_poly_clear(F);
    _arb_vec_clear(T, len);
    _arb_vec_clear(U, len);
}
