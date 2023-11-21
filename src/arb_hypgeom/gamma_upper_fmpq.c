/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb_hypgeom.h"

static void
mag_div_fmpq(mag_t res, const mag_t x, const fmpq_t a)
{
    mag_t t;
    mag_init(t);
    mag_set_fmpz_lower(t, fmpq_numref(a));
    mag_div(res, x, t);
    mag_set_fmpz(t, fmpq_denref(a));
    mag_mul(res, res, t);
    mag_clear(t);
}

static void
mag_pow_fmpq_fast(mag_t res, const mag_t x, const fmpq_t e)
{
    fmpz_t b;
    fmpz_init(b);

    if (mag_cmp_2exp_si(x, 0) >= 0)
    {
        fmpz_cdiv_q(b, fmpq_numref(e), fmpq_denref(e));
        mag_pow_fmpz(res, x, b);
    }
    else
    {
        fmpz_fdiv_q(b, fmpq_numref(e), fmpq_denref(e));
        mag_pow_fmpz(res, x, b);
    }

    fmpz_clear(b);
}

slong
_arb_hypgeom_gamma_upper_fmpq_inf_choose_N_1(mag_t err, const fmpq_t a, const arb_t z, int prefactor, const mag_t abs_tol)
{
    slong N, aa, ab;
    fmpz_t az1, az2;
    fmpq_t a1;
    mag_t t, u;

    fmpz_init(az1);
    fmpz_init(az2);
    fmpq_init(a1);
    mag_init(t);
    mag_init(u);

    fmpz_fdiv_q(az1, fmpq_numref(a), fmpq_denref(a));
    fmpz_cdiv_q(az2, fmpq_numref(a), fmpq_denref(a));

    if (!fmpz_fits_si(az1) || !fmpz_fits_si(az2))
    {
        mag_inf(err);
        N = -1;
    }
    else
    {
        aa = fmpz_get_si(az1);
        ab = fmpz_get_si(az2);

        /* prefactor z^(a-1) * exp(-z) */
        if (prefactor)
        {
            arb_get_mag_lower(t, z);
            mag_expinv(t, t);
            fmpq_sub_ui(a1, a, 1);
            arb_get_mag(u, z);
            mag_pow_fmpq_fast(u, u, a1);
            mag_mul(err, t, u);
        }
        else
        {
            mag_one(err);
        }

        if (mag_is_inf(err))
        {
            N = -1;
        }
        else
        {
            arb_get_mag_lower(t, z);
            mag_inv(t, t);

            for (N = 1; ; N++)
            {
                mag_mul_ui(u, err, FLINT_MAX(FLINT_ABS(aa - N), FLINT_ABS(ab - N)));
                mag_mul(u, u, t);

                if (N >= ab - 1 && mag_cmp(u, abs_tol) < 0)
                {
                    mag_swap(err, u);
                    break;
                }

                /* Stop if terms are increasing, unless a is a positive integer in
                   which case the series will terminate eventually. */
                if (mag_cmp(u, err) > 0 && !(aa == ab && aa >= 1))
                {
                    mag_inf(err);
                    N = -1;
                    break;
                }

                mag_swap(err, u);
            }
        }
    }

    fmpz_clear(az1);
    fmpz_clear(az2);
    mag_clear(t);
    mag_clear(u);
    fmpq_clear(a1);

    return N;
}

slong
_arb_hypgeom_gamma_upper_fmpq_inf_choose_N(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol)
{
    return _arb_hypgeom_gamma_upper_fmpq_inf_choose_N_1(err, a, z, 1, abs_tol);
}

slong
_arb_hypgeom_gamma_upper_fmpq_inf_choose_N_rel(mag_t err, const fmpq_t a, const arb_t z, slong prec)
{
    mag_t tol;
    slong N;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);
    N = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N_1(err, a, z, 0, tol);
    mag_clear(tol);
    return N;
}

static void
upper_bsplit(arb_t M, arb_t S, arb_t Q, const fmpz_t ap, const fmpz_t aq, const arb_t z, slong na, slong nb, int cont, slong prec)
{
    if (nb - na == 0)
    {
        arb_zero(M);
        arb_zero(S);
        arb_one(Q);
    }
    else if (nb - na == 1)
    {
        fmpz_t t;
        fmpz_init_set(t, ap);
        fmpz_submul_ui(t, aq, na + 1);
        fmpz_neg(t, t);
        arb_set_fmpz(M, t);
        arb_mul_fmpz(S, z, aq, prec);
        arb_neg(S, S);
        arb_set(Q, S);
        fmpz_clear(t);
    }
    else
    {
        slong m;
        arb_t M2, S2, Q2;

        m = na + (nb - na) / 2;

        arb_init(M2);
        arb_init(S2);
        arb_init(Q2);

        upper_bsplit(M, S, Q, ap, aq, z, na, m, 1, prec);
        upper_bsplit(M2, S2, Q2, ap, aq, z, m, nb, cont, prec);

        /* todo: squaring opt; power table */
        arb_mul(S, S, Q2, prec);
        arb_addmul(S, M, S2, prec);

        if (cont)
            arb_mul(M, M, M2, prec);

        arb_mul(Q, Q, Q2, prec);

        arb_clear(M2);
        arb_clear(S2);
        arb_clear(Q2);
    }
}

void
_arb_hypgeom_gamma_upper_fmpq_inf_bsplit(arb_t res, const fmpq_t a, const arb_t z, slong N, slong prec)
{
    arb_t M, S, Q;
    fmpq_t a1;

    arb_init(M);
    arb_init(S);
    arb_init(Q);
    fmpq_init(a1);

    N = FLINT_MAX(N, 0);

    upper_bsplit(M, S, Q, fmpq_numref(a), fmpq_denref(a), z, 0, N, 0, prec);
    arb_div(S, S, Q, prec);

    fmpq_sub_ui(a1, a, 1);
    arb_pow_fmpq(M, z, a1, prec);
    arb_mul(S, S, M, prec);

    arb_neg(M, z);
    arb_exp(M, M, prec);
    arb_mul(res, S, M, prec);

    arb_clear(M);
    arb_clear(S);
    arb_clear(Q);
    fmpq_clear(a1);
}

/* select N and bound error for z**a * exp(-z) / a * sum(z**n / rf(a+1, n)) */
slong
_arb_hypgeom_gamma_lower_fmpq_0_choose_N(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol)
{
    slong N, aa, ab, c;
    fmpz_t az1, az2;
    mag_t t, u;

    fmpz_init(az1);
    fmpz_init(az2);
    mag_init(t);
    mag_init(u);

    fmpz_fdiv_q(az1, fmpq_numref(a), fmpq_denref(a));
    fmpz_cdiv_q(az2, fmpq_numref(a), fmpq_denref(a));

    if (!fmpz_fits_si(az1) || !fmpz_fits_si(az2))
    {
        mag_inf(err);
        N = -1;
    }
    else
    {
        aa = fmpz_get_si(az1);
        ab = fmpz_get_si(az2);

        /* prefactor z^a * exp(-z) / a */
        arb_get_mag_lower(t, z);
        mag_expinv(t, t);
        arb_get_mag(u, z);
        mag_pow_fmpq_fast(u, u, a);
        mag_mul(err, t, u);
        mag_div_fmpq(err, err, a);

        arb_get_mag(t, z);

        for (N = 1; ; N++)
        {
            c = FLINT_MIN(FLINT_ABS(aa + N), FLINT_ABS(ab + N));

            if (c == 0)
            {
                fmpq_t q;
                fmpq_init(q);
                fmpq_add_ui(q, a, N);
                mag_div_fmpq(err, err, q);
                fmpq_clear(q);
            }
            else
            {
                mag_div_ui(err, err, c);
            }

            mag_mul(err, err, t);

            /* todo: condition can be relaxed */
            /* todo: faster check (compare t) */
            if ((aa + N) > 0 && mag_cmp(err, abs_tol) < 0)
            {
                mag_div_ui(u, t, aa + N);
                mag_geom_series(u, u, 0);

                mag_mul(u, err, u);

                if (mag_cmp(u, abs_tol) < 0)
                {
                    mag_swap(err, u);
                    break;
                }
            }
        }
    }

    fmpz_clear(az1);
    fmpz_clear(az2);
    mag_clear(t);
    mag_clear(u);

    return N;
}


/* todo: squaring opt; power table */
static void
lower_bsplit(arb_t M, arb_t S, arb_t Q, const fmpz_t ap, const fmpz_t aq, const arb_t z, slong na, slong nb, int cont, slong prec)
{
    if (nb - na == 0)
    {
        arb_zero(M);
        arb_zero(S);
        arb_one(Q);
    }
    else if (nb - na == 1)
    {
        fmpz_t t;
        fmpz_init_set(t, ap);
        fmpz_addmul_ui(t, aq, na + 1);
        arb_set_fmpz(S, t);
        arb_set(Q, S);
        arb_mul_fmpz(M, z, aq, prec);
        fmpz_clear(t);
    }
    else
    {
        slong m;
        arb_t M2, S2, Q2;

        m = na + (nb - na) / 2;

        arb_init(M2);
        arb_init(S2);
        arb_init(Q2);

        lower_bsplit(M, S, Q, ap, aq, z, na, m, 1, prec);
        lower_bsplit(M2, S2, Q2, ap, aq, z, m, nb, cont, prec);

        arb_mul(S, S, Q2, prec);
        arb_addmul(S, M, S2, prec);

        if (cont)
            arb_mul(M, M, M2, prec);

        arb_mul(Q, Q, Q2, prec);

        arb_clear(M2);
        arb_clear(S2);
        arb_clear(Q2);
    }
}

void
_arb_hypgeom_gamma_lower_fmpq_0_bsplit(arb_t res, const fmpq_t a, const arb_t z, slong N, slong prec)
{
    arb_t M, S, Q;

    arb_init(M);
    arb_init(S);
    arb_init(Q);

    N = FLINT_MAX(N, 0);

    lower_bsplit(M, S, Q, fmpq_numref(a), fmpq_denref(a), z, 0, N, 0, prec);
    arb_div(S, S, Q, prec);

    arb_pow_fmpq(M, z, a, prec);
    arb_mul(S, S, M, prec);

    arb_neg(M, z);
    arb_exp(M, M, prec);
    arb_mul(S, S, M, prec);

    arb_div_fmpz(S, S, fmpq_numref(a), prec);
    arb_mul_fmpz(res, S, fmpq_denref(a), prec);

    arb_clear(M);
    arb_clear(S);
    arb_clear(Q);
}

/* bounded by 1/z^n * sum z^k / k! */
slong
_arb_hypgeom_gamma_upper_singular_si_choose_N(mag_t err, slong n, const arb_t z, const mag_t abs_tol)
{
    slong k;
    mag_t t, u, zm;

    mag_init(t);
    mag_init(u);
    mag_init(zm);

    arb_get_mag(zm, z);

    arb_get_mag_lower(t, z);
    mag_inv(t, t);
    mag_pow_ui(t, t, n);

    for (k = 1; ; k++)
    {
        mag_mul(t, t, zm);
        mag_div_ui(t, t, k);

        if (mag_cmp(t, abs_tol) < 0)
        {
            mag_div_ui(u, zm, k);
            mag_geom_series(u, u, 0);
            mag_mul(u, t, u);

            if (mag_cmp(u, abs_tol) < 0)
            {
                mag_swap(err, t);
                break;
            }
        }
    }

    mag_clear(t);
    mag_clear(u);
    mag_clear(zm);

    return k;
}

static void
singular_bsplit(arb_t M, arb_t S, arb_t Q, slong n, const arb_t z, slong na, slong nb, int cont, slong prec)
{
    if (nb - na == 0)
    {
        arb_zero(M);
        arb_zero(S);
        arb_one(Q);
    }
    else if (nb - na == 1)
    {
        fmpz_t t;
        slong k;

        k = na;

        if (k == n)
            arb_neg(M, z);
        else
            arb_mul_si(M, z, n - k, prec);

        arb_set_si(S, (k != n) ? (k + 1) : 0);

        fmpz_init_set_si(t, k + 1);
        if (k != n)
            fmpz_mul_si(t, t, k - n);
        arb_set_fmpz(Q, t);
        fmpz_clear(t);
    }
    else
    {
        slong m;
        arb_t M2, S2, Q2;

        m = na + (nb - na) / 2;

        arb_init(M2);
        arb_init(S2);
        arb_init(Q2);

        singular_bsplit(M, S, Q, n, z, na, m, 1, prec);
        singular_bsplit(M2, S2, Q2, n, z, m, nb, cont, prec);

        arb_mul(S, S, Q2, prec);
        arb_addmul(S, M, S2, prec);

        if (cont)
            arb_mul(M, M, M2, prec);

        arb_mul(Q, Q, Q2, prec);

        arb_clear(M2);
        arb_clear(S2);
        arb_clear(Q2);
    }
}

void
_arb_hypgeom_gamma_upper_singular_si_bsplit(arb_t res, slong n, const arb_t z, slong N, slong prec)
{
    arb_t M, S, Q;

    arb_init(M);
    arb_init(S);
    arb_init(Q);

    N = FLINT_MAX(N, 0);

    singular_bsplit(M, S, Q, n, z, 0, N, 0, prec);

    /* (-1)**n/fac(n) * (digamma(n+1) - ln(z)) - (S/Q)/z**n */
    arb_pow_ui(M, z, n, prec);
    arb_mul(Q, Q, M, prec);
    arb_div(S, S, Q, prec);

    arb_set_ui(M, n + 1);
    arb_digamma(M, M, prec);
    arb_log(Q, z, prec);
    arb_sub(M, M, Q, prec);
    arb_fac_ui(Q, n, prec);
    arb_div(M, M, Q, prec);
    if (n & 1)
        arb_neg(M, M);

    arb_sub(res, M, S, prec);

    arb_clear(M);
    arb_clear(S);
    arb_clear(Q);
}
