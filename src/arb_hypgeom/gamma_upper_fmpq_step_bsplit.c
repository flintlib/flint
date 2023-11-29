/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpq.h"
#include "arb_mat.h"
#include "arb_hypgeom.h"

#ifdef __GNUC__
# define ceil __builtin_ceil
#else
# include <math.h>
#endif

static void
taylor_M(mag_t M, const arb_t a, const arb_t z, const mag_t x, slong Rexp)
{
    arb_t t, u;
    arb_init(t);
    arb_init(u);

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

    arb_clear(t);
    arb_clear(u);
}

/* choose N such that M * C^N / (1 - C) <= tol */
/* todo: fix */
static slong
mag_geom_choose_N(const mag_t M, const mag_t C, const mag_t tol)
{
    mag_t t, u;
    slong N;

    if (mag_is_finite(M) && mag_is_zero(C))
        return 1;

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

    N = ceil(mag_get_d(t));
    N = FLINT_MAX(N, 1);

    mag_clear(t);
    mag_clear(u);

    return N;
}

static void
taylor_bound(mag_t err, const arb_t a, const arb_t z, const mag_t x, slong Rexp, slong N)
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
        taylor_M(M, a, z, x, Rexp);
        mag_mul(err, err, M);

        mag_mul_2exp_si(err, err, Rexp);
        mag_div_ui(err, err, N);
    }

    mag_clear(C);
    mag_clear(M);
}

static slong
taylor_N(const arb_t a, const arb_t z, const mag_t x, slong Rexp, const mag_t abs_tol)
{
    mag_t C, M;
    slong N;

    mag_init(C);
    mag_init(M);

    /* C = x / R */
    mag_mul_2exp_si(C, x, -Rexp);

    if (mag_cmp_2exp_si(C, 0) < 0)
    {
        taylor_M(M, a, z, x, Rexp);
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

static void
arb_hypgeom_gamma_upper_taylor_choose(slong * res_N, mag_t err, const arb_t a, const arb_t z, const mag_t x, const mag_t abs_tol)
{
    slong N, New;
    mag_t zlow;
    slong Rexp;

    mag_init(zlow);
    arb_get_mag_lower(zlow, z);

    Rexp = 0;
    while (mag_cmp_2exp_si(zlow, Rexp + 1) < 0)
        Rexp--;

    N = taylor_N(a, z, x, Rexp, abs_tol);

    while (N > 1 && mag_cmp_2exp_si(x, Rexp - 1) < 0)
    {
        New = taylor_N(a, z, x, Rexp - 1, abs_tol);

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
            New = taylor_N(a, z, x, Rexp + 1, abs_tol);

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

    *res_N = N;
    taylor_bound(err, a, z, x, Rexp, N);

    if (mag_cmp(err, abs_tol) > 0)
    {
        /* TODO: Find a nice way to print this using flint_throw */
        printf("err = "); mag_printd(err, 10);
        printf("\nabs_tol = "); mag_printd(abs_tol, 10);
        printf("\na = "); arb_printd(a, 10);
        printf("\nz = "); arb_printd(z, 10);
        printf("\nx = "); mag_printd(x, 10);
        flint_abort();
    }

    mag_clear(zlow);
}

static void
gamma_upper_taylor_bsplit(arb_mat_t M, arb_t Q,
    const fmpz_t ap, const fmpz_t aq, const arb_t z0, const arb_t x, const arb_t x2, slong a, slong b, int cont, slong prec)
{
    if (b - a == 0)
    {
        arb_mat_one(M);
    }
    else if (b - a == 1)
    {
        slong n;
        fmpz_t t;

        n = a;
        fmpz_init(t);

        /* Q = -z0*(n+1)*(n+2)*aq */
        fmpz_mul2_uiui(t, aq, n + 1, n + 2);
        arb_mul_fmpz(Q, z0, t, prec);
        arb_neg(Q, Q);

        /* x Q */
        arb_mul(arb_mat_entry(M, 0, 1), Q, x, prec);

        /* aq n x */
        fmpz_mul_ui(t, aq, n);
        arb_mul_fmpz(arb_mat_entry(M, 1, 0), x, t, prec);

        /* x*(-ap + aq*(n + z0 + 1))*(n + 1) */
        arb_add_ui(arb_mat_entry(M, 1, 1), z0, n + 1, prec);
        arb_mul_fmpz(arb_mat_entry(M, 1, 1), arb_mat_entry(M, 1, 1), aq, prec);
        arb_sub_fmpz(arb_mat_entry(M, 1, 1), arb_mat_entry(M, 1, 1), ap, prec);
        arb_mul_ui(arb_mat_entry(M, 1, 1), arb_mat_entry(M, 1, 1), n + 1, prec);
        arb_mul(arb_mat_entry(M, 1, 1), arb_mat_entry(M, 1, 1), x, prec);

        arb_set(arb_mat_entry(M, 2, 0), Q);
        arb_set(arb_mat_entry(M, 2, 2), Q);

        fmpz_clear(t);
    }
    else
    {
        arb_mat_t M1, M2;
        arb_t Q2;
        slong m;

        arb_mat_init(M1, 3, 3);
        arb_mat_init(M2, 3, 3);
        arb_init(Q2);

        m = a + (b - a) / 2;

        gamma_upper_taylor_bsplit(M1, Q, ap, aq, z0, x, x2, a, m, 1, prec);
        gamma_upper_taylor_bsplit(M2, Q2, ap, aq, z0, x, x2, m, b, cont, prec);

        if (cont)
        {
            /* todo: parallel slowdown; arb_mat_mul needs retuning */
            arb_mat_mul_classical(M, M2, M1, prec);
        }
        else
        {
            arb_mat_transpose(M1, M1);

            arb_dot(arb_mat_entry(M, 2, 0), NULL, 0, arb_mat_entry(M1, 0, 0), 1, arb_mat_entry(M2, 2, 0), 1, 3, prec);
            arb_dot(arb_mat_entry(M, 2, 1), NULL, 0, arb_mat_entry(M1, 1, 0), 1, arb_mat_entry(M2, 2, 0), 1, 3, prec);
            arb_dot(arb_mat_entry(M, 2, 2), NULL, 0, arb_mat_entry(M1, 2, 0), 1, arb_mat_entry(M2, 2, 0), 1, 3, prec);
        }

        arb_mul(Q, Q2, Q, prec);

        arb_mat_clear(M1);
        arb_mat_clear(M2);
        arb_clear(Q2);
    }
}

/*
Given Gz0 = Gamma(a, z0) and expmz0 = exp(-z0), compute Gz1 = Gamma(a, z1)
*/
void
_arb_gamma_upper_fmpq_step_bsplit(arb_t Gz1, const fmpq_t a, const arb_t z0, const arb_t z1, const arb_t Gz0, const arb_t expmz0, const mag_t abs_tol, slong prec)
{
    arb_t x, Q, a_real;
    arb_mat_t M;
    mag_t xmag, err;
    slong N;
    fmpq_t a1;

    if (arb_is_zero(z0))
    {
        mag_init(err);
        arb_init(x);
        N = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(err, a, z1, abs_tol);
        _arb_hypgeom_gamma_lower_fmpq_0_bsplit(Gz1, a, z1, N, prec);
        arb_add_error_mag(Gz1, err);
        arb_sub(Gz1, Gz0, Gz1, prec);
        arb_clear(x);
        mag_clear(err);
        return;
    }

    mag_init(xmag);
    mag_init(err);
    arb_init(x);
    arb_init(Q);
    arb_init(a_real);
    fmpq_init(a1);
    arb_mat_init(M, 3, 3);

    arb_sub(x, z1, z0, prec);
    arb_get_mag(xmag, x);

    arb_set_fmpq(a_real, a, 53);
    arb_hypgeom_gamma_upper_taylor_choose(&N, err, a_real, z0, xmag, abs_tol);

    gamma_upper_taylor_bsplit(M, Q, fmpq_numref(a), fmpq_denref(a), z0, x, NULL, 0, N, 0, prec);

    arb_mul(arb_mat_entry(M, 2, 0), arb_mat_entry(M, 2, 0), Gz0, prec);

    fmpq_sub_ui(a1, a, 1);
    arb_pow_fmpq(arb_mat_entry(M, 0, 0), z0, a1, prec);
    arb_mul(arb_mat_entry(M, 0, 0), arb_mat_entry(M, 0, 0), expmz0, prec);
    arb_submul(arb_mat_entry(M, 2, 0), arb_mat_entry(M, 2, 1), arb_mat_entry(M, 0, 0), prec);

    arb_div(Gz1, arb_mat_entry(M, 2, 0), Q, prec);

    arb_add_error_mag(Gz1, err);

    mag_clear(xmag);
    mag_clear(err);
    arb_clear(x);
    arb_clear(Q);
    arb_clear(a_real);
    fmpq_clear(a1);
    arb_mat_clear(M);
}
