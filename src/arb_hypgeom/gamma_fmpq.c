/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "fmpz_poly.h"
#include "fmpq.h"
#include "arb_hypgeom.h"
#include "hypgeom.h"

void
arb_gamma_const_1_3_eval(arb_t s, slong prec)
{
    hypgeom_t series;
    arb_t t, u;
    slong wp = prec + 4 + 2 * FLINT_BIT_COUNT(prec);

    arb_init(t);
    arb_init(u);

    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "1  1");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "4  5 -46 108 -72");
    fmpz_poly_set_str(series->Q, "4  0 0 0 512000");

    prec += FLINT_CLOG2(prec);

    arb_hypgeom_infsum(s, t, series, wp, wp);

    arb_sqrt_ui(u, 10, wp);
    arb_mul(t, t, u, wp);

    arb_const_pi(u, wp);
    arb_pow_ui(u, u, 4, wp);
    arb_mul_ui(u, u, 12, wp);
    arb_mul(s, s, u, wp);

    arb_div(s, s, t, wp);
    arb_root_ui(s, s, 2, wp);
    arb_root_ui(s, s, 3, prec);

    hypgeom_clear(series);
    arb_clear(t);
    arb_clear(u);
}

ARB_DEF_CACHED_CONSTANT(arb_gamma_const_1_3, arb_gamma_const_1_3_eval)

void
arb_gamma_const_1_4_eval(arb_t x, slong prec)
{
    arb_t t, u;
    slong wp = prec + 4 + 2 * FLINT_BIT_COUNT(prec);

    arb_init(t);
    arb_init(u);

    arb_one(t);
    arb_sqrt_ui(u, 2, wp);
    arb_agm(x, t, u, wp);

    arb_const_pi(t, wp);
    arb_mul_2exp_si(t, t, 1);
    arb_sqrt(u, t, wp);
    arb_mul(u, u, t, wp);

    arb_div(x, u, x, wp);
    arb_sqrt(x, x, wp);

    arb_clear(t);
    arb_clear(u);
}

ARB_DEF_CACHED_CONSTANT(arb_gamma_const_1_4, arb_gamma_const_1_4_eval)


void arb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n, const arb_t x, int use_reflect, int digamma, slong prec);
void arb_hypgeom_gamma_stirling_inner(arb_t s, const arb_t z, slong N, slong prec);

void
arb_hypgeom_gamma_fmpq_stirling(arb_t y, const fmpq_t a, slong prec)
{
    int reflect;
    slong r, n, wp;
    arb_t x, t, u, v;

    wp = (slong) fmpz_bits(fmpq_numref(a)) - (slong) fmpz_bits(fmpq_denref(a));
    wp = FLINT_MAX(wp, 0);
    wp += prec + FLINT_BIT_COUNT(prec);

    arb_init(x);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    arb_set_fmpq(x, a, wp);
    arb_hypgeom_gamma_stirling_choose_param(&reflect, &r, &n, x, 1, 0, wp);

    if (reflect)
    {
        /* gamma(x) = (rf(1-x, r) * pi) / (gamma(1-x+r) sin(pi x)) */
        fmpq_t b;
        fmpq_init(b);
        fmpz_sub(fmpq_numref(b), fmpq_denref(a), fmpq_numref(a));
        fmpz_set(fmpq_denref(b), fmpq_denref(a));
        arb_rising_fmpq_ui(u, b, r, wp);
        fmpq_clear(b);
        arb_const_pi(v, wp);
        arb_mul(u, u, v, wp);
        arb_sub_ui(t, x, 1, wp);
        arb_neg(t, t);
        arb_add_ui(t, t, r, wp);
        arb_hypgeom_gamma_stirling_inner(v, t, n, wp);
        arb_exp(v, v, wp);
        arb_sin_pi_fmpq(t, a, wp);
        arb_mul(v, v, t, wp);
    }
    else
    {
        /* gamma(x) = gamma(x+r) / rf(x,r) */
        arb_add_ui(t, x, r, wp);
        arb_hypgeom_gamma_stirling_inner(u, t, n, wp);
        arb_exp(u, u, prec);
        arb_rising_fmpq_ui(v, a, r, wp);
    }

    arb_div(y, u, v, prec);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    arb_clear(x);
}

void
arb_hypgeom_gamma_small_frac(arb_t y, unsigned int p, unsigned int q, slong prec)
{
    slong wp = prec + 4;

    if (q == 1)
    {
        arb_one(y);
    }
    else if (q == 2)  /* p = 1 */
    {
        arb_const_sqrt_pi(y, prec);
    }
    else if (q == 3)
    {
        if (p == 1)
        {
            arb_gamma_const_1_3(y, prec);
        }
        else  /* p = 2 */
        {
            arb_t t;
            arb_init(t);
            arb_gamma_const_1_3(y, wp);
            arb_sqrt_ui(t, 3, wp);
            arb_mul(y, y, t, wp);
            arb_const_pi(t, wp);
            arb_div(y, t, y, prec);
            arb_mul_2exp_si(y, y, 1);
            arb_clear(t);
        }
    }
    else if (q == 4)
    {
        if (p == 1)
        {
            arb_gamma_const_1_4(y, prec);
        }
        else  /* p = 3 */
        {
            arb_t t;
            arb_init(t);
            arb_sqrt_ui(y, 2, wp);
            arb_const_pi(t, wp);
            arb_mul(y, y, t, wp);
            arb_gamma_const_1_4(t, wp);
            arb_div(y, y, t, prec);
            arb_clear(t);
        }
    }
    else if (q == 6)
    {
        arb_t t;
        arb_init(t);
        arb_const_pi(t, wp);
        arb_div_ui(t, t, 3, wp);
        arb_sqrt(t, t, wp);
        arb_set_ui(y, 2);
        arb_root_ui(y, y, 3, wp);
        arb_mul(t, t, y, wp);
        arb_gamma_const_1_3(y, wp);
        arb_mul(y, y, y, prec);

        if (p == 1)
        {
            arb_div(y, y, t, prec);
        }
        else  /* p = 5 */
        {
            arb_div(y, t, y, wp);
            arb_const_pi(t, wp);
            arb_mul(y, y, t, prec);
            arb_mul_2exp_si(y, y, 1);
        }

        arb_clear(t);
    }
    else
    {
        flint_throw(FLINT_ERROR, "small fraction not implemented!\n");
    }
}

slong _arb_compute_bs_exponents(slong * tab, slong n);
slong _arb_get_exp_pos(const slong * tab, slong step);

static void
bsplit2(arb_t P, arb_t Q, const fmpz_t zp, const fmpz_t zq,
    const slong * xexp, arb_srcptr xpow,
    ulong N, slong a, slong b, int cont, slong prec)
{
    if (b - a == 1)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set(t, zp);
        fmpz_addmul_ui(t, zq, a + 1);
        arb_set_fmpz(P, t);
        arb_set(Q, P);
        fmpz_clear(t);
    }
    else
    {
        arb_t Pb, Qb;
        slong step, i, m;

        arb_init(Pb);
        arb_init(Qb);

        step = (b - a) / 2;
        m = a + step;

        bsplit2(P, Q, zp, zq, xexp, xpow, N, a, m, 1, prec);
        bsplit2(Pb, Qb, zp, zq, xexp, xpow, N, m, b, 1, prec);

        arb_mul(P, P, Pb, prec);
        arb_mul(Q, Q, Pb, prec);

        i = (step == 1) ? 0 : _arb_get_exp_pos(xexp, step);
        arb_addmul(Q, Qb, xpow + i, prec);

        arb_clear(Pb);
        arb_clear(Qb);
    }
}

static void
bsplit3(arb_t P, arb_t Q, const fmpz_t zp, const fmpz_t zq,
    const slong * xexp, arb_srcptr xpow,
    ulong N, slong a, slong b, int cont, slong prec)
{
    if (b - a == 1)
    {
        fmpz_t t;
        fmpz_init(t);
        arb_set(P, xpow + 0);  /* N zq */
        fmpz_set(t, zp);
        fmpz_submul_ui(t, zq, a + 1);  /* zp - (a + 1) zq */
        arb_set_fmpz(Q, t);
        fmpz_clear(t);
    }
    else
    {
        arb_t Pb, Qb;
        slong step, i, m;

        arb_init(Pb);
        arb_init(Qb);

        step = (b - a) / 2;
        m = a + step;

        bsplit3(P, Q, zp, zq, xexp, xpow, N, a, m, 1, prec);
        bsplit3(Pb, Qb, zp, zq, xexp, xpow, N, m, b, 1, prec);

        i = _arb_get_exp_pos(xexp, m - a);

        arb_mul(P, P, xpow + i, prec);
        if (b - m != m - a)
            arb_mul(P, P, xpow + 0, prec);

        arb_addmul(P, Q, Pb, prec);

        if (cont)
        {
            arb_mul(Q, Q, Qb, prec);
        }
        else
        {
            i = _arb_get_exp_pos(xexp, m - a);

            arb_mul(Q, xpow + i, xpow + i, prec);
            if (b - m != m - a)
                arb_mul(Q, Q, xpow + 0, prec);
        }

        arb_clear(Pb);
        arb_clear(Qb);
    }
}

double d_lambertw_branch1(double x);

static ulong
more_trailing_zeros(ulong N)
{
    ulong bc, N2;

    bc = FLINT_BIT_COUNT(N);

    if (bc >= 8)
    {
        N2 = (N >> (bc - 5)) << (bc - 5);
        N2 += ((N2 != N) << (bc - 5));
        return N2;
    }
    else
    {
        return N;
    }
}

#define C_LOG2 0.69314718055994530942
#define C_EXP1 2.7182818284590452354

static void
build_bsplit_power_table(arb_ptr xpow, const slong * xexp, slong len, slong prec)
{
    slong i;

    for (i = 1; i < len; i++)
    {
        if (xexp[i] == 2 * xexp[i-1])
        {
            arb_mul(xpow + i, xpow + i - 1, xpow + i - 1, prec);
        }
        else if (xexp[i] == 2 * xexp[i-2]) /* prefer squaring if possible */
        {
            arb_mul(xpow + i, xpow + i - 2, xpow + i - 2, prec);
        }
        else if (xexp[i] == 2 * xexp[i-1] + 1)
        {
            arb_mul(xpow + i, xpow + i - 1, xpow + i - 1, prec);
            arb_mul(xpow + i, xpow + i, xpow, prec);
        }
        else if (xexp[i] == 2 * xexp[i-2] + 1)
        {
            arb_mul(xpow + i, xpow + i - 2, xpow + i - 2, prec);
            arb_mul(xpow + i, xpow + i, xpow, prec);
        }
        else
        {
            flint_throw(FLINT_ERROR, "power table has the wrong structure!\n");
        }
    }
}

/* assumes z in [1, 2] */
static void
arb_hypgeom_gamma_fmpq_general_off1(arb_t res, const fmpq_t z, slong prec)
{
    slong wp, N, n, n2, length, length2, wp2;
    double alpha;
    arb_t P, Q;
    slong *xexp, *xexp2;
    arb_ptr xpow;
    mag_t err, err2;

    wp = prec + 30;

    alpha = 0.52;  /* tuning parameter between 0.5 and 1 */

    N = alpha * C_LOG2 * wp;
    N = more_trailing_zeros(N);
    alpha = N / (C_LOG2 * wp);

    /* Terms in convergent series */
    n = (1 - alpha) / d_lambertw((1 - alpha) / (alpha * C_EXP1)) * C_LOG2 * wp;

    /* Precision and terms in asymptotic series */
    wp2 = wp * (1 - alpha);
    wp2 = FLINT_MAX(wp2, 30);
    n2 = (alpha - 1) / d_lambertw_branch1((alpha - 1) / (alpha * C_EXP1)) * C_LOG2 * wp;
    n2 = FLINT_MAX(n2, 2);  /* binary splitting correctness */

    mag_init(err);
    mag_init(err2);
    arb_init(P);
    arb_init(Q);

    /* compute the powers of x = N*zq that will appear (at least x^1) */
    xexp = flint_calloc(2 * FLINT_BITS, sizeof(slong));
    xexp2 = flint_calloc(2 * FLINT_BITS, sizeof(slong));

    length = _arb_compute_bs_exponents(xexp, n);
    length2 = _arb_compute_bs_exponents(xexp2, n2);

    xpow = _arb_vec_init(FLINT_MAX(length, length2));

    arb_set_fmpz(xpow + 0, fmpq_denref(z));
    arb_mul_ui(xpow + 0, xpow + 0, N, wp);

    build_bsplit_power_table(xpow, xexp, length, wp);

    /* 1F1(1, 1+z, N) */
    bsplit2(P, Q, fmpq_numref(z), fmpq_denref(z), xexp, xpow, N, 0, n, 0, wp);
    arb_div(P, Q, P, wp);

    /* Convergent series error bound: N^n / n! (1 + (N/n) + ...) */
    mag_set_ui(err, N);
    mag_pow_ui(err, err, n);
    mag_rfac_ui(err2, n);
    mag_mul(err, err, err2);
    mag_set_ui(err2, N);
    mag_div_ui(err2, err2, n);
    mag_geom_series(err2, err2, 0);
    mag_mul(err, err, err2);
    arb_add_error_mag(P, err);

    /* divide 1F1 by z */
    arb_mul_fmpz(P, P, fmpq_denref(z), wp);
    arb_div_fmpz(P, P, fmpq_numref(z), wp);
    arb_swap(res, P);

    build_bsplit_power_table(xpow, xexp2, length2, wp2);

    bsplit3(P, Q, fmpq_numref(z), fmpq_denref(z), xexp2, xpow, N, 0, n2, 0, wp2);
    arb_div(P, P, Q, wp2);

    /* 2F0 error bound (bounded by first omitted term) */
    mag_fac_ui(err, n2);
    mag_set_ui_lower(err2, N);
    mag_pow_ui_lower(err2, err2, n2);
    mag_div(err, err, err2);
    arb_add_error_mag(P, err);

    /* N^z * exp(-N) * (1F1/z + 2F0/N)  */
    arb_div_ui(P, P, N, wp2);

    arb_add(res, res, P, wp);
    arb_set_ui(Q, N);
    arb_pow_fmpq(Q, Q, z, wp);
    arb_mul(res, res, Q, wp);

    arb_set_si(Q, -N);
    arb_exp(Q, Q, wp);
    arb_mul(res, res, Q, wp);

    _arb_vec_clear(xpow, FLINT_MAX(length, length2));
    flint_free(xexp);
    flint_free(xexp2);

    arb_clear(P);
    arb_clear(Q);
    mag_clear(err);
    mag_clear(err2);
}

/* assumes z in (0, 1] */
void
arb_hypgeom_gamma_fmpq_hyp(arb_t res, const fmpq_t z, slong prec)
{
    fmpq_t t;
    fmpq_init(t);
    fmpq_add_ui(t, z, 1);
    arb_hypgeom_gamma_fmpq_general_off1(res, t, prec);
    arb_mul_fmpz(res, res, fmpq_denref(z), prec + 30);
    arb_div_fmpz(res, res, fmpq_numref(z), prec);
    fmpq_clear(t);
}

void
arb_hypgeom_gamma_fmpq_outward(arb_t y, const fmpq_t x, slong prec)
{
    fmpq_t a;
    fmpz_t n;
    fmpz p, q;
    slong m;
    arb_t t, u;

    fmpq_init(a);
    fmpz_init(n);
    arb_init(t);
    arb_init(u);

    /* write x = a + n with 0 < a <= 1 */
    if (fmpz_is_one(fmpq_denref(x)))
    {
        fmpq_one(a);
        fmpz_sub_ui(n, fmpq_numref(x), 1);
    }
    else
    {
        fmpz_fdiv_qr(n, fmpq_numref(a), fmpq_numref(x), fmpq_denref(x));
        fmpz_set(fmpq_denref(a), fmpq_denref(x));
    }

    if (!fmpz_fits_si(n))
    {
        flint_throw(FLINT_ERROR, "gamma: too large fmpq to reduce to 0!\n");
    }

    m = fmpz_get_si(n);

    /* evaluate gamma(a) */
    p = *fmpq_numref(a);
    q = *fmpq_denref(a);

    if (q == 1 || q == 2 || q == 3 || q == 4 || q == 6)
    {
        arb_hypgeom_gamma_small_frac(t, p, q, prec + 4 * (m != 0));
    }
    else
    {
        arb_hypgeom_gamma_fmpq_hyp(t, a, prec + 4 * (m != 0));
    }

    /* argument reduction */
    if (m >= 0)
    {
        arb_rising_fmpq_ui(u, a, m, prec + 4);
        arb_mul(y, t, u, prec);
    }
    else
    {
        arb_rising_fmpq_ui(u, x, -m, prec + 4);
        arb_div(y, t, u, prec);
    }

    fmpq_clear(a);
    fmpz_clear(n);
    arb_clear(t);
    arb_clear(u);
}

int
arb_hypgeom_gamma_fmpq_taylor(arb_t y, const fmpq_t x, slong prec)
{
    fmpq_t a;
    fmpz_t n;
    slong m;
    arb_t t;
    int success;

    fmpq_init(a);
    fmpz_init(n);
    arb_init(t);

    success = 1;

    /* write x = a + n with 0 < a <= 1 */
    if (fmpz_is_one(fmpq_denref(x)))
    {
        fmpq_one(a);
        fmpz_sub_ui(n, fmpq_numref(x), 1);
    }
    else
    {
        fmpz_fdiv_qr(n, fmpq_numref(a), fmpq_numref(x), fmpq_denref(x));
        fmpz_set(fmpq_denref(a), fmpq_denref(x));
    }

    if (!fmpz_fits_si(n))
    {
        success = 0;
        goto cleanup;
    }

    m = fmpz_get_si(n);

    if (m < -(40 + (prec - 40) / 4))
    {
        success = 0;
        goto cleanup;
    }

    if (m > 70 + (prec - 40) / 8)
    {
        success = 0;
        goto cleanup;
    }

    arb_set_fmpq(t, a, prec + 4);
    success = arb_hypgeom_gamma_taylor(t, t, 0, prec + 4);

    if (success)
    {
        arb_t u;
        arb_init(u);

        if (m >= 0)
        {
            arb_rising_fmpq_ui(u, a, m, prec + 4);
            arb_mul(y, t, u, prec);
        }
        else
        {
            arb_rising_fmpq_ui(u, x, -m, prec + 4);
            arb_div(y, t, u, prec);
        }

        arb_clear(u);
    }

cleanup:
    fmpq_clear(a);
    fmpz_clear(n);
    arb_clear(t);

    return success;
}

void
arb_hypgeom_gamma_fmpq(arb_t y, const fmpq_t x, slong prec)
{
    fmpz p, q;

    p = *fmpq_numref(x);
    q = *fmpq_denref(x);

    if ((q == 1 || q == 2 || q == 3 || q == 4 || q == 6) && !COEFF_IS_MPZ(p))
    {
        if (q == 1)
        {
            if (p <= 0)
            {
                arb_indeterminate(y);
                return;
            }

            if (p < 1200 || 1.44265 * (p*log(p) - p) < 15.0 * prec)
            {
                fmpz_t t;
                fmpz_init(t);
                fmpz_fac_ui(t, p - 1);
                arb_set_round_fmpz(y, t, prec);
                fmpz_clear(t);
                return;
            }
        }

        p = FLINT_ABS(p);

        if (p < q * 500.0 || p < q * (500.0 + 0.1 * prec * sqrt(prec)))
        {
            arb_hypgeom_gamma_fmpq_outward(y, x, prec);
            return;
        }
    }

    if (q != 1 && prec > 7000 + 300 * fmpz_bits(fmpq_denref(x)) &&
        (slong) fmpz_bits(fmpq_numref(x)) - (slong) fmpz_bits(fmpq_denref(x)) < FLINT_BITS &&
        fabs(fmpq_get_d(x)) < 0.03 * prec * sqrt(prec))
    {
        arb_hypgeom_gamma_fmpq_outward(y, x, prec);
        return;
    }

    if (fmpz_bits(fmpq_denref(x)) > 0.1 * prec || fmpz_bits(fmpq_numref(x)) > 0.1 * prec)
    {
        slong wp;

        wp = (slong) fmpz_bits(fmpq_numref(x)) - (slong) fmpz_bits(fmpq_denref(x));
        wp = FLINT_MAX(wp, 0);
        wp = FLINT_MIN(wp, 10 * prec);
        wp += prec + 4;

        arb_set_fmpq(y, x, wp);

        if (!arb_hypgeom_gamma_taylor(y, y, 0, prec))
            arb_hypgeom_gamma_stirling(y, y, 0, prec);

        return;
    }

    if (arb_hypgeom_gamma_fmpq_taylor(y, x, prec))
        return;

    arb_hypgeom_gamma_fmpq_stirling(y, x, prec);
}

void
arb_hypgeom_gamma_fmpz(arb_t y, const fmpz_t x, slong prec)
{
    fmpq_t t;
    *fmpq_numref(t) = *x;
    *fmpq_denref(t) = WORD(1);
    arb_hypgeom_gamma_fmpq(y, t, prec);
}

