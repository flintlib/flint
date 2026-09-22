/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpz_vec.h"
#include "arb.h"
#include "fixed.h"
#include "fmpzi.h"
#include "acb.h"
#include "arb/impl.h"

#define TERMINATOR -32768


/* The relation tables, primes and weights live in the fixed module
   (fixed/rel_tab.c, fixed_rel_table); the 13-prime tables there are
   the ones this file used to carry. */

static void
_arb_log_reduce_fixed(slong * rel, const short * d, const double * epsilon, const double * epsilon_inv,
    const fmpz * alpha, const float * weights,
    slong num_alpha, const fmpz_t x, slong prec, double max_weight)
{
    slong i, j, n;
    slong * new_rel;
    const short * d_row;
    double dalpha;
    double weight, dx;
    fmpz_t t;

    new_rel = flint_malloc(num_alpha * sizeof(slong));
    fmpz_init(t);

    for (i = 0; i < num_alpha; i++)
        rel[i] = 0;

    dx = fmpz_get_d(x);
    dx = ldexp(dx, -prec);

    /* Reduce by the first alpha, which is assumed to be free */
    dalpha = fmpz_get_d(alpha + 0);
    dalpha = ldexp(dalpha, -prec);
    n = floor(dx / dalpha + 0.5);
    dx -= dalpha * n;
    rel[0] = n;

    /* Recompute accurately if there is significant cancellation */
    if (FLINT_ABS(n) > 10)
    {
        fmpz_set(t, x);
        fmpz_submul_si(t, alpha + 0, rel[0]);
        dx = fmpz_get_d(t);
        dx = ldexp(dx, -prec);
    }

    for (i = 0; ; i++)
    {
        d_row = d + i * num_alpha;

        if (d_row[0] == TERMINATOR)
            break;

        for (j = 0; j < num_alpha; j++)
            new_rel[j] = d_row[j];

        n = floor(dx * epsilon_inv[i] + 0.5);

        if (n != 0)
        {
            weight = 0.0;
            for (j = 0; j < num_alpha; j++)
            {
                new_rel[j] = rel[j] + n * new_rel[j];
                if (j != 0)
                    weight += FLINT_ABS(new_rel[j]) * weights[j] * 1.442695;
            }

            if (weight > max_weight)
                break;

            for (j = 0; j < num_alpha; j++)
                rel[j] = new_rel[j];

            dx -= n * epsilon[i];
        }

        if (i % 8 == 7)
        {
            fmpz_set(t, x);
            for (j = 0; j < num_alpha; j++)
                fmpz_submul_si(t, alpha + j, rel[j]);

            dx = fmpz_get_d(t);
            dx = ldexp(dx, -prec);
        }
    }

    fmpz_clear(t);
    flint_free(new_rel);
}

static void
rel_product(fmpz_t p, fmpz_t q, const ulong * primes, const slong * rel, slong len)
{
    slong i;

    if (len <= 4)
    {
        fmpz_t r;
        fmpz_init(r);

        for (i = 0; i < len; i++)
        {
            fmpz_ui_pow_ui(r, primes[i], FLINT_ABS(rel[i]));

            if (rel[i] >= 0)
                fmpz_mul(p, p, r);
            else
                fmpz_mul(q, q, r);
        }

        fmpz_clear(r);
    }
    else
    {
        fmpz_t p2, q2;

        fmpz_init_set_ui(p2, 1);
        fmpz_init_set_ui(q2, 1);

        rel_product(p, q, primes, rel, len / 2);
        rel_product(p2, q2, primes + len / 2, rel + len / 2, len - len / 2);

        fmpz_mul(p, p, p2);
        fmpz_mul(q, q, q2);

        fmpz_clear(p2);
        fmpz_clear(q2);
    }
}

/* todo: error propagation */
static void
_arb_exp_arf_precomp(arb_t res, const arf_t x, slong prec, int minus_one,
    slong num_logs, arb_srcptr logs, const ulong * primes,
    const float * weights,
    const short * log_rel_d,
    const double * epsilon, const double * epsilon_inv, double max_weight)
{
    arb_t t;
    fmpz_t p, r, q;
    slong wp;
    slong * rel;
    slong i;
    fmpz * alpha;
    slong mag;
    mag_t err, err2;

    mag = arf_abs_bound_lt_2exp_si(x);

    arb_init(t);

    rel = flint_malloc(num_logs * sizeof(slong));

    alpha = _fmpz_vec_init(num_logs);
    fmpz_init(r);

    if (prec <= 10000)
        wp = 256;
    else if (prec <= 100000)
        wp = 512;
    else
        wp = 768;

    for (i = 0; i < num_logs; i++)
        arf_get_fmpz_fixed_si(alpha + i, arb_midref(logs + i), -wp);

    arf_get_fmpz_fixed_si(r, x, -wp);
    _arb_log_reduce_fixed(rel, log_rel_d, epsilon, epsilon_inv,
        alpha, weights, num_logs, r, wp, max_weight);

    fmpz_clear(r);
    _fmpz_vec_clear(alpha, num_logs);

    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);
    if (minus_one && mag < 0)
        wp += (-mag);
    else if (mag > 0)
        wp += mag;

    arb_set_arf(t, x);
    arb_dot_si(t, t, 1, logs, 1, rel, 1, num_logs, wp);
    arb_exp_arf_generic(res, arb_midref(t), wp, 0);

    /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) */
    mag_init(err);
    mag_init(err2);
    arb_get_mag(err, res);
    mag_expm1(err2, arb_radref(t));
    mag_mul(arb_radref(res), err, err2);
    mag_clear(err);
    mag_clear(err2);

    fmpz_init(p);
    fmpz_init(q);

    fmpz_one(p);
    fmpz_one(q);
    rel_product(p, q, primes + 1, rel + 1, num_logs - 1);

    arb_mul_fmpz(res, res, p, wp);
    arb_div_fmpz(res, res, q, wp);
    arb_mul_2exp_si(res, res, rel[0]);

    if (minus_one)
        arb_sub_ui(res, res, 1, prec);
    else
        arb_set_round(res, res, prec);

    flint_free(rel);

    fmpz_clear(p);
    fmpz_clear(q);
    arb_clear(t);
}

void
arb_exp_arf_log_reduction(arb_t res, const arf_t x, slong prec, int minus_one)
{
    slong wp;
    slong mag;

    mag = arf_abs_bound_lt_2exp_si(x);

    if (mag < -prec / 16 || mag < -768 || arf_bits(x) < prec / 128)
    {
        arb_exp_arf_generic(res, x, prec, minus_one);
        return;
    }

    /* multiprecision log(2) reduction not implemented here */
    if ((FLINT_BITS == 32 && mag > 20) || (FLINT_BITS == 64 && mag > 40))
    {
        arb_exp_arf_huge(res, x, mag, prec, minus_one);
        return;
    }

    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);
    wp += FLINT_BITS;
    if (minus_one && mag < 0)
        wp += (-mag);
    else if (mag > 0)
        wp += mag;

    _arb_log_p_ensure_cached(wp);

    {
        const fixed_rel_struct * tab = fixed_rel_table(0, ARB_LOG_PRIME_CACHE_NUM);
        _arb_exp_arf_precomp(res, x, prec, minus_one,
            ARB_LOG_PRIME_CACHE_NUM,
            _arb_log_p_cache_vec(),
            tab->primes, tab->weights,
            tab->d, tab->epsilon, tab->epsilon_inv, prec);
    }
}

static void
gaussian_rel_product(fmpzi_t p, fmpzi_t q, const signed char * primes, const slong * rel, slong len)
{
    slong i;

    if (len <= 4)
    {
        fmpzi_t r;
        fmpzi_init(r);

        for (i = 0; i < len; i++)
        {
            fmpzi_set_si_si(r, primes[2 * i], primes[2 * i + 1]);
            fmpzi_pow_ui(r, r, FLINT_ABS(rel[i]));

            if (rel[i] >= 0)
                fmpzi_mul(p, p, r);
            else
                fmpzi_mul(q, q, r);
        }

        fmpzi_clear(r);
    }
    else
    {
        fmpzi_t p2, q2;

        fmpzi_init(p2);
        fmpzi_init(q2);

        fmpzi_one(p2);
        fmpzi_one(q2);

        gaussian_rel_product(p, q, primes, rel, len / 2);
        gaussian_rel_product(p2, q2, primes + 2 * (len / 2), rel + len / 2, len - len / 2);

        fmpzi_mul(p, p, p2);
        fmpzi_mul(q, q, q2);

        fmpzi_clear(p2);
        fmpzi_clear(q2);
    }
}

static void
_arb_sin_cos_arf_precomp(arb_t res1, arb_t res2, const arf_t x, slong prec,
    slong num_logs, arb_srcptr logs, const signed char * primes,
    const float * weights,
    const short * log_rel_d,
    const double * epsilon, const double * epsilon_inv, double max_weight)
{
    arb_t t;
    acb_t u, v;
    slong wp;
    slong * rel;
    slong i;
    fmpz * alpha;
    fmpz_t r;
    slong mag;

    arb_init(t);

    rel = flint_malloc(num_logs * sizeof(slong));

    alpha = _fmpz_vec_init(num_logs);
    fmpz_init(r);

    if (prec <= 10000)
        wp = 256;
    else if (prec <= 100000)
        wp = 512;
    else
        wp = 768;

    for (i = 0; i < num_logs; i++)
        arf_get_fmpz_fixed_si(alpha + i, arb_midref(logs + i), -wp);

    arf_get_fmpz_fixed_si(r, x, -wp);
    _arb_log_reduce_fixed(rel, log_rel_d, epsilon, epsilon_inv,
        alpha, weights, num_logs, r, wp, max_weight);

/*
    {
        slong i;
        for (i = 0; i < num_logs; i++)
            printf("%ld ", rel[i]);
        printf("\n");
    }
*/

    fmpz_clear(r);
    _fmpz_vec_clear(alpha, num_logs);

    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);

    mag = arf_abs_bound_lt_2exp_si(x);
    mag = FLINT_MAX(mag, 0);
    wp += mag;

    arb_set_arf(t, x);
    arb_dot_si(t, t, 1, logs, 1, rel, 1, num_logs, wp);

    acb_init(u);
    acb_init(v);

    arb_sin_cos_arf_generic(acb_imagref(u), acb_realref(u), arb_midref(t), wp);

    arb_add_error_mag(acb_imagref(u), arb_radref(t));
    arb_add_error_mag(acb_realref(u), arb_radref(t));

    {
        fmpzi_t p, r, s, q;

        fmpzi_init(p);
        fmpzi_init(q);
        fmpzi_init(r);
        fmpzi_init(s);

        fmpzi_one(p);
        fmpzi_one(q);

        gaussian_rel_product(p, q, primes + 2, rel + 1, num_logs - 1);

        fmpzi_conj(r, p);
        fmpzi_conj(s, q);

        fmpzi_mul(p, p, s);
        fmpzi_mul(q, q, r);

        /* printf("bits %ld %ld %ld\n", prec, _fmpz_vec_max_bits(p, 2), _fmpz_vec_max_bits(q, 2)); */

        arb_set_fmpz(acb_realref(v), fmpzi_realref(p));
        arb_set_fmpz(acb_imagref(v), fmpzi_imagref(p));
        acb_mul(u, u, v, wp);

        arb_set_fmpz(acb_realref(v), fmpzi_realref(q));
        arb_set_fmpz(acb_imagref(v), fmpzi_imagref(q));
        acb_div(u, u, v, wp);

        if ((rel[0] & 3) == 1)
            acb_mul_onei(u, u);
        else if ((rel[0] & 3) == 2)
            acb_neg(u, u);
        else if ((rel[0] & 3) == 3)
            acb_div_onei(u, u);

        fmpzi_clear(p);
        fmpzi_clear(q);
        fmpzi_clear(r);
        fmpzi_clear(s);
    }

    if (res1 != NULL)
        arb_set_round(res1, acb_imagref(u), prec);

    if (res2 != NULL)
        arb_set_round(res2, acb_realref(u), prec);

    flint_free(rel);

    arb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
arb_sin_cos_arf_atan_reduction(arb_t res1, arb_t res2, const arf_t x, slong prec)
{
    slong wp;
    slong mag;

    mag = arf_abs_bound_lt_2exp_si(x);

    if (mag < -prec / 16 || mag < -768 || arf_bits(x) < prec / 128)
    {
        arb_sin_cos_arf_generic(res1, res2, x, prec);
        return;
    }

    /* multiprecision pi reduction not implemented here */
    if ((FLINT_BITS == 32 && mag > 20) || (FLINT_BITS == 64 && mag > 40))
    {
        arb_sin_cos_arf_generic(res1, res2, x, prec);
        return;
    }

    wp = prec + 5 + 2 * FLINT_BIT_COUNT(prec);
    if (mag > 0)
        wp += mag;

    _arb_atan_gauss_p_ensure_cached(wp);

    {
        const fixed_rel_struct * tab = fixed_rel_table(1, ARB_ATAN_GAUSS_PRIME_CACHE_NUM);
        _arb_sin_cos_arf_precomp(res1, res2, x, prec,
            ARB_ATAN_GAUSS_PRIME_CACHE_NUM,
            _arb_atan_gauss_p_cache_vec(),
            _fixed_gaussian_primes, tab->weights,
            tab->d, tab->epsilon, tab->epsilon_inv, 0.5 * prec);
    }
}
