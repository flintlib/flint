/*
    Copyright (C) 2013-2014 Fredrik Johansson
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "thread_support.h"
#include "ulong_extras.h"
#include "arb.h"
#include "fmpz_vec.h"
#include "fixed.h"
#include "arb/impl.h"

/* The Machin-type sets live in the fixed module (fixed/machin_tab.c,
   fixed_machin_table); the Gaussian primes in fixed/atan_gauss.c. */

typedef struct
{
    const fixed_machin_struct * tab;
    arb_ptr res;
    slong prec;
    int hyperbolic;
}
atan_work;

static void
parallel_atan_worker(slong i, atan_work * work)
{
    fmpz_t p, q;

    fmpz_init(p);
    fmpz_init(q);

    fmpz_one(p);
    fixed_machin_get_x(q, work->tab, i);

    arb_atan_frac_bsplit(work->res + i, p, q, work->hyperbolic, work->prec);

    fmpz_clear(p);
    fmpz_clear(q);
}

void
arb_log_primes_vec_bsplit(arb_ptr res, slong n, slong prec)
{
    ulong den, prime;
    ulong * primes;
    fmpz * crow;
    const fixed_machin_struct * mt;
    slong i, j, k, wp, ln;
    arb_ptr y;
    arb_t t;
    fmpz_t p, q;
    n_primes_t iter;

    wp = prec + 64;

    mt = fixed_machin_table(0, n);
    ln = mt->num; den = mt->den;
    crow = _fmpz_vec_init(ln);

    y = _arb_vec_init(ln);
    arb_init(t);
    fmpz_init(p);
    fmpz_init(q);

    primes = flint_malloc(sizeof(ulong) * n);
    n_primes_init(iter);
    for (i = 0; i < n; i++)
        primes[i] = n_primes_next(iter);
    n_primes_clear(iter);

    {
        atan_work work;

        work.tab = mt;
        work.res = y;
        work.prec = wp;
        work.hyperbolic = 1;
        flint_parallel_do((do_func_t) parallel_atan_worker, &work, ln, -1, FLINT_PARALLEL_STRIDED);
    }

    for (i = 0; i < FLINT_MIN(n, ln); i++)
    {
        fixed_machin_get_c_row(crow, mt, i);
        arb_dot_fmpz(res + i, NULL, 0, y, 1, crow, 1, ln, wp);
        if (den == 1)
            arb_set_round(res + i, res + i, prec);
        else
            arb_div_ui(res + i, res + i, den, prec);
    }

    /* todo: sieving instead of factoring */
    /* todo: parallel comp here as well */
    for (i = ln; i < n; i++)
    {
        n_factor_t fac;

        prime = primes[i];

        fmpz_one(p);
        fmpz_set_ui(q, 2 * prime * prime - 1);

        arb_atan_frac_bsplit(res + i, p, q, 1, wp);
        arb_mul_2exp_si(res + i, res + i, 1);

        n_factor_init(&fac);
        n_factor(&fac, (prime - 1) / 2, 1);

        for (j = 0; j < fac.num; j++)
            for (k = 0; k < i; k++)
                if (fac.p[j] == primes[k])
                    arb_addmul_ui(res + i, res + k, fac.exp[j], wp);

        n_factor_init(&fac);
        n_factor(&fac, (primes[i] + 1) / 2, 1);

        for (j = 0; j < fac.num; j++)
            for (k = 0; k < i; k++)
                if (fac.p[j] == primes[k])
                    arb_addmul_ui(res + i, res + k, fac.exp[j], wp);

        arb_mul_2exp_si(res + i, res + i, -1);

        arb_add(res + i, res + i, res + 0, prec);
    }

    _arb_vec_clear(y, ln);
    _fmpz_vec_clear(crow, ln);
    arb_clear(t);
    fmpz_clear(p);
    fmpz_clear(q);
    flint_free(primes);
}

FLINT_TLS_PREFIX arb_struct _arb_log_p_cache[ARB_LOG_PRIME_CACHE_NUM];
FLINT_TLS_PREFIX slong _arb_log_p_cache_prec = 0;

static void _arb_log_p_cleanup(void)
{
    slong i;
    for (i = 0; i < ARB_LOG_PRIME_CACHE_NUM; i++)
        arb_clear(_arb_log_p_cache + i);
    _arb_log_p_cache_prec = 0;
}

arb_srcptr _arb_log_p_cache_vec(void)
{
    return _arb_log_p_cache;
}


void _arb_log_p_ensure_cached(slong prec)
{
    slong i, wp;

    if (_arb_log_p_cache_prec < prec)
    {
        if (_arb_log_p_cache_prec == 0)
        {
            for (i = 0; i < ARB_LOG_PRIME_CACHE_NUM; i++)
                arb_init(_arb_log_p_cache + i);

            flint_register_cleanup_function(_arb_log_p_cleanup);
        }

        wp = prec + 32;

        if (wp <= ARB_LOG_TAB2_PREC - 16)
        {
            for (i = 0; i < ARB_LOG_PRIME_CACHE_NUM; i++)
            {
                slong exp, exp_fix;
                slong n;
                arb_ptr res = _arb_log_p_cache + i;

                n = ARB_LOG_TAB2_PREC / FLINT_BITS;

                /* exponent of log(prime(i+1)) */
                exp = (i >= 1) + (i >= 4) + (i >= 16) + (i >= 429);

                /* just reading the table is known to give the correct rounding */
                _arf_set_round_mpn(arb_midref(res), &exp_fix, arb_log_p_tab[i], n, 0, wp, ARF_RND_NEAR);
                exp += exp_fix;
                _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp);

                /* 1/2 ulp error */
                _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), exp - wp);
                MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
            }
        }
        else
        {
            prec = FLINT_MAX(prec, _arb_log_p_cache_prec * 1.25);

            arb_log_primes_vec_bsplit(_arb_log_p_cache, ARB_LOG_PRIME_CACHE_NUM, prec + 32);
        }

        _arb_log_p_cache_prec = prec;
    }
}

static int factor_smooth(ulong * c, ulong n)
{
    slong i;

    for (i = 0; i < ARB_LOG_PRIME_CACHE_NUM; i++)
        c[i] = 0;

    /* Hardcoded so that the compiler can remove divisions */
    while (n != 1 && n % 2 == 0) { n /= 2; c[0]++; }
    while (n != 1 && n % 3 == 0) { n /= 3; c[1]++; }
    while (n != 1 && n % 5 == 0) { n /= 5; c[2]++; }
    while (n != 1 && n % 7 == 0) { n /= 7; c[3]++; }
    while (n != 1 && n % 11 == 0) { n /= 11; c[4]++; }
    while (n != 1 && n % 13 == 0) { n /= 13; c[5]++; }
    while (n != 1 && n % 17 == 0) { n /= 17; c[6]++; }
    while (n != 1 && n % 19 == 0) { n /= 19; c[7]++; }
    while (n != 1 && n % 23 == 0) { n /= 23; c[8]++; }
    while (n != 1 && n % 29 == 0) { n /= 29; c[9]++; }
    while (n != 1 && n % 31 == 0) { n /= 31; c[10]++; }
    while (n != 1 && n % 37 == 0) { n /= 37; c[11]++; }
    while (n != 1 && n % 41 == 0) { n /= 41; c[12]++; }

    return n == 1;
}

/* todo: use in log_ui in appropriates ranges */
int
_arb_log_ui_smooth(arb_t res, ulong n, slong prec)
{
    ulong c[ARB_LOG_PRIME_CACHE_NUM];

    if (factor_smooth(c, n))
    {
        _arb_log_p_ensure_cached(prec);
        arb_dot_ui(res, NULL, 0, _arb_log_p_cache, 1, c, 1, ARB_LOG_PRIME_CACHE_NUM, prec);
        return 1;
    }
    else
    {
        return 0;
    }
}

void
arb_atan_gauss_primes_vec_bsplit(arb_ptr res, slong n, slong prec)
{
    fmpz * crow;
    const fixed_machin_struct * mt;
    slong i, j, wp, ln;
    arb_ptr y;
    arb_t t;
    fmpz_t p, q;
    ulong den;

    /* not implemented */
    if (n > 64)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    wp = prec + 64;

    mt = fixed_machin_table(1, n);
    ln = mt->num; den = mt->den;
    crow = _fmpz_vec_init(ln);

    y = _arb_vec_init(ln);
    arb_init(t);
    fmpz_init(p);
    fmpz_init(q);

    {
        atan_work work;

        work.tab = mt;
        work.res = y;
        work.prec = wp;
        work.hyperbolic = 0;
        flint_parallel_do((do_func_t) parallel_atan_worker, &work, ln, -1, FLINT_PARALLEL_STRIDED);
    }

    for (i = 0; i < FLINT_MIN(n, ln); i++)
    {
        fixed_machin_get_c_row(crow, mt, i);
        arb_dot_fmpz(res + i, NULL, 0, y, 1, crow, 1, ln, wp);
        if (den == 1)
            arb_set_round(res + i, res + i, prec);
        else
            arb_div_ui(res + i, res + i, den, prec);
    }

    for (i = ln; i < n; i++)
    {
        double best = 100, t;
        slong xa, xb, ya, yb;
        slong best_j = 0;

        xa = _fixed_gaussian_primes[2 * i];
        xb = _fixed_gaussian_primes[2 * i + 1];

        for (j = 0; j < FLINT_MIN(i, 100); j++)
        {
            ya = _fixed_gaussian_primes[2 * j];
            yb = _fixed_gaussian_primes[2 * j + 1];

            t = (xb*ya - xa*yb) / (double) (xa*ya + xb*yb);

            if (fabs(t) < best)
            {
                best = fabs(t);
                best_j = j;
            }
        }

        ya = _fixed_gaussian_primes[2 * best_j];
        yb = _fixed_gaussian_primes[2 * best_j + 1];

        fmpz_set_si(p, xb*ya - xa*yb);
        fmpz_set_si(q, xa*ya + xb*yb);

        arb_atan_frac_bsplit(res + i, p, q, 0, wp);
        arb_add(res + i, res + i, res + best_j, prec);
    }

    _arb_vec_clear(y, ln);
    _fmpz_vec_clear(crow, ln);
    arb_clear(t);
    fmpz_clear(p);
    fmpz_clear(q);
}

FLINT_TLS_PREFIX arb_struct _arb_atan_gauss_p_cache[ARB_ATAN_GAUSS_PRIME_CACHE_NUM];
FLINT_TLS_PREFIX slong _arb_atan_gauss_p_cache_prec = 0;

static void _arb_atan_gauss_p_cleanup(void)
{
    slong i;
    for (i = 0; i < ARB_ATAN_GAUSS_PRIME_CACHE_NUM; i++)
        arb_clear(_arb_atan_gauss_p_cache + i);
    _arb_atan_gauss_p_cache_prec = 0;
}

arb_srcptr _arb_atan_gauss_p_cache_vec(void)
{
    return _arb_atan_gauss_p_cache;
}

void _arb_atan_gauss_p_ensure_cached(slong prec)
{
    slong i, wp;

    if (_arb_atan_gauss_p_cache_prec < prec)
    {
        if (_arb_atan_gauss_p_cache_prec == 0)
        {
            for (i = 0; i < ARB_ATAN_GAUSS_PRIME_CACHE_NUM; i++)
                arb_init(_arb_atan_gauss_p_cache + i);

            flint_register_cleanup_function(_arb_atan_gauss_p_cleanup);
        }

        wp = prec + 32;

        /* todo */
        if (wp <= ARB_ATAN_TAB2_PREC - 16)
        {
            for (i = 0; i < ARB_ATAN_GAUSS_PRIME_CACHE_NUM; i++)
            {
                slong exp, exp_fix;
                slong n;
                static const char exponents[24] = {0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1};
                arb_ptr res = _arb_atan_gauss_p_cache + i;

                n = ARB_LOG_TAB2_PREC / FLINT_BITS;

                if (i >= 24)
                    flint_throw(FLINT_ERROR, "(%s)\n", __func__);
                /* exponent of 2*atan(x) */
                exp = exponents[i] + 1;

                /* just reading the table is known to give the correct rounding */
                _arf_set_round_mpn(arb_midref(res), &exp_fix, arb_atan_gauss_tab[i], n, 0, wp, ARF_RND_NEAR);
                exp += exp_fix;
                _fmpz_set_si_small(ARF_EXPREF(arb_midref(res)), exp);

                /* 1/2 ulp error */
                _fmpz_set_si_small(MAG_EXPREF(arb_radref(res)), exp - wp);
                MAG_MAN(arb_radref(res)) = MAG_ONE_HALF;
            }
        }
        else
        {
            prec = FLINT_MAX(prec, _arb_atan_gauss_p_cache_prec * 1.25);

            arb_atan_gauss_primes_vec_bsplit(_arb_atan_gauss_p_cache, ARB_ATAN_GAUSS_PRIME_CACHE_NUM, prec + 32);
            _arb_vec_scalar_mul_2exp_si(_arb_atan_gauss_p_cache, _arb_atan_gauss_p_cache, ARB_ATAN_GAUSS_PRIME_CACHE_NUM, 1);
        }

        _arb_atan_gauss_p_cache_prec = prec;
    }
}
