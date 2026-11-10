/*
    Copyright (C) 2013-2014, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"

/* The logarithms of primes and the arguments of Gaussian primes are
   computed by mp_real (machin_bsplit.c: the Machin-type sets by binary
   splitting, Zuniga's series), which also keeps the shared caches of
   their floors for its argument reductions. */

void
arb_log_primes_vec_bsplit(arb_ptr res, slong n, slong prec)
{
    ulong * primes;
    n_primes_t iter;
    fmpz_t one, q;
    slong i, j, k, wp;

    if (n <= 0)
        return;

    if (n <= MP_REAL_LOG_PRIMES_MAX)
    {
        _mp_real_log_primes_vec_arb(res, n, prec);
        return;
    }

    /* beyond the tabulated sets, from the earlier primes: log p = 2
       atanh(1 / (2 p^2 - 1)) / 2 + (log((p - 1) / 2) + log((p + 1) /
       2)) / 2 + log 2, with (p +- 1) / 2 factored over the earlier
       primes */
    wp = prec + 64;
    _mp_real_log_primes_vec_arb(res, MP_REAL_LOG_PRIMES_MAX, wp);

    primes = flint_malloc(sizeof(ulong) * n);
    n_primes_init(iter);
    for (i = 0; i < n; i++)
        primes[i] = n_primes_next(iter);
    n_primes_clear(iter);

    fmpz_init_set_ui(one, 1);
    fmpz_init(q);

    for (i = MP_REAL_LOG_PRIMES_MAX; i < n; i++)
    {
        ulong pm[2];
        int s;

        fmpz_set_ui(q, primes[i]);
        fmpz_mul_ui(q, q, primes[i]);
        fmpz_mul_2exp(q, q, 1);
        fmpz_sub_ui(q, q, 1);
        arb_atan_frac_bsplit(res + i, one, q, 1, wp);
        arb_mul_2exp_si(res + i, res + i, 1);

        pm[0] = (primes[i] - 1) / 2;
        pm[1] = (primes[i] + 1) / 2;
        for (s = 0; s < 2; s++)
        {
            n_factor_t fac;
            n_factor_init(&fac);
            n_factor(&fac, pm[s], 1);
            for (j = 0; j < fac.num; j++)
                for (k = 0; k < i; k++)
                    if (fac.p[j] == primes[k])
                        arb_addmul_ui(res + i, res + k, fac.exp[j], wp);
        }

        arb_mul_2exp_si(res + i, res + i, -1);
        arb_add(res + i, res + i, res + 0, wp);
    }

    for (i = 0; i < n; i++)
        arb_set_round(res + i, res + i, prec);

    fmpz_clear(one);
    fmpz_clear(q);
    flint_free(primes);
}

void
arb_atan_gauss_primes_vec_bsplit(arb_ptr res, slong n, slong prec)
{
    if (n <= 0)
        return;
    if (n > MP_REAL_ATAN_GAUSS_MAX)
        flint_throw(FLINT_ERROR, "(%s): at most %d Gaussian primes\n", __func__, MP_REAL_ATAN_GAUSS_MAX);
    /* (mp_real gives 2 arg(pi_j)) */
    _mp_real_atan_gauss_vec_arb(res, n, prec + 1);
    _arb_vec_scalar_mul_2exp_si(res, res, n, -1);
}
