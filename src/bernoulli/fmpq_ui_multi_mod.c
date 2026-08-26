/*
    Copyright (C) 2021, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "bernoulli.h"
#include "bernoulli/impl.h"
#include "arb.h"

#define TIMING 0
#define DEBUG 0

typedef struct
{
    ulong n;
    nn_ptr primes;
    nn_ptr residues;
}
mod_p_param_t;

static void
mod_p_worker(slong i, void * param)
{
    mod_p_param_t * p = (mod_p_param_t *) param;

    p->residues[i] = bernoulli_mod_p_harvey(p->n, p->primes[i]);
}

void
_bernoulli_fmpq_ui_multi_mod(fmpz_t num, fmpz_t den, ulong n, double alpha)
{
    n_primes_t prime_iter;
    slong i, bits, mod_bits, zeta_bits, num_primes;
    ulong p;
    nn_ptr primes, residues;
    mag_t primes_product;
    fmpz_t M;
#if TIMING
    double t1, t2;
#endif

    if (n < 10 || n % 2 != 0)
    {
        _bernoulli_fmpq_ui_zeta(num, den, n);
        return;
    }

    if (alpha < 0)
    {
        if (n < 18000)
            alpha = 0.0;
        else if (n < 60000)
            alpha = 0.005 + 3.6e-6 * n;
        else
            alpha = FLINT_MIN(0.18 + 0.5e-6 * n, 0.28);
    }

#if TIMING
    t1 = clock();
#endif

    arith_bernoulli_number_denom(den, n);

    bits = arith_bernoulli_number_size(n) + fmpz_bits(den) + 2;
    mod_bits = bits * alpha;
    zeta_bits = bits - mod_bits;

    num_primes = 0;
    mag_init(primes_product);
    mag_one(primes_product);

    n_primes_init(prime_iter);

    p = 5;
    n_primes_jump_after(prime_iter, 5);

    for ( ; mag_cmp_2exp_si(primes_product, mod_bits) < 0; p = n_primes_next(prime_iter))
    {
        if (n % (p - 1) != 0)
        {
            mag_mul_ui_lower(primes_product, primes_product, p);
            num_primes++;
        }
    }

#if DEBUG
    flint_printf("\nn = %lu, bits = %lu, num_primes = %ld\n", n, bits, num_primes);
#endif

    primes = flint_malloc(sizeof(ulong) * num_primes);
    residues = flint_malloc(sizeof(ulong) * num_primes);

    p = 5;
    n_primes_jump_after(prime_iter, 5);

    for (i = 0; i < num_primes; p = n_primes_next(prime_iter))
    {
        if (n % (p - 1) != 0)
        {
            primes[i] = p;
            i++;
        }
    }

    n_primes_clear(prime_iter);

#if TIMING
    t2 = clock();
    flint_printf("init time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    flint_printf("num_primes = %ld\n", num_primes);
#endif

    {
        mod_p_param_t param;
        param.n = n;
        param.primes = primes;
        param.residues = residues;

        flint_parallel_do(mod_p_worker, &param, num_primes, 0, FLINT_PARALLEL_STRIDED /* | FLINT_PARALLEL_VERBOSE */);
    }

#if TIMING
    t2 = clock();
    flint_printf("mod time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    flint_printf("start CRT\n");
    t1 = clock();
#endif

    fmpz_init(M);
    fmpz_multi_CRT_ui_once(num, M, residues, primes, num_primes, 0);
    fmpz_mul(num, num, den);
    fmpz_mod(num, num, M);

    if (n % 4 == 0)
    {
        fmpz_sub(num, M, num);
        fmpz_neg(num, num);
    }

#if TIMING
    flint_printf("end CRT\n");
    t2 = clock();
    flint_printf("CRT time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
    t1 = clock();
#endif

    if (zeta_bits > 0)
    {
        slong prec;
        arb_t b;
        fmpz_t t;

        arb_init(b);
        fmpz_init(t);

        for (prec = zeta_bits + 10; ; prec += 32)
        {
            arb_bernoulli_ui_zeta(b, n, prec);
            arb_mul_fmpz(b, b, den, prec);
            arb_sub_fmpz(b, b, num, prec);
            arb_div_fmpz(b, b, M, prec);

            if (arb_get_unique_fmpz(t, b))
            {
                fmpz_addmul(num, t, M);
                break;
            }

            flint_printf("bernoulli: n = %wu, bits = %wd, mod = %wd, zeta = %wd: get_unique_fmpz failed!\n", n, bits, mod_bits, zeta_bits);
        }

        arb_clear(b);
        fmpz_clear(t);
    }

#if TIMING
    flint_printf("end zeta\n");
    t2 = clock();
    flint_printf("zeta time = %f\n", (t2 - t1) / (double) CLOCKS_PER_SEC);
#endif

    flint_free(primes);
    flint_free(residues);
    fmpz_clear(M);
    mag_clear(primes_product);
}
