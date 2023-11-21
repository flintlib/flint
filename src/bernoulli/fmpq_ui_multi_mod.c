/*
    Copyright (C) 2021, 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "bernoulli.h"
#include "arb.h"

#define TIMING 0
#define DEBUG 0

typedef struct
{
    fmpz r;
    fmpz m;
}
crt_res_t;

typedef struct
{
    mp_srcptr residues;
    mp_srcptr primes;
}
crt_args_t;

static void
crt_init(crt_res_t * x, crt_args_t * args)
{
    fmpz_init(&x->r);
    fmpz_init(&x->m);
}

static void
crt_clear(crt_res_t * x, crt_args_t * args)
{
    fmpz_clear(&x->r);
    fmpz_clear(&x->m);
}

static void
_fmpz_crt_combine(fmpz_t r1r2, fmpz_t m1m2, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2)
{
    fmpz_invmod(m1m2, m1, m2);
    fmpz_mul(m1m2, m1m2, m1);
    fmpz_sub(r1r2, r2, r1);
    fmpz_mul(r1r2, r1r2, m1m2);
    fmpz_add(r1r2, r1r2, r1);
    fmpz_mul(m1m2, m1, m2);
    fmpz_mod(r1r2, r1r2, m1m2);
}

static void
crt_combine(crt_res_t * res, crt_res_t * left, crt_res_t * right, crt_args_t * args)
{
    _fmpz_crt_combine(&res->r, &res->m, &left->r, &left->m, &right->r, &right->m);
}

static void
crt_basecase(crt_res_t * res, slong a, slong b, crt_args_t * args)
{
    if (b - a == 0)
    {
        fmpz_zero(&res->r);
        fmpz_one(&res->m);
    }
    else if (b - a == 1)
    {
        fmpz_set_ui(&res->r, args->residues[a]);
        fmpz_set_ui(&res->m, args->primes[a]);
    }
    else
    {
        crt_res_t left, right;
        slong m = a + (b - a) / 2;

        crt_init(&left, args);
        crt_init(&right, args);

        crt_basecase(&left, a, m, args);
        crt_basecase(&right, m, b, args);
        crt_combine(res, &left, &right, args);

        crt_clear(&left, args);
        crt_clear(&right, args);
    }
}

/* todo: optimize basecase and move to flint */
void
_arb_tree_crt(fmpz_t r, fmpz_t m, mp_srcptr residues, mp_srcptr primes, slong len)
{
    crt_res_t res;
    crt_args_t args;

    res.r = *r;
    res.m = *m;

    args.residues = residues;
    args.primes = primes;

    flint_parallel_binary_splitting(&res,
        (bsplit_basecase_func_t) crt_basecase,
        (bsplit_merge_func_t) crt_combine,
        sizeof(crt_res_t),
        (bsplit_init_func_t) crt_init,
        (bsplit_clear_func_t) crt_clear,
        &args, 0, len, 20, -1, 0);

    *r = res.r;
    *m = res.m;

    return;
}

typedef struct
{
    ulong n;
    mp_ptr primes;
    mp_ptr residues;
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
    mp_ptr primes, residues;
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

    primes = flint_malloc(sizeof(mp_limb_t) * num_primes);
    residues = flint_malloc(sizeof(mp_limb_t) * num_primes);

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
    _arb_tree_crt(num, M, residues, primes, num_primes);
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

