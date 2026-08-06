/*
    Copyright (C) 2011, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "nmod.h"
#include "fmpz.h"
#include "arith.h"
#include "thread_support.h"

static ulong
arith_bell_number_nmod2(const unsigned int * divtab, ulong n, ulong p)
{
    ulong s, t, u, v, s2, s1, s0, t1, t0;
    ulong qq[3];
    slong i;
    nn_ptr facs, pows;
    ulong one_red, i_red;
    nmod_t mod;
    nmod_redc_ctx_t ctx;

    nmod_init(&mod, p);
    nmod_redc_ctx_init_nmod(ctx, mod);

    TMP_INIT;
    TMP_START;

    facs = TMP_ALLOC(2 * (n + 1) * sizeof(ulong));
    pows = facs + (n + 1);

    /* Compute inverse factorials */
    /* We actually compute (n! / i!) and divide out (n!)^2 at the end */
    one_red = nmod_redc_set_nmod(1, ctx);

    facs[n] = one_red;
    i_red = nmod_redc_set_nmod(n, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        facs[i] = nmod_redc_fast_mul(facs[i + 1], i_red, ctx);
        i_red = nmod_redc_fast_sub(i_red, one_red, ctx);
    }

    /* Compute powers */
    pows[0] = 0;
    pows[1] = i_red = one_red;
    for (i = 2; i <= n; i++)
    {
        i_red = nmod_redc_fast_add(i_red, one_red, ctx);
        if (divtab[2 * i] == 1)
            pows[i] = _nmod_redc_fast_pow_ui(i_red, n, ctx);
        else
            pows[i] = nmod_redc_fast_mul(pows[divtab[2 * i]], pows[divtab[2 * i + 1]], ctx);
    }

    s2 = s1 = s0 = 0;

    for (t = i = 0; i <= n; i++)
    {
        if (i % 2 == 0)
            t = nmod_redc_fast_add(t, facs[i], ctx);
        else
            t = nmod_redc_fast_sub(t, facs[i], ctx);

        u = pows[n - i];
        v = facs[n - i];
        u = nmod_redc_fast_mul(u, v, ctx);

        // u = nmod_redc_fast_mul(u, t, ctx);
        // s = nmod_redc_fast_add(s, u, ctx);
        // Faster: introduce an extra factor R: */
        umul_ppmm(t1, t0, u, t);
        add_sssaaaaaa(s2, s1, s0, s2, s1, s0, 0, t1, t0);
    }

    qq[2] = s2;
    qq[1] = s1;
    qq[0] = s0;
    s = mpn_mod_1(qq, 3, mod.n);  // could use a preinv method

    /* Twice to remove the extra factor R introduced in the main sum. */
    s = nmod_redc_get_nmod(s, ctx);
    s = nmod_redc_get_nmod(s, ctx);

    /* Remove (n!)^2 */
    u = nmod_redc_get_nmod(facs[0], ctx);
    u = nmod_inv(u, mod);
    u = nmod_mul(u, u, mod);
    s = nmod_mul(s, u, mod);

    TMP_END;

    return s;
}

static void
divisor_table(unsigned int * tab, slong len)
{
    slong i, j;

    for (i = 0; i < len; i++)
    {
        tab[2 * i] = 1;
        tab[2 * i + 1] = i;
    }

    for (i = 2; i < len; i++)
    {
        for (j = 2; j <= i && i * j < len; j++)
        {
            tab[2 * i * j] = j;
            tab[2 * i * j + 1] = i;
        }
    }
}

typedef struct
{
    unsigned int * divtab;
    nn_ptr primes;
    nn_ptr residues;
    ulong n;
}
work_t;

static void
worker(slong k, void * _work)
{
    work_t * work = (work_t *) _work;
    work->residues[k] = arith_bell_number_nmod2(work->divtab, work->n, work->primes[k]);
}

void
arith_bell_number_multi_mod(fmpz_t res, ulong n)
{
    fmpz_comb_temp_t temp;
    fmpz_comb_t comb;
    slong k, num_primes;
    flint_bitcnt_t size, prime_bits;
    work_t work;

    if (n <= 1)
    {
        fmpz_one(res);
        return;
    }

    size = arith_bell_number_size(n) + 1;
    prime_bits = FLINT_BITS - 3;
    num_primes = (size + prime_bits - 1) / prime_bits;

    work.n = n;
    work.primes = flint_malloc(num_primes * sizeof(ulong));
    work.residues = flint_malloc(num_primes * sizeof(ulong));
    work.divtab = flint_malloc(2 * sizeof(unsigned int) * (n + 1));

    divisor_table(work.divtab, n + 1);

    work.primes[0] = n_nextprime(UWORD(1) << prime_bits, 0);
    for (k = 1; k < num_primes; k++)
        work.primes[k] = n_nextprime(work.primes[k-1], 0);

    flint_parallel_do((do_func_t) worker, &work, num_primes, -1, FLINT_PARALLEL_UNIFORM);

    fmpz_comb_init(comb, work.primes, num_primes);
    fmpz_comb_temp_init(temp, comb);
    fmpz_multi_CRT_ui(res, work.residues, comb, temp, 0);
    fmpz_comb_clear(comb);
    fmpz_comb_temp_clear(temp);

    flint_free(work.primes);
    flint_free(work.residues);
    flint_free(work.divtab);
}
