/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_vec.h"
#include "fmpz.h"
#include <math.h>


static int nmod_discrete_log_pohlig_hellman_table_entry_struct_cmp(
    const nmod_discrete_log_pohlig_hellman_table_entry_struct * lhs,
    const nmod_discrete_log_pohlig_hellman_table_entry_struct * rhs)
{
    return (lhs->gammapow < rhs->gammapow) ? -1 : (lhs->gammapow > rhs->gammapow);
}

void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)
{
    L->num_factors = 0;
    L->entries = NULL;
    nmod_init(&L->mod, 2);
}

void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)
{
    slong i;
    nmod_discrete_log_pohlig_hellman_entry_struct * Li;

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        flint_free(Li->table);
    }

    if (L->entries)
    {
        flint_free(L->entries);
    }
}


static slong _pow_ui_cost(ulong pow)
{
    slong cost = -2;
    while (pow != 0)
    {
        cost += 1 + (pow&1);
        pow = pow/2;
    }
    return FLINT_MAX(cost, 0);
}

/*
    Assume that p is prime, don't check. Return an estimate on the number of
    multiplications need for one run.
*/
double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, mp_limb_t p)
{
    slong i;
    ulong c;
    nmod_discrete_log_pohlig_hellman_entry_struct * Li;
    n_factor_t factors;
    double total_cost;

    /* just free everything and allocate again for now */
    nmod_discrete_log_pohlig_hellman_clear(L);

    n_factor_init(&factors);
    n_factor(&factors, p - 1, 1);

    nmod_init(&L->mod, p);
    L->entries = NULL;
    L->num_factors = factors.num;
    if (L->num_factors > 0)
    {
        L->entries = (nmod_discrete_log_pohlig_hellman_entry_struct*) flint_malloc(
                     L->num_factors*sizeof(nmod_discrete_log_pohlig_hellman_entry_struct));
    }

    for (i = 0; i < L->num_factors; i++)
    {
        fmpz_t pipow, pm1, temp, recp;

        Li = L->entries + i;

        Li->exp = factors.exp[i];
        Li->prime = factors.p[i];

        fmpz_init(recp);
        fmpz_init(temp);
        fmpz_init_set_ui(pipow, Li->prime);
        fmpz_pow_ui(pipow, pipow, Li->exp);
        fmpz_init_set_ui(pm1, p - 1);
        fmpz_divexact(recp, pm1, pipow);
        fmpz_invmod(temp, recp, pipow);
        fmpz_mul(temp, temp, recp);

        Li->idem = fmpz_fdiv_ui(temp, p - 1);

        Li->co = fmpz_get_ui(recp);
        Li->startinge = fmpz_get_ui(pipow)/Li->prime;

        fmpz_clear(pipow);
        fmpz_clear(pm1);
        fmpz_clear(temp);
        fmpz_clear(recp);
    }

    /* alpha will be a primitive root */
    L->alpha = 0;
try_alpha:
    L->alpha++;
    if (L->alpha >= p)
    {
        /* L is corrupted */
        /* factors appears to need no cleanup */
        flint_throw(FLINT_ERROR, "Exception in nmod_discrete_log_pohlig"
                  "_hellman_precompute_prime: Could not find primitive root.");
    }
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        Li->gamma = nmod_pow_ui(L->alpha, (p - 1) / Li->prime, L->mod);
        if (Li->gamma == 1)
        {
            goto try_alpha;
        }
    }

    L->alphainv = nmod_inv(L->alpha, L->mod);

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        Li->gammainv = nmod_inv(Li->gamma, L->mod);

        Li->startingbeta = nmod_pow_ui(L->alphainv, Li->co, L->mod);

        Li->dbound = ceil(sqrt((double) Li->prime));
        Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        while (Li->cbound > 100)
        {
            Li->dbound *= 2;
            Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        }

        FLINT_ASSERT(Li->dbound > 0);
        FLINT_ASSERT(Li->cbound > 0);
        Li->table = (nmod_discrete_log_pohlig_hellman_table_entry_struct *) flint_malloc(
                   Li->cbound*sizeof(nmod_discrete_log_pohlig_hellman_table_entry_struct));

        for (c = 0; c < Li->cbound; c++)
        {
            Li->table[c].cm = c*Li->dbound;
            Li->table[c].gammapow = nmod_pow_ui(Li->gamma, Li->table[c].cm, L->mod);
        }
        qsort(Li->table, Li->cbound,
                sizeof(nmod_discrete_log_pohlig_hellman_table_entry_struct),
               (int(*)(const void*, const void*))
                                  nmod_discrete_log_pohlig_hellman_table_entry_struct_cmp);
        for (c = 1; c < Li->cbound; c++)
        {
            FLINT_ASSERT(Li->table[c - 1].gammapow < Li->table[c].gammapow);
            FLINT_ASSERT(Li->table[c].gammapow
                          == nmod_pow_ui(Li->gamma, Li->table[c].cm, L->mod));
        }
    }

    total_cost = 0;
    for (i = 0; i < L->num_factors; i++)
    {
        double this_cost = 0;
        ulong e;
        slong j;
        Li = L->entries + i;
        this_cost += _pow_ui_cost(Li->co);
        e = Li->startinge;
        j = 0;
        do {
            this_cost += _pow_ui_cost(e);
            this_cost += Li->dbound*(1 + log(Li->cbound)); /* bsgs search */
            this_cost += 2*log(Li->prime); /* some power < Li->prime */
            e = e / Li->prime;
        } while (++j < Li->exp);
        total_cost += this_cost;
    }

    return total_cost;  
}

/* return x such that y = alpha^x mod p, alpha is the p.r. L->alpha*/
ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, mp_limb_t y)
{
    slong i, j;
    ulong x, q, r, e, x0 = 0, x1 = 0, x2 = 0, pp0, pp1, acc, g, pipow;
    ulong lo, mid, hi, d;
    mp_limb_t beta, z, w;
    nmod_discrete_log_pohlig_hellman_entry_struct * Li;

    FLINT_ASSERT(y != 0);
    FLINT_ASSERT(y < L->mod.n);

    i = 0;
    if (i < L->num_factors && L->entries[i].prime == 2)
    {
        Li = L->entries + i;
        FLINT_ASSERT(Li->prime == 2);

        z = nmod_pow_ui(y, Li->co, L->mod);
        beta = Li->startingbeta;
        e = Li->startinge;
        j = 0;
        pipow = 1; /* Li->prime^j */
        acc = 0;
        do {
            w = nmod_pow_ui(z, e, L->mod);
            /* solve Li->gamma ^ g == w mod p */

            if (w == 1)
            {
                g = 0;
            }
            else
            {
                if (w != Li->gamma)
                {
                    goto cleanup_and_throw;
                }
                g = 1;
                z = nmod_mul(z, beta, L->mod);
            }
            beta = nmod_mul(beta, beta, L->mod);
            acc += g*pipow;
            pipow = pipow*2;
            e = e / 2;
        } while (++j < Li->exp);

        umul_ppmm(pp1, pp0, acc, Li->idem);
        add_sssaaaaaa(x2, x1, x0, x2, x1, x0, WORD(0), pp1, pp0);
        i = 1;
    }

    for (; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        z = nmod_pow_ui(y, Li->co, L->mod);
        beta = Li->startingbeta;
        e = Li->startinge;
        j = 0;
        pipow = 1; /* Li->prime^j */
        acc = 0;
        do {
            w = nmod_pow_ui(z, e, L->mod);
            /* solve Li->gamma ^ g == w mod p */
            d = 0;
            while (1)
            {
                lo = 0; hi = Li->cbound;
                while (hi - lo > 4)
                {
                    mid = lo + (hi - lo)/2;
                    if (Li->table[mid].gammapow == w)
                    {
                        g = Li->table[mid].cm + d;
                        goto found_g;
                    }
                    if (Li->table[mid].gammapow > w)
                        hi = mid;
                    else
                        lo = mid;
                }
                while (lo < hi)
                {
                    if (Li->table[lo].gammapow == w)
                    {
                        g = Li->table[lo].cm + d;
                        goto found_g;
                    }
                    lo++;
                }
                w = nmod_mul(w, Li->gammainv, L->mod);
                d++;
                if (d >= Li->dbound)
                {
                    goto cleanup_and_throw;
                }
            }
        found_g:
            FLINT_ASSERT(g < Li->prime);
            z = nmod_mul(z, nmod_pow_ui(beta, g, L->mod), L->mod);
            beta = nmod_pow_ui(beta, Li->prime, L->mod);
            acc += g*pipow;
            pipow = pipow*Li->prime;
            e = e / Li->prime;
        } while (++j < Li->exp);

        umul_ppmm(pp1, pp0, acc, Li->idem);
        add_sssaaaaaa(x2, x1, x0, x2, x1, x0, WORD(0), pp1, pp0);
    }

    udiv_qrnnd(q, r, x2, x1, L->mod.n - 1);
    udiv_qrnnd(q, x, r, x0, L->mod.n - 1);
    return x;

cleanup_and_throw:

    /* nothing currently to cleanup */
    flint_throw(FLINT_ERROR, "Exception in nmod_discrete_log_pohlig"
                                         "_hellman_run: Could not find log.");
    return 0;
}
