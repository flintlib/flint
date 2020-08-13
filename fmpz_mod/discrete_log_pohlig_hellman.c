/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"
#include <math.h>

/*
    Assumption on fmpz_mod_discrete_log_pohlig_hellman_t:
        p is prime.
        The prime factors of p - 1 all fit a ulong.
    If the prime factors of p - 1 do not fit a ulong you do not want to calculate dlog mod p.

    so we have
    p - 1 = p1^e1 * ... * pn^en for ulong pi and ei

    The assumption p is prime could be removed, but then phi(p) needs to be calculated by someone somewhere.
*/

static int fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct_cmp(
    const fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct * lhs,
    const fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct * rhs)
{
    return fmpz_cmp(lhs->gammapow, rhs->gammapow);
}

void fmpz_mod_discrete_log_pohlig_hellman_init(fmpz_mod_discrete_log_pohlig_hellman_t L)
{
    fmpz_t two;
    fmpz_init_set_ui(two, 2);

    L->num_factors = 0;
    L->entries = NULL;
    fmpz_init(L->alpha);
    fmpz_init(L->alphainv);
    fmpz_init(L->pm1);
    fmpz_mod_ctx_init(L->fpctx, two);

    fmpz_clear(two);
    return;
}


void fmpz_mod_discrete_log_pohlig_hellman_clear(fmpz_mod_discrete_log_pohlig_hellman_t L)
{
    slong i;
    ulong c;
    fmpz_mod_discrete_log_pohlig_hellman_entry_struct * Li;

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        fmpz_clear(Li->idem);
        fmpz_clear(Li->co);
        fmpz_clear(Li->startinge);
        fmpz_clear(Li->startingbeta);
        fmpz_clear(Li->gamma);
        fmpz_clear(Li->gammainv);
        for (c = 0; c < Li->cbound; c++)
        {
            fmpz_clear(Li->table[c].gammapow);
        }
        flint_free(Li->table);
    }

    if (L->entries)
    {
        flint_free(L->entries);
    }

    fmpz_clear(L->alpha);
    fmpz_clear(L->alphainv);
    fmpz_clear(L->pm1);

    fmpz_mod_ctx_clear(L->fpctx);

    return;
}

static slong _pow_fmpz_cost(const fmpz_t pow)
{
    slong cost = fmpz_bits(pow) + fmpz_popcnt(pow) - 2;
    return FLINT_MAX(cost, 0);    
}

/*
    Assume that p is prime, don't check. Return an estimate on the number of
    multiplications need for one run.
*/
double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(
    fmpz_mod_discrete_log_pohlig_hellman_t L,
    const fmpz_t p)
{
    slong i;
    ulong c;
    fmpz_mod_discrete_log_pohlig_hellman_entry_struct * Li;
    fmpz_factor_t factors;
    fmpz_t temp;
    double total_cost;

    /* just clear it and allocate everything again */
    fmpz_mod_discrete_log_pohlig_hellman_clear(L);

    fmpz_init(L->alpha);
    fmpz_init(L->alphainv);
    fmpz_init(L->pm1);
    fmpz_mod_ctx_init(L->fpctx, p);

    fmpz_init(temp);

    fmpz_factor_init(factors);
    fmpz_sub_ui(L->pm1, p, 1);
    fmpz_factor(factors, L->pm1);
    L->num_factors = factors->num;
    L->entries = NULL;
    if (L->num_factors > 0)
    {
        L->entries = (fmpz_mod_discrete_log_pohlig_hellman_entry_struct*) flint_malloc(
                 L->num_factors*sizeof(fmpz_mod_discrete_log_pohlig_hellman_entry_struct));
    }
    for (i = 0; i < L->num_factors; i++)
    {
        fmpz_t pipow, recp;

        Li = L->entries + i;

        fmpz_init(Li->idem);
        fmpz_init(Li->co);
        fmpz_init(Li->startinge);
        fmpz_init(Li->startingbeta);
        fmpz_init(Li->gamma);
        fmpz_init(Li->gammainv);

        if (!fmpz_abs_fits_ui(factors->p + i))
        {
            /* L is corrupted */
            fmpz_clear(temp);
            fmpz_factor_clear(factors);
            flint_throw(FLINT_ERROR, "Exception in fmpz_mod_discrete_log_"
                  "pohlig_hellman_precompute_prime: Prime factor is large.\n");
        }
        Li->exp = factors->exp[i];
        Li->prime = fmpz_get_ui(factors->p + i);

        fmpz_init(recp);
        fmpz_init_set_ui(pipow, Li->prime);
        fmpz_pow_ui(pipow, pipow, Li->exp);
        fmpz_divexact(recp, L->pm1, pipow);
        fmpz_invmod(temp, recp, pipow);
        fmpz_mul(temp, temp, recp);

        fmpz_mod(Li->idem, temp, L->pm1);

        fmpz_set(Li->co, recp);
        fmpz_divexact_ui(Li->startinge, pipow, Li->prime);

        fmpz_clear(pipow);
        fmpz_clear(recp);
    }
    fmpz_factor_clear(factors);

    /* alpha will be a primitive root */
    fmpz_zero(L->alpha);
try_alpha:
    fmpz_add_ui(L->alpha, L->alpha, 1);
    if (fmpz_cmp(L->alpha, p) >= 0)
    {
        /* L is corrupted */
        fmpz_clear(temp);
        /* factors was already cleared */
        flint_throw(FLINT_ERROR, "Exception in fmpz_mod_discrete_log_pohlig_"
                   "hellman_precompute_prime: Could not find primitive root.");
    }
    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;
        fmpz_divexact_ui(temp, L->pm1, Li->prime);

        fmpz_mod_pow_fmpz(Li->gamma, L->alpha, temp, L->fpctx);

        if (fmpz_is_one(Li->gamma))
        {
            goto try_alpha;
        }
    }

    fmpz_mod_inv(L->alphainv, L->alpha, L->fpctx);

    for (i = 0; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        fmpz_mod_inv(Li->gammainv, Li->gamma, L->fpctx);
        fmpz_mod_pow_fmpz(Li->startingbeta, L->alphainv, Li->co, L->fpctx);

        Li->dbound = ceil(sqrt((double) Li->prime));
        Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        while (Li->cbound > 100)
        {
            Li->dbound *= 2;
            Li->cbound = (Li->prime + Li->dbound - 1)/Li->dbound;
        }

        FLINT_ASSERT(Li->dbound > 0);
        FLINT_ASSERT(Li->cbound > 0);
        Li->table = (fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct *)
                 flint_malloc(Li->cbound*sizeof(fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct));

        for (c = 0; c < Li->cbound; c++)
        {
            Li->table[c].cm = c*Li->dbound;
            fmpz_init(Li->table[c].gammapow);
            fmpz_mod_pow_ui(Li->table[c].gammapow, Li->gamma, Li->table[c].cm, L->fpctx);
        }
        qsort(Li->table, Li->cbound,
                sizeof(fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct),
               (int(*)(const void*, const void*)) fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct_cmp);
        for (c = 1; c < Li->cbound; c++)
        {
            FLINT_ASSERT(fmpz_cmp(Li->table[c - 1].gammapow,
                                    Li->table[c].gammapow) < 0);
        }
    }

    fmpz_clear(temp);

    total_cost = 0;
    for (i = 0; i < L->num_factors; i++)
    {
        double this_cost = 0;
        fmpz_t e;
        slong j;

        Li = L->entries + i;
        this_cost += _pow_fmpz_cost(Li->co);
        fmpz_init_set(e, Li->startinge);

        j = 0;
        do {
            this_cost += _pow_fmpz_cost(e);
            this_cost += Li->dbound*(1 + log(Li->cbound)); /* bsgs search */
            this_cost += 2*log(Li->prime); /* some power < Li->prime */
            fmpz_divexact_ui(e, e, Li->prime);
        } while (++j < Li->exp);
        total_cost += this_cost;
        fmpz_clear(e);
    }

    return total_cost;  

}

/* return with xx such that xx = alpha^y mod p, alpha is the p.r. L->alpha*/
void fmpz_mod_discrete_log_pohlig_hellman_run(
    fmpz_t xx,
    const fmpz_mod_discrete_log_pohlig_hellman_t L,
    const fmpz_t y)
{
    slong i, j;
    ulong g;
    fmpz_t x;
    fmpz_t pipow, e, acc;
    ulong lo, mid, hi, d;
    fmpz_t beta, z, w, temp;
    fmpz_mod_discrete_log_pohlig_hellman_entry_struct * Li;

    fmpz_init(x);
    fmpz_init(acc);
    fmpz_init(pipow);
    fmpz_init(e);
    fmpz_init(beta);
    fmpz_init(z);
    fmpz_init(w);
    fmpz_init(temp);

    FLINT_ASSERT(!fmpz_is_zero(y));
    FLINT_ASSERT(fmpz_mod_is_canonical(y, L->fpctx));

    i = 0;
    if (i < L->num_factors && L->entries[i].prime == 2)
    {
        Li = L->entries + i;
        FLINT_ASSERT(Li->prime == 2);

        fmpz_mod_pow_fmpz(z, y, Li->co, L->fpctx);
        fmpz_set(beta, Li->startingbeta);
        fmpz_set(e, Li->startinge);
        j = 0;
        fmpz_one(pipow); /* Li->prime^j */
        fmpz_zero(acc);
        do {
            fmpz_mod_pow_fmpz(w, z, e, L->fpctx);
            /* solve Li->gamma ^ g == w mod p */
            if (fmpz_is_one(w))
            {
                g = 0;
            }
            else
            {
                if (!fmpz_equal(w, Li->gamma))
                {
                    goto cleanup_and_throw;
                }
                g = 1;
                fmpz_mod_mul(z, z, beta, L->fpctx);
            }
            fmpz_mod_mul(beta, beta, beta, L->fpctx);
            fmpz_addmul_ui(acc, pipow, g);
            fmpz_mul_2exp(pipow, pipow, 1);
            fmpz_tdiv_q_2exp(e, e, 1);
        } while (++j < Li->exp);

        fmpz_addmul(x, acc, Li->idem);
        i = 1;
    }

    for (; i < L->num_factors; i++)
    {
        Li = L->entries + i;

        fmpz_mod_pow_fmpz(z, y, Li->co, L->fpctx);
        fmpz_set(beta, Li->startingbeta);
        fmpz_set(e, Li->startinge);
        j = 0;
        fmpz_one(pipow); /* Li->prime^j */
        fmpz_zero(acc);
        do {
            fmpz_mod_pow_fmpz(w, z, e, L->fpctx);
            /* solve Li->gamma ^ g == w mod p */
            d = 0;
            while (1)
            {
                lo = 0; hi = Li->cbound;
                while (hi - lo > 4)
                {
                    int cmp;
                    mid = lo + (hi - lo)/2;
                    cmp = fmpz_cmp(Li->table[mid].gammapow, w);
                    if (cmp == 0)
                    {
                        g = Li->table[mid].cm + d;
                        goto found_g;
                    }
                    else if (cmp > 0)
                    {
                        hi = mid;
                    }
                    else
                    {
                        lo = mid;
                    }
                }
                while (lo < hi)
                {
                    if (fmpz_equal(Li->table[lo].gammapow, w))
                    {
                        g = Li->table[lo].cm + d;
                        goto found_g;
                    }
                    lo++;
                }
                fmpz_mod_mul(w, w, Li->gammainv, L->fpctx);
                d++;
                if (d >= Li->dbound)
                {
                    goto cleanup_and_throw;
                }
            }
        found_g:
            FLINT_ASSERT(g < Li->prime);
            fmpz_mod_pow_ui(temp, beta, g, L->fpctx);
            fmpz_mod_mul(z, z, temp, L->fpctx);
            fmpz_mod_pow_ui(beta, beta, Li->prime, L->fpctx);
            fmpz_addmul_ui(acc, pipow, g);
            fmpz_mul_ui(pipow, pipow, Li->prime);
            fmpz_divexact_ui(e, e, Li->prime);
        } while (++j < Li->exp);

        fmpz_addmul(x, acc, Li->idem);
    }

    fmpz_mod(xx, x, L->pm1);
    fmpz_clear(acc);
    fmpz_clear(pipow);
    fmpz_clear(e);
    fmpz_clear(beta);
    fmpz_clear(z);
    fmpz_clear(w);
    fmpz_clear(temp);
    fmpz_clear(x);
    return;

cleanup_and_throw:

    fmpz_clear(acc);
    fmpz_clear(pipow);
    fmpz_clear(e);
    fmpz_clear(beta);
    fmpz_clear(z);
    fmpz_clear(w);
    fmpz_clear(temp);
    fmpz_clear(x);
    flint_throw(FLINT_ERROR, "Exception in fmpz_mod_discrete_log_pohlig_"
                                          "hellman_run: Could not find log.");
    return;


}
