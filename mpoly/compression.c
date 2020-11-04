/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "long_extras.h"

/*
    input is S and the sizes n, l
    outputs S', V, D, deg satisfy

    det(V) = +- 1
    V.S'[i] + D = S[i]
    deg[j] = max_j S[i][j]
    0 = min_j S[i][j]

    where S'[i] / S[i] is the i^th row of S' / S
*/
slong _mpoly_compress_exps(
    slong * V,      /* n*n matrix */
    slong * D,      /* n vector */
    slong * deg,    /* n vector */
    slong * S,      /* lxn matrix */
    slong n,
    slong l)
{
    slong mind, maxd, tot;
    slong * minp = FLINT_ARRAY_ALLOC(n, slong);
    slong * maxp = FLINT_ARRAY_ALLOC(n, slong);
    slong * minn = FLINT_ARRAY_ALLOC(n, slong);
    slong * maxn = FLINT_ARRAY_ALLOC(n, slong);
    slong * perm = FLINT_ARRAY_ALLOC(n, slong);
    slong * tmp = maxn;
    slong i, j, k, m;
    slong best_loc_i, best_loc_j, best_min, best_deg;
    fmpz_t this_prod, best_prod;

    fmpz_init(this_prod);
    fmpz_init(best_prod);

    FLINT_ASSERT(n > 0);
    FLINT_ASSERT(l > 1);

    for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
        V[i*n + j] = (i == j);

    for (j = 0; j < n; j++)
        minp[j] = maxp[j] = S[0*n + j];
    for (i = 1; i < l; i++)
    for (j = 0; j < n; j++)
    {
        minp[j] = FLINT_MIN(minp[j], S[i*n + j]);
        maxp[j] = FLINT_MAX(maxp[j], S[i*n + j]);
    }

    for (j = 0; j < n; j++)
    {
        D[j] = minp[j];
        deg[j] = 1 + maxp[j] - minp[j];
    }

    for (i = 0; i < l; i++)
    for (j = 0; j < n; j++)
        S[i*n + j] -= D[j];

    tot = 0;
    for (j = 0; j < n; j++)
    {
        if (z_add_checked(&tot, tot, deg[j]))
            goto done;
    }

again:

    fmpz_one(best_prod);
    for (j = 0; j < n; j++)
        fmpz_mul_si(best_prod, best_prod, deg[j]);

    best_loc_i = -1;
    best_loc_j = -1;
    best_min = 0;
    for (i = 0; i < n; i++)
    {
        slong this_best_j, this_best_min = WORD_MAX, this_best_deg;

        mind = WORD_MAX;
        maxd = WORD_MIN;
        for (j = 0; j < n; j++)
        {
            minp[j] = minn[j] = WORD_MAX;
            maxp[j] = maxn[j] = WORD_MIN;
        }
        for (k = 0; k < l; k++)
        {
            slong * Sk = S + k*n;
            tot = 0;
            for (j = 0; j < n; j++)
            {
                tot += Sk[j];
                minp[j] = FLINT_MIN(minp[j], Sk[i] + Sk[j]);
                maxp[j] = FLINT_MAX(maxp[j], Sk[i] + Sk[j]);
                minn[j] = FLINT_MIN(minn[j], Sk[i] - Sk[j]);
                maxn[j] = FLINT_MAX(maxn[j], Sk[i] - Sk[j]);
            }
            mind = FLINT_MIN(mind, tot);
            maxd = FLINT_MAX(maxd, tot);
        }

        this_best_deg = deg[i];
        this_best_j = n + 1; /* something > n */
        if (1 + maxd - mind < this_best_deg)
        {
            this_best_j = 0;
            this_best_min = mind;
            this_best_deg = 1 + maxd - mind;
        }
        for (j = 0; j < n; j++)
        {
            if (j == i)
                continue;
            if (1 + maxp[j] - minp[j] < this_best_deg)
            {
                this_best_j = 1 + j;
                this_best_min = minp[j];
                this_best_deg = 1 + maxp[j] - minp[j];
            }
            if (1 + maxn[j] - minn[j] < this_best_deg)
            {
                this_best_j = -1 - j;
                this_best_min = minn[j];
                this_best_deg = 1 + maxn[j] - minn[j];
            }
        }

        if (this_best_j > n)
            continue;

        fmpz_set_si(this_prod, this_best_deg);
        for (j = 0; j < n; j++)
            if (j != i)
                fmpz_mul_si(this_prod, this_prod, deg[j]);

        if (fmpz_cmp(this_prod, best_prod) < 0)
        {
            fmpz_swap(best_prod, this_prod);
            best_loc_i = i;
            best_loc_j = this_best_j;
            best_min = this_best_min;
            best_deg = this_best_deg;
        }
    }

    if (best_loc_i >= 0)
    {
        i = best_loc_i;
        j = best_loc_j;
        deg[i] = best_deg;
        if (j < 0)
        {
            j = -j - 1;
            for (k = 0; k < l; k++)
                S[k*n + i] += -S[k*n + j] - best_min;
            for (k = 0; k < n; k++)
            {
                D[k] += best_min*V[k*n + i];
                V[k*n + j] += V[k*n + i];
            }
        }
        else if (j > 0)
        {
            j = j - 1;
            for (k = 0; k < l; k++)
                S[k*n + i] += S[k*n + j] - best_min;
            for (k = 0; k < n; k++)
            {
                D[k] += best_min*V[k*n + i];
                V[k*n + j] += -V[k*n + i];
            }         
        }
        else
        {
            for (k = 0; k < l; k++)
            {
                slong tot = 0;
                for (j = 0; j < n; j++)
                    tot += S[k*n + j];
                S[k*n + i] = tot - best_min;
            }
            for (k = 0; k < n; k++)
            {
                D[k] += best_min*V[k*n + i];
                for (j = 0; j < n; j++)
                    if (j != i)
                        V[k*n + j] += -V[k*n + i];
            }
        }

        goto again;
    }

done:

    for (i = 0; i < n; i++)
        tmp[i] = i;
    for (i = 1; i < n; i++)
        for (j = i; j > 0 && deg[tmp[j]] < deg[tmp[j - 1]]; j--)
            SLONG_SWAP(tmp[j], tmp[j - 1]);
    m = 1;
    while (m < n && deg[tmp[n - (m + 1)]] > 1)
        m++;
    for (i = 0; i < n; i++)
        perm[i] = (i < m) ? tmp[n - m + i] : tmp[i - m];

    for (i = 0; i < n; i++)
        tmp[i] = deg[perm[i]] - 1;
    for (i = 0; i < n; i++)
        deg[i] = tmp[i];

    for (k = 0; k < l; k++)
    {
        slong * Sk = S + k*n;
        for (i = 0; i < n; i++)
            tmp[i] = Sk[perm[i]];
        for (i = 0; i < n; i++)
            Sk[i] = tmp[i];
    }

    for (k = 0; k < n; k++)
    {
        slong * Vk = V + k*n;
        for (i = 0; i < n; i++)
            tmp[i] = Vk[perm[i]];
        for (i = 0; i < n; i++)
            Vk[i] = tmp[i];
    }

    flint_free(minp);
    flint_free(maxp);
    flint_free(minn);
    flint_free(maxn);
    flint_free(perm);
    fmpz_clear(this_prod);
    fmpz_clear(best_prod);

    return m;
}


void mpoly_compression_clear(mpoly_compression_t M)
{
    flint_free(M->umat);
    flint_free(M->deltas);
    flint_free(M->degs);
    flint_free(M->exps);
}


void mpoly_compression_init(
    mpoly_compression_t M,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mpoly_ctx_t mctx)
{
    int same;
    slong i, j;
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong nvars = mctx->nvars;

    M->nvars = nvars;
    M->umat = FLINT_ARRAY_ALLOC(nvars*nvars, slong);
    M->deltas = FLINT_ARRAY_ALLOC(nvars, slong);
    M->degs = FLINT_ARRAY_ALLOC(nvars, slong);
    M->exps = FLINT_ARRAY_ALLOC(Alen*nvars, slong);
    M->exps_alloc = Alen*nvars;
    for (i = 0; i < Alen; i++)
        mpoly_get_monomial_ui_sp((ulong *)M->exps + nvars*i, Aexps + N*i, Abits, mctx);

    M->mvars = _mpoly_compress_exps(M->umat, M->deltas, M->degs, M->exps, nvars, Alen);

    FLINT_ASSERT(M->mvars > 0);

    same = (M->mvars == nvars) && (mctx->ord == ORD_LEX);
    if (same)
    {
        for (i = 0; i < nvars; i++)
        for (j = 0; j < nvars; j++)
            same &= (M->umat[i*nvars + j] == (i == j));
    }

    M->is_trivial = same;
}

