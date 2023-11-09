/*
    Copyright (C) 2018-2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpoly.h"

void mpoly_gcd_info_init(mpoly_gcd_info_t I, slong nvars)
{
    char * d;

    FLINT_ASSERT(nvars > 0);

    d = (char *) flint_malloc(nvars*(10*sizeof(ulong) + 12*sizeof(slong)));

    I->data = d;

    I->Amax_exp         = (ulong *) d; d += nvars*sizeof(ulong);
    I->Amin_exp         = (ulong *) d; d += nvars*sizeof(ulong);
    I->Astride          = (ulong *) d; d += nvars*sizeof(ulong);
    I->Adeflate_deg     = (slong *) d; d += nvars*sizeof(slong);
    I->Alead_count      = (slong *) d; d += nvars*sizeof(slong);
    I->Atail_count      = (slong *) d; d += nvars*sizeof(slong);

    I->Bmax_exp         = (ulong *) d; d += nvars*sizeof(ulong);
    I->Bmin_exp         = (ulong *) d; d += nvars*sizeof(ulong);
    I->Bstride          = (ulong *) d; d += nvars*sizeof(ulong);
    I->Bdeflate_deg     = (slong *) d; d += nvars*sizeof(slong);
    I->Blead_count      = (slong *) d; d += nvars*sizeof(slong);
    I->Btail_count      = (slong *) d; d += nvars*sizeof(slong);

    I->Gmin_exp           = (ulong *) d; d += nvars*sizeof(ulong);
    I->Abarmin_exp        = (ulong *) d; d += nvars*sizeof(ulong);
    I->Bbarmin_exp        = (ulong *) d; d += nvars*sizeof(ulong);
    I->Gstride            = (ulong *) d; d += nvars*sizeof(ulong);
    I->Gdeflate_deg_bound = (slong *) d; d += nvars*sizeof(slong);
    I->Gterm_count_est    = (slong *) d; d += nvars*sizeof(slong);

    I->hensel_perm   = (slong *) d; d += nvars*sizeof(slong);
    I->brown_perm   = (slong *) d; d += nvars*sizeof(slong);
    I->zippel_perm  = (slong *) d; d += nvars*sizeof(slong);
    I->zippel2_perm = (slong *) d; d += nvars*sizeof(slong);
}

void mpoly_gcd_info_clear(mpoly_gcd_info_t I)
{
    flint_free(I->data);
}

/*
    Scan A and fill in the min and max exponents of each variable along
    with the count of terms attached to each.
*/
void mpoly_gcd_info_limits(ulong * Amax_exp, ulong * Amin_exp,
                       slong * Amax_exp_count, slong * Amin_exp_count,
                       const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                        const mpoly_ctx_t mctx)
{
    ulong * exps;
    slong i, j, N;
    slong nvars = mctx->nvars;
    TMP_INIT;

    FLINT_ASSERT(Alength > 0);
    FLINT_ASSERT(Abits <= FLINT_BITS);

    TMP_START;

    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    N = mpoly_words_per_exp(Abits, mctx);

    i = 0;
    mpoly_get_monomial_ui(exps, Aexps + N*i, Abits, mctx);
    for (j = 0; j < nvars; j++)
    {
        Amin_exp[j] = exps[j];
        Amax_exp[j] = exps[j];
        Amin_exp_count[j] = 1;
        Amax_exp_count[j] = 1;
    }
    for (i = 1; i < Alength; i++)
    {
        mpoly_get_monomial_ui(exps, Aexps + N*i, Abits, mctx);

        for (j = 0; j < nvars; j++)
        {
            if (Amin_exp[j] > exps[j])
            {
                Amin_exp[j] = exps[j];
                Amin_exp_count[j] = 1;
            }
            else if (Amin_exp[j] == exps[j])
            {
                Amin_exp_count[j] += 1;
            }

            if (Amax_exp[j] < exps[j])
            {
                Amax_exp[j] = exps[j];
                Amax_exp_count[j] = 1;
            }
            else if (Amax_exp[j] == exps[j])
            {
                Amax_exp_count[j] += 1;
            }
        }
    }

    TMP_END;
}


/*
    For each variable v, let SA[v] be the set of exponents of variable v in A.
    Ditto for SB[v]. The function computes
        strides[v] = GCD(SA[v] - min(SA[v]), SB[v] - min(SB[v]))
    It is assumed that {A|B}{max|min}_exp are correct.
*/
void mpoly_gcd_info_stride(ulong * strides,
          const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                             const ulong * Amax_exp, const ulong * Amin_exp,
          const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength,
                             const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                        const mpoly_ctx_t mctx)
{
    slong i, j, NA, NB;
    slong nvars = mctx->nvars;
    ulong mask;
    ulong * exps;
    TMP_INIT;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Bbits <= FLINT_BITS);

    for (j = 0; j < nvars; j++)
    {
        strides[j] = n_gcd(Amax_exp[j] - Amin_exp[j],
                           Bmax_exp[j] - Bmin_exp[j]);
    }

    TMP_START;
    exps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    NA = mpoly_words_per_exp(Abits, mctx);

    for (i = 0; i < Alength; i++)
    {
        mpoly_get_monomial_ui(exps, Aexps + NA*i, Abits, mctx);
        mask = 0;
        for (j = 0; j < nvars; j++)
        {
            strides[j] = n_gcd(strides[j], exps[j] - Amin_exp[j]);
            mask |= strides[j];
        }
        if (mask < UWORD(2))
        {
            goto cleanup;
        }
    }

    NB = mpoly_words_per_exp(Bbits, mctx);

    for (i = 0; i < Blength; i++)
    {
        mpoly_get_monomial_ui(exps, Bexps + NB*i, Bbits, mctx);
        mask = 0;
        for (j = 0; j < nvars; j++)
        {
            strides[j] = n_gcd(strides[j], exps[j] - Bmin_exp[j]);
            mask |= strides[j];
        }
        if (mask < UWORD(2))
        {
            goto cleanup;
        }
    }

cleanup:

    TMP_END;

    return;
}

void mpoly_gcd_info_set_perm(
    mpoly_gcd_info_t I,
    slong Alength,
    slong Blength,
    const mpoly_ctx_t mctx)
{
    slong j, m;

    I->Adensity = Alength;
    I->Bdensity = Blength;

    m = 0;
    for (j = 0; j < mctx->nvars; j++)
    {
        if (I->Amax_exp[j] > I->Amin_exp[j])
        {
            FLINT_ASSERT(I->Gstride[j] != UWORD(0));
            FLINT_ASSERT((I->Amax_exp[j] - I->Amin_exp[j]) % I->Gstride[j] == 0);
            FLINT_ASSERT((I->Bmax_exp[j] - I->Bmin_exp[j]) % I->Gstride[j] == 0);

            I->Adensity /= UWORD(1) + (ulong)(I->Adeflate_deg[j]);
            I->Bdensity /= UWORD(1) + (ulong)(I->Bdeflate_deg[j]);

            I->hensel_perm[m] = j;
            I->brown_perm[m] = j;
            I->zippel_perm[m] = j;
            I->zippel2_perm[m] = j;
            m++;
        }
        else
        {
            FLINT_ASSERT(I->Amax_exp[j] == I->Amin_exp[j]);
            FLINT_ASSERT(I->Bmax_exp[j] == I->Bmin_exp[j]);
        }
    }

    I->mvars = m;

    I->can_use = 0;
}


void mpoly_gcd_info_measure_hensel(
    mpoly_gcd_info_t I,
    slong Alength,
    slong Blength,
    const mpoly_ctx_t mctx)
{
    slong i, k;
    slong m = I->mvars;
    slong * perm = I->hensel_perm;
    flint_bitcnt_t abits, bbits;
    double te, tg, ta, tb;
    double stgab, mtgab, iblend, eblend;

    /* need at least 2 variables */
    if (m < 2)
        return;

    abits = FLINT_BIT_COUNT(Alength);
    bbits = FLINT_BIT_COUNT(Blength);

    te = tg = ta = tb = 1;
    for (i = 0; i < m; i++)
    {
        double x;

        k = perm[i];
        if (abits + FLINT_BIT_COUNT(I->Adeflate_deg[k]) > FLINT_BITS ||
            bbits + FLINT_BIT_COUNT(I->Bdeflate_deg[k]) > FLINT_BITS)
        {
            /* each variable is eventually converted to dense storage */
            return;
        }

        te *= 1 + FLINT_MAX(I->Adeflate_deg[k], I->Bdeflate_deg[k]);

        x = I->Gdeflate_deg_bound[k];
        tg *= 1 + x + 0.005*x*x;

        x = FLINT_MAX(0, I->Adeflate_deg[k] - I->Gdeflate_deg_bound[k]);
        ta *= 1 + x + 0.005*x*x;

        x = FLINT_MAX(0, I->Bdeflate_deg[k] - I->Gdeflate_deg_bound[k]);
        tb *= 1 + x + 0.005*x*x;
    }

    iblend = 1;
    eblend = 1;

    stgab = tg + ta + tb;
    mtgab = FLINT_MIN(tg, ta);
    mtgab = FLINT_MIN(mtgab, tb);

    I->can_use |= MPOLY_GCD_USE_HENSEL;
    I->hensel_time = 0.005*te*(I->Adensity + I->Bdensity)*eblend +
                     0.004*(iblend*stgab + (1 - iblend)*mtgab);
}

/*
    limit past which we should not test divisibility
*/
slong mpoly_gcd_info_get_brown_upper_limit(
    const mpoly_gcd_info_t I,
    slong var,
    slong bound)
{
    if (I == NULL || !I->Gdeflate_deg_bounds_are_nice)
    {
        return 0;
    }
    else
    {
        slong k, max;
        slong density;
        slong limit;

        FLINT_ASSERT(var < I->mvars);
        k = I->brown_perm[var];
        max = FLINT_MAX(I->Adeflate_deg[k], I->Bdeflate_deg[k]);
        bound = FLINT_MAX(bound, 1 + max);
        density = 0.5*(I->Adensity + I->Bdensity);
        limit = bound*((1.125 - density)*(1.125 - density))*0.375;
        return FLINT_MIN(limit, bound/2);
    }
}

void mpoly_gcd_info_measure_brown(
    mpoly_gcd_info_t I,
    slong Alength,
    slong Blength,
    const mpoly_ctx_t mctx)
{
    slong i, k;
    slong m = I->mvars;
    slong * perm = I->brown_perm;
    flint_bitcnt_t abits, bbits;
    double te, tg, ta, tb;
    double stgab, mtgab, iblend, eblend;

    /* need at least 2 variables */
    if (m < 2)
        return;

    abits = FLINT_BIT_COUNT(Alength);
    bbits = FLINT_BIT_COUNT(Blength);

    te = tg = ta = tb = 1;
    for (i = 0; i < m; i++)
    {
        double x;

        k = perm[i];
        if (abits + FLINT_BIT_COUNT(I->Adeflate_deg[k]) > FLINT_BITS ||
            bbits + FLINT_BIT_COUNT(I->Bdeflate_deg[k]) > FLINT_BITS)
        {
            /* each variable is eventually converted to dense storage */
            return;
        }

        te *= 1 + FLINT_MAX(I->Adeflate_deg[k], I->Bdeflate_deg[k]);

        x = I->Gdeflate_deg_bound[k];
        tg *= 1 + x + 0.005*x*x;

        x = FLINT_MAX(0, I->Adeflate_deg[k] - I->Gdeflate_deg_bound[k]);
        ta *= 1 + x + 0.005*x*x;

        x = FLINT_MAX(0, I->Bdeflate_deg[k] - I->Gdeflate_deg_bound[k]);
        tb *= 1 + x + 0.005*x*x;
    }

    iblend = 1;
    eblend = 1;
    if (I->Gdeflate_deg_bounds_are_nice)
    {
        slong k = perm[m - 1];
        slong limit = mpoly_gcd_info_get_brown_upper_limit(I, m - 1, 0);
        slong expected_stab;

        expected_stab = FLINT_MIN(I->Adeflate_deg[k], I->Bdeflate_deg[k]);
        expected_stab = expected_stab - I->Gdeflate_deg_bound[k];
        expected_stab = FLINT_MIN(expected_stab, I->Gdeflate_deg_bound[k]);

        if (expected_stab < limit)
        {
            slong bound = 1 + FLINT_MAX(I->Adeflate_deg[k], I->Bdeflate_deg[k]);
            iblend = (I->Adensity + I->Bdensity);
            iblend = FLINT_MIN(iblend, 1);
            iblend = FLINT_MAX(iblend, 0.01);
            eblend = 0.25 + 0.75*(double)(expected_stab)/(double)(bound);
        }
    }

    stgab = tg + ta + tb;
    mtgab = FLINT_MIN(tg, ta);
    mtgab = FLINT_MIN(mtgab, tb);

    I->can_use |= MPOLY_GCD_USE_BROWN;
    I->brown_time = 0.005*te*(I->Adensity + I->Bdensity)*eblend +
                    0.004*(iblend*stgab + (1 - iblend)*mtgab);
}

void mpoly_gcd_info_measure_bma(
    mpoly_gcd_info_t I,
    slong Alength,
    slong Blength,
    const mpoly_ctx_t mctx)
{
    slong i, j, k;
    slong m = I->mvars;
    slong * perm = I->zippel2_perm;
    slong max_main_degree;
    double Glength, Glead_count_X, Gtail_count_X, Glead_count_Y, Gtail_count_Y;
    double evals, bivar, reconstruct;

    /* need at least 3 variables */
    if (m < 3)
        return;

    /* figure out the two main variables y_0, y_1 */
    for (k = 0; k < 2; k++)
    {
        slong main_var;
        ulong count, deg, new_count, new_deg;

        main_var = k;
        j = perm[main_var];
        count = FLINT_MIN(I->Alead_count[j], I->Blead_count[j]);
        deg = FLINT_MAX(I->Adeflate_deg[j], I->Bdeflate_deg[j]);
        for (i = k + 1; i < m; i++)
        {
            j = perm[i];
            new_count = FLINT_MIN(I->Alead_count[j], I->Blead_count[j]);
            new_deg = FLINT_MAX(I->Adeflate_deg[j], I->Bdeflate_deg[j]);
            if (new_deg + new_count/256 < deg + count/256)
            {
                count = new_count;
                deg = new_deg;
                main_var = i;
            }
        }

        if (main_var != k)
        {
            slong t = perm[main_var];
            perm[main_var] = perm[k];
            perm[k] = t;
        }
    }

    max_main_degree = 0;
    for (i = 0; i < 2; i++)
    {
        k = perm[i];
        max_main_degree = FLINT_MAX(max_main_degree, I->Adeflate_deg[k]);
        max_main_degree = FLINT_MAX(max_main_degree, I->Bdeflate_deg[k]);
    }

    /* two main variables must be packed into bits = FLINT_BITS/2 */
    if (FLINT_BIT_COUNT(max_main_degree) >= FLINT_BITS/2)
        return;

    /* estimate length of gcd */
    Glength = 0.5*(I->Adensity + I->Bdensity);
    for (i = 0; i < m; i++)
    {
        k = perm[i];
        Glength *= UWORD(1) + (ulong)(I->Gdeflate_deg_bound[k]);
    }

    /* estimate number of lead/tail terms of G wrt X,Y */
    {
        double a, b;
        double Alead_density_X, Atail_density_X;
        double Blead_density_X, Btail_density_X;
        double Alead_density_Y, Atail_density_Y;
        double Blead_density_Y, Btail_density_Y;

        k = perm[0];
        a = I->Adensity*(UWORD(1) + (ulong)(I->Adeflate_deg[k]))/Alength;
        b = I->Bdensity*(UWORD(1) + (ulong)(I->Bdeflate_deg[k]))/Blength;
        Alead_density_X = a*I->Alead_count[k];
        Atail_density_X = a*I->Atail_count[k];
        Blead_density_X = b*I->Blead_count[k];
        Btail_density_X = b*I->Btail_count[k];

        k = perm[1];
        a = I->Adensity*(UWORD(1) + (ulong)(I->Adeflate_deg[k]))/Alength;
        b = I->Bdensity*(UWORD(1) + (ulong)(I->Bdeflate_deg[k]))/Blength;
        Alead_density_Y = a*I->Alead_count[k];
        Atail_density_Y = a*I->Atail_count[k];
        Blead_density_Y = b*I->Blead_count[k];
        Btail_density_Y = b*I->Btail_count[k];

        Glead_count_X = 0.5*(Alead_density_X + Blead_density_X);
        Gtail_count_X = 0.5*(Atail_density_X + Btail_density_X);
        Glead_count_Y = 0.5*(Alead_density_Y + Blead_density_Y);
        Gtail_count_Y = 0.5*(Atail_density_Y + Btail_density_Y);
        for (i = 0; i < m; i++)
        {
            k = perm[i];
            if (i != 0)
            {
                Glead_count_X *= UWORD(1) + (ulong)(I->Gdeflate_deg_bound[k]);
                Gtail_count_X *= UWORD(1) + (ulong)(I->Gdeflate_deg_bound[k]);
            }

            if (i != 1)
            {
                Glead_count_Y *= UWORD(1) + (ulong)(I->Gdeflate_deg_bound[k]);
                Gtail_count_Y *= UWORD(1) + (ulong)(I->Gdeflate_deg_bound[k]);
            }
        }
    }

    /* evaluations needed is the max length of the coefficients of G wrt X,Y */
    {
        double Gmax_terms_X, Gmax_terms_Y;

        k = perm[0];
        Gmax_terms_X = Glength/(UWORD(1) + (ulong)(I->Gterm_count_est[k]));
        Gmax_terms_X = FLINT_MAX(Gmax_terms_X, Glead_count_X);
        Gmax_terms_X = FLINT_MAX(Gmax_terms_X, Gtail_count_X);
        Gmax_terms_X = FLINT_MAX(Gmax_terms_X, 1);

        k = perm[1];
        Gmax_terms_Y = Glength/(UWORD(1) + (ulong)(I->Gterm_count_est[k]));
        Gmax_terms_Y = FLINT_MAX(Gmax_terms_Y, Glead_count_Y);
        Gmax_terms_Y = FLINT_MAX(Gmax_terms_Y, Gtail_count_Y);
        Gmax_terms_Y = FLINT_MAX(Gmax_terms_Y, 1);

        evals = Gmax_terms_X*Gmax_terms_Y/(1 + Glength);
    }

    /* time for bivar gcd */
    {
        double te, tg, ta, tb;
        te = tg = ta = tb = 1;
        for (i = 0; i < 2; i++)
        {
            k = perm[i];
            /* already checked reasonable degrees with max_main_degree above */
            te *= 1 + FLINT_MAX(I->Adeflate_deg[k], I->Bdeflate_deg[k]);
            tg *= 1 + I->Gdeflate_deg_bound[k];
            ta *= 1 + FLINT_MAX(0, I->Adeflate_deg[k] - I->Gdeflate_deg_bound[k]);
            tb *= 1 + FLINT_MAX(0, I->Bdeflate_deg[k] - I->Gdeflate_deg_bound[k]);
        }
        bivar = te + (tg + ta + tb)*0.1;
    }

    /* time for reconstruction */
    {
        reconstruct = I->Gterm_count_est[perm[0]];
        reconstruct += I->Gterm_count_est[perm[1]];
        reconstruct = Glength*Glength/(1 + reconstruct);
    }

    I->can_use |= MPOLY_GCD_USE_ZIPPEL2;
    I->zippel2_time = 0.00000002*bivar*evals*(Alength + Blength) +
                      0.0003*reconstruct;
}

void mpoly_gcd_info_measure_zippel(
    mpoly_gcd_info_t I,
    slong Alength,
    slong Blength,
    const mpoly_ctx_t mctx)
{
    slong i, j, k;
    slong m = I->mvars;
    slong * perm = I->zippel_perm;

    /* need at least two variables */
    if (m < 2)
        return;

    /* figure out a main variable y_0 */
    {
        slong main_var;
        ulong count, deg, new_count, new_deg;

        main_var = 0;
        j = I->zippel_perm[main_var];
        count = FLINT_MIN(I->Atail_count[j], I->Alead_count[j]);
        count = FLINT_MIN(count, I->Btail_count[j]);
        count = FLINT_MIN(count, I->Blead_count[j]);
        deg = FLINT_MAX(I->Adeflate_deg[j], I->Bdeflate_deg[j]);
        for (i = 1; i < m; i++)
        {
            j = perm[i];
            new_count = FLINT_MIN(I->Atail_count[j], I->Alead_count[j]);
            new_count = FLINT_MIN(new_count, I->Btail_count[j]);
            new_count = FLINT_MIN(new_count, I->Blead_count[j]);
            new_deg = FLINT_MAX(I->Adeflate_deg[j], I->Bdeflate_deg[j]);
            if (new_count < count || (new_count == count && new_deg < deg))
            {
                count = new_count;
                deg = new_deg;
                main_var = i;
            }
        }

        if (main_var != 0)
        {
            slong t = perm[main_var];
            perm[main_var] = perm[0];
            perm[0] = t;
        }
    }

    /* sort with hope that ddeg(G,y_1) >= ddeg(G,y_2) ... >= ddeg(G,y_m) */
    for (k = 1; k + 1 < m; k++)
    {
        slong var;
        ulong deg, new_deg;

        var = k;
        j = perm[var];
        deg = FLINT_MIN(I->Adeflate_deg[j], I->Bdeflate_deg[j]);
        for (i = k + 1; i < m; i++)
        {
            j = perm[i];
            new_deg = FLINT_MIN(I->Adeflate_deg[j], I->Bdeflate_deg[j]);
            if (new_deg > deg)
            {
                deg = new_deg;
                var = i;
            }
        }
        if (var != k)
        {
            slong t = I->zippel_perm[var];
            perm[var] = perm[k];
            perm[k] = t;
        }
    }

    I->can_use |= MPOLY_GCD_USE_ZIPPEL;
    I->zippel_time = 0.3456;
}

void mpoly_gcd_info_measure_zippel2(
    mpoly_gcd_info_t I,
    slong Alength,
    slong Blength,
    const mpoly_ctx_t mctx)
{
    slong i, j, k;
    slong m = I->mvars;
    slong * perm = I->zippel2_perm;
    slong max_main_degree;

    /* need at least 3 variables */
    if (m < 3)
        return;

    /* figure out the two main variables y_0, y_1 */

#define NEEDS_SWAP                                                            \
    FLINT_MIN(I->Adeflate_deg[perm[j]], I->Bdeflate_deg[perm[j]]) <           \
    FLINT_MIN(I->Adeflate_deg[perm[j-1]], I->Bdeflate_deg[perm[j-1]])         \

    for (i = 1; i < m; i++)
        for (j = i; j > 0 && NEEDS_SWAP; j--)
            FLINT_SWAP(slong, perm[j], perm[j - 1]);


#define NEEDS_SWAP2                                                           \
    FLINT_MIN(I->Adeflate_deg[perm[j]], I->Bdeflate_deg[perm[j]]) >           \
    FLINT_MIN(I->Adeflate_deg[perm[j-1]], I->Bdeflate_deg[perm[j-1]])         \

    for (i = 3; i < m; i++)
        for (j = i; j > 2 && NEEDS_SWAP; j--)
            FLINT_SWAP(slong, perm[j], perm[j - 1]);


    max_main_degree = 0;
    for (i = 0; i < 2; i++)
    {
        k = perm[i];
        max_main_degree = FLINT_MAX(max_main_degree, I->Adeflate_deg[k]);
        max_main_degree = FLINT_MAX(max_main_degree, I->Bdeflate_deg[k]);
    }

    /* two main variables must be packed into bits = FLINT_BITS/2 */
    if (FLINT_BIT_COUNT(max_main_degree) >= FLINT_BITS/2)
        return;

    I->can_use |= MPOLY_GCD_USE_ZIPPEL2;
    I->zippel2_time = 0.243;
}
