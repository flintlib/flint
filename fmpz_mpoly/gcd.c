/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"
#include "ulong_extras.h"


/*
    Assume when B is converted to univar format, its length would be one.
    Gcd is gcd of coefficients of univar(A) and B (modulo some shifts).
*/
static int _try_missing_var(fmpz_mpoly_t G, mp_bitcnt_t Gbits, slong var,
                                         const fmpz_mpoly_t A, ulong Ashift,
                                         const fmpz_mpoly_t B, ulong Bshift,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fmpz_mpoly_t tG;
    fmpz_mpoly_univar_t Ax;

    fmpz_mpoly_init(tG, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);

    success = fmpz_mpoly_to_univar(Ax, A, var, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Ax->length > 0);
    success = _fmpz_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx);

    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _fmpz_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }

    fmpz_mpoly_swap(G, tG, ctx);
    _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

cleanup:

    fmpz_mpoly_clear(tG, ctx);
    fmpz_mpoly_univar_clear(Ax, ctx);

    return success;
}

/*
    return 1 for success or 0 for failure
*/
static int _try_zippel(fmpz_mpoly_t G, mp_bitcnt_t Gbits, ulong * Gstride,
         const fmpz_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
                   const slong * Amax_exp_count, const slong * Amin_exp_count,
         const fmpz_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                   const slong * Bmax_exp_count, const slong * Bmin_exp_count,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong n, m;
    int success;
    ulong * Gshift, * Addeg, * Bddeg;
    mpoly_zipinfo_t zinfo;
    mp_bitcnt_t ABbits;
    flint_rand_t randstate;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu;
    fmpz_mpoly_t Acontent, Bcontent;
    fmpz_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    flint_randinit(randstate);

    /*
        Let the variables in A and B be
            x_0, x_1, ..., x_{n-1}
        where n = ctx->minfo->nvars, and x_0 is most significant

        Recall from _fmpz_mpoly_gcd that all variables are either
            missing from both ess(A) and ess(B) (non-essential), or
            present in both ess(A) and ess(B) (essential)
        and that there are at least 2 variables in the essential case.

        Let y_0, y_1, ..., y_{m} with m >= 1 denote the variables present
        in both ess(A) and ess(B). Each y_i is one of the x_j and the variables
        are ordered as y_0 > ... > y_{m} with LEX order.
        The Zippel algorithm will operate in Z[y_0][y_1,...,y_m] and it
        only operates with Z[y_1,...,y_m] in LEX.

        When converting to the mpolyu format via
        fmpz_mpoly_to_mpolyu_perm_deflate, the non-essential variables
        will be immediately striped out and the remaining variables will be
        mapped according to the permutation in zinfo->perm as

            y_k = x_perm[k] ^ Gstride[perm[k]]

        When converting out of the mpolyu format via
        fmpz_mpoly_from_mpolyu_perm_inflate, the contribution of the
        non-essential variables will be put back in.
    */

    n = ctx->minfo->nvars;
    m = 0;
    for (j = 0; j < n; j++)
    {
        if (Amax_exp[j] > Amin_exp[j])
        {
            FLINT_ASSERT(Bmax_exp[j] > Bmin_exp[j]);
            m++;
        }
        else
        {
            FLINT_ASSERT(Bmax_exp[j] == Bmin_exp[j]);
        }
    }

    /* need at least two variables */
    if (m < 2)
        return 0;

    /* interpolation will continue in m variables */
    mpoly_zipinfo_init(zinfo, m);

    /* uctx is context for Z[y_1,...,y_{m-1}]*/
    fmpz_mpoly_ctx_init(uctx, m - 1, ORD_LEX);

    Gshift = (ulong *) flint_malloc(n*sizeof(ulong));
    /* degrees after deflation */
    Addeg = (ulong *) flint_malloc(n*sizeof(ulong));
    Bddeg = (ulong *) flint_malloc(n*sizeof(ulong));

    /* fill in a valid zinfo->perm and set deflation values used when converting */
    i = 0;
    for (j = 0; j < n; j++)
    {
        Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);
        if (Amax_exp[j] > Amin_exp[j])
        {
            FLINT_ASSERT(Gstride[j] != UWORD(0));
            FLINT_ASSERT((Amax_exp[j] - Amin_exp[j]) % Gstride[j] == UWORD(0));
            FLINT_ASSERT((Bmax_exp[j] - Bmin_exp[j]) % Gstride[j] == UWORD(0));

            Addeg[j] = (Amax_exp[j] - Amin_exp[j]) / Gstride[j];
            Bddeg[j] = (Bmax_exp[j] - Bmin_exp[j]) / Gstride[j];

            zinfo->perm[i] = j;
            i++;
        }
        else
        {
            Addeg[j] = 0;
            Bddeg[j] = 0;
        }
    }
    FLINT_ASSERT(i == m);

    /* work out a favourable permation to zinfo->perm */

    /* figure out a main variable y_0 */
    {
        slong main_var;
        ulong count, deg, new_count, new_deg;

        main_var = 0;
        j = zinfo->perm[main_var];
        count = FLINT_MIN(FLINT_MIN(Amin_exp_count[j], Amax_exp_count[j]),
                          FLINT_MIN(Bmin_exp_count[j], Bmax_exp_count[j]));
        deg = FLINT_MAX(Addeg[j], Bddeg[j]);
        for (i = 1; i < m; i++)
        {
            j = zinfo->perm[i];
            new_count = FLINT_MIN(FLINT_MIN(Amin_exp_count[j], Amax_exp_count[j]),
                                  FLINT_MIN(Bmin_exp_count[j], Bmax_exp_count[j]));
            new_deg = FLINT_MAX(Addeg[j], Bddeg[j]);

            if (new_count < count || (new_count == count && new_deg < deg))
            {
                count = new_count;
                deg = new_deg;
                main_var = i;
            }
        }

        if (main_var != 0)
        {
            slong t = zinfo->perm[main_var];
            zinfo->perm[main_var] = zinfo->perm[0];
            zinfo->perm[0] = t;
        }
    }

    /* sort with hope that ddeg(G,y_1) >= ddeg(G,y_2) ... >= ddeg(G,y_m) */
    for (k = 1; k + 1 < m; k++)
    {
        slong var;
        ulong deg, new_deg;

        var = k;
        j = zinfo->perm[var];
        deg = FLINT_MIN(Addeg[j], Bddeg[j]);
        for (i = k + 1; i < m; i++)
        {
            j = zinfo->perm[i];
            new_deg = FLINT_MIN(Addeg[j], Bddeg[j]);
            if (new_deg > deg)
            {
                deg = new_deg;
                var = i;
            }
        }
        if (var != k)
        {
            slong t = zinfo->perm[var];
            zinfo->perm[var] = zinfo->perm[k];
            zinfo->perm[k] = t;
        }
    }

    /* fill in the degrees of the y_i */
    for (i = 0; i < m; i++)
    {
        j = zinfo->perm[i];
        FLINT_ASSERT(Addeg[j] != 0);
        FLINT_ASSERT(Bddeg[j] != 0);
        zinfo->Adegs[i] = Addeg[j];
        zinfo->Bdegs[i] = Bddeg[j];
    }

    ABbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpolyu_init(Au, ABbits, uctx);
    fmpz_mpolyu_init(Bu, ABbits, uctx);
    fmpz_mpolyu_init(Gu, ABbits, uctx);

    fmpz_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx,
                                               zinfo->perm, Amin_exp, Gstride);
    fmpz_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx,
                                               zinfo->perm, Bmin_exp, Gstride);

    FLINT_ASSERT(Au->bits == ABbits);
    FLINT_ASSERT(Au->bits == ABbits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    fmpz_mpoly_init3(Acontent, 0, ABbits, uctx);
    fmpz_mpoly_init3(Bcontent, 0, ABbits, uctx);
    fmpz_mpolyu_init(Abar, ABbits, uctx);
    fmpz_mpolyu_init(Bbar, ABbits, uctx);
    fmpz_mpolyu_init(Gbar, ABbits, uctx);

    /* remove content from A and B */
    success = fmpz_mpolyu_content_mpoly(Acontent, Au, uctx);
    success = success && fmpz_mpolyu_content_mpoly(Bcontent, Bu, uctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);
    FLINT_ASSERT(Bcontent->bits == ABbits);
    fmpz_mpolyu_divexact_mpoly(Abar, Au, Acontent, uctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, Bu, Bcontent, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */

    /* compute GCD */
    success = fmpz_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);

    fmpz_mpolyu_mul_mpoly(Gu, Gbar, Acontent, uctx);
    fmpz_mpoly_from_mpolyu_perm_inflate(G, Gbits, ctx, Gu, uctx,
                                                 zinfo->perm, Gshift, Gstride);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    fmpz_mpolyu_clear(Abar, uctx);
    fmpz_mpolyu_clear(Bbar, uctx);
    fmpz_mpolyu_clear(Gbar, uctx);
    fmpz_mpoly_clear(Acontent, uctx);
    fmpz_mpoly_clear(Bcontent, uctx);

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    flint_free(Gshift);
    flint_free(Addeg);
    flint_free(Bddeg);

    return success;
}


/*
    return 1 for success or 0 for failure
*/
static int _try_berlekamp_massey(fmpz_mpoly_t G, mp_bitcnt_t Gbits, ulong * Gstride,
         const fmpz_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
                   const slong * Amax_exp_count, const slong * Amin_exp_count,
         const fmpz_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                   const slong * Bmax_exp_count, const slong * Bmin_exp_count,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong n, m;
    int success;
    ulong * Gshift, * Addeg, * Bddeg;
    slong * perm;
    mp_bitcnt_t ABbits;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Auu, Buu, Guu;
    fmpz_mpoly_t Acontent, Bcontent, Gamma;
    fmpz_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    n = ctx->minfo->nvars;
    m = 0;
    for (j = 0; j < n; j++)
    {
        if (Amax_exp[j] > Amin_exp[j])
        {
            FLINT_ASSERT(Bmax_exp[j] > Bmin_exp[j]);
            m++;
        }
        else
        {
            FLINT_ASSERT(Amax_exp[j] == Amin_exp[j]);
            FLINT_ASSERT(Bmax_exp[j] == Bmin_exp[j]);
        }
    }

    /* need at least 3 variables, but lets require 4 for now */
    if (m < 4)
        return 0;

    /* uctx is context for Z[y_2,...,y_{m - 1}]*/
    fmpz_mpoly_ctx_init(uctx, m - 2, ORD_LEX);

    perm = (slong *) flint_malloc(m*sizeof(slong));
    Gshift = (ulong *) flint_malloc(n*sizeof(ulong));
    /* degrees after deflation */
    Addeg = (ulong *) flint_malloc(n*sizeof(ulong));
    Bddeg = (ulong *) flint_malloc(n*sizeof(ulong));

    /* fill in a valid perm and set deflation values used when converting */
    i = 0;
    for (j = 0; j < n; j++)
    {
        Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);
        if (Amax_exp[j] > Amin_exp[j])
        {
            FLINT_ASSERT(Gstride[j] != UWORD(0));
            FLINT_ASSERT((Amax_exp[j] - Amin_exp[j]) % Gstride[j] == UWORD(0));
            FLINT_ASSERT((Bmax_exp[j] - Bmin_exp[j]) % Gstride[j] == UWORD(0));

            Addeg[j] = (Amax_exp[j] - Amin_exp[j]) / Gstride[j];
            Bddeg[j] = (Bmax_exp[j] - Bmin_exp[j]) / Gstride[j];

            perm[i] = j;
            i++;
        }
        else
        {
            Addeg[j] = 0;
            Bddeg[j] = 0;
        }
    }
    FLINT_ASSERT(i == m);

    /* figure out the two main variables y_0, y_1 */
    for (k = 0; k < 2; k++)
    {
        slong main_var;
        ulong count, deg, new_count, new_deg;

        main_var = k;
        j = perm[main_var];
        count = FLINT_MIN(Amax_exp_count[j], Bmax_exp_count[j]);
        deg = FLINT_MAX(Addeg[j], Bddeg[j]);
        for (i = k + 1; i < m; i++)
        {
            j = perm[i];
            new_count = FLINT_MIN(Amax_exp_count[j], Bmax_exp_count[j]);
            new_deg = FLINT_MAX(Addeg[j], Bddeg[j]);

            if (new_deg + new_count/2 < deg + count/2)
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

    ABbits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpolyu_init(Auu, ABbits, uctx);
    fmpz_mpolyu_init(Buu, ABbits, uctx);
    fmpz_mpolyu_init(Guu, ABbits, uctx);

    fmpz_mpoly_to_mpolyuu_perm_deflate(Auu, uctx, A, ctx,
                                                      perm, Amin_exp, Gstride);
    fmpz_mpoly_to_mpolyuu_perm_deflate(Buu, uctx, B, ctx,
                                                      perm, Bmin_exp, Gstride);

    FLINT_ASSERT(Auu->bits == ABbits);
    FLINT_ASSERT(Buu->bits == ABbits);
    FLINT_ASSERT(Auu->length > 1);
    FLINT_ASSERT(Buu->length > 1);

    fmpz_mpoly_init3(Acontent, 0, ABbits, uctx);
    fmpz_mpoly_init3(Bcontent, 0, ABbits, uctx);
    fmpz_mpoly_init3(Gamma, 0, ABbits, uctx);
    fmpz_mpolyu_init(Abar, ABbits, uctx);
    fmpz_mpolyu_init(Bbar, ABbits, uctx);
    fmpz_mpolyu_init(Gbar, ABbits, uctx);

    /* remove content from A and B */
    success = fmpz_mpolyu_content_mpoly(Acontent, Auu, uctx);
    success = success && fmpz_mpolyu_content_mpoly(Bcontent, Buu, uctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);
    FLINT_ASSERT(Bcontent->bits == ABbits);
    fmpz_mpolyu_divexact_mpoly(Abar, Auu, Acontent, uctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, Buu, Bcontent, uctx);

    /* compute GCD of leading coefficients */
    FLINT_ASSERT(A->length > 0 && B->length > 0);
    _fmpz_mpoly_gcd(Gamma, ABbits, Abar->coeffs + 0, Bbar->coeffs + 0, uctx);
    if (!success)
        goto cleanup;

    success = fmpz_mpolyuu_gcd_berlekamp_massey(Gbar, Abar, Bbar, Gamma, uctx);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _fmpz_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;

    fmpz_mpolyu_mul_mpoly(Guu, Gbar, Acontent, uctx);

    fmpz_mpoly_from_mpolyuu_perm_inflate(G, Gbits, ctx, Guu, uctx,
                                                        perm, Gshift, Gstride);
    if (fmpz_sgn(G->coeffs + 0) < 0)
        fmpz_mpoly_neg(G, G, ctx);

    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    fmpz_mpolyu_clear(Abar, uctx);
    fmpz_mpolyu_clear(Bbar, uctx);
    fmpz_mpolyu_clear(Gbar, uctx);
    fmpz_mpoly_clear(Acontent, uctx);
    fmpz_mpoly_clear(Bcontent, uctx);
    fmpz_mpoly_clear(Gamma, uctx);

    fmpz_mpolyu_clear(Auu, uctx);
    fmpz_mpolyu_clear(Buu, uctx);
    fmpz_mpolyu_clear(Guu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    flint_free(Gshift);
    flint_free(Addeg);
    flint_free(Bddeg);
    flint_free(perm);

    return success;
}



/*
    return 1 for success or 0 for failure
*/
static int _try_brown(fmpz_mpoly_t G, mp_bitcnt_t Gbits, ulong * Gstride,
         const fmpz_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
         const fmpz_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong j;
    slong denseAsize, denseBsize;
    slong m, n = ctx->minfo->nvars;
    slong * perm;
    ulong * Gshift;
    mp_bitcnt_t ABbits;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu, Abaru, Bbaru;
    fmpz_t cA, cB, cG;
    TMP_INIT;

    /* first see if a dense algo is a good idea */
    if (n > 7)
        return 0;

    denseAsize = WORD(1);
    denseBsize = WORD(1);
    for (j = 0; j < n; j++)
    {
        if (Amax_exp[j] > Amin_exp[j])
        {
            ulong hi;

            FLINT_ASSERT(Gstride[j] != UWORD(0));
            FLINT_ASSERT(Gstride[j] != UWORD(0));
            FLINT_ASSERT((Amax_exp[j] - Amin_exp[j]) % Gstride[j] == UWORD(0));
            FLINT_ASSERT((Bmax_exp[j] - Bmin_exp[j]) % Gstride[j] == UWORD(0));

            umul_ppmm(hi, denseAsize, denseAsize,
                                   (Amax_exp[j] - Amin_exp[j])/Gstride[j] + 1);
            if (hi != 0 || denseAsize <= 0
                        || denseAsize > 100000000)
            {
                return 0;
            }

            umul_ppmm(hi, denseBsize, denseBsize,
                                   (Bmax_exp[j] - Bmin_exp[j])/Gstride[j] + 1);
            if (hi != 0 || denseBsize <= 0
                        || denseBsize > 100000000)
            {
                return 0;
            }
        }
    }

    if (   denseAsize/A->length > 4
        || denseBsize/B->length > 4)
    {
        return 0;
    }

    TMP_START;
    Gshift = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    perm = (slong *) TMP_ALLOC(n*sizeof(slong)); /* only first m entries used */

    /* fill in perm and set shift of GCD */
    m = 0;
    for (j = 0; j < n; j++)
    {
        Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);
        if (Amax_exp[j] > Amin_exp[j])
        {
            FLINT_ASSERT(Bmax_exp[j] > Bmin_exp[j]);
            perm[m] = j;
            m++;
        }
        else
        {
            FLINT_ASSERT(Amax_exp[j] == Amin_exp[j]);
            FLINT_ASSERT(Bmax_exp[j] == Bmin_exp[j]);
        }
    }

    /* need at least 2 variables */
    if (m < 2)
    {
        success = 0;
        goto cleanup_tmp;
    }

    /* TODO: find a favourable permutation of perm */

    ABbits = FLINT_MAX(A->bits, B->bits);

    /* uctx is for Z[y_1, ..., y_{m-1}] */
    fmpz_mpoly_ctx_init(uctx, m - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, ABbits, uctx);
    fmpz_mpolyu_init(Bu, ABbits, uctx);
    fmpz_mpolyu_init(Gu, ABbits, uctx);
    fmpz_mpolyu_init(Abaru, ABbits, uctx);
    fmpz_mpolyu_init(Bbaru, ABbits, uctx);

    fmpz_init(cA);
    fmpz_init(cB);
    fmpz_init(cG);

    fmpz_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx, perm, Amin_exp, Gstride);
    fmpz_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx, perm, Bmin_exp, Gstride);

    fmpz_mpolyu_content_fmpz(cA, Au, uctx);
    fmpz_mpolyu_content_fmpz(cB, Bu, uctx);
    fmpz_gcd(cG, cA, cB);
    fmpz_mpolyu_divexact_fmpz(Au, Au, cA, uctx);
    fmpz_mpolyu_divexact_fmpz(Bu, Bu, cB, uctx);
    success = fmpz_mpolyu_gcd_brown(Gu, Abaru, Bbaru, Au, Bu, uctx);
    if (success)
    {
        fmpz_mpoly_from_mpolyu_perm_inflate(G, Gbits, ctx, Gu, uctx, perm, Gshift, Gstride);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_neg(cG, cG);
        fmpz_mpoly_scalar_mul_fmpz(G, G, cG, ctx);
    }

    fmpz_clear(cA);
    fmpz_clear(cB);
    fmpz_clear(cG);

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpolyu_clear(Abaru, uctx);
    fmpz_mpolyu_clear(Bbaru, uctx);
    fmpz_mpoly_ctx_clear(uctx);

cleanup_tmp:

    TMP_END;

    return success;
}


/*
    return 1 for success or 0 for failure
*/
static int _try_prs(fmpz_mpoly_t G, mp_bitcnt_t Gbits,
         const fmpz_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
                   const slong * Amax_exp_count, const slong * Amin_exp_count,
         const fmpz_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                   const slong * Bmax_exp_count, const slong * Bmin_exp_count,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success = 1;
    slong i, j, var, n = ctx->minfo->nvars;
    ulong work, newwork;
    fmpz_mpoly_t ac, bc, gc, gabc, g;
    fmpz_mpoly_univar_t ax, bx, gx;

    FLINT_ASSERT(Gbits <= FLINT_BITS);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    /* find a variable with monomial leading or trailing term */
    var = -WORD(1);
    work = -UWORD(1);
    for (j = 0; j < n; j++)
    {
        if (Amax_exp[j] > Amin_exp[j])
        {
            FLINT_ASSERT(Bmax_exp[j] > Bmin_exp[j]);
            if (   Amax_exp_count[j] == 1
                || Amin_exp_count[j] == 1
                || Bmax_exp_count[j] == 1
                || Bmin_exp_count[j] == 1
               )
            {
                newwork = FLINT_MAX(Amax_exp[j] - Amin_exp[j],
                                    Bmax_exp[j] - Bmin_exp[j]);
                if (var < 0)
                {
                    var = j;
                    work = newwork;
                }
                else if (newwork < work)
                {
                    var = j;
                    work = newwork;
                }
            }
        }
    }

    /* stringent condition ensures that we only try this for low degrees */
    if (var < 0 || work + n > 5)
    {
        return 0;
    }

    fmpz_mpoly_init(ac, ctx);
    fmpz_mpoly_init(bc, ctx);
    fmpz_mpoly_init(gc, ctx);
    fmpz_mpoly_init(gabc, ctx);
    fmpz_mpoly_init(g, ctx);

    fmpz_mpoly_univar_init(ax, ctx);
    fmpz_mpoly_univar_init(bx, ctx);
    fmpz_mpoly_univar_init(gx, ctx);

    /*
        note: we have no control over the intermediate bits thoughout this
        algorithm. We just pack into Gbits at the very and.
    */

    success = fmpz_mpoly_to_univar(ax, A, var, ctx);
    if (!success)
        goto cleanup;

    success = fmpz_mpoly_to_univar(bx, B, var, ctx);
    if (!success)
        goto cleanup;

    gx->var = var;

    /* remove content from Ax */
    FLINT_ASSERT(ax->length > 1);
    FLINT_ASSERT(Amin_exp[var] == ax->exps[ax->length - 1]);
    success = fmpz_mpoly_gcd(ac, ax->coeffs + 0, ax->coeffs + 1, ctx);
    if (!success)
        goto cleanup;
    for (i = 2; i < ax->length; i++)
    {
        success = fmpz_mpoly_gcd(ac, ac, ax->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }
    for (i = 0; i < ax->length; i++)
    {
        if (!fmpz_mpoly_divides(ax->coeffs + i, ax->coeffs + i, ac, ctx))
            FLINT_ASSERT(0 && "no divides");
        ax->exps[i] -= Amin_exp[var];
    }

    /* remove content from Bx */
    FLINT_ASSERT(bx->length > 1);
    FLINT_ASSERT(Bmin_exp[var] == bx->exps[bx->length - 1]);
    success = fmpz_mpoly_gcd(bc, bx->coeffs + 0, bx->coeffs + 1, ctx);
    if (!success)
        goto cleanup;
    for (i = 2; i < bx->length; i++)
    {
        success = fmpz_mpoly_gcd(bc, bc, bx->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }
    for (i = 0; i < bx->length; i++)
    {
        if (!fmpz_mpoly_divides(bx->coeffs + i, bx->coeffs + i, bc, ctx))
            FLINT_ASSERT(0 && "no divides");
        bx->exps[i] -= Bmin_exp[var];
    }

    /* compute pseudo gcd of contentless polynomials */
    FLINT_ASSERT(ax->exps[0] > 0);
    FLINT_ASSERT(bx->exps[0] > 0);
    if (ax->exps[0] >= bx->exps[0])
    {
        _fmpz_mpoly_univar_pgcd(gx, ax, bx, ctx);
    }
    else
    {
        _fmpz_mpoly_univar_pgcd(gx, bx, ax, ctx);
    }
    FLINT_ASSERT(gx->length > 0);

    /* try to divide out easy content from gcd */
    FLINT_ASSERT(gx->length > 0);
    if ((gx->coeffs + 0)->length != 1
        && (gx->coeffs + gx->length - 1)->length != 1)
    {
        if (   (ax->coeffs + 0)->length == 1
            || (bx->coeffs + 0)->length == 1)
        {
            fmpz_mpoly_term_content(gc, gx->coeffs + 0, ctx);
            if (!fmpz_mpoly_divides(gabc, gx->coeffs + 0, gc, ctx))
                FLINT_ASSERT(0 && "no divides");
            for (i = 0; i < gx->length; i++)
            {
                if (!fmpz_mpoly_divides(gx->coeffs + i, gx->coeffs + i, gabc, ctx))
                    FLINT_ASSERT(0 && "not lead div");
            }
        }
        else if (   (ax->coeffs + ax->length - 1)->length == 1
                 || (bx->coeffs + bx->length - 1)->length == 1)
        {
            fmpz_mpoly_term_content(gc, gx->coeffs + gx->length - 1, ctx);
            if (!fmpz_mpoly_divides(gabc, gx->coeffs + gx->length - 1, gc, ctx))
                FLINT_ASSERT(0 && "no divides");
            for (i = 0; i < gx->length; i++)
            {
                if (!fmpz_mpoly_divides(gx->coeffs + i, gx->coeffs + i, gabc, ctx))
                    FLINT_ASSERT(0 && "not trail div");
            }
        }
        else
        {
            success = fmpz_mpoly_gcd(gc, ax->coeffs + 0, bx->coeffs + 0, ctx);
            if (!success)
                goto cleanup;
            if (gc->length == 1)
            {
                fmpz_mpoly_term_content(gc, gx->coeffs + 0, ctx);
                if (!fmpz_mpoly_divides(gabc, gx->coeffs + 0, gc, ctx))
                    FLINT_ASSERT(0 && "no divides");
                for (i = 0; i < gx->length; i++)
                {
                    if (!fmpz_mpoly_divides(gx->coeffs + i, gx->coeffs + i, gabc, ctx))
                        FLINT_ASSERT(0 && "not lead div");
                }
            }
        }
    }

    /* divide out monomial content from gcd */
    fmpz_mpoly_term_content(gc, gx->coeffs + 0, ctx);
    for (i = 1; i < gx->length; i++)
    {
        success = fmpz_mpoly_gcd(gc, gc, gx->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }
    for (i = 0; i < gx->length; i++)
    {
        if (!fmpz_mpoly_divides(gx->coeffs + i, gx->coeffs + i, gc, ctx))
            FLINT_ASSERT(0 && "no divides");
    }

    /* divide out content from gcd */
    fmpz_mpoly_set(gc, gx->coeffs + 0, ctx);
    for (i = 1; i < gx->length; i++)
    {
        success = fmpz_mpoly_gcd(gc, gc, gx->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }
    for (i = 0; i < gx->length; i++)
    {
        if (!fmpz_mpoly_divides(gx->coeffs + i, gx->coeffs + i, gc, ctx))
            FLINT_ASSERT(0 && "no divides");
    }

    /* put back real content */
    fmpz_mpoly_gcd(gabc, ac, bc, ctx);
    for (i = 0; i < gx->length; i++)
    {
        fmpz_mpoly_mul(gx->coeffs + i, gx->coeffs + i, gabc, ctx);
        gx->exps[i] += FLINT_MIN(Amin_exp[var], Bmin_exp[var]);
    }

    FLINT_ASSERT(Gbits <= FLINT_BITS);

    fmpz_mpoly_from_univar_bits(G, Gbits, gx, ctx);

    FLINT_ASSERT(G->length > 0);
    if (fmpz_sgn(G->coeffs + 0) < 0)
    {
        fmpz_mpoly_neg(G, G, ctx);
    }

    success = 1;

cleanup:

    fmpz_mpoly_clear(ac, ctx);
    fmpz_mpoly_clear(bc, ctx);
    fmpz_mpoly_clear(gc, ctx);
    fmpz_mpoly_clear(gabc, ctx);
    fmpz_mpoly_clear(g, ctx);

    fmpz_mpoly_univar_clear(ax, ctx);
    fmpz_mpoly_univar_clear(bx, ctx);
    fmpz_mpoly_univar_clear(gx, ctx);

    return success;
}


/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _fmpz_mpoly_gcd(fmpz_mpoly_t G, mp_bitcnt_t Gbits,
                               const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int success;
    slong v_in_both;
    slong v_in_either;
    slong v_in_A_only;
    slong v_in_B_only;
    slong j;
    slong nvars = ctx->minfo->nvars;
    ulong * Amax_exp, * Amin_exp;
    ulong * Bmax_exp, * Bmin_exp;
    slong * Amax_exp_count, * Amin_exp_count;
    slong * Bmax_exp_count, * Bmin_exp_count;
    ulong * Gstride, * Gshift;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(Gbits <= FLINT_BITS);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    if (A->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }

    Amax_exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
    Amin_exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
    Bmax_exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
    Bmin_exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
    Amax_exp_count = (slong *) flint_malloc(nvars*sizeof(slong));
    Amin_exp_count = (slong *) flint_malloc(nvars*sizeof(slong));
    Bmax_exp_count = (slong *) flint_malloc(nvars*sizeof(slong));
    Bmin_exp_count = (slong *) flint_malloc(nvars*sizeof(slong));
    Gstride = (ulong *) flint_malloc(nvars*sizeof(ulong));
    Gshift = (ulong *) flint_malloc(nvars*sizeof(ulong));

    mpoly_gcd_info_limits(Amax_exp, Amin_exp, Amax_exp_count, Amin_exp_count,
                                      A->exps, A->bits, A->length, ctx->minfo);
    mpoly_gcd_info_limits(Bmax_exp, Bmin_exp, Bmax_exp_count, Bmin_exp_count,
                                      B->exps, B->bits, B->length, ctx->minfo);

    /* set ess(p) := p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    if (A->length == B->length)
    {
        success = _fmpz_mpoly_gcd_monomial_cofactors_sp(G, Gbits,
                                                   A, Amax_exp, Amin_exp,
                                                   B, Bmax_exp, Bmin_exp, ctx);
        if (success)
        {
            goto cleanup;
        }
    }

    /* check if ess(A) and ess(B) have a variable v_in_both in common */
    v_in_both = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (Amax_exp[j] > Amin_exp[j] && Bmax_exp[j] > Bmin_exp[j])
        {
            v_in_both = j;
            break;
        }
    }
    if (v_in_both == -WORD(1))
    {
        /*
            The variables in ess(A) and ess(B) are disjoint.
            gcd is trivial to compute.
        */
        fmpz_t gA, gB;

        fmpz_init(gA);
        fmpz_init(gB);
        _fmpz_vec_content(gA, A->coeffs, A->length);
        _fmpz_vec_content(gB, B->coeffs, B->length);

        for (j = 0; j < nvars; j++)
            Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        fmpz_mpoly_fit_length(G, 1, ctx);
        fmpz_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, Gshift, Gbits, ctx->minfo);
        fmpz_gcd(G->coeffs + 0, gA, gB);
        _fmpz_mpoly_set_length(G, 1, ctx);

        fmpz_clear(gA);
        fmpz_clear(gB);

        success = 1;
        goto cleanup;
    }

    /* check if ess(A) and ess(B) depend on another variable v_in_either */
    FLINT_ASSERT(0 <= v_in_both);
    FLINT_ASSERT(v_in_both < nvars);

    v_in_either = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (j == v_in_both)
            continue;

        if (Amax_exp[j] > Amin_exp[j] || Bmax_exp[j] > Bmin_exp[j])
        {
            v_in_either = j;
            break;
        }
    }

    if (v_in_either == -WORD(1))
    {
        /*
            The ess(A) and ess(B) depend on only one variable v_in_both
            Calculate gcd using univariates
        */
        fmpz_poly_t a, b, g;

        for (j = 0; j < nvars; j++)
            Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        _fmpz_mpoly_to_fmpz_poly_deflate(a, A, v_in_both, Amin_exp, Gstride, ctx);
        _fmpz_mpoly_to_fmpz_poly_deflate(b, B, v_in_both, Bmin_exp, Gstride, ctx);
        fmpz_poly_gcd(g, a, b);
        _fmpz_mpoly_from_fmpz_poly_inflate(G, Gbits, g, v_in_both,
                                                          Gshift, Gstride, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);

        success = 1;
        goto cleanup;
    }

    /* check if there is a variable in ess(A) that is not in ess(B) */
    v_in_A_only = -WORD(1);
    v_in_B_only = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (Amax_exp[j] > Amin_exp[j] && Bmax_exp[j] == Bmin_exp[j])
        {
            v_in_A_only = j;
            break;
        }
        if (Bmax_exp[j] > Bmin_exp[j] && Amax_exp[j] == Amin_exp[j])
        {
            v_in_B_only = j;
            break;
        }
    }
    if (v_in_A_only != -WORD(1))
    {
        success = _try_missing_var(G, Gbits, v_in_A_only,
                                                A, Amin_exp[v_in_A_only],
                                                B, Bmin_exp[v_in_A_only], ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, Gbits, v_in_B_only,
                                                B, Bmin_exp[v_in_B_only],
                                                A, Amin_exp[v_in_B_only], ctx);
        goto cleanup;
    }

    /* all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
        and there are at least two in the latter case
    */

    success = _try_prs(G, Gbits,
                   A, Amax_exp, Amin_exp, Amax_exp_count, Amin_exp_count,
                   B, Bmax_exp, Bmin_exp, Bmax_exp_count, Bmin_exp_count, ctx);
    if (success)
        goto cleanup;

    mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

    success = _try_brown(G, Gbits, Gstride, A, Amax_exp, Amin_exp,
                                            B, Bmax_exp, Bmin_exp, ctx);
    if (success)
        goto cleanup;

    success = _try_berlekamp_massey(G, Gbits, Gstride,
                   A, Amax_exp, Amin_exp, Amax_exp_count, Amin_exp_count,
                   B, Bmax_exp, Bmin_exp, Bmax_exp_count, Bmin_exp_count, ctx);
    if (success)
        goto cleanup;

    success = _try_zippel(G, Gbits, Gstride,
                   A, Amax_exp, Amin_exp, Amax_exp_count, Amin_exp_count,
                   B, Bmax_exp, Bmin_exp, Bmax_exp_count, Bmin_exp_count, ctx);

cleanup:

    flint_free(Amax_exp);
    flint_free(Amin_exp);
    flint_free(Bmax_exp);
    flint_free(Bmin_exp);
    flint_free(Amax_exp_count);
    flint_free(Amin_exp_count);
    flint_free(Bmax_exp_count);
    flint_free(Bmin_exp_count);
    flint_free(Gstride);
    flint_free(Gshift);

    return success;
}


int fmpz_mpoly_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    mp_bitcnt_t Gbits;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (B->length == 0)
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        }
        else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        /* usual gcd's go right down here */
        return _fmpz_mpoly_gcd(G, Gbits, A, B, ctx);
    }

    if (A->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _fmpz_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }
    else if (_fmpz_mpoly_gcd_monomial_cofactors(G, A, B, ctx))
    {
        return 1;
    }
    else
    {
        /*
            The gcd calculation is unusual.
            First see if both inputs fit into FLINT_BITS.
            Then, try deflation as a last resort.
        */

        int success;
        int useAnew = 0;
        int useBnew = 0;
        slong k;
        fmpz * Ashift, * Astride;
        fmpz * Bshift, * Bstride;
        fmpz * Gshift, * Gstride;
        fmpz_mpoly_t Anew;
        fmpz_mpoly_t Bnew;

        fmpz_mpoly_init(Anew, ctx);
        fmpz_mpoly_init(Bnew, ctx);

        if (A->bits > FLINT_BITS)
        {
            useAnew = fmpz_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
                goto could_not_repack;
        }

        if (B->bits > FLINT_BITS)
        {
            useBnew = fmpz_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
                goto could_not_repack;
        }

        success = _fmpz_mpoly_gcd(G, FLINT_BITS, useAnew ? Anew : A,
                                                 useBnew ? Bnew : B, ctx);
        goto cleanup;

could_not_repack:

        /*
            One of A or B could not be repacked into FLINT_BITS. See if
            they both fit into FLINT_BITS after deflation.
        */

        Ashift  = _fmpz_vec_init(ctx->minfo->nvars);
        Astride = _fmpz_vec_init(ctx->minfo->nvars);
        Bshift  = _fmpz_vec_init(ctx->minfo->nvars);
        Bstride = _fmpz_vec_init(ctx->minfo->nvars);
        Gshift  = _fmpz_vec_init(ctx->minfo->nvars);
        Gstride = _fmpz_vec_init(ctx->minfo->nvars);

        fmpz_mpoly_deflation(Ashift, Astride, A, ctx);
        fmpz_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        fmpz_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!fmpz_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        fmpz_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!fmpz_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        success = _fmpz_mpoly_gcd(G, FLINT_BITS, Anew, Bnew, ctx);

        if (success)
        {
            fmpz_mpoly_inflate(G, G, Gshift, Gstride, ctx);

            /* inflation may have changed the lc */
            FLINT_ASSERT(G->length > 0);
            if (fmpz_sgn(G->coeffs + 0) < 0)
            {
                fmpz_mpoly_neg(G, G, ctx);
            }
        }

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

        fmpz_mpoly_clear(Anew, ctx);
        fmpz_mpoly_clear(Bnew, ctx);

        return success;
    }
}
