/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

/*
    Assume when B is converted to univar format, its length would be one.
    Gcd is gcd of coefficients of univar(A) and B (modulo some shifts).
*/
static int _try_missing_var(fq_nmod_mpoly_t G, flint_bitcnt_t Gbits, slong var,
                                        const fq_nmod_mpoly_t A, ulong Ashift,
                                        const fq_nmod_mpoly_t B, ulong Bshift,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    fq_nmod_mpoly_t tG;
    fq_nmod_mpoly_univar_t Ax;

    fq_nmod_mpoly_init(tG, ctx);
    fq_nmod_mpoly_univar_init(Ax, ctx);

    success = fq_nmod_mpoly_to_univar(Ax, A, var, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Ax->length > 0);
    success = _fq_nmod_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _fq_nmod_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }

    fq_nmod_mpoly_swap(G, tG, ctx);
    _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

cleanup:

    fq_nmod_mpoly_clear(tG, ctx);
    fq_nmod_mpoly_univar_clear(Ax, ctx);

    return success;
}


/*
    return 1 for success or 0 for failure
*/
static int _try_zippel(fq_nmod_mpoly_t G, flint_bitcnt_t Gbits, ulong * Gstride,
      const fq_nmod_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
                   const slong * Amax_exp_count, const slong * Amin_exp_count,
      const fq_nmod_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                   const slong * Bmax_exp_count, const slong * Bmin_exp_count,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong n, m;
    int success;
    ulong * Gshift, * Addeg, * Bddeg;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t ABbits;
    flint_rand_t randstate;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyu_t Au, Bu, Gu;
    fq_nmod_mpoly_t Acontent, Bcontent;
    fq_nmod_mpolyu_t Abar, Bbar, Gbar;

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

        Let y_m, y_0, ..., y_{m-1} with m >= 1 denote the variables present
        in both ess(A) and ess(B). Each y_i is one of the x_j and the variables
        are ordered as y_{m} > y_0 > ... > y_{m-1} with LEX order.
        The Zippel algorithm will operate in Zp[y_0,...,y_{m-1}][y_m] and it
        only operates with Zp[y_0,...,y_{m-1}] in LEX.

        When converting to the mpolyu format via
        nmod_mpoly_to_mpolyu_perm_deflate, the non-essential variables
        will be immediately striped out and the remaining variables will be
        mapped according to the permutation in zinfo->perm as

            y_k = x_perm[k] ^ Gstride[perm[k]]

        When converting out of the mpolyu format via
        nmod_mpoly_from_mpolyu_perm_inflate, the contribution of the
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

    /* uctx is context for Fq[y_1,...,y_{m-1}]*/
    fq_nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->fqctx);

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
        for (i = 0; i < m; i++)
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

    /* sort with hope that ddeg(G,y_1) >= ddeg(G,y_2) ... >= ddeg(G,y_{m-1}) */
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

    fq_nmod_mpolyu_init(Au, ABbits, uctx);
    fq_nmod_mpolyu_init(Bu, ABbits, uctx);
    fq_nmod_mpolyu_init(Gu, ABbits, uctx);

    _fq_nmod_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx,
                                               zinfo->perm, Amin_exp, Gstride);
    _fq_nmod_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx,
                                               zinfo->perm, Bmin_exp, Gstride);

    FLINT_ASSERT(Au->bits == ABbits);
    FLINT_ASSERT(Bu->bits == ABbits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    fq_nmod_mpoly_init(Acontent, uctx);
    fq_nmod_mpoly_init(Bcontent, uctx);
    fq_nmod_mpolyu_init(Abar, ABbits, uctx);
    fq_nmod_mpolyu_init(Bbar, ABbits, uctx);
    fq_nmod_mpolyu_init(Gbar, ABbits, uctx);

    /* compute content of A */
    success = _fq_nmod_mpoly_gcd(Acontent, ABbits, Au->coeffs + 0,
                                                  Au->coeffs + 1, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == ABbits);

    for (i = 2; i < Au->length; i++)
    {
        success = _fq_nmod_mpoly_gcd(Acontent, ABbits, Acontent,
                                                      Au->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Acontent->bits == ABbits);
    }

    /* compute content of B */
    success = _fq_nmod_mpoly_gcd(Bcontent, ABbits, Bu->coeffs + 0,
                                                  Bu->coeffs + 1, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Bcontent->bits == ABbits);
    for (i = 2; i < Bu->length; i++)
    {
        success = _fq_nmod_mpoly_gcd(Bcontent, ABbits, Bcontent,
                                                      Bu->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Bcontent->bits == ABbits);
    }

    /* remove content from A and B */
    fq_nmod_mpolyu_divexact_mpoly(Abar, Au, Acontent, uctx);
    fq_nmod_mpolyu_divexact_mpoly(Bbar, Bu, Bcontent, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */

    /* compute GCD */
    success = fq_nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _fq_nmod_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);

    fq_nmod_mpolyu_mul_mpoly(Gu, Gbar, Acontent, uctx);
    _fq_nmod_mpoly_from_mpolyu_perm_inflate(G, Gbits, ctx, Gu, uctx,
                                                 zinfo->perm, Gshift, Gstride);
    fq_nmod_mpoly_make_monic(G, G, ctx);

    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    fq_nmod_mpolyu_clear(Abar, uctx);
    fq_nmod_mpolyu_clear(Bbar, uctx);
    fq_nmod_mpolyu_clear(Gbar, uctx);
    fq_nmod_mpoly_clear(Acontent, uctx);
    fq_nmod_mpoly_clear(Bcontent, uctx);

    fq_nmod_mpolyu_clear(Au, uctx);
    fq_nmod_mpolyu_clear(Bu, uctx);
    fq_nmod_mpolyu_clear(Gu, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);

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
static int _try_brown(fq_nmod_mpoly_t G, flint_bitcnt_t Gbits, ulong * Gstride,
      const fq_nmod_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
      const fq_nmod_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int success;
    slong j;
    slong denseAsize, denseBsize;
    slong m, n = ctx->minfo->nvars;
    slong * perm;
    ulong * Gshift;
    flint_bitcnt_t ABbits;
    fq_nmod_mpoly_ctx_t uctx;
    fq_nmod_mpolyun_t An, Bn, Gn, Abarn, Bbarn;

    /* first see if a dense algorithm is a good idea */

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

    Gshift = (ulong *) flint_malloc(n*sizeof(ulong));
    perm = (slong *) flint_malloc(n*sizeof(slong)); /* only first m entries used */

    /* fill in perm and set shift of GCD */
    m = 0;
    for (j = 0; j < n; j++)
    {
        Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);
        if (Amax_exp[j] > Amin_exp[j])
        {
            perm[m] = j;
            m++;
        }
    }

    /* need at least 2 variables */
    if (m < 2)
    {
        success = 0;
        goto cleanup1;
    }

    /* TODO: find a favourable permutation of perm */

    ABbits = FLINT_MAX(A->bits, B->bits);

    fq_nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->fqctx);

    fq_nmod_mpolyun_init(An, ABbits, uctx);
    fq_nmod_mpolyun_init(Bn, ABbits, uctx);
    fq_nmod_mpolyun_init(Gn, ABbits, uctx);
    fq_nmod_mpolyun_init(Abarn, ABbits, uctx);
    fq_nmod_mpolyun_init(Bbarn, ABbits, uctx);

    _fq_nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                                      perm, Amin_exp, Gstride);
    _fq_nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx,
                                                      perm, Bmin_exp, Gstride);
    success = fq_nmod_mpolyun_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 2, uctx);
    if (!success)
    {
        _fq_nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                                      perm, Amin_exp, Gstride);
        _fq_nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx,
                                                      perm, Bmin_exp, Gstride);
        success = fq_nmod_mpolyun_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 2, uctx);
    }
    if (success)
    {
        _fq_nmod_mpoly_from_mpolyun_perm_inflate(G, Gbits, ctx, Gn, uctx,
                                                        perm, Gshift, Gstride);
        fq_nmod_mpoly_make_monic(G, G, ctx);
        FLINT_ASSERT(G->bits == Gbits);
    }

    fq_nmod_mpolyun_clear(An, uctx);
    fq_nmod_mpolyun_clear(Bn, uctx);
    fq_nmod_mpolyun_clear(Gn, uctx);
    fq_nmod_mpolyun_clear(Abarn, uctx);
    fq_nmod_mpolyun_clear(Bbarn, uctx);
    fq_nmod_mpoly_ctx_clear(uctx);

cleanup1:

    flint_free(Gshift);
    flint_free(perm);

    return success;
}



/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, flint_bitcnt_t Gbits,
                            const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
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
    ulong * Gstride;
    TMP_INIT;

    if (A->length == 1)
    {
        return _fq_nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _fq_nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }

    TMP_START;

    Amax_exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    Amin_exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    Bmax_exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    Bmin_exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    Amax_exp_count = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Amin_exp_count = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bmax_exp_count = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bmin_exp_count = (slong *) TMP_ALLOC(nvars*sizeof(slong));

    mpoly_gcd_info_limits(Amax_exp, Amin_exp, Amax_exp_count, Amin_exp_count,
                                      A->exps, A->bits, A->length, ctx->minfo);
    mpoly_gcd_info_limits(Bmax_exp, Bmin_exp, Bmax_exp_count, Bmin_exp_count,
                                      B->exps, B->bits, B->length, ctx->minfo);

    /* set ess(p) := p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    if (A->length == B->length)
    {
        success = _fq_nmod_mpoly_gcd_monomial_cofactors_sp(G, Gbits,
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
        ulong * minexps;

        minexps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
            minexps[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        fq_nmod_mpoly_fit_length(G, 1, ctx);
        fq_nmod_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, minexps, Gbits, ctx->minfo);
        fq_nmod_one(G->coeffs + 0, ctx->fqctx);
        _fq_nmod_mpoly_set_length(G, 1, ctx);

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
        ulong * Gshift;
        fq_nmod_poly_t a, b, g;

        Gshift = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
            Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        Gstride = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

        fq_nmod_poly_init(a, ctx->fqctx);
        fq_nmod_poly_init(b, ctx->fqctx);
        fq_nmod_poly_init(g, ctx->fqctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(a, A, v_in_both, Amin_exp, Gstride, ctx);
        _fq_nmod_mpoly_to_fq_nmod_poly_deflate(b, B, v_in_both, Bmin_exp, Gstride, ctx);
        fq_nmod_poly_gcd(g, a, b, ctx->fqctx);
        _fq_nmod_mpoly_from_fq_nmod_poly_inflate(G, Gbits, g, v_in_both,
                                                         Gshift, Gstride, ctx);
        fq_nmod_poly_clear(a, ctx->fqctx);
        fq_nmod_poly_clear(b, ctx->fqctx);
        fq_nmod_poly_clear(g, ctx->fqctx);

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

    Gstride = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
    mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

    success = _try_brown(G, Gbits, Gstride, A, Amax_exp, Amin_exp,
                                            B, Bmax_exp, Bmin_exp, ctx);
    if (success)
        goto cleanup;

    success = _try_zippel(G, Gbits, Gstride,
                   A, Amax_exp, Amin_exp, Amax_exp_count, Amin_exp_count,
                   B, Bmax_exp, Bmin_exp, Bmax_exp_count, Bmin_exp_count, ctx);
cleanup:

    TMP_END;

    return success;
}


int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A,
                       const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t Gbits;

    if (fq_nmod_mpoly_is_zero(A, ctx))
    {
        if (fq_nmod_mpoly_is_zero(B, ctx))
            fq_nmod_mpoly_zero(G, ctx);
        else
            fq_nmod_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (fq_nmod_mpoly_is_zero(B, ctx))
    {
        fq_nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        /* usual gcd's go right down here */
        return _fq_nmod_mpoly_gcd(G, Gbits, A, B, ctx);
    }

    if (A->length == 1)
    {
        return _fq_nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _fq_nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }
    else if (_fq_nmod_mpoly_gcd_monomial_cofactors(G, A, B, ctx))
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
        fq_nmod_mpoly_t Anew;
        fq_nmod_mpoly_t Bnew;

        fq_nmod_mpoly_init(Anew, ctx);
        fq_nmod_mpoly_init(Bnew, ctx);

        if (A->bits > FLINT_BITS)
        {
            useAnew = fq_nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
                goto could_not_repack;
        }

        if (B->bits > FLINT_BITS)
        {
            useBnew = fq_nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
                goto could_not_repack;
        }

        success = _fq_nmod_mpoly_gcd(G, FLINT_BITS, useAnew ? Anew : A,
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

        fq_nmod_mpoly_deflation(Ashift, Astride, A, ctx);
        fq_nmod_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        fq_nmod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!fq_nmod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        fq_nmod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!fq_nmod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        success = _fq_nmod_mpoly_gcd(G, FLINT_BITS, Anew, Bnew, ctx);

        if (success)
        {
            fq_nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
            fq_nmod_mpoly_make_monic(G, G, ctx);
        }

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

        fq_nmod_mpoly_clear(Anew, ctx);
        fq_nmod_mpoly_clear(Bnew, ctx);

        return success;
    }
}
