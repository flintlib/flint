/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    Assume when B is converted to univar format, its length would be one.
    Gcd is gcd of coefficients of univar(A) and B (modulo some shifts).
*/
static int _try_missing_var(nmod_mpoly_t G, mp_bitcnt_t Gbits, slong var,
                                         const nmod_mpoly_t A, ulong Ashift,
                                         const nmod_mpoly_t B, ulong Bshift,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i;
    nmod_mpoly_t tG;
    nmod_mpoly_univar_t Ax;

    nmod_mpoly_init(tG, ctx);
    nmod_mpoly_univar_init(Ax, ctx);

    success = nmod_mpoly_to_univar(Ax, A, var, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Ax->length > 0);
    success = _nmod_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _nmod_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }

    nmod_mpoly_swap(G, tG, ctx);
    _mpoly_gen_shift_left(G->exps, G->bits, G->length,
                                   var, FLINT_MIN(Ashift, Bshift), ctx->minfo);

cleanup:

    nmod_mpoly_clear(tG, ctx);
    nmod_mpoly_univar_clear(Ax, ctx);

    return success;
}


/*
    return 1 for success or 0 for failure
*/
static int _try_zippel(nmod_mpoly_t G, mp_bitcnt_t Gbits, ulong * Gstride,
         const nmod_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
                   const slong * Amax_exp_count, const slong * Amin_exp_count,
         const nmod_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                   const slong * Bmax_exp_count, const slong * Bmin_exp_count,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong n, m;
    int success;
    ulong * Gshift, * Addeg, * Bddeg;
    mpoly_zipinfo_t zinfo;
    mp_bitcnt_t new_bits;
    flint_rand_t randstate;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    nmod_mpoly_t Acontent, Bcontent;
    nmod_mpolyu_t Abar, Bbar, Gbar;
    TMP_INIT;

    TMP_START;

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
    m = -WORD(1);
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
    FLINT_ASSERT(m >= 1);

    /* interpolation will continue in m + 1 variables */
    mpoly_zipinfo_init(zinfo, m + 1);

    /* uctx is context for Zp[y_0,...,y_{m-1}]*/
    nmod_mpoly_ctx_init(uctx, m, ORD_LEX, ctx->ffinfo->mod.n);

    Gshift = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    /* degrees after deflation */
    Addeg = (ulong *) TMP_ALLOC(n*sizeof(ulong));
    Bddeg = (ulong *) TMP_ALLOC(n*sizeof(ulong));

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
    FLINT_ASSERT(i == m + 1);

    /* work out a favourable permation to zinfo->perm */

    /* figure out a main variable y_m */
    {
        slong main_var;
        ulong count, deg, new_count, new_deg;

        main_var = m;
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

        if (main_var != m)
        {
            slong t = zinfo->perm[main_var];
            zinfo->perm[main_var] = zinfo->perm[m];
            zinfo->perm[m] = t;
        }
    }

    /* sort with hope that ddeg(G,y_0) >= ddeg(G,y_1) ... >= ddeg(G,y_{m-1}) */
    for (k = 0; k + 1 < m; k++)
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
    for (i = 0; i <= m; i++)
    {
        j = zinfo->perm[i];
        FLINT_ASSERT(Addeg[j] != 0);
        FLINT_ASSERT(Bddeg[j] != 0);
        zinfo->Adegs[i] = Addeg[j];
        zinfo->Bdegs[i] = Bddeg[j];
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    nmod_mpolyu_init(Au, new_bits, uctx);
    nmod_mpolyu_init(Bu, new_bits, uctx);
    nmod_mpolyu_init(Gu, new_bits, uctx);

    nmod_mpoly_to_mpolyu_perm_deflate(Au, A,
                                   zinfo->perm, Amin_exp, Gstride, uctx, ctx);
    nmod_mpoly_to_mpolyu_perm_deflate(Bu, B,
                                   zinfo->perm, Bmin_exp, Gstride, uctx, ctx);


    FLINT_ASSERT(Au->bits == Bu->bits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    nmod_mpoly_init(Acontent, uctx);
    nmod_mpoly_init(Bcontent, uctx);
    nmod_mpolyu_init(Abar, Au->bits, uctx);
    nmod_mpolyu_init(Bbar, Au->bits, uctx);
    nmod_mpolyu_init(Gbar, Au->bits, uctx);

    /* compute content of A */
    success = _nmod_mpoly_gcd(Acontent, new_bits, Au->coeffs + 0,
                                                  Au->coeffs + 1, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == new_bits);

    for (i = 2; i < Au->length; i++)
    {
        success = _nmod_mpoly_gcd(Acontent, new_bits, Acontent,
                                                      Au->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Acontent->bits == new_bits);
    }

    /* compute content of B */
    success = _nmod_mpoly_gcd(Bcontent, new_bits, Bu->coeffs + 0,
                                                  Bu->coeffs + 1, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Bcontent->bits == new_bits);
    for (i = 2; i < Bu->length; i++)
    {
        success = _nmod_mpoly_gcd(Bcontent, new_bits, Bcontent,
                                                      Bu->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Bcontent->bits == new_bits);
    }

    /* remove content from A and B */
    nmod_mpolyu_divexact_mpoly(Abar, Au, Acontent, uctx);
    nmod_mpolyu_divexact_mpoly(Bbar, Bu, Bcontent, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */

    /* compute GCD */
    success = nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _nmod_mpoly_gcd(Acontent, new_bits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == new_bits);

    nmod_mpolyu_mul_mpoly(Gu, Gbar, Acontent, uctx);
    nmod_mpoly_from_mpolyu_perm_inflate(G, Gbits, Gu,
                                      zinfo->perm, Gshift, Gstride, uctx, ctx);
    nmod_mpoly_make_monic(G, G, ctx);
    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    nmod_mpolyu_clear(Abar, uctx);
    nmod_mpolyu_clear(Bbar, uctx);
    nmod_mpolyu_clear(Gbar, uctx);
    nmod_mpoly_clear(Acontent, uctx);
    nmod_mpoly_clear(Bcontent, uctx);

    nmod_mpolyu_clear(Au, uctx);
    nmod_mpolyu_clear(Bu, uctx);
    nmod_mpolyu_clear(Gu, uctx);
    nmod_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    TMP_END;

    return success;
}

/*
    return 1 for success or 0 for failure
*/
static int _try_brown(nmod_mpoly_t G, mp_bitcnt_t Gbits, ulong * Gstride,
         const nmod_mpoly_t A, const ulong * Amax_exp, const ulong * Amin_exp,
         const nmod_mpoly_t B, const ulong * Bmax_exp, const ulong * Bmin_exp,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong j;
    slong denseAsize, denseBsize;
    slong m, n = ctx->minfo->nvars;
    slong * perm;
    ulong * Gshift;
    nmod_mpolyd_t Ad, Bd, Gd, Abar, Bbar;
    TMP_INIT;

    /*
        Let the variables in A and B be
            x_0, x_1, ..., x_{n-1}
        where n = ctx->minfo->nvars, and x_0 is most significant

        Recall from _fmpz_mpoly_gcd that all variables are either
            missing from both ess(A) and ess(B) (non-essential), or
            present in both ess(A) and ess(B) (essential)
        and that there are at least 2 variables in the essential case.

        Let y_0, ..., y_{m-1} with m >= 2 denote the variables present
        in both ess(A) and ess(B). Each y_i is one of the x_j and the variables
        are ordered as y_0 > ... > y_{m-1} with LEX order.
        The Brown algorithm will operate in Z[y_0,...,y_{m-1}] and it
        only operates with Z[y_0,...,y_{m-1}] in LEX because the coefficients
        are stored in a dense format.

        When converting to the mpolyd format via
        fmpz_mpoly_to_mpolyd_perm_deflate, the non-essential variables
        will be immediately striped out and the remaining variables will be
        mapped according to the permutation in perm as

            y_k = x_perm[k] ^ Gstride[perm[k]]

        When converting out of the mpolyd format via
        fmpz_mpoly_from_mpolyd_perm_inflate, the contribution of the
        non-essential variables will be put back in.
    */

    /* first see if dense representation is a good idea */
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
            perm[m] = j;
            m++;
        }
    }
    FLINT_ASSERT(m >= 2);

    /* TODO: see how performance depends on perm and try reorder */

    nmod_mpolyd_init(Ad, m);
    nmod_mpolyd_init(Bd, m);
    nmod_mpolyd_init(Gd, m);
    nmod_mpolyd_init(Abar, m);
    nmod_mpolyd_init(Bbar, m);

    nmod_mpoly_to_nmod_mpolyd_perm_deflate(Ad, m, A,
                                       perm, Amin_exp, Gstride, Amax_exp, ctx);
    nmod_mpoly_to_nmod_mpolyd_perm_deflate(Bd, m, B,
                                       perm, Bmin_exp, Gstride, Bmax_exp, ctx);

    success = nmod_mpolyd_gcd_brown_smprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
    if (!success)
    {
        nmod_mpoly_to_nmod_mpolyd_perm_deflate(Ad, m, A,
                                       perm, Amin_exp, Gstride, Amax_exp, ctx);
        nmod_mpoly_to_nmod_mpolyd_perm_deflate(Bd, m, B,
                                       perm, Bmin_exp, Gstride, Bmax_exp, ctx);
        success = nmod_mpolyd_gcd_brown_lgprime(Gd, Abar, Bbar, Ad, Bd, ctx->ffinfo);
    }
    if (!success)
        goto cleanup;

    nmod_mpoly_from_nmod_mpolyd_perm_inflate(G, Gbits, ctx, Gd,
                                                        perm, Gshift, Gstride);
    nmod_mpoly_make_monic(G, G, ctx);
    FLINT_ASSERT(G->bits == Gbits);

    success = 1;

cleanup:

    nmod_mpolyd_clear(Bbar);
    nmod_mpolyd_clear(Abar);
    nmod_mpolyd_clear(Gd);
    nmod_mpolyd_clear(Bd);
    nmod_mpolyd_clear(Ad);

    TMP_END;

    return success;
}



/*
    The function must pack its answer into bits = Gbits <= FLINT_BITS
    Both A and B have to be packed into bits <= FLINT_BITS

    return is 1 for success, 0 for failure.
*/
int _nmod_mpoly_gcd(nmod_mpoly_t G, mp_bitcnt_t Gbits,
                               const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
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
        return _nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
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
        success = _nmod_mpoly_gcd_monomial_cofactors_sp(G, Gbits,
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

        nmod_mpoly_fit_length(G, 1, ctx);
        nmod_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, minexps, Gbits, ctx->minfo);
        G->coeffs[0] = UWORD(1);
        _nmod_mpoly_set_length(G, 1, ctx);

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
        nmod_poly_t a, b, g;

        Gshift = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
            Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        Gstride = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

        nmod_poly_init_preinv(a, ctx->ffinfo->mod.n, ctx->ffinfo->mod.ninv);
        nmod_poly_init_preinv(b, ctx->ffinfo->mod.n, ctx->ffinfo->mod.ninv);
        nmod_poly_init_preinv(g, ctx->ffinfo->mod.n, ctx->ffinfo->mod.ninv);
        _nmod_mpoly_to_nmod_poly_deflate(a, A, v_in_both, Amin_exp, Gstride, ctx);
        _nmod_mpoly_to_nmod_poly_deflate(b, B, v_in_both, Bmin_exp, Gstride, ctx);
        nmod_poly_gcd(g, a, b);
        _nmod_mpoly_from_nmod_poly_inflate(G, Gbits, g, v_in_both,
                                                          Gshift, Gstride, ctx);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(g);

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


int nmod_mpoly_gcd(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    mp_bitcnt_t Gbits;

    if (nmod_mpoly_is_zero(A, ctx))
    {
        if (nmod_mpoly_is_zero(B, ctx))
            nmod_mpoly_zero(G, ctx);
        else
            nmod_mpoly_make_monic(G, B, ctx);
        return 1;
    }

    if (nmod_mpoly_is_zero(B, ctx))
    {
        nmod_mpoly_make_monic(G, A, ctx);
        return 1;
    }

    Gbits = FLINT_MIN(A->bits, B->bits);

    if (A->bits <= FLINT_BITS && B->bits <= FLINT_BITS)
    {
        /* usual gcd's go right down here */
        return _nmod_mpoly_gcd(G, Gbits, A, B, ctx);
    }

    if (A->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
    }
    else if (_nmod_mpoly_gcd_monomial_cofactors(G, A, B, ctx))
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
        nmod_mpoly_t Anew;
        nmod_mpoly_t Bnew;

        nmod_mpoly_init(Anew, ctx);
        nmod_mpoly_init(Bnew, ctx);

        if (A->bits > FLINT_BITS)
        {
            useAnew = nmod_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
                goto could_not_repack;
        }

        if (B->bits > FLINT_BITS)
        {
            useBnew = nmod_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
                goto could_not_repack;
        }

        success = _nmod_mpoly_gcd(G, FLINT_BITS, useAnew ? Anew : A,
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

        nmod_mpoly_deflation(Ashift, Astride, A, ctx);
        nmod_mpoly_deflation(Bshift, Bstride, B, ctx);
        _fmpz_vec_min(Gshift, Ashift, Bshift, ctx->minfo->nvars);
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_gcd(Gstride + k, Astride + k, Bstride + k);
        }

        success = 0;

        nmod_mpoly_deflate(Anew, A, Ashift, Gstride, ctx);
        if (Anew->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Anew, Anew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        nmod_mpoly_deflate(Bnew, B, Bshift, Gstride, ctx);
        if (Bnew->bits > FLINT_BITS)
        {
            if (!nmod_mpoly_repack_bits(Bnew, Bnew, FLINT_BITS, ctx))
                goto deflate_cleanup;
        }

        success = _nmod_mpoly_gcd(G, FLINT_BITS, Anew, Bnew, ctx);

        if (success)
        {
            nmod_mpoly_inflate(G, G, Gshift, Gstride, ctx);
            nmod_mpoly_make_monic(G, G, ctx);
        }

deflate_cleanup:

        _fmpz_vec_clear(Ashift, ctx->minfo->nvars);
        _fmpz_vec_clear(Astride, ctx->minfo->nvars);
        _fmpz_vec_clear(Bshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Bstride, ctx->minfo->nvars);
        _fmpz_vec_clear(Gshift, ctx->minfo->nvars);
        _fmpz_vec_clear(Gstride, ctx->minfo->nvars);

cleanup:

        nmod_mpoly_clear(Anew, ctx);
        nmod_mpoly_clear(Bnew, ctx);

        return success;
    }
}
