/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    Scan A and fill in the min and max exponents of each variable along
    with the count of terms attached to each.
*/
static void _init_info(mpoly_gcd_var_info_struct * info, const fmpz_mpoly_t A,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    ulong * Aexp;
    slong i, j, N;
    slong nvars = ctx->minfo->nvars;
    TMP_INIT;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    TMP_START;

    Aexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    i = 0;
    mpoly_get_monomial_ui(Aexp, A->exps + N*i, A->bits, ctx->minfo);
    for (j = 0; j < nvars; j++)
    {
        info[j].min_exp = Aexp[j];
        info[j].max_exp = Aexp[j];
        info[j].min_exp_term_count = 1;
        info[j].max_exp_term_count = 1;
    }
    for (i = 1; i < A->length; i++)
    {
        mpoly_get_monomial_ui(Aexp, A->exps + N*i, A->bits, ctx->minfo);

        for (j = 0; j < nvars; j++)
        {
            if (info[j].min_exp > Aexp[j])
            {
                info[j].min_exp = Aexp[j];
                info[j].min_exp_term_count = 1;            
            }
            else if (info[j].min_exp == Aexp[j])
            {
                info[j].min_exp_term_count += 1;
            }

            if (info[j].max_exp < Aexp[j])
            {
                info[j].max_exp = Aexp[j];
                info[j].max_exp_term_count = 1;            
            }
            else if (info[j].max_exp == Aexp[j])
            {
                info[j].max_exp_term_count += 1;
            }
        }
    }

    TMP_END;
}


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
    ulong Gshift;
    fmpz_mpoly_t tG;
    fmpz_mpoly_univar_t Ax;

    fmpz_mpoly_init(tG, ctx);
    fmpz_mpoly_univar_init(Ax, ctx);

    success = fmpz_mpoly_to_univar(Ax, A, var, ctx);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Ax->length > 0);

    if (Bshift == 0)
    {
        success = _fmpz_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx);
    }
    else
    {
        fmpz_mpoly_t Bs;
        fmpz_mpoly_init(Bs, ctx);
        fmpz_mpoly_set(Bs, B, ctx);
        _fmpz_mpoly_gen_shift_right(Bs, var, Bshift, ctx);
        success = _fmpz_mpoly_gcd(tG, Gbits, Bs, Ax->coeffs + 0, ctx);
        fmpz_mpoly_clear(Bs, ctx);
    }

    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _fmpz_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx);
        if (!success)
            goto cleanup;
    }

    fmpz_mpoly_swap(G, tG, ctx);

    Gshift = FLINT_MIN(Ashift, Bshift);
    if (Gshift != 0)
    {
        _fmpz_mpoly_gen_shift_left(G, var, Gshift, ctx);
    }

cleanup:

    fmpz_mpoly_clear(tG, ctx);
    fmpz_mpoly_univar_clear(Ax, ctx);

    return success;
}


static int _try_zippel(fmpz_mpoly_t G, mp_bitcnt_t Gbits,
                 const fmpz_mpoly_t A, const mpoly_gcd_var_info_struct * Ainfo,
                 const fmpz_mpoly_t B, const mpoly_gcd_var_info_struct * Binfo,
                                                   const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    slong n, m;
    int success;
    ulong * defaults;
    mpoly_zipinfo_t zinfo;
    mp_bitcnt_t new_bits;
    flint_rand_t randstate;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu;
    slong ABminshift;
    fmpz_mpoly_t Acontent, Bcontent;
    fmpz_mpolyu_t Abar, Bbar, Gbar;
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
        The Zippel algorithm will operate in Z[y_0,...,y_{m-1}][y_m] and it
        only operate with Z[y_0,...,y_{m-1}] in LEX.

        When converting to the mpolyu format, the non-essential variables
        will be immediately striped out and the remaining variables will be
        mapped according to the permutation in zinfo->perm as
            x_i = y_{perm[i]}

        When converting out of the mpolyu format, the contribution of the
        non-essential variables needs to be put back in.
    */
    n = ctx->minfo->nvars;
    m = -WORD(1);
    for (j = 0; j < n; j++)
    {
        if (Ainfo[j].max_exp > Ainfo[j].min_exp)
        {
            FLINT_ASSERT(Binfo[j].max_exp > Binfo[j].min_exp);
            m++;
        }
        else
        {
            FLINT_ASSERT(Binfo[j].max_exp == Binfo[j].min_exp);
        }
    }
    FLINT_ASSERT(m >= 1);

    /* interpolation will continue in m + 1 variables */
    mpoly_zipinfo_init(zinfo, m + 1);

    /* uctx is context for Z[y_0,...,y_{m-1}]*/
    fmpz_mpoly_ctx_init(uctx, m, ORD_LEX);

    defaults = (ulong *) TMP_ALLOC(n*sizeof(ulong));

    /* fill in zinfo->perm and set default values to use when converting */
    i = 0;
    for (j = 0; j < n; j++)
    {
        if (Ainfo[j].max_exp > Ainfo[j].min_exp)
        {
            zinfo->perm[i] = j;
            zinfo->Adegs[i] = Ainfo[j].max_exp - Ainfo[j].min_exp;
            zinfo->Bdegs[i] = Binfo[j].max_exp - Binfo[j].min_exp;
            i++;
            /* default[j] should not matter in this case - set to absurd */
            defaults[j] = -UWORD(1);
        }
        else
        {
            defaults[j] = FLINT_MIN(Ainfo[j].min_exp, Binfo[j].min_exp);
        }
    }
    FLINT_ASSERT(i == m + 1);

    /* TODO: see how performance depends on zinfo->perm and try reorder */

    new_bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);

    fmpz_mpoly_to_mpolyu_perm_new(Au, A, zinfo->perm, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm_new(Bu, B, zinfo->perm, uctx, ctx);

    FLINT_ASSERT(Au->bits == Bu->bits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    fmpz_mpoly_init(Acontent, uctx);
    fmpz_mpoly_init(Bcontent, uctx);
    fmpz_mpolyu_init(Abar, Au->bits, uctx);
    fmpz_mpolyu_init(Bbar, Au->bits, uctx);
    fmpz_mpolyu_init(Gbar, Au->bits, uctx);

    /* compute content of A */
    success = _fmpz_mpoly_gcd(Acontent, new_bits, Au->coeffs + 0,
                                                  Au->coeffs + 1, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == new_bits);

    for (i = 2; i < Au->length; i++)
    {
        success = _fmpz_mpoly_gcd(Acontent, new_bits, Acontent,
                                                      Au->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Acontent->bits == new_bits);
    }

    /* compute content of B */
    success = _fmpz_mpoly_gcd(Bcontent, new_bits, Bu->coeffs + 0,
                                                  Bu->coeffs + 1, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Bcontent->bits == new_bits);
    for (i = 2; i < Bu->length; i++)
    {
        success = _fmpz_mpoly_gcd(Bcontent, new_bits, Bcontent,
                                                      Bu->coeffs + i, uctx);
        if (!success)
            goto cleanup;
        FLINT_ASSERT(Bcontent->bits == new_bits);
    }

    /* remove content from A and B */
    fmpz_mpolyu_divexact_mpoly(Abar, Au, Acontent, uctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, Bu, Bcontent, uctx);
    ABminshift = FLINT_MIN(Abar->exps[Abar->length - 1],
                           Bbar->exps[Bbar->length - 1]);
    fmpz_mpolyu_shift_right(Abar, Abar->exps[Abar->length - 1]);
    fmpz_mpolyu_shift_right(Bbar, Bbar->exps[Bbar->length - 1]);

    /* compute GCD */
    success = fmpz_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _fmpz_mpoly_gcd(Acontent, new_bits, Acontent, Bcontent, uctx);
    if (!success)
        goto cleanup;
    FLINT_ASSERT(Acontent->bits == new_bits);

    fmpz_mpolyu_shift_left(Gbar, ABminshift);
    fmpz_mpolyu_mul_mpoly(Gu, Gbar, Acontent, uctx);
    fmpz_mpoly_from_mpolyu_perm_new(G, Gbits, Gu, zinfo->perm, defaults, uctx, ctx);
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

    TMP_END;

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
    mpoly_gcd_var_info_struct * Ainfo, * Binfo;
    TMP_INIT;

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

    TMP_START;

    Ainfo = (mpoly_gcd_var_info_struct *) TMP_ALLOC(nvars
                                           *sizeof(mpoly_gcd_var_info_struct));
    Binfo = (mpoly_gcd_var_info_struct *) TMP_ALLOC(nvars
                                           *sizeof(mpoly_gcd_var_info_struct));
    _init_info(Ainfo, A, ctx);
    _init_info(Binfo, B, ctx);

    /* set ess(p) = p/term_content(p) */

    /* check if the cofactors could be monomials, i.e. ess(A) == ess(B) */
    if (A->length == B->length)
    {
        success = _fmpz_mpoly_gcd_monomial_cofactors_sp(G, Gbits,
                                                      A, Ainfo, B, Binfo, ctx);
        if (success)
        {
            goto cleanup;
        }
    }

    /* check if ess(A) and ess(B) have a variable v_in_both in common */
    v_in_both = -WORD(1);
    for (j = 0; j < nvars; j++)
    {
        if (   Ainfo[j].max_exp > Ainfo[j].min_exp
            && Binfo[j].max_exp > Binfo[j].min_exp
           )
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
        ulong * minexps;

        fmpz_init(gA);
        fmpz_init(gB);
        _fmpz_vec_content(gA, A->coeffs, A->length);
        _fmpz_vec_content(gB, B->coeffs, B->length);

        minexps = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
        {
            minexps[j] = FLINT_MIN(Ainfo[j].min_exp, Binfo[j].min_exp);
        }

        fmpz_mpoly_fit_length(G, 1, ctx);
        fmpz_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, minexps, Gbits, ctx->minfo);
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

        if (   Ainfo[j].max_exp > Ainfo[j].min_exp
            || Binfo[j].max_exp > Binfo[j].min_exp
           )
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
        ulong * minA, * minB, * minG;
        fmpz_poly_t a, b, g;

        minA = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        minB = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        minG = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));
        for (j = 0; j < nvars; j++)
        {
            minA[j] = Ainfo[j].min_exp;
            minB[j] = Binfo[j].min_exp;
            minG[j] = FLINT_MIN(Ainfo[j].min_exp, Binfo[j].min_exp);
        }
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        _fmpz_mpoly_to_fmpz_poly_shifts(a, A, minA, v_in_both, ctx);
        _fmpz_mpoly_to_fmpz_poly_shifts(b, B, minB, v_in_both, ctx);
        fmpz_poly_gcd(g, a, b);
        _fmpz_mpoly_from_fmpz_poly_shifts(G, Gbits, g, minG, v_in_both, ctx);
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
        if (   Ainfo[j].max_exp > Ainfo[j].min_exp
            && Binfo[j].max_exp == Binfo[j].min_exp
           )
        {
            v_in_A_only = j;
            break;
        }
        if (   Ainfo[j].max_exp == Ainfo[j].min_exp
            && Binfo[j].max_exp > Binfo[j].min_exp
           )
        {
            v_in_B_only = j;
            break;
        }
    }
    if (v_in_A_only != -WORD(1))
    {
        success = _try_missing_var(G, Gbits, v_in_A_only,
                                            A, Ainfo[v_in_A_only].min_exp,
                                            B, Binfo[v_in_A_only].min_exp, ctx);
        goto cleanup;
    }
    if (v_in_B_only != -WORD(1))
    {
        success = _try_missing_var(G, Gbits, v_in_B_only,
                                            B, Binfo[v_in_B_only].min_exp,
                                            A, Ainfo[v_in_B_only].min_exp, ctx);
        goto cleanup;
    }

    /* all variable are now either
            missing from both ess(A) and ess(B), or
            present in both ess(A) and ess(B)
    */

    {
        success = _try_zippel(G, Gbits, A, Ainfo, B, Binfo, ctx);
    }

cleanup:

    TMP_END;

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
        Since we do not have truly sparse GCD algorithms,
        if either of A or B cannot be packed into FLINT_BITS and neither is a
        monomial, assume that the gcd computation is hopeless and fail.
    */
        int success;
        int useAnew = 0;
        int useBnew = 0;
        fmpz_mpoly_t Anew;
        fmpz_mpoly_t Bnew;

        if (A->bits > FLINT_BITS)
        {
            fmpz_mpoly_init(Anew, ctx);
            useAnew = fmpz_mpoly_repack_bits(Anew, A, FLINT_BITS, ctx);
            if (!useAnew)
            {
                useAnew = 1;
                success = 0;
                goto cleanup;
            }
        }

        if (B->bits > FLINT_BITS)
        {
            fmpz_mpoly_init(Bnew, ctx);
            useBnew = fmpz_mpoly_repack_bits(Bnew, B, FLINT_BITS, ctx);
            if (!useBnew)
            {
                useBnew = 1;
                success = 0;
                goto cleanup;
            }
        }

        success = _fmpz_mpoly_gcd(G, FLINT_BITS, useAnew ? Anew : A,
                                                 useBnew ? Bnew : B, ctx);

cleanup:
        if (useAnew)
            fmpz_mpoly_clear(Anew, ctx);
        if (useBnew)
            fmpz_mpoly_clear(Bnew, ctx);

        return success;
    }
}
