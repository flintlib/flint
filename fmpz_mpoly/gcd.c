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
/*
flint_printf("_fmpz_mpoly_gcd called Gbits = %wd\n", Gbits);
printf("A: "); fmpz_mpoly_print_pretty(A, NULL, ctx); printf("\n");
printf("B: "); fmpz_mpoly_print_pretty(B, NULL, ctx); printf("\n");
*/
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

    success = fmpz_mpoly_gcd_brown(G, A, B, ctx);

cleanup:

    TMP_END;
/*
flint_printf("_fmpz_mpoly_gcd returning %d\n", success);
printf("G: "); fmpz_mpoly_print_pretty(G, NULL, ctx); printf("\n");
*/
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
