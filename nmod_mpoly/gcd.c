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

    success = nmod_mpoly_gcd_brown(G, A, B, ctx);

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
        return 0;
    }
}
