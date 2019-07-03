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
static int _try_missing_var(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    slong var,
    const nmod_mpoly_t A,
    ulong Ashift,
    const nmod_mpoly_t B,
    ulong Bshift,
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
    success = _nmod_mpoly_gcd(tG, Gbits, B, Ax->coeffs + 0, ctx, NULL, 0);
    if (!success)
        goto cleanup;

    for (i = 1; i < Ax->length; i++)
    {
        success = _nmod_mpoly_gcd(tG, Gbits, tG, Ax->coeffs + i, ctx, NULL, 0);
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
static int _try_zippel(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    ulong * Gstride,
    const nmod_mpoly_t A,
    const ulong * Amax_exp,
    const ulong * Amin_exp,
    const slong * Amax_exp_count,
    const slong * Amin_exp_count,
    const nmod_mpoly_t B,
    const ulong * Bmax_exp,
    const ulong * Bmin_exp,
    const slong * Bmax_exp_count,
    const slong * Bmin_exp_count,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j, k;
    slong n, m;
    int success;
    ulong * Gshift, * Addeg, * Bddeg;
    mpoly_zipinfo_t zinfo;
    flint_bitcnt_t ABbits;
    flint_rand_t randstate;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyu_t Au, Bu, Gu;
    nmod_mpoly_t Acontent, Bcontent;
    nmod_mpolyu_t Abar, Bbar, Gbar;

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
        The Zippel algorithm will operate in Zp[y_0][y_1,...,y_m] and it
        only operates with Zp[y_1,...,y_m] in LEX.

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

    /* interpolation will continue in m  variables */
    mpoly_zipinfo_init(zinfo, m);

    /* uctx is context for Zp[y_1,...,y_{m-1}]*/
    nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->ffinfo->mod.n);

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
        count = FLINT_MIN(Amin_exp_count[j], Amax_exp_count[j]);
        count = FLINT_MIN(count, Bmin_exp_count[j]);
        count = FLINT_MIN(count, Bmax_exp_count[j]);
        deg = FLINT_MAX(Addeg[j], Bddeg[j]);
        for (i = 0; i < m; i++)
        {
            j = zinfo->perm[i];
            new_count = FLINT_MIN(Amin_exp_count[j], Amax_exp_count[j]);
            new_count = FLINT_MIN(new_count, Bmin_exp_count[j]);
            new_count = FLINT_MIN(new_count, Bmax_exp_count[j]);
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

    nmod_mpolyu_init(Au, ABbits, uctx);
    nmod_mpolyu_init(Bu, ABbits, uctx);
    nmod_mpolyu_init(Gu, ABbits, uctx);

    nmod_mpoly_to_mpolyu_perm_deflate(Au, uctx, A, ctx,
                                      zinfo->perm, Amin_exp, Gstride, NULL, 0);
    nmod_mpoly_to_mpolyu_perm_deflate(Bu, uctx, B, ctx,
                                      zinfo->perm, Bmin_exp, Gstride, NULL, 0);

    FLINT_ASSERT(Au->bits == ABbits);
    FLINT_ASSERT(Bu->bits == ABbits);
    FLINT_ASSERT(Au->length > 1);
    FLINT_ASSERT(Bu->length > 1);

    nmod_mpoly_init3(Acontent, 0, ABbits, uctx);
    nmod_mpoly_init3(Bcontent, 0, ABbits, uctx);
    nmod_mpolyu_init(Abar, ABbits, uctx);
    nmod_mpolyu_init(Bbar, ABbits, uctx);
    nmod_mpolyu_init(Gbar, ABbits, uctx);

    /* remove content from A and B */
    success = nmod_mpolyu_content_mpoly(Acontent, Au, uctx, NULL, 0);
    success = success
           && nmod_mpolyu_content_mpoly(Bcontent, Bu, uctx, NULL, 0);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);
    FLINT_ASSERT(Bcontent->bits == ABbits);
    nmod_mpolyu_divexact_mpoly(Abar, Au, Acontent, uctx);
    nmod_mpolyu_divexact_mpoly(Bbar, Bu, Bcontent, uctx);

    /* after removing content, degree bounds in zinfo are still valid bounds */

    /* compute GCD */
    success = nmod_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, uctx, zinfo, randstate);
    if (!success)
        goto cleanup;

    /* put back content */
    success = _nmod_mpoly_gcd(Acontent, ABbits, Acontent, Bcontent, uctx,
                                                                      NULL, 0);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Acontent->bits == ABbits);
    nmod_mpolyu_mul_mpoly(Gu, Gbar, Acontent, uctx);
    nmod_mpoly_from_mpolyu_perm_inflate(G, Gbits, ctx, Gu, uctx,
                                                 zinfo->perm, Gshift, Gstride);
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

    flint_free(Gshift);
    flint_free(Addeg);
    flint_free(Bddeg);

    return success;
}



typedef struct
{
    nmod_mpolyun_struct * Pn;
    const nmod_mpoly_ctx_struct * uctx;
    const nmod_mpoly_struct * P;
    const nmod_mpoly_ctx_struct * ctx;
    const slong * perm;
    const ulong * shift, * stride;
    const thread_pool_handle * handles;
    slong num_handles;
}
_convertn_arg_struct;

typedef _convertn_arg_struct _convertn_arg_t[1];

static void _worker_convertn(void * varg)
{
    _convertn_arg_struct * arg = (_convertn_arg_struct *) varg;

    nmod_mpoly_to_mpolyun_perm_deflate(arg->Pn, arg->uctx, arg->P, arg->ctx,
           arg->perm, arg->shift, arg->stride, arg->handles, arg->num_handles);
}

/*
    return 1 for success or 0 for failure
*/
static int _try_brown(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    ulong * Gstride,
    const nmod_mpoly_t A,
    const ulong * Amax_exp,
    const ulong * Amin_exp,
    const nmod_mpoly_t B,
    const ulong * Bmax_exp,
    const ulong * Bmin_exp,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
{
    int success;
    slong j;
    slong denseAsize, denseBsize;
    slong m, n = ctx->minfo->nvars;
    slong * perm;
    ulong * Gshift;
    flint_bitcnt_t ABbits;
    nmod_mpoly_ctx_t uctx;
    nmod_mpolyun_t An, Bn, Gn, Abarn, Bbarn;
    nmod_poly_stack_t Sp;

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

    nmod_mpoly_ctx_init(uctx, m - 1, ORD_LEX, ctx->ffinfo->mod.n);
    nmod_poly_stack_init(Sp, ABbits, uctx);
    nmod_mpolyun_init(An, ABbits, uctx);
    nmod_mpolyun_init(Bn, ABbits, uctx);
    nmod_mpolyun_init(Gn, ABbits, uctx);
    nmod_mpolyun_init(Abarn, ABbits, uctx);
    nmod_mpolyun_init(Bbarn, ABbits, uctx);

    if (num_handles > 0)
    {
        slong m = mpoly_divide_threads(num_handles, A->length, B->length);
        _convertn_arg_t arg;

        FLINT_ASSERT(m >= 0);
        FLINT_ASSERT(m < num_handles);

        arg->Pn = Bn;
        arg->uctx = uctx;
        arg->P = B;
        arg->ctx = ctx;
        arg->perm = perm;
        arg->shift = Bmin_exp;
        arg->stride = Gstride;
        arg->handles = handles + (m + 1);
        arg->num_handles = num_handles - (m + 1);

        thread_pool_wake(global_thread_pool, handles[m], _worker_convertn, arg);

        nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                      perm, Amin_exp, Gstride, handles + 0, m);

        thread_pool_wait(global_thread_pool, handles[m]);
    }
    else
    {
        nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                             perm, Amin_exp, Gstride, NULL, 0);
        nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx,
                                             perm, Bmin_exp, Gstride, NULL, 0);
    }

    FLINT_ASSERT(An->bits == ABbits);
    FLINT_ASSERT(Bn->bits == ABbits);
    FLINT_ASSERT(An->length > 1);
    FLINT_ASSERT(Bn->length > 1);

    success = (num_handles > 0)
        ? nmod_mpolyun_gcd_brown_smprime_threaded(Gn, Abarn, Bbarn, An, Bn,
                                            m - 2, uctx, handles, num_handles)
        : nmod_mpolyun_gcd_brown_smprime(Gn, Abarn, Bbarn, An, Bn,
                                                              m - 2, uctx, Sp);

    if (!success)
    {
        nmod_mpoly_to_mpolyun_perm_deflate(An, uctx, A, ctx,
                                             perm, Amin_exp, Gstride, NULL, 0);
        nmod_mpoly_to_mpolyun_perm_deflate(Bn, uctx, B, ctx,
                                             perm, Bmin_exp, Gstride, NULL, 0);
        success = nmod_mpolyun_gcd_brown_lgprime(Gn, Abarn, Bbarn, An, Bn,
                                                                  m - 2, uctx);
    }

    if (success)
    {
        nmod_mpoly_from_mpolyun_perm_inflate(G, Gbits, ctx, Gn, uctx,
                                                        perm, Gshift, Gstride);
        nmod_mpoly_make_monic(G, G, ctx);
        FLINT_ASSERT(G->bits == Gbits);
    }

    nmod_mpolyun_clear(An, uctx);
    nmod_mpolyun_clear(Bn, uctx);
    nmod_mpolyun_clear(Gn, uctx);
    nmod_mpolyun_clear(Abarn, uctx);
    nmod_mpolyun_clear(Bbarn, uctx);
    nmod_poly_stack_clear(Sp);
    nmod_mpoly_ctx_clear(uctx);

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
int _nmod_mpoly_gcd(
    nmod_mpoly_t G,
    flint_bitcnt_t Gbits,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    const thread_pool_handle * handles,
    slong num_handles)
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

    if (A->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, B, A, ctx);
    }
    else if (B->length == 1)
    {
        return _nmod_mpoly_gcd_monomial(G, Gbits, A, B, ctx);
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
        for (j = 0; j < nvars; j++)
            Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        nmod_mpoly_fit_length(G, 1, ctx);
        nmod_mpoly_fit_bits(G, Gbits, ctx);
        G->bits = Gbits;
        mpoly_set_monomial_ui(G->exps, Gshift, Gbits, ctx->minfo);
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
        nmod_poly_t a, b, g;

        for (j = 0; j < nvars; j++)
            Gshift[j] = FLINT_MIN(Amin_exp[j], Bmin_exp[j]);

        mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

        nmod_poly_init_mod(a, ctx->ffinfo->mod);
        nmod_poly_init_mod(b, ctx->ffinfo->mod);
        nmod_poly_init_mod(g, ctx->ffinfo->mod);
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
    mpoly_gcd_info_stride(Gstride,
                  A->exps, A->bits, A->length, Amax_exp, Amin_exp,
                  B->exps, B->bits, B->length, Bmax_exp, Bmin_exp, ctx->minfo);

    success = _try_brown(G, Gbits, Gstride, A, Amax_exp, Amin_exp,
                                            B, Bmax_exp, Bmin_exp, ctx,
                                                         handles, num_handles);
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


int nmod_mpoly_gcd_threaded(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx,
    slong thread_limit)
{
    slong i;
    flint_bitcnt_t Gbits;
    int success;
    thread_pool_handle * handles;
    slong num_handles;

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

        /* get workers */
        handles = NULL;
        num_handles = 0;
        if (global_thread_pool_initialized)
        {
            slong max_num_handles = thread_pool_get_size(global_thread_pool);
            max_num_handles = FLINT_MIN(thread_limit - 1, max_num_handles);
            if (max_num_handles > 0)
            {
                handles = (thread_pool_handle *) flint_malloc(
                                   max_num_handles*sizeof(thread_pool_handle));
                num_handles = thread_pool_request(global_thread_pool,
                                                     handles, max_num_handles);
            }
        }

        success = _nmod_mpoly_gcd(G, Gbits, A, B, ctx, handles, num_handles);

        for (i = 0; i < num_handles; i++)
        {
            thread_pool_give_back(global_thread_pool, handles[i]);
        }
        if (handles)
        {
            flint_free(handles);
        }

        return success;
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
                                             useBnew ? Bnew : B, ctx, NULL, 0);
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

        success = _nmod_mpoly_gcd(G, FLINT_BITS, Anew, Bnew, ctx, NULL, 0);

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

int nmod_mpoly_gcd(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    return nmod_mpoly_gcd_threaded(G, A, B, ctx, MPOLY_DEFAULT_THREAD_LIMIT);
}
