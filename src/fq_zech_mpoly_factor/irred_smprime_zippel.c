/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_poly_factor.h"
#include "fq_zech_mpoly_factor.h"

static void fq_zech_mpoly_delete_duplicate_terms(
    fq_zech_mpoly_t A,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(fq_zech_equal(A->coeffs + j, A->coeffs + i, ctx->fqctx));
            continue;
        }
        j++;
        fq_zech_set(A->coeffs + j, A->coeffs + i, ctx->fqctx);
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}


slong fq_zech_mpolyu_find_term(const fq_zech_mpolyu_t A, ulong e)
{
    slong i;
    for (i = 0; i < A->length; i++)
        if (A->exps[i] == e)
            return i;
    return -1;
}


void fq_zech_poly_product_roots(fq_zech_poly_t P, fq_zech_struct * r,
                                            slong n, const fq_zech_ctx_t fqctx)
{
    slong i;
    fq_zech_poly_t B;
    fq_zech_t a;
    fq_zech_init(a, fqctx);
    fq_zech_poly_init(B, fqctx);
    fq_zech_poly_one(P, fqctx);
    fq_zech_poly_gen(B, fqctx);
    for (i = 0; i < n; i++)
    {
        fq_zech_neg(a, r + i, fqctx);
        fq_zech_poly_set_coeff(B, 0, a, fqctx);
        fq_zech_poly_mul(P, P, B, fqctx);
    }
    fq_zech_clear(a, fqctx);
    fq_zech_poly_clear(B, fqctx);
}


fq_zech_mpoly_struct * _fq_zech_mpolyu_get_coeff(fq_zech_mpolyu_t A,
                                     ulong pow, const fq_zech_mpoly_ctx_t uctx);

void _fq_zech_mpoly_monomial_evals(
    fq_zech_struct * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_zech_struct * alpha,
    slong vstart,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    fq_zech_struct * LUTvalue;
    slong LUTlen;
    fq_zech_t xpoweval;
    ulong * inputexpmask;

    inputexpmask = FLINT_ARRAY_ALLOC(N, ulong);
    LUToffset = FLINT_ARRAY_ALLOC(N*FLINT_BITS, slong);
    LUTmask   = FLINT_ARRAY_ALLOC(N*FLINT_BITS, ulong);
    LUTvalue  = FLINT_ARRAY_ALLOC(N*FLINT_BITS, fq_zech_struct);
    for (i = 0; i < N*FLINT_BITS; i++)
        fq_zech_init(LUTvalue + i, ctx->fqctx);
    fq_zech_init(xpoweval, ctx->fqctx);

    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < N; j++)
        {
            inputexpmask[j] |= (Aexps + N*i)[j];
        }
    }

    LUTlen = 0;
    for (j = nvars - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, ctx->minfo);

        fq_zech_set(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fq_zech_set(LUTvalue + LUTlen, xpoweval, ctx->fqctx);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            fq_zech_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        fq_zech_one(xpoweval, ctx->fqctx);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fq_zech_mul(xpoweval, xpoweval, LUTvalue + j, ctx->fqctx);
            }
        }
        fq_zech_set(E + i, xpoweval, ctx->fqctx);
    }

    flint_free(inputexpmask);
    flint_free(LUToffset);
    flint_free(LUTmask);
    flint_free(LUTvalue);
}


static void _fq_zech_mpoly_monomial_evals_indirect(
    fq_zech_struct * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong * Aind,
    slong Alen,
    const fq_zech_struct * alpha,
    slong vstart,
    slong vstop,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * LUToffset;
    ulong * LUTmask;
    fq_zech_struct * LUTvalue;
    slong LUTlen;
    fq_zech_t xpoweval;
    ulong * inputexpmask;
    const ulong * thisAexp;

    FLINT_ASSERT(0 <= vstart);
    FLINT_ASSERT(vstart < vstop);
    FLINT_ASSERT(vstop <= ctx->minfo->nvars);

    inputexpmask = FLINT_ARRAY_ALLOC(N, ulong);
    LUToffset = FLINT_ARRAY_ALLOC(N*FLINT_BITS, slong);
    LUTmask   = FLINT_ARRAY_ALLOC(N*FLINT_BITS, ulong);
    LUTvalue  = FLINT_ARRAY_ALLOC(N*FLINT_BITS, fq_zech_struct);
    for (i = 0; i < N*FLINT_BITS; i++)
        fq_zech_init(LUTvalue + i, ctx->fqctx);
    fq_zech_init(xpoweval, ctx->fqctx);

    mpoly_monomial_zero(inputexpmask, N);
    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        for (j = 0; j < N; j++)
            inputexpmask[j] |= thisAexp[j];
    }

    LUTlen = 0;
    for (j = vstop - 1; j >= vstart; j--)
    {
        mpoly_gen_offset_shift_sp(&offset, &shift, j, Abits, ctx->minfo);

        fq_zech_set(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            fq_zech_set(LUTvalue + LUTlen, xpoweval, ctx->fqctx);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            fq_zech_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        fq_zech_one(xpoweval, ctx->fqctx);
        for (j = 0; j < LUTlen; j++)
        {
            if ((thisAexp[LUToffset[j]] & LUTmask[j]) != 0)
            {
                fq_zech_mul(xpoweval, xpoweval, LUTvalue + j, ctx->fqctx);
            }
        }
        fq_zech_set(E + i, xpoweval, ctx->fqctx);
    }

    flint_free(inputexpmask);
    flint_free(LUToffset);
    flint_free(LUTmask);
    flint_free(LUTvalue);
}


int fq_zech_zip_find_coeffs_new(
    fq_zech_struct * coeffs,             /* length mlength */
    const fq_zech_struct * monomials,    /* length mlength */
    slong mlength,
    const fq_zech_struct * evals,        /* length elength */
    slong elength,
    const fq_zech_struct * master,       /* length (mlength + 1) */
    fq_zech_struct * temp,               /* length mlength */
    const fq_zech_ctx_t ctx)
{
    int success;
    slong i, j;
    fq_zech_t V, V0, T, S, r, p0;

    fq_zech_init(V, ctx);
    fq_zech_init(V0, ctx);
    fq_zech_init(T, ctx);
    fq_zech_init(S, ctx);
    fq_zech_init(r, ctx);
    fq_zech_init(p0, ctx);

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        fq_zech_zero(V0, ctx);
        fq_zech_zero(T, ctx);
        fq_zech_zero(S, ctx);
        fq_zech_set(r, monomials + i, ctx);
        for (j = mlength; j > 0; j--)
        {
            fq_zech_mul(T, r, T, ctx);
            fq_zech_add(T, T, master + j, ctx);

            fq_zech_mul(S, r, S, ctx);
            fq_zech_add(S, S, T, ctx);

            fq_zech_mul(p0, evals + (j - 1), T, ctx);
            fq_zech_add(V0, V0, p0, ctx);
        }
        /* roots[i] should be a root of master */
#ifdef FLINT_WANT_ASSERT
        fq_zech_mul(p0, r, T, ctx);
        fq_zech_add(p0, p0, master + 0, ctx);
        FLINT_ASSERT(fq_zech_is_zero(p0, ctx));
#endif
        fq_zech_set(V, V0, ctx);
        fq_zech_mul(S, S, r, ctx);
        if (fq_zech_is_zero(S, ctx))
        {
            success = -1;
            goto cleanup;
        }

        fq_zech_inv(p0, S, ctx);
        fq_zech_mul(coeffs + i, V, p0, ctx);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        fq_zech_pow_ui(temp + j, monomials + j, mlength, ctx);

    for (i = mlength; i < elength; i++)
    {
        fq_zech_zero(V0, ctx);
        fq_zech_zero(S, ctx);
        for (j = 0; j < mlength; j++)
        {
            fq_zech_mul(temp + j, temp + j, monomials + j, ctx);
            fq_zech_mul(p0, coeffs + j, temp + j, ctx);
            fq_zech_add(V0, V0, p0, ctx);
        }
        fq_zech_set(V, V0, ctx);
        if (!fq_zech_equal(V, evals + i, ctx))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    fq_zech_clear(V, ctx);
    fq_zech_clear(V0, ctx);
    fq_zech_clear(T, ctx);
    fq_zech_clear(S, ctx);
    fq_zech_clear(r, ctx);
    fq_zech_clear(p0, ctx);

    return success;
}



static void _clearit(
    n_polyun_t W,
    mpoly_rbtree_ui_t T,
    slong idx)
{
    mpoly_rbnode_ui_struct * nodes = T->nodes + 2;

    FLINT_ASSERT(0 <= idx && idx < T->length);

    if (nodes[idx].right >= 0)
        _clearit(W, T, nodes[idx].right);

    FLINT_ASSERT(W->length < W->alloc);
    W->exps[W->length] = nodes[idx].key;
    W->coeffs[W->length] = ((n_poly_struct *) T->data)[idx];
    W->length++;

    if (nodes[idx].left >= 0)
        _clearit(W, T, nodes[idx].left);
}


static void fq_zech_mpoly_set_eval_helper3(
    fq_zech_polyun_t EH,
    const fq_zech_mpoly_t A,
    slong yvar,
    const fq_zech_struct * alpha,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    slong yoff, xoff, zoff;
    slong yshift, xshift, zshift;
    fq_zech_struct * p;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    const fq_zech_struct * Acoeffs = A->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    mpoly_rbtree_ui_t W;

    n_polyun_init(T);

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    mpoly_rbtree_ui_init(W, sizeof(n_poly_struct));
    for (i = 0; i < Alen; i++)
    {
        n_poly_struct * Wc;
        int its_new;

        y = (Aexps[N*i + yoff] >> yshift) & mask;
        x = (Aexps[N*i + xoff] >> xshift) & mask;
        z = (Aexps[N*i + zoff] >> zshift) & mask;
        Wc = mpoly_rbtree_ui_lookup(W, &its_new, pack_exp3(y, x, z));
        if (its_new)
        {
            n_poly_init2(Wc, 4);
            Wc->coeffs[0] = i;
            Wc->length = 1;
        }
        else
        {
            n_poly_fit_length(Wc, Wc->length + 1);
            Wc->coeffs[Wc->length] = i;
            Wc->length++;
        }
    }

    FLINT_ASSERT(W->length > 0);

    T->coeffs = FLINT_ARRAY_ALLOC(W->length, n_poly_struct);
    T->exps   = FLINT_ARRAY_ALLOC(W->length, ulong);
    T->alloc = W->length;
    T->length = 0;
    _clearit(T, W, W->nodes[2 - 1].left);
    mpoly_rbtree_ui_clear(W);

    fq_zech_polyun_fit_length(EH, T->length, ctx->fqctx);
    EH->length = T->length;

    for (i = 0; i < T->length; i++)
    {
        EH->exps[i] = T->exps[i];
        n = T->coeffs[i].length;
        fq_zech_poly_fit_length(EH->coeffs + i, 3*n, ctx->fqctx);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;
        _fq_zech_mpoly_monomial_evals_indirect(p, Aexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        for (j = n - 1; j >= 0; j--)
        {
            fq_zech_set(p + 3*j + 2, p + j, ctx->fqctx);
            fq_zech_set(p + 3*j + 0, p + 3*j + 2, ctx->fqctx);
            fq_zech_set(p + 3*j + 1, Acoeffs + ind[j], ctx->fqctx);
        }
    }

    n_polyun_clear(T);
}


/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets coming from different y's (shouldn't happen)
*/
static slong fq_zech_mpoly_set_eval_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    fq_zech_polyun_t EH,
    fq_zech_mpolyu_t H,
    const fq_zech_mpoly_t B,
    const fq_zech_struct * alpha,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const fq_zech_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    fq_zech_struct * p;
    fq_zech_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const fq_zech_struct * Bcoeffs = B->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    ulong deg;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == H->bits);
    FLINT_ASSERT(Blen > 0);

    /* init T */
    {
        mpoly_rbtree_ui_t W;
        n_poly_struct * Wc;
        slong yoff, xoff, zoff;
        slong yshift, xshift, zshift;
        int its_new;

        mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
        mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
        mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

        deg = (Bexps[N*0 + xoff] >> xshift) & mask;

        mpoly_rbtree_ui_init(W, sizeof(n_poly_struct));
        for (i = 0; i < Blen; i++)
        {
            y = (Bexps[N*i + yoff] >> yshift) & mask;
            x = (Bexps[N*i + xoff] >> xshift) & mask;
            z = (Bexps[N*i + zoff] >> zshift) & mask;

            FLINT_ASSERT(x <= deg);

            Wc = mpoly_rbtree_ui_lookup(W, &its_new, pack_exp3(y, x, z));
            if (its_new)
            {
                n_poly_init2(Wc, 4);
                Wc->coeffs[0] = i;
                Wc->length = 1;
            }
            else
            {
                n_poly_fit_length(Wc, Wc->length + 1);
                Wc->coeffs[Wc->length] = i;
                Wc->length++;
            }
        }

        FLINT_ASSERT(W->length > 0);

        T->coeffs = FLINT_ARRAY_ALLOC(W->length, n_poly_struct);
        T->exps   = FLINT_ARRAY_ALLOC(W->length, ulong);
        T->alloc = W->length;
        T->length = 0;
        _clearit(T, W, W->nodes[2 - 1].left);
        mpoly_rbtree_ui_clear(W);
    }

    fq_zech_polyun_fit_length(EH, T->length, ctx->fqctx);
    EH->length = T->length;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EH->exps[i] = T->exps[i];
        y = extract_exp(EH->exps[i], 2, 3);
        x = extract_exp(EH->exps[i], 1, 3);
        z = extract_exp(EH->exps[i], 0, 3);
        n = T->coeffs[i].length;
        fq_zech_poly_fit_length(EH->coeffs + i, 3*n, ctx->fqctx);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;
        _fq_zech_mpoly_monomial_evals_indirect(p, Bexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fq_zech_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            fq_zech_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            for (j = 0; j < n; j++)
            {
                fq_zech_set(Hc->coeffs + old_len + j, p + j, ctx->fqctx);
                mpoly_monomial_set(Hc->exps + N*(old_len + j),
                                   Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                fq_zech_mpoly_sort_terms(Hc, ctx);
                fq_zech_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            fq_zech_set(p + 3*j + 2, p + j, ctx->fqctx);
            fq_zech_set(p + 3*j + 0, p + 3*j + 2, ctx->fqctx);
            fq_zech_set(p + 3*j + 1, Bcoeffs + ind[j], ctx->fqctx);
        }
    }

    n_polyun_clear(T);

    *deg_ = deg;

    return zip_length;
}


static void fq_zech_poly_eval_step(
    fq_zech_t res,
    fq_zech_poly_t A,
    const fq_zech_ctx_t ctx)
{
    slong i, Alen = A->length;
    fq_zech_struct * Acoeffs = A->coeffs;
    fq_zech_t t;

    FLINT_ASSERT(3*Alen <= A->alloc);

    if (Alen < 1)
    {
        fq_zech_zero(res, ctx);
        return;
    }

    fq_zech_init(t, ctx);

    i = 0;
    fq_zech_mul(res, Acoeffs + (3*i + 0), Acoeffs + (3*i + 1), ctx);
    fq_zech_mul(Acoeffs + (3*i + 0), Acoeffs + (3*i + 0), Acoeffs + (3*i + 2), ctx);
    for (i = 1; i < Alen; i++)
    {
        fq_zech_mul(t, Acoeffs + (3*i + 0), Acoeffs + (3*i + 1), ctx);
        fq_zech_add(res, res, t, ctx);
        fq_zech_mul(Acoeffs + (3*i + 0), Acoeffs + (3*i + 0), Acoeffs + (3*i + 2), ctx);
    }

    fq_zech_clear(t, ctx);
}


void fq_zech_polyu_eval_step(
    fq_zech_polyu_t E,
    fq_zech_polyun_t A,
    const fq_zech_ctx_t ctx)
{
    slong Ai, Ei;

    fq_zech_polyu_fit_length(E, A->length, ctx);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = A->exps[Ai];
        fq_zech_poly_eval_step(E->coeffs + Ei, A->coeffs + Ai, ctx);
        Ei += !fq_zech_is_zero(E->coeffs + Ei, ctx);
    }
    E->length = Ei;
}


void fq_zech_polyu3_add_zip_limit1(
    fq_zech_polyun_t Z,
    const fq_zech_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length,
    const fq_zech_ctx_t ctx)
{
    const fq_zech_poly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    fq_zech_poly_struct * Zcoeffs = Z->coeffs;
    ulong * Zexps = Z->exps;
    slong Ai, ai, Zi, j;

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->coeffs[Zi].length == cur_length);
    }

    Ai = -1;
    ai = -1;
    do {
        Ai++;
    } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
    if (Ai < A->length)
        ai = fq_zech_poly_degree(Acoeffs + Ai, ctx);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        if (Aexps[Ai] + ai > Zexps[Zi])
        {
            /* missing from Z */
            fq_zech_polyun_fit_length(Z, Z->length + 1, ctx);
            Zcoeffs = Z->coeffs;
            Zexps = Z->exps;

            for (j = Z->length; j > Zi; j--)
            {
                fq_zech_poly_swap(Zcoeffs + j, Zcoeffs + j - 1, ctx);
                FLINT_SWAP(ulong, Zexps[j], Zexps[j - 1]);
            }

            Z->length++;

            Zexps[Zi] = Aexps[Ai] + ai;
            fq_zech_poly_fit_length(Zcoeffs + Zi, fit_length, ctx);
            Zcoeffs[Zi].length = cur_length;

            for (j = 0; j < cur_length; j++)
                fq_zech_zero(Zcoeffs[Zi].coeffs + j, ctx);

            goto in_both;
        }
        else if (Aexps[Ai] + ai < Zexps[Zi])
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
            fq_zech_zero(Zcoeffs[Zi].coeffs + cur_length, ctx);
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
        }
        else
        {
    in_both:
            FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
            FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
            fq_zech_set(Zcoeffs[Zi].coeffs + cur_length, Acoeffs[Ai].coeffs + ai, ctx);
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && fq_zech_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
            if (ai < 0)
            {
                do {
                    Ai++;
                } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = fq_zech_poly_degree(Acoeffs + Ai, ctx);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;

        fq_zech_polyun_fit_length(Z, Zi + A->length - Ai, ctx);
        Zcoeffs = Z->coeffs;
        Zexps = Z->exps;

        Zexps[Zi] = Aexps[Ai] + ai;
        fq_zech_poly_fit_length(Zcoeffs + Zi, fit_length, ctx);
        Zcoeffs[Zi].length = cur_length;

        for (j = 0; j < cur_length; j++)
            fq_zech_zero(Zcoeffs[Zi].coeffs + j, ctx);

        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
        fq_zech_set(Zcoeffs[Zi].coeffs + cur_length, Acoeffs[Ai].coeffs + ai, ctx);
        Zcoeffs[Zi].length = cur_length + 1;

        Z->length = ++Zi;

        do {
            ai--;
        } while (ai >= 0 && fq_zech_is_zero(Acoeffs[Ai].coeffs + ai, ctx));
        if (ai < 0)
        {
            do {
                Ai++;
            } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
            if (Ai < A->length)
                ai = fq_zech_poly_degree(Acoeffs + Ai, ctx);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
        fq_zech_zero(Zcoeffs[Zi].coeffs + cur_length, ctx);
        Zcoeffs[Zi].length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->coeffs[Zi].length == cur_length + 1);
    }
}


static int fq_zech_mpoly_from_zip(
    fq_zech_mpoly_t B,
    const fq_zech_polyun_t Z,
    fq_zech_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const fq_zech_mpoly_ctx_t ctx,
    fq_zech_polyun_t M)
{
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    fq_zech_struct * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    fq_zech_mpoly_struct * Hc;
    slong Hlen = H->length;

    FLINT_ASSERT(bits == H->bits);

    fq_zech_polyun_fit_length(M, Hlen + 1, ctx->fqctx);
    for (i = 0; i <= Hlen; i++)
        M->coeffs[i].length = 0;

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    /* x is most significant in ctx, so keeping the lc_x in B is easy */
    FLINT_ASSERT(xvar == 0);

    for (Bi = 0; Bi < B->length; Bi++)
    {
        x = (((B->exps + N*Bi)[xoff] >> xshift) & mask);
        FLINT_ASSERT(x <= deg);
        if (x != deg)
            break;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        y = extract_exp(Z->exps[Zi], 2, 3);
        x = extract_exp(Z->exps[Zi], 1, 3);
        z = extract_exp(Z->exps[Zi], 0, 3);
        FLINT_ASSERT(x < deg);
        Hi = fq_zech_mpolyu_find_term(H, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        fq_zech_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->coeffs[Hi].length < 1)
        {
            fq_zech_poly_product_roots(M->coeffs + Hi,
                                           Hc->coeffs, Hc->length, ctx->fqctx);
        }

        fq_zech_poly_fit_length(M->coeffs + Hlen, Hc->length, ctx->fqctx);

        success = fq_zech_zip_find_coeffs_new(Bcoeffs + Bi, Hc->coeffs, Hc->length,
                                Z->coeffs[Zi].coeffs, Z->coeffs[Zi].length,
                                M->coeffs[Hi].coeffs, M->coeffs[Hlen].coeffs,
                                                                   ctx->fqctx);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (fq_zech_is_zero(Bcoeffs + j, ctx->fqctx))
                continue;

            FLINT_ASSERT(Bi < B->alloc);

            fq_zech_set(Bcoeffs + Bi, Bcoeffs + j, ctx->fqctx);
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    fq_zech_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(fq_zech_mpoly_is_canonical(B, ctx));

    return 1;
}

/* bit counts of all degrees should be < FLINT_BITS/3 */
int fq_zech_mpoly_hlift_zippel(
    slong m,
    fq_zech_mpoly_struct * B,
    slong r,
    const fq_zech_struct * alpha,
    const fq_zech_mpoly_t A,
    const slong * degs,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    slong i;
    slong zip_fails_remaining;
    slong req_zip_images, cur_zip_image;
    fq_zech_mpolyu_struct * H;
    fq_zech_polyun_struct M[1], Aeh[1], * Beh, * BBeval, * Z;
    fq_zech_polyu_struct Aeval[1], * Beval;
    fq_zech_struct * beta;
    flint_bitcnt_t bits = A->bits;
    fq_zech_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

#ifdef FLINT_WANT_ASSERT
    {
        fq_zech_mpoly_t T;
        slong j, * check_degs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);

        fq_zech_mpoly_init(T, ctx);

        fq_zech_mpoly_degrees_si(check_degs, A, ctx);
        for (j = 0; j < ctx->minfo->nvars; j++)
            FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);

        fq_zech_mpoly_one(T, ctx);
        for (i = 0; i < r; i++)
        {
            fq_zech_mpoly_degrees_si(check_degs, B + i, ctx);
            for (j = 0; j < ctx->minfo->nvars; j++)
                FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);
            fq_zech_mpoly_mul(T, T, B + i, ctx);
        }
        fq_zech_mpoly_sub(T, A, T, ctx);

        fq_zech_mpoly_evaluate_one_fq_zech(T, T, m, alpha + m - 1, ctx);
        FLINT_ASSERT(fq_zech_mpoly_is_zero(T, ctx));

        fq_zech_mpoly_clear(T, ctx);
        flint_free(check_degs);
    }
#endif

    beta = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, fq_zech_struct);
    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_zech_init(beta + i, ctx->fqctx);

    Bdegs = FLINT_ARRAY_ALLOC(r, ulong);
    H = FLINT_ARRAY_ALLOC(r, fq_zech_mpolyu_struct);
    Beh = FLINT_ARRAY_ALLOC(r, fq_zech_polyun_struct);
    Beval = FLINT_ARRAY_ALLOC(r, fq_zech_polyu_struct);
    BBeval = FLINT_ARRAY_ALLOC(r, fq_zech_polyun_struct);
    Z = FLINT_ARRAY_ALLOC(r, fq_zech_polyun_struct);

    fq_zech_polyun_init(Aeh, ctx->fqctx);
    fq_zech_polyu_init(Aeval, ctx->fqctx);
    fq_zech_polyun_init(M, ctx->fqctx);
    for (i = 0; i < r; i++)
    {
        fq_zech_mpolyu_init(H + i, bits, ctx);
        fq_zech_polyun_init(Beh + i, ctx->fqctx);
        fq_zech_polyu_init(Beval + i, ctx->fqctx);
        fq_zech_polyun_init(BBeval + i, ctx->fqctx);
        fq_zech_polyun_init(Z + i, ctx->fqctx);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = fq_zech_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fq_zech_rand(beta + i, state, ctx->fqctx);
        if (fq_zech_is_zero(beta + i, ctx->fqctx))
            fq_zech_one(beta + i, ctx->fqctx);
    }

    fq_zech_mpoly_set_eval_helper3(Aeh, A, m, beta, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_zip_images;
        this_zip_images = fq_zech_mpoly_set_eval_helper_and_zip_form3(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_zip_images);
        FLINT_ASSERT(Bdegs[i] > 0);

        Z[i].length = 0;
    }

    cur_zip_image = 0;

next_zip_image:

    fq_zech_polyu_eval_step(Aeval, Aeh, ctx->fqctx);

    for (i = 0; i < r; i++)
        fq_zech_polyu_eval_step(Beval + i, Beh + i, ctx->fqctx);

    success = fq_zech_polyu3_hlift(r, BBeval, Aeval, Beval,
                                             alpha + m - 1, degs0, ctx->fqctx);
    if (success < 1)
    {
        if (--zip_fails_remaining >= 0)
            goto choose_betas;

        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        fq_zech_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                    cur_zip_image, req_zip_images, ctx->fqctx);
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = fq_zech_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx, M);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    fq_zech_mpoly_init3(T1, A->length, bits, ctx);
    fq_zech_mpoly_init3(T2, A->length, bits, ctx);
    fq_zech_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fq_zech_mpoly_mul(T2, T1, B + i, ctx);
        fq_zech_mpoly_swap(T1, T2, ctx);
    }

    success = fq_zech_mpoly_equal(T1, A, ctx);
    fq_zech_mpoly_clear(T1, ctx);
    fq_zech_mpoly_clear(T2, ctx);

cleanup:

    fq_zech_polyun_clear(Aeh, ctx->fqctx);
    fq_zech_polyu_clear(Aeval, ctx->fqctx);
    fq_zech_polyun_clear(M, ctx->fqctx);
    for (i = 0; i < r; i++)
    {
        fq_zech_mpolyu_clear(H + i, ctx);
        fq_zech_polyun_clear(Beh + i, ctx->fqctx);
        fq_zech_polyu_clear(Beval + i, ctx->fqctx);
        fq_zech_polyun_clear(BBeval + i, ctx->fqctx);
        fq_zech_polyun_clear(Z + i, ctx->fqctx);
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_zech_clear(beta + i, ctx->fqctx);
    flint_free(beta);

    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

    return success;
}


/*
    return 1: success
           0: failed
          -1: exception
*/
int fq_zech_mpoly_factor_irred_smprime_zippel(
    fq_zech_mpolyv_t fac,
    const fq_zech_mpoly_t A,
    const fq_zech_mpoly_factor_t lcAfac,
    const fq_zech_mpoly_t lcA,
    const fq_zech_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fq_zech_struct * alpha;
    fq_zech_poly_struct * alphabetas;
    fq_zech_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fq_zech_mpolyv_t tfac;
    fq_zech_mpoly_t t, Acopy;
    fq_zech_mpoly_struct * newA;
    fq_zech_poly_t Abfc;
    fq_zech_bpoly_t Ab;
    fq_zech_tpoly_t Abfp;
    fq_zech_mpoly_t m, mpow;
    fq_zech_mpolyv_t new_lcs, lc_divs;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(fq_zech_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    if (fq_zech_ctx_degree(ctx->fqctx) <= n_clog(A->length, fq_zech_ctx_mod(ctx->fqctx).n))
        return 0;

    fq_zech_mpoly_init(Acopy, ctx);
    fq_zech_mpoly_init(m, ctx);
    fq_zech_mpoly_init(mpow, ctx);

    fq_zech_mpolyv_init(new_lcs, ctx);
    fq_zech_mpolyv_init(lc_divs, ctx);

    fq_zech_poly_init(Abfc, ctx->fqctx);
    fq_zech_tpoly_init(Abfp, ctx->fqctx);
    fq_zech_bpoly_init(Ab, ctx->fqctx);

    degs    = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
	alpha   = FLINT_ARRAY_ALLOC(n, fq_zech_struct);
    alphabetas = FLINT_ARRAY_ALLOC(n, fq_zech_poly_struct);
    Aevals  = FLINT_ARRAY_ALLOC(n, fq_zech_mpoly_struct);
	for (i = 0; i < n; i++)
    {
        fq_zech_init(alpha + i, ctx->fqctx);
        fq_zech_poly_init(alphabetas + i, ctx->fqctx);
		fq_zech_mpoly_init(Aevals + i, ctx);
    }
    fq_zech_mpolyv_init(tfac, ctx);
	fq_zech_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
	fq_zech_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
    {
        fq_zech_rand(alpha + i, state, ctx->fqctx);
        if (fq_zech_is_zero(alpha + i, ctx->fqctx))
            fq_zech_one(alpha + i, ctx->fqctx);
    }

	for (i = n - 1; i >= 0; i--)
	{
        fq_zech_mpoly_evaluate_one_fq_zech(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
		fq_zech_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != degs[j])
				goto next_alpha;
	}

    fq_zech_mpoly_get_fq_zech_poly(Abfc, Aevals + 0, 0, ctx);
    if (!fq_zech_poly_is_squarefree(Abfc, ctx->fqctx))
        goto next_alpha;

    alphabetas_tries_remaining = 2 + alphabetas_length;

next_alphabetas:

    if (--alphabetas_tries_remaining < 0)
    {
        if (++alphabetas_length > 5)
        {
            success = 0;
            goto cleanup;
        }
        goto next_alpha;
    }

    for (i = 0; i < n; i++)
    {
        fq_zech_poly_fit_length(alphabetas + i, alphabetas_length, ctx->fqctx);
        fq_zech_set(alphabetas[i].coeffs + 0, alpha + i, ctx->fqctx);
        for (j = 1; j < alphabetas_length; j++)
            fq_zech_rand(alphabetas[i].coeffs + j, state, ctx->fqctx);
        alphabetas[i].length = alphabetas_length;
        _fq_zech_poly_normalise(alphabetas + i, ctx->fqctx);
    }

    _fq_zech_mpoly_eval_to_bpoly(Ab, A, alphabetas, ctx);
    success = fq_zech_bpoly_factor_smprime(Abfc, Abfp, Ab, 0, ctx->fqctx);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        fq_zech_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fq_zech_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    fq_zech_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = fq_zech_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            fq_zech_mpoly_one(lc_divs->coeffs + i, ctx);
    }

    success = fq_zech_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fq_zech_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fq_zech_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fq_zech_mpoly_is_one(mpow, ctx))
    {
        newA = (fq_zech_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fq_zech_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fq_zech_mpoly_degrees_si(degs, newA, ctx);

    for (i = 0; i < n + 1; i++)
    {
        if (FLINT_BIT_COUNT(degs[i]) >= FLINT_BITS/3)
        {
            success = -1;
            goto cleanup;
        }
    }

    fq_zech_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fq_zech_mpoly_evaluate_one_fq_zech(t, mpow, i + 1, alpha + i, ctx);
        fq_zech_mpoly_swap(t, mpow, ctx);
        fq_zech_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fq_zech_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fq_zech_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_zech_mpoly_evaluate_one_fq_zech(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fq_zech_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        fq_zech_t q;
        fq_zech_init(q, ctx->fqctx);
        FLINT_ASSERT(fq_zech_mpoly_is_fq_zech(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fq_zech_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fq_zech_mpoly_set_fq_zech_bpoly_var1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        fq_zech_inv(q, fac->coeffs[i].coeffs + 0, ctx->fqctx);
        fq_zech_mul(q, q, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->fqctx);
        fq_zech_mpoly_scalar_mul_fq_zech(fac->coeffs + i, fac->coeffs + i, q, ctx);
        fq_zech_clear(q, ctx->fqctx);
    }

    fq_zech_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fq_zech_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        if (k > 2)
        {
            success = fq_zech_mpoly_hlift_zippel(k, tfac->coeffs, r, alpha,
                                  k < n ? Aevals + k : newA, degs, ctx, state);
        }
        else
        {
            success = fq_zech_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);
        }

        if (!success)
            goto next_alphabetas;

        fq_zech_mpolyv_swap(tfac, fac, ctx);
    }

    if (!fq_zech_mpoly_is_fq_zech(m, ctx))
    {
        fq_zech_mpoly_univar_t u;
        fq_zech_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            fq_zech_mpoly_to_univar(u, fac->coeffs + i, 0, ctx);
            success = fq_zech_mpoly_univar_content_mpoly(t, u, ctx);
            if (!success)
            {
                fq_zech_mpoly_univar_clear(u, ctx);
                success = -1;
                goto cleanup;
            }
            success = fq_zech_mpoly_divides(fac->coeffs + i,
                                            fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
        fq_zech_mpoly_univar_clear(u, ctx);
    }

    for (i = 0; i < r; i++)
        fq_zech_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);

    success = 1;

cleanup:

    fq_zech_mpolyv_clear(new_lcs, ctx);
    fq_zech_mpolyv_clear(lc_divs, ctx);

    fq_zech_bpoly_clear(Ab, ctx->fqctx);
    fq_zech_poly_clear(Abfc, ctx->fqctx);
    fq_zech_tpoly_clear(Abfp, ctx->fqctx);

	for (i = 0; i < n; i++)
    {
		fq_zech_mpoly_clear(Aevals + i, ctx);
        fq_zech_poly_clear(alphabetas + i, ctx->fqctx);
        fq_zech_clear(alpha + i, ctx->fqctx);
    }
    flint_free(alphabetas);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(degs);
    flint_free(degeval);

    fq_zech_mpolyv_clear(tfac, ctx);
    fq_zech_mpoly_clear(t, ctx);

    fq_zech_mpoly_clear(Acopy, ctx);
    fq_zech_mpoly_clear(m, ctx);
    fq_zech_mpoly_clear(mpow, ctx);

#ifdef FLINT_WANT_ASSERT
    if (success)
    {
        fq_zech_mpoly_t prod;
        fq_zech_mpoly_init(prod, ctx);
        fq_zech_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fq_zech_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fq_zech_mpoly_equal(prod, A, ctx));
        fq_zech_mpoly_clear(prod, ctx);
    }
#endif

	return success;
}
