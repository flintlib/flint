/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly_factor.h"


static void fq_nmod_mpoly_delete_duplicate_terms(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(fq_nmod_equal(A->coeffs + j, A->coeffs + i, ctx->fqctx));
            continue;
        }
        j++;
        fq_nmod_set(A->coeffs + j, A->coeffs + i, ctx->fqctx);
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}


static slong fq_nmod_mpolyu_find_term(const fq_nmod_mpolyu_t A, ulong e)
{
    slong i;
    for (i = 0; i < A->length; i++)
        if (A->exps[i] == e)
            return i;
    return -1;
}


fq_nmod_mpoly_struct * _fq_nmod_mpolyu_get_coeff(fq_nmod_mpolyu_t A,
                                    ulong pow, const fq_nmod_mpoly_ctx_t uctx);

void fq_nmod_poly_product_roots(
    fq_nmod_poly_t master,
    const fq_nmod_struct * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx);

void n_poly_fq_product_roots_n_fq(
    n_poly_t master,
    const mp_limb_t * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    fq_nmod_poly_t p;
    fq_nmod_struct * m = FLINT_ARRAY_ALLOC(mlength, fq_nmod_struct);

    fq_nmod_poly_init(p, ctx);
    for (i = 0; i < mlength; i++)
    {
        fq_nmod_init(m + i, ctx);
        n_fq_get_fq_nmod(m + i, monomials + d*i, ctx);
    }

    fq_nmod_poly_product_roots(p, m, mlength, ctx);

    n_poly_fq_set_fq_nmod_poly(master, p, ctx);

    fq_nmod_poly_clear(p, ctx);
    for (i = 0; i < mlength; i++)
        fq_nmod_clear(m + i, ctx);
    flint_free(m);
}

void n_poly_fq_product_roots_fq_nmod(
    n_poly_t master,
    fq_nmod_struct * monomials,
    slong mlength,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_poly_t p;
    fq_nmod_poly_init(p, ctx);
    fq_nmod_poly_product_roots(p, monomials, mlength, ctx);
    n_poly_fq_set_fq_nmod_poly(master, p, ctx);
    fq_nmod_poly_clear(p, ctx);
}


void _fq_nmod_mpoly_monomial_evals(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * alpha,
    slong vstart,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t * tmp, * xpoweval;
    ulong * inputexpmask;

    inputexpmask = FLINT_ARRAY_ALLOC(N, ulong);
    LUToffset = FLINT_ARRAY_ALLOC(N*FLINT_BITS, slong);
    LUTmask   = FLINT_ARRAY_ALLOC(N*FLINT_BITS, ulong);
    LUTvalue  = FLINT_ARRAY_ALLOC(d*N*FLINT_BITS, mp_limb_t);
    tmp = FLINT_ARRAY_ALLOC(12*d, mp_limb_t);
    xpoweval = tmp + 6*d;

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

        n_fq_set_fq_nmod(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            _n_fq_set(LUTvalue + d*LUTlen, xpoweval, d);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            _n_fq_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx, tmp);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        _n_fq_one(xpoweval, d);
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                _n_fq_mul(xpoweval, xpoweval, LUTvalue + d*j, ctx->fqctx, tmp);
            }
        }
        _n_fq_set(E + d*i, xpoweval, d);
    }

    flint_free(inputexpmask);
    flint_free(LUToffset);
    flint_free(LUTmask);
    flint_free(LUTvalue);
    flint_free(tmp);
}


static void _fq_nmod_mpoly_monomial_evals_indirect(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong * Aind,
    slong Alen,
    const fq_nmod_struct * alpha,
    slong vstart,
    slong vstop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t * tmp, * xpoweval;
    ulong * inputexpmask;
    const ulong * thisAexp;

    FLINT_ASSERT(0 <= vstart);
    FLINT_ASSERT(vstart < vstop);
    FLINT_ASSERT(vstop <= ctx->minfo->nvars);

    inputexpmask = FLINT_ARRAY_ALLOC(N, ulong);
    LUToffset = FLINT_ARRAY_ALLOC(N*FLINT_BITS, slong);
    LUTmask   = FLINT_ARRAY_ALLOC(N*FLINT_BITS, ulong);
    LUTvalue  = FLINT_ARRAY_ALLOC(d*N*FLINT_BITS, mp_limb_t);
    tmp = FLINT_ARRAY_ALLOC(12*d, mp_limb_t);
    xpoweval = tmp + 6*d;

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

        n_fq_set_fq_nmod(xpoweval, alpha + j, ctx->fqctx); /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            _n_fq_set(LUTvalue + d*LUTlen, xpoweval, d);
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            _n_fq_mul(xpoweval, xpoweval, xpoweval, ctx->fqctx, tmp);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        _n_fq_one(xpoweval, d);
        for (j = 0; j < LUTlen; j++)
        {
            if ((thisAexp[LUToffset[j]] & LUTmask[j]) != 0)
            {
                _n_fq_mul(xpoweval, xpoweval, LUTvalue + d*j, ctx->fqctx, tmp);
            }
        }
        _n_fq_set(E + d*i, xpoweval, d);
    }

    flint_free(inputexpmask);
    flint_free(LUToffset);
    flint_free(LUTmask);
    flint_free(LUTvalue);
    flint_free(tmp);
}

int fq_nmod_zip_find_coeffs_new_fq_nmod(
    fq_nmod_struct * coeffs,             /* length mlength */
    const fq_nmod_struct * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length d*elength */
    slong elength,
    const mp_limb_t * master,       /* length d*(mlength + 1) */
    mp_limb_t * temp,               /* length d*mlength */
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int success;
    slong i, j;
    mp_limb_t * tmp = FLINT_ARRAY_ALLOC(12*d, mp_limb_t);
    mp_limb_t * V = tmp + 6*d;
    mp_limb_t * V0 = V + d;
    mp_limb_t * T = V0 + d;
    mp_limb_t * S = T + d;
    mp_limb_t * r = S + d;
    mp_limb_t * p0 = r + d;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        _n_fq_zero(V0, d);
        _n_fq_zero(T, d);
        _n_fq_zero(S, d);
        n_fq_set_fq_nmod(r, monomials + i, ctx);
        for (j = mlength; j > 0; j--)
        {
            _n_fq_mul(T, r, T, ctx, tmp);
            _n_fq_add(T, T, master + d*j, ctx);

            _n_fq_mul(S, r, S, ctx, tmp);
            _n_fq_add(S, S, T, ctx);

            _n_fq_mul(p0, evals + d*(j - 1), T, ctx, tmp);
            _n_fq_add(V0, V0, p0, ctx);
        }
        /* roots[i] should be a root of master */
#if WANT_ASSERT
        _n_fq_mul(p0, r, T, ctx, tmp);
        _n_fq_add(p0, p0, master + d*0, ctx);
        FLINT_ASSERT(_n_fq_is_zero(p0, d));
#endif
        _n_fq_set(V, V0, d);
        _n_fq_mul(S, S, r, ctx, tmp);
        if (_n_fq_is_zero(S, d))
        {
            success = -1;
            goto cleanup;
        }

        _n_fq_inv(p0, S, ctx, tmp);
        _n_fq_mul(p0, V, p0, ctx, tmp);
        n_fq_get_fq_nmod(coeffs + i, p0, ctx);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
    {
        n_fq_set_fq_nmod(p0, monomials + j, ctx);
        _n_fq_pow_ui(temp + d*j, p0, mlength, ctx);
    }

    for (i = mlength; i < elength; i++)
    {
        _n_fq_zero(V0, d);
        _n_fq_zero(S, d);
        for (j = 0; j < mlength; j++)
        {
            n_fq_set_fq_nmod(p0, monomials + j, ctx);
            _n_fq_mul(temp + d*j, temp + d*j, p0, ctx, tmp);
            n_fq_set_fq_nmod(p0, coeffs + j, ctx);
            _n_fq_mul(p0, p0, temp + d*j, ctx, tmp);
            _n_fq_add(V0, V0, p0, ctx);
        }
        _n_fq_set(V, V0, d);
        if (!_n_fq_equal(V, evals + d*i, d))
        {
            success = 0;
            goto cleanup;
        }
    }

    success = 1;

cleanup:

    flint_free(tmp);

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
    W->terms[W->length].exp = nodes[idx].key;
    W->terms[W->length].coeff[0] = ((n_poly_struct *) T->data)[idx];
    W->length++;

    if (nodes[idx].left >= 0)
        _clearit(W, T, nodes[idx].left);
}


static void fq_nmod_mpoly_set_eval_helper3(
    n_polyun_t EH,
    const fq_nmod_mpoly_t A,
    slong yvar,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    slong yoff, xoff, zoff;
    slong yshift, xshift, zshift;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    const fq_nmod_struct * Acoeffs = A->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    mpoly_rbtree_ui_t W;

    n_polyun_init(T);

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    mpoly_rbtree_ui_init(W);
    for (i = 0; i < Alen; i++)
    {
        n_poly_struct * Wc;
        int its_new;

        y = (Aexps[N*i + yoff] >> yshift) & mask;
        x = (Aexps[N*i + xoff] >> xshift) & mask;
        z = (Aexps[N*i + zoff] >> zshift) & mask;
        Wc = mpoly_rbtree_ui_lookup(W, &its_new, pack_exp3(y, x, z),
                                                        sizeof(n_poly_struct));
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

    T->terms = flint_malloc(W->length*sizeof(n_polyun_term_struct));
    T->alloc = W->length;
    T->length = 0;
    _clearit(T, W, W->nodes[2 - 1].left);
    mpoly_rbtree_ui_clear(W);

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;
    EHterms = EH->terms;

    for (i = 0; i < T->length; i++)
    {
        EHterms[i].exp = T->terms[i].exp;
        n = T->terms[i].coeff->length;
        n_poly_fit_length(EHterms[i].coeff, d*3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _fq_nmod_mpoly_monomial_evals_indirect(p, Aexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        for (j = n - 1; j >= 0; j--)
        {
            _n_fq_set(p + d*(3*j + 2), p + d*j, d);
            _n_fq_set(p + d*(3*j + 0), p + d*(3*j + 2), d);
            n_fq_set_fq_nmod(p + d*(3*j + 1), Acoeffs + ind[j], ctx->fqctx);
        }
    }

    n_polyun_clear(T);
}


/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets comming from different y's (shouldn't happen)
*/
static slong fq_nmod_mpoly_set_eval_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    fq_nmod_mpolyu_t H,
    const fq_nmod_mpoly_t B,
    const fq_nmod_struct * alpha,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    n_polyun_term_struct * EHterms;
    mp_limb_t * p;
    fq_nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const fq_nmod_struct * Bcoeffs = B->coeffs;
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

        mpoly_rbtree_ui_init(W);
        for (i = 0; i < Blen; i++)
        {
            y = (Bexps[N*i + yoff] >> yshift) & mask;
            x = (Bexps[N*i + xoff] >> xshift) & mask;
            z = (Bexps[N*i + zoff] >> zshift) & mask;

            FLINT_ASSERT(x <= deg);

            Wc = mpoly_rbtree_ui_lookup(W, &its_new, pack_exp3(y, x, z),
                                                        sizeof(n_poly_struct));
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

        T->terms = flint_malloc(W->length*sizeof(n_polyun_term_struct));
        T->alloc = W->length;
        T->length = 0;
        _clearit(T, W, W->nodes[2 - 1].left);
        mpoly_rbtree_ui_clear(W);
    }

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;
    EHterms = EH->terms;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EHterms[i].exp = T->terms[i].exp;
        y = extract_exp(EHterms[i].exp, 2, 3);
        x = extract_exp(EHterms[i].exp, 1, 3);
        z = extract_exp(EHterms[i].exp, 0, 3);
        n = T->terms[i].coeff->length;
        n_poly_fit_length(EHterms[i].coeff, d*3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _fq_nmod_mpoly_monomial_evals_indirect(p, Bexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fq_nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            fq_nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            for (j = 0; j < n; j++)
            {
                n_fq_get_fq_nmod(Hc->coeffs + old_len + j, p + d*j, ctx->fqctx);
                mpoly_monomial_set(Hc->exps + N*(old_len + j),
                                   Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                fq_nmod_mpoly_sort_terms(Hc, ctx);
                fq_nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            _n_fq_set(p + d*(3*j + 2), p + d*j, d);
            _n_fq_set(p + d*(3*j + 0), p + d*(3*j + 2), d);
            n_fq_set_fq_nmod(p + d*(3*j + 1), Bcoeffs + ind[j], ctx->fqctx);
        }
    }

    n_polyun_clear(T);

    *deg_ = deg;

    return zip_length;
}


static void fq_nmod_poly_eval_step(
    mp_limb_t * res,
    n_poly_t A,
    const fq_nmod_ctx_t ctx,
    mp_limb_t * tmp)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    mp_limb_t * sum = tmp + 4*d;

    FLINT_ASSERT(d*3*Alen <= A->alloc);

    if (Alen < 1)
    {
        _n_fq_zero(res, d);
        return;
    }

    i = 0;
    _n_fq_mul2(sum, Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 1), ctx);
    _n_fq_mul(Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 2), ctx, tmp);
    for (i = 1; i < Alen; i++)
    {
        _n_fq_madd2(sum, Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 1), ctx, tmp);
        _n_fq_mul(Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 0), Acoeffs + d*(3*i + 2), ctx, tmp);
    }
    _n_fq_reduce2(res, sum, ctx, tmp);
}


void fq_nmod_polyu_eval_step(
    n_polyu_t E,
    n_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Ai, Ei;
    n_polyun_term_struct * Aterms = A->terms;
    mp_limb_t * tmp = FLINT_ARRAY_ALLOC(8*d, mp_limb_t);

    n_polyu_fit_length(E, d*A->length);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = Aterms[Ai].exp;
        fq_nmod_poly_eval_step(E->coeffs + d*Ei, Aterms[Ai].coeff, ctx, tmp);
        Ei += !_n_fq_is_zero(E->coeffs + d*Ei, d);
    }
    E->length = Ei;

    flint_free(tmp);
}


void fq_nmod_polyu3_add_zip_limit1(
    n_polyun_t Z,
    const n_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    const n_polyun_term_struct * At = A->terms;
    const n_polyun_term_struct * Ait;
    n_polyun_term_struct * Zit;
    slong Ai, ai, Zi, j;

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->terms[Zi].coeff->length == cur_length);
    }


    Ai = -1;
    ai = -1;
    do {
        Ai++;
        Ait = At + Ai;
    } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
    if (Ai < A->length)
        ai = n_poly_degree(Ait->coeff);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        Ait = At + Ai;
        if (Ait->exp + ai > Zit->exp)
        {
            /* missing from Z */
            n_polyun_fit_length(Z, Z->length + 1);
            for (j = Z->length; j > Zi; j--)
                n_polyun_term_swap(Z->terms + j, Z->terms + j - 1);
            Z->length++;
            Zit = Z->terms + Zi;
            Zit->exp = Ait->exp + ai;
            n_poly_fit_length(Zit->coeff, d*fit_length);
            Zit->coeff->length = cur_length;
            flint_mpn_zero(Zit->coeff->coeffs, d*cur_length);
            goto in_both;            
        }
        else if (Ait->exp + ai < Zit->exp)
        {
            /* missing from A */
            FLINT_ASSERT(d*(cur_length + 1) <= Zit->coeff->alloc);
            _n_fq_zero(Zit->coeff->coeffs + d*cur_length, d);
            Zit->coeff->length = cur_length + 1;
            Zi++;
        }
        else
        {
in_both:
            FLINT_ASSERT(cur_length == Zit->coeff->length);
            FLINT_ASSERT(d*(cur_length + 1) <= Zit->coeff->alloc);
            _n_fq_set(Zit->coeff->coeffs + d*cur_length, Ait->coeff->coeffs + d*ai, d);
            Zit->coeff->length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Ait->coeff->coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai++;
                    Ait = At + Ai;
                } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = n_poly_degree(Ait->coeff);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;
        n_polyun_fit_length(Z, Zi + A->length - Ai);
        Zit = Z->terms + Zi;
        Zit->exp = Ait->exp + ai;
        n_poly_fit_length(Zit->coeff, d*fit_length);
        Zit->coeff->length = cur_length;
        flint_mpn_zero(Zit->coeff->coeffs, d*cur_length);
        Z->length = ++Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(d*(cur_length + 1) <= Zit->coeff->alloc);
        _n_fq_set(Zit->coeff->coeffs + d*cur_length, Ait->coeff->coeffs + d*ai, d);
        Zit->coeff->length = cur_length + 1;
        do {
            ai--;
        } while (ai >= 0 && _n_fq_is_zero(Ait->coeff->coeffs + d*ai, d));
        if (ai < 0)
        {
            do {
                Ai++;
                Ait = At + Ai;
            } while (Ai < A->length && extract_exp(Ait->exp, 1, 3) >= deg1);
            if (Ai < A->length)
                ai = n_poly_degree(Ait->coeff);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        Zit = Z->terms + Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(d*(cur_length + 1) <= Zit->coeff->alloc);
        _n_fq_zero(Zit->coeff->coeffs + d*cur_length, d);
        Zit->coeff->length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->terms[Zi].coeff->length == cur_length + 1);
    }
}

/*
    for each Y^y*X^x*Z^z in B with x = deg,
        keep the Y^y*X^x*Z^z*poly(x1,...) in B
    for each Y^y*X^x*Z^z in Z,
        assert that x < deg
        if there is no Y^0*X^x*Z^y in H, fail
        find coefficients of poly using this entry in H
        output Y^y*X^x*Z^z*poly(x1,...) to A
    sort A

    return
        -1: singular vandermonde matrix encountered
        0:  inconsistent system encountered
        1:  success
*/
static int fq_nmod_mpoly_from_zip(
    fq_nmod_mpoly_t B,
    const n_polyun_t Z,
    fq_nmod_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const fq_nmod_mpoly_ctx_t ctx,
    n_polyun_t M)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    fq_nmod_struct * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    n_polyun_term_struct * Zt = Z->terms;
    fq_nmod_mpoly_struct * Hc;
    slong Hlen = H->length;

    FLINT_ASSERT(bits == H->bits);

    n_polyun_fit_length(M, Hlen + 1);
    for (i = 0; i <= Hlen; i++)
        M->terms[i].coeff->length = 0;

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
        y = extract_exp(Zt[Zi].exp, 2, 3);
        x = extract_exp(Zt[Zi].exp, 1, 3);
        z = extract_exp(Zt[Zi].exp, 0, 3);
        FLINT_ASSERT(x < deg);
        Hi = fq_nmod_mpolyu_find_term(H, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        fq_nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->terms[Hi].coeff->length < 1)
        {
            n_poly_fq_product_roots_fq_nmod(M->terms[Hi].coeff,
                                           Hc->coeffs, Hc->length, ctx->fqctx);
        }

        n_poly_fit_length(M->terms[Hlen].coeff, d*Hc->length);

        success = fq_nmod_zip_find_coeffs_new_fq_nmod(Bcoeffs + Bi, Hc->coeffs,
                    Hc->length, Zt[Zi].coeff->coeffs, Zt[Zi].coeff->length,
                    M->terms[Hi].coeff->coeffs, M->terms[Hlen].coeff->coeffs,
                                                             ctx->fqctx);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (fq_nmod_is_zero(Bcoeffs + j, ctx->fqctx))
                continue;
            fq_nmod_set(Bcoeffs + Bi, Bcoeffs + j, ctx->fqctx);
            FLINT_ASSERT(Bi < B->alloc);
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    fq_nmod_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(fq_nmod_mpoly_is_canonical(B, ctx));

    return 1;
}


/* bit counts of all degrees should be < FLINT_BITS/3 */
int fq_nmod_mpoly_hlift_zippel(
    slong m,
    fq_nmod_mpoly_struct * B,
    slong r,
    const fq_nmod_struct * alpha,
    const fq_nmod_mpoly_t A,
    const slong * degs,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    int success;
    slong i;
    slong zip_fails_remaining;
    slong req_zip_images, cur_zip_image;
    fq_nmod_mpolyu_struct * H;
    n_polyun_struct M[1], Aeh[1], * Beh, * BBeval, * Z;
    n_polyu_struct Aeval[1], * Beval;
    fq_nmod_struct * beta;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

#if WANT_ASSERT
    {
        fq_nmod_mpoly_t T;
        slong j, * check_degs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);

        fq_nmod_mpoly_init(T, ctx);

        fq_nmod_mpoly_degrees_si(check_degs, A, ctx);
        for (j = 0; j < ctx->minfo->nvars; j++)
            FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);

        fq_nmod_mpoly_one(T, ctx);
        for (i = 0; i < r; i++)
        {
            fq_nmod_mpoly_degrees_si(check_degs, B + i, ctx);
            for (j = 0; j < ctx->minfo->nvars; j++)
                FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);
            fq_nmod_mpoly_mul(T, T, B + i, ctx);
        }
        fq_nmod_mpoly_sub(T, A, T, ctx);

        fq_nmod_mpoly_evaluate_one_fq_nmod(T, T, m, alpha + m - 1, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_is_zero(T, ctx));

        fq_nmod_mpoly_clear(T, ctx);
        flint_free(check_degs);
    }
#endif


    beta = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, fq_nmod_struct);
    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_nmod_init(beta + i, ctx->fqctx);

    Bdegs = FLINT_ARRAY_ALLOC(r, ulong);
    H = FLINT_ARRAY_ALLOC(r, fq_nmod_mpolyu_struct);
    Beh = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    Beval = FLINT_ARRAY_ALLOC(r, n_polyu_struct);
    BBeval = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    Z = FLINT_ARRAY_ALLOC(r, n_polyun_struct);

    n_polyun_init(Aeh);
    n_polyu_init(Aeval);
    n_polyun_init(M);
    for (i = 0; i < r; i++)
    {
        fq_nmod_mpolyu_init(H + i, bits, ctx);
        n_polyun_init(Beh + i);
        n_polyu_init(Beval + i);
        n_polyun_init(BBeval + i);
        n_polyun_init(Z + i);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = fq_nmod_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fq_nmod_rand(beta + i, state, ctx->fqctx);
        if (fq_nmod_is_zero(beta + i, ctx->fqctx))
            fq_nmod_one(beta + i, ctx->fqctx);
    }

    fq_nmod_mpoly_set_eval_helper3(Aeh, A, m, beta, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_zip_images;
        this_zip_images = fq_nmod_mpoly_set_eval_helper_and_zip_form3(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_zip_images);
        FLINT_ASSERT(Bdegs[i] > 0);

        Z[i].length = 0;
    }

    cur_zip_image = 0;

next_zip_image:

    fq_nmod_polyu_eval_step(Aeval, Aeh, ctx->fqctx);

    for (i = 0; i < r; i++)
        fq_nmod_polyu_eval_step(Beval + i, Beh + i, ctx->fqctx);

    success = n_polyu3_fq_hlift(r, BBeval, Aeval, Beval,
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
        fq_nmod_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                    cur_zip_image, req_zip_images, ctx->fqctx);
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = fq_nmod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx, M);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    fq_nmod_mpoly_init3(T1, A->length, bits, ctx);
    fq_nmod_mpoly_init3(T2, A->length, bits, ctx);
    fq_nmod_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fq_nmod_mpoly_mul(T2, T1, B + i, ctx);
        fq_nmod_mpoly_swap(T1, T2, ctx);
    }

    success = fq_nmod_mpoly_equal(T1, A, ctx);
    fq_nmod_mpoly_clear(T1, ctx);
    fq_nmod_mpoly_clear(T2, ctx);

cleanup:

    n_polyun_clear(Aeh);
    n_polyu_clear(Aeval);
    n_polyun_clear(M);
    for (i = 0; i < r; i++)
    {
        fq_nmod_mpolyu_clear(H + i, ctx);
        n_polyun_clear(Beh + i);
        n_polyu_clear(Beval + i);
        n_polyun_clear(BBeval + i);
        n_polyun_clear(Z + i);
    }

    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_nmod_clear(beta + i, ctx->fqctx);
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
int fq_nmod_mpoly_factor_irred_smprime_zippel(
    fq_nmod_mpolyv_t fac,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_factor_t lcAfac,
    const fq_nmod_mpoly_t lcA,
    const fq_nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    int alphas_tries_remaining, alphabetas_tries_remaining, alphabetas_length;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k, r;
    fq_nmod_struct * alpha;
    n_poly_struct * alphabetas;
    fq_nmod_mpoly_struct * Aevals;
    slong * degs, * degeval;
    fq_nmod_mpolyv_t tfac;
    fq_nmod_mpoly_t t, Acopy;
    fq_nmod_mpoly_struct * newA;
    n_poly_t Abfc;
    n_bpoly_t Ab;
    n_tpoly_t Abfp;
    fq_nmod_mpoly_t m, mpow;
    fq_nmod_mpolyv_t new_lcs, lc_divs;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(A->length > 1);
    FLINT_ASSERT(fq_nmod_is_one(A->coeffs + 0, ctx->fqctx));
    FLINT_ASSERT(A->bits <= FLINT_BITS);

    if (ctx->fqctx->modulus->length < n_clog(A->length, ctx->fqctx->modulus->mod.n))
        return 0;

    fq_nmod_mpoly_init(Acopy, ctx);
    fq_nmod_mpoly_init(m, ctx);
    fq_nmod_mpoly_init(mpow, ctx);

    fq_nmod_mpolyv_init(new_lcs, ctx);
    fq_nmod_mpolyv_init(lc_divs, ctx);

    n_poly_init(Abfc);
    n_tpoly_init(Abfp);
    n_bpoly_init(Ab);

    degs    = FLINT_ARRAY_ALLOC(n + 1, slong);
    degeval = FLINT_ARRAY_ALLOC(n + 1, slong);
	alpha   = FLINT_ARRAY_ALLOC(n, fq_nmod_struct);
    alphabetas = FLINT_ARRAY_ALLOC(n, n_poly_struct);
    Aevals  = FLINT_ARRAY_ALLOC(n, fq_nmod_mpoly_struct);
	for (i = 0; i < n; i++)
    {
        fq_nmod_init(alpha + i, ctx->fqctx);
        n_poly_init(alphabetas + i);
		fq_nmod_mpoly_init(Aevals + i, ctx);
    }
    fq_nmod_mpolyv_init(tfac, ctx);
	fq_nmod_mpoly_init(t, ctx);

    /* init done */

    alphabetas_length = 2;
    alphas_tries_remaining = 10;
	fq_nmod_mpoly_degrees_si(degs, A, ctx);

next_alpha:

    if (--alphas_tries_remaining < 0)
	{
		success = 0;
        goto cleanup;
	}

    for (i = 0; i < n; i++)
    {
        fq_nmod_rand(alpha + i, state, ctx->fqctx);
        if (fq_nmod_is_zero(alpha + i, ctx->fqctx))
            fq_nmod_one(alpha + i, ctx->fqctx);
    }

    /* ensure degrees do not drop under evaluation */
	for (i = n - 1; i >= 0; i--)
	{
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
		fq_nmod_mpoly_degrees_si(degeval, Aevals + i, ctx);
		for (j = 0; j <= i; j++)
			if (degeval[j] != degs[j])
				goto next_alpha;
	}

    /* make sure univar is squarefree */
	fq_nmod_mpoly_derivative(t, Aevals + 0, 0, ctx);
	fq_nmod_mpoly_gcd(t, t, Aevals + 0, ctx);
	if (!fq_nmod_mpoly_is_one(t, ctx))
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
        n_poly_fit_length(alphabetas + i, d*alphabetas_length);
        n_fq_set_fq_nmod(alphabetas[i].coeffs + d*0, alpha + i, ctx->fqctx);
        for (j = d; j < d*alphabetas_length; j++)
            alphabetas[i].coeffs[j] = n_urandint(state, ctx->fqctx->mod.n);
        alphabetas[i].length = alphabetas_length;
        _n_poly_fq_normalise(alphabetas + i, d);
    }

    _fq_nmod_mpoly_eval_rest_to_n_bpoly_fq(Ab, A, alphabetas, ctx);
    success = n_bpoly_fq_factor_smprime(Abfc, Abfp, Ab, 0, ctx->fqctx);
    if (!success)
    {
        FLINT_ASSERT(0 && "this should not happen");
        goto next_alpha;
    }

    r = Abfp->length;

    if (r < 2)
    {
        fq_nmod_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fq_nmod_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = fq_nmod_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                       Abfc, Abfp->coeffs, r, alphabetas, ctx);
        if (!success)
            goto next_alphabetas;
    }
    else
    {
        for (i = 0; i < r; i++)
            fq_nmod_mpoly_one(lc_divs->coeffs + i, ctx);
    }

    success = fq_nmod_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx);
    FLINT_ASSERT(success);
    for (i = 1; i < r; i++)
    {
        success = fq_nmod_mpoly_divides(m, m, lc_divs->coeffs + i, ctx);
        FLINT_ASSERT(success);
    }

    fq_nmod_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fq_nmod_mpoly_is_one(mpow, ctx))
    {
        newA = (fq_nmod_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fq_nmod_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_degrees_si(degs, newA, ctx);

    for (i = 0; i < n + 1; i++)
    {
        if (FLINT_BIT_COUNT(degs[i]) >= FLINT_BITS/3)
        {
            success = -1;
            goto cleanup;
        }
    }

    fq_nmod_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(t, mpow, i + 1, alpha + i, ctx);
        fq_nmod_mpoly_swap(t, mpow, ctx);
        fq_nmod_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fq_nmod_mpolyv_fit_length(new_lcs, (n + 1)*r, ctx);
    i = n;
    for (j = 0; j < r; j++)
    {
        fq_nmod_mpoly_mul(new_lcs->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fq_nmod_mpoly_evaluate_one_fq_nmod(new_lcs->coeffs + i*r + j,
                       new_lcs->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fq_nmod_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        fq_nmod_t q;
        fq_nmod_init(q, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fq_nmod_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fq_nmod_mpoly_set_n_bpoly_fq_var1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        fq_nmod_inv(q, fac->coeffs[i].coeffs + 0, ctx->fqctx);
        fq_nmod_mul(q, q, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->fqctx);
        fq_nmod_mpoly_scalar_mul_fq_nmod(fac->coeffs + i, fac->coeffs + i, q, ctx);
        fq_nmod_clear(q, ctx->fqctx);
    }

    fq_nmod_mpolyv_fit_length(tfac, r, ctx);
    tfac->length = r;
    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _fq_nmod_mpoly_set_lead0(tfac->coeffs + i, fac->coeffs + i,
                                               new_lcs->coeffs + k*r + i, ctx);
        }

        if (k > 2)
        {
            success = fq_nmod_mpoly_hlift_zippel(k, tfac->coeffs, r, alpha,
                                  k < n ? Aevals + k : newA, degs, ctx, state);
        }
        else
        {
            success = fq_nmod_mpoly_hlift(k, tfac->coeffs, r, alpha,
                                         k < n ? Aevals + k : newA, degs, ctx);
        }

        if (!success)
            goto next_alphabetas;

        fq_nmod_mpolyv_swap(tfac, fac, ctx);
    }

    if (!fq_nmod_mpoly_is_fq_nmod(m, ctx))
    {
        fq_nmod_mpoly_univar_t u;
        fq_nmod_mpoly_univar_init(u, ctx);
        for (i = 0; i < r; i++)
        {
            fq_nmod_mpoly_to_univar(u, fac->coeffs + i, 0, ctx);
            success = _fq_nmod_mpoly_vec_content_mpoly(t, u->coeffs, u->length, ctx);
            if (!success)
            {
                fq_nmod_mpoly_univar_clear(u, ctx);
                success = -1;
                goto cleanup;
            }
            success = fq_nmod_mpoly_divides(fac->coeffs + i,
                                            fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
        fq_nmod_mpoly_univar_clear(u, ctx);
    }

    for (i = 0; i < r; i++)
        fq_nmod_mpoly_make_monic(fac->coeffs + i, fac->coeffs + i, ctx);

    success = 1;

cleanup:

    fq_nmod_mpolyv_clear(new_lcs, ctx);
    fq_nmod_mpolyv_clear(lc_divs, ctx);

    n_poly_clear(Abfc);
    n_tpoly_clear(Abfp);
    n_bpoly_clear(Ab);

	for (i = 0; i < n; i++)
    {
		fq_nmod_mpoly_clear(Aevals + i, ctx);
        n_poly_clear(alphabetas + i);
        fq_nmod_clear(alpha + i, ctx->fqctx);
    }
    flint_free(alphabetas);
    flint_free(alpha);
    flint_free(Aevals);
    flint_free(degs);
    flint_free(degeval);

    fq_nmod_mpolyv_clear(tfac, ctx);
    fq_nmod_mpoly_clear(t, ctx);

    fq_nmod_mpoly_clear(Acopy, ctx);
    fq_nmod_mpoly_clear(m, ctx);
    fq_nmod_mpoly_clear(mpow, ctx);

#if WANT_ASSERT
    if (success)
    {
        fq_nmod_mpoly_t prod;
        fq_nmod_mpoly_init(prod, ctx);
        fq_nmod_mpoly_one(prod, ctx);
        for (i = 0; i < fac->length; i++)
            fq_nmod_mpoly_mul(prod, prod, fac->coeffs + i, ctx);
        FLINT_ASSERT(fq_nmod_mpoly_equal(prod, A, ctx));
        fq_nmod_mpoly_clear(prod, ctx);
    }
#endif

	return success;
}
