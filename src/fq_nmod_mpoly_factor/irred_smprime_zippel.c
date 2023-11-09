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


static void _sort_and_delete_duplicates(
    fq_nmod_mpoly_t A,
    slong d,
    const mpoly_ctx_t mctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, mctx);

    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && mpoly_monomial_lt_nomask(A->exps + N*(j - 1),
                                                  A->exps + N*j, N); j--)
    {
        mpoly_monomial_swap(A->exps + N*(j - 1), A->exps + N*j, N);
        _n_fq_swap(A->coeffs + d*(j - 1), A->coeffs + d*(j), d);
    }

    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(_n_fq_equal(A->coeffs + d*j, A->coeffs + d*i, d));
            continue;
        }
        j++;
        _n_fq_set(A->coeffs + d*j, A->coeffs + d*i, d);
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
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


static void fq_nmod_mpoly_set_eval_helper3(
    n_polyun_t EH,
    const fq_nmod_mpoly_t A,
    slong yvar,
    n_poly_struct * caches,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    ulong y, x, z;
    slong yoff, xoff, zoff, * off;
    slong yshift, xshift, zshift, * shift;
    mp_limb_t * p;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    const mp_limb_t * Acoeffs = A->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    mpoly_rbtree_ui_t W;
    TMP_INIT;

    TMP_START;

    n_polyun_init(T);

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    off = (slong *) TMP_ALLOC(2*yvar*sizeof(slong));
    shift = off + yvar;
    for (i = 2; i < yvar; i++)
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, bits, ctx->minfo);

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

    T->exps = FLINT_ARRAY_ALLOC(W->length, ulong);
    T->coeffs = FLINT_ARRAY_ALLOC(W->length, n_poly_struct);
    T->alloc = W->length;
    T->length = 0;
    _clearit(T, W, W->nodes[2 - 1].left);
    mpoly_rbtree_ui_clear(W);

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;

    for (i = 0; i < T->length; i++)
    {
        EH->exps[i] = T->exps[i];
        n = T->coeffs[i].length;
        n_poly_fit_length(EH->coeffs + i, d*3*n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Ai = ind[j];

            /* set cur = monomial eval */
            _n_fq_one(p + d*j, d);
            for (k = 2; k < yvar; k++)
            {
                ulong ei = (Aexps[N*Ai + off[k]] >> shift[k]) & mask;
                n_fq_pow_cache_mulpow_ui(p + d*j, p + d*j, ei, caches + 3*k + 0,
                               caches + 3*k + 1, caches + 3*k + 2, ctx->fqctx);
            }

            /* copy cur to inc */
            _n_fq_set(p + d*(1*n + j), p + d*(0*n + j), d);

            /* copy coeff */
            _n_fq_set(p + d*(2*n + j), Acoeffs + d*Ai, d);
        }
    }

    n_polyun_clear(T);

    TMP_END;
}

static void fq_nmod_mpoly_set_evalp_helper3(
    n_polyun_t EH,
    const fq_nmod_mpoly_t A,
    slong yvar,
    n_poly_struct * caches,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    ulong y, x, z;
    slong yoff, xoff, zoff, * off;
    slong yshift, xshift, zshift, * shift;
    mp_limb_t * p;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    const mp_limb_t * Acoeffs = A->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    mpoly_rbtree_ui_t W;
    TMP_INIT;

    TMP_START;

    n_polyun_init(T);

    mpoly_gen_offset_shift_sp(&yoff, &yshift, yvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&xoff, &xshift, xvar, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&zoff, &zshift, zvar, bits, ctx->minfo);

    off = (slong *) TMP_ALLOC(2*yvar*sizeof(slong));
    shift = off + yvar;
    for (i = 2; i < yvar; i++)
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, bits, ctx->minfo);

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

    T->exps = FLINT_ARRAY_ALLOC(W->length, ulong);
    T->coeffs = FLINT_ARRAY_ALLOC(W->length, n_poly_struct);
    T->alloc = W->length;
    T->length = 0;
    _clearit(T, W, W->nodes[2 - 1].left);
    mpoly_rbtree_ui_clear(W);

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;

    for (i = 0; i < T->length; i++)
    {
        EH->exps[i] = T->exps[i];
        n = T->coeffs[i].length;
        n_poly_fit_length(EH->coeffs + i, (d + 2)*n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Ai = ind[j];

            /* set cur = monomial eval */
            p[j] = 1;
            for (k = 2; k < yvar; k++)
            {
                ulong ei = (Aexps[N*Ai + off[k]] >> shift[k]) & mask;
                p[j] = nmod_pow_cache_mulpow_ui(p[j], ei, caches + 3*k + 0,
                          caches + 3*k + 1, caches + 3*k + 2, ctx->fqctx->mod);
            }

            /* copy cur to inc */
            p[j + n] = p[j];

            /* copy coeff */
            _n_fq_set(p + 2*n + d*j, Acoeffs + d*Ai, d);
        }
    }

    TMP_END;

    n_polyun_clear(T);
}


/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets coming from different y's (shouldn't happen)
*/
static slong fq_nmod_mpoly_set_eval_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    fq_nmod_mpolyu_t H,
    const fq_nmod_mpoly_t B,
    n_poly_struct * caches,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    slong * off, * shift;
    ulong y, x, z;
    mp_limb_t * p;
    fq_nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const mp_limb_t * Bcoeffs = B->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    ulong deg;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == H->bits);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    off = (slong *) TMP_ALLOC(2*yvar*sizeof(slong));
    shift = off + yvar;
    for (i = 2; i < yvar; i++)
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, bits, ctx->minfo);

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

        T->exps = FLINT_ARRAY_ALLOC(W->length, ulong);
        T->coeffs = FLINT_ARRAY_ALLOC(W->length, n_poly_struct);
        T->alloc = W->length;
        T->length = 0;
        _clearit(T, W, W->nodes[2 - 1].left);
        mpoly_rbtree_ui_clear(W);
    }

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EH->exps[i] = T->exps[i];
        y = extract_exp(EH->exps[i], 2, 3);
        x = extract_exp(EH->exps[i], 1, 3);
        z = extract_exp(EH->exps[i], 0, 3);
        n = T->coeffs[i].length;
        n_poly_fit_length(EH->coeffs + i, d*3*n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Bi = ind[j];

            /* set cur = monomial eval */
            _n_fq_one(p + d*j, d);
            for (k = 2; k < yvar; k++)
            {
                ulong ei = (Bexps[N*Bi + off[k]] >> shift[k]) & mask;
                n_fq_pow_cache_mulpow_ui(p + d*j, p + d*j, ei, caches + 3*k + 0,
                               caches + 3*k + 1, caches + 3*k + 2, ctx->fqctx);
            }

            /* copy cur to inc */
            _n_fq_set(p + d*(1*n + j), p + d*(0*n + j), d);

            /* copy coeff */
            _n_fq_set(p + d*(2*n + j), Bcoeffs + d*Bi, d);
        }

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fq_nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            fq_nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            for (j = 0; j < n; j++)
            {
                _n_fq_set(Hc->coeffs + d*(old_len + j), p + d*j, d);
                mpoly_monomial_set(Hc->exps + N*(old_len + j), Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                _sort_and_delete_duplicates(Hc, d, ctx->minfo);
            }
        }
    }

    n_polyun_clear(T);

    TMP_END;

    *deg_ = deg;

    return zip_length;
}


static slong fq_nmod_mpoly_set_evalp_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    fq_nmod_mpolyu_t H,
    const fq_nmod_mpoly_t B,
    n_poly_struct * caches,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    slong * off, * shift;
    ulong y, x, z;
    mp_limb_t * p;
    fq_nmod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const mp_limb_t * Bcoeffs = B->coeffs;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    ulong * ind;
    n_polyun_t T;
    ulong deg;
    TMP_INIT;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == H->bits);
    FLINT_ASSERT(Blen > 0);

    TMP_START;

    off = (slong *) TMP_ALLOC(2*yvar*sizeof(slong));
    shift = off + yvar;
    for (i = 2; i < yvar; i++)
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, bits, ctx->minfo);

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

        T->exps = FLINT_ARRAY_ALLOC(W->length, ulong);
        T->coeffs = FLINT_ARRAY_ALLOC(W->length, n_poly_struct);
        T->alloc = W->length;
        T->length = 0;
        _clearit(T, W, W->nodes[2 - 1].left);
        mpoly_rbtree_ui_clear(W);
    }

    n_polyun_fit_length(EH, T->length);
    EH->length = T->length;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EH->exps[i] = T->exps[i];
        y = extract_exp(EH->exps[i], 2, 3);
        x = extract_exp(EH->exps[i], 1, 3);
        z = extract_exp(EH->exps[i], 0, 3);
        n = T->coeffs[i].length;
        n_poly_fit_length(EH->coeffs + i, (d + 2)*n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Bi = ind[j];

            /* set cur = monomial eval */
            p[j] = 1;
            for (k = 2; k < yvar; k++)
            {
                ulong ei = (Bexps[N*Bi + off[k]] >> shift[k]) & mask;
                p[j] = nmod_pow_cache_mulpow_ui(p[j], ei, caches + 3*k + 0,
                          caches + 3*k + 1, caches + 3*k + 2, ctx->fqctx->mod);
            }

            /* copy cur to inc */
            p[j + n] = p[j];

            /* copy coeff */
            _n_fq_set(p + 2*n + d*j, Bcoeffs + d*Bi, d);
        }

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fq_nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            _fq_nmod_mpoly_fit_length(&Hc->coeffs, &Hc->coeffs_alloc, 1,
                                      &Hc->exps, &Hc->exps_alloc, N, n);
            old_len = Hc->length;
            for (j = 0; j < n; j++)
            {
                Hc->coeffs[old_len + j] = p[j];
                mpoly_monomial_set(Hc->exps + N*(old_len + j), Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                _sort_and_delete_duplicates(Hc, 1, ctx->minfo);
            }
        }
    }

    n_polyun_clear(T);

    TMP_END;

    *deg_ = deg;

    return zip_length;
}


static void fq_nmod_polyu_eval_step(
    n_polyu_t E,
    n_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Ai, Ei, n;
    mp_limb_t * p;

    n_polyu_fit_length(E, d*A->length);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = A->exps[Ai];

        n = A->coeffs[Ai].length;
        p = A->coeffs[Ai].coeffs;
        FLINT_ASSERT(A->coeffs[Ai].alloc >= d*3*n);
        _n_fq_zip_eval_step(E->coeffs + d*Ei,
                                p + d*(0*n), p + d*(1*n), p + d*(2*n), n, ctx);

        Ei += !_n_fq_is_zero(E->coeffs + d*Ei, d);
    }
    E->length = Ei;
}

static void fq_nmod_polyu_evalp_step(
    n_polyu_t E,
    n_polyun_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong Ai, Ei, n;
    mp_limb_t * p;

    n_polyu_fit_length(E, d*A->length);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = A->exps[Ai];

        n = A->coeffs[Ai].length;
        p = A->coeffs[Ai].coeffs;
        FLINT_ASSERT(A->coeffs[Ai].alloc >= (d + 2)*n);
        _n_fqp_zip_eval_step(E->coeffs + d*Ei,
                                    p + 0*n, p + 1*n, p + 2*n, n, d, ctx->mod);

        Ei += !_n_fq_is_zero(E->coeffs + d*Ei, d);
    }
    E->length = Ei;
}


static void fq_nmod_polyu3_add_zip_limit1(
    n_polyun_t Z,
    const n_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    const n_fq_poly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    n_fq_poly_struct * Zcoeffs = Z->coeffs;
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
        ai = n_poly_degree(Acoeffs + Ai);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        if (Aexps[Ai] + ai > Zexps[Zi])
        {
            /* missing from Z */
            n_polyun_fit_length(Z, Z->length + 1);
            Zcoeffs = Z->coeffs;
            Zexps = Z->exps;

            for (j = Z->length; j > Zi; j--)
            {
                n_poly_swap(Zcoeffs + j, Zcoeffs + j - 1);
                FLINT_SWAP(ulong, Zexps[j], Zexps[j - 1]);
            }

            Z->length++;

            Zexps[Zi] = Aexps[Ai] + ai;
            n_poly_fit_length(Zcoeffs + Zi, d*fit_length);
            Zcoeffs[Zi].length = cur_length;
            flint_mpn_zero(Zcoeffs[Zi].coeffs, d*cur_length);
            goto in_both;
        }
        else if (Aexps[Ai] + ai < Zexps[Zi])
        {
            /* missing from A */
            FLINT_ASSERT(d*(cur_length + 1) <= Zcoeffs[Zi].alloc);
            _n_fq_zero(Zcoeffs[Zi].coeffs + d*cur_length, d);
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
        }
        else
        {
    in_both:
            FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
            FLINT_ASSERT(d*(cur_length + 1) <= Zcoeffs[Zi].alloc);
            _n_fq_set(Zcoeffs[Zi].coeffs + d*cur_length, Acoeffs[Ai].coeffs + d*ai, d);
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai++;
                } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;

        n_polyun_fit_length(Z, Zi + A->length - Ai);
        Zcoeffs = Z->coeffs;
        Zexps = Z->exps;

        Zexps[Zi] = Aexps[Ai] + ai;
        n_poly_fit_length(Zcoeffs + Zi, d*fit_length);
        Zcoeffs[Zi].length = cur_length;
        flint_mpn_zero(Zcoeffs[Zi].coeffs, d*cur_length);
        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(d*(cur_length + 1) <= Zcoeffs[Zi].alloc);
        _n_fq_set(Zcoeffs[Zi].coeffs + d*cur_length, Acoeffs[Ai].coeffs + d*ai, d);
        Zcoeffs[Zi].length = cur_length + 1;

        Z->length = ++Zi;

        do {
            ai--;
        } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
        if (ai < 0)
        {
            do {
                Ai++;
            } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
            if (Ai < A->length)
                ai = n_poly_degree(Acoeffs + Ai);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(d*(cur_length + 1) <= Zcoeffs[Zi].alloc);
        _n_fq_zero(Zcoeffs[Zi].coeffs + d*cur_length, d);
        Zcoeffs[Zi].length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->coeffs[Zi].length == cur_length + 1);
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
    n_polyun_t M,   /* temp */
    n_poly_stack_t St)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    mp_limb_t * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    fq_nmod_mpoly_struct * Hc;
    slong Hlen = H->length;

    FLINT_ASSERT(bits == H->bits);

    n_polyun_fit_length(M, Hlen + 1);
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
        Hi = mpoly_monomial_index1_nomask(H->exps, H->length, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        fq_nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->coeffs[Hi].length < 1)
            n_fq_poly_product_roots_n_fq(M->coeffs + Hi,
                                       Hc->coeffs, Hc->length, ctx->fqctx, St);

        n_poly_fit_length(M->coeffs + Hlen, d*Hc->length);

        success = _n_fq_zip_vand_solve(Bcoeffs + d*Bi, Hc->coeffs, Hc->length,
                                Z->coeffs[Zi].coeffs, Z->coeffs[Zi].length,
                                M->coeffs[Hi].coeffs, M->coeffs[Hlen].coeffs,
                                                                   ctx->fqctx);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (_n_fq_is_zero(Bcoeffs + d*j, d))
                continue;

            FLINT_ASSERT(N*Bi < B->exps_alloc);
            FLINT_ASSERT(d*Bi < B->coeffs_alloc);

            _n_fq_set(Bcoeffs + d*Bi, Bcoeffs + d*j, d);
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

static int fq_nmod_mpoly_from_zipp(
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
    mp_limb_t * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    fq_nmod_mpoly_struct * Hc;
    slong Hlen = H->length;

    FLINT_ASSERT(bits == H->bits);

    n_polyun_fit_length(M, Hlen + 1);
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
        Hi = mpoly_monomial_index1_nomask(H->exps, H->length, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        fq_nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->coeffs[Hi].length < 1)
            n_poly_mod_product_roots_nmod_vec(M->coeffs + Hi,
                                      Hc->coeffs, Hc->length, ctx->fqctx->mod);

        n_poly_fit_length(M->coeffs + Hlen, Hc->length);

        success = _n_fqp_zip_vand_solve(Bcoeffs + d*Bi, Hc->coeffs, Hc->length,
                                Z->coeffs[Zi].coeffs, Z->coeffs[Zi].length,
                                M->coeffs[Hi].coeffs, M->coeffs[Hlen].coeffs,
                                                                   ctx->fqctx);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (_n_fq_is_zero(Bcoeffs + d*j, d))
                continue;

            FLINT_ASSERT(d*Bi < B->coeffs_alloc);
            FLINT_ASSERT(N*Bi < B->exps_alloc);

            _n_fq_set(Bcoeffs + d*Bi, Bcoeffs + d*j, d);
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
    int success, betas_in_fp;
    slong i;
    slong zip_fails_remaining;
    slong req_zip_images, cur_zip_image;
    fq_nmod_mpolyu_struct * H;
    n_polyun_struct M[1], Aeh[1], * Beh, * BBeval, * Z;
    n_polyu_struct Aeval[1], * Beval;
    fq_nmod_struct * beta;
    n_poly_struct * caches;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_mpoly_t T1, T2;
    n_poly_bpoly_stack_t St;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

#ifdef FLINT_WANT_ASSERT
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

    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    beta = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, fq_nmod_struct);
    for (i = 0; i < ctx->minfo->nvars; i++)
        fq_nmod_init(beta + i, ctx->fqctx);

    /* caches for powers of the betas */
    caches = FLINT_ARRAY_ALLOC(3*ctx->minfo->nvars, n_poly_struct);
    for (i = 0; i < 3*ctx->minfo->nvars; i++)
        n_poly_init(caches + i);

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
    betas_in_fp = (ctx->fqctx->mod.norm < FLINT_BITS/4);
    if (betas_in_fp)
    {
        for (i = 2; i < m; i++)
        {
            ulong bb = n_urandint(state, ctx->fqctx->mod.n - 3) + 2;
            fq_nmod_set_ui(beta + i, bb, ctx->fqctx);
            nmod_pow_cache_start(bb, caches + 3*i + 0,
                                     caches + 3*i + 1, caches + 3*i + 2);
        }

        fq_nmod_mpoly_set_evalp_helper3(Aeh, A, m, caches, ctx);
    }
    else
    {
        for (i = 2; i < m; i++)
        {
            fq_nmod_rand_not_zero(beta + i, state, ctx->fqctx);
            n_fq_pow_cache_start_fq_nmod(beta + i, caches + 3*i + 0,
                               caches + 3*i + 1, caches + 3*i + 2, ctx->fqctx);
        }

        fq_nmod_mpoly_set_eval_helper3(Aeh, A, m, caches, ctx);
    }

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_images = betas_in_fp ?
                         fq_nmod_mpoly_set_evalp_helper_and_zip_form3(
                            Bdegs + i, Beh + i, H + i, B + i, caches, m, ctx) :
                         fq_nmod_mpoly_set_eval_helper_and_zip_form3(
                            Bdegs + i, Beh + i, H + i, B + i, caches, m, ctx);

        req_zip_images = FLINT_MAX(req_zip_images, this_images);
        FLINT_ASSERT(Bdegs[i] > 0);
        Z[i].length = 0;
    }

    for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
    {
        if (betas_in_fp)
        {
            fq_nmod_polyu_evalp_step(Aeval, Aeh, ctx->fqctx);
            for (i = 0; i < r; i++)
                fq_nmod_polyu_evalp_step(Beval + i, Beh + i, ctx->fqctx);
        }
        else
        {
            fq_nmod_polyu_eval_step(Aeval, Aeh, ctx->fqctx);
            for (i = 0; i < r; i++)
                fq_nmod_polyu_eval_step(Beval + i, Beh + i, ctx->fqctx);
        }

        success = n_fq_polyu3_hlift(r, BBeval, Aeval, Beval,
                                         alpha + m - 1, degs0, ctx->fqctx, St);
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
    }

    for (i = 0; i < r; i++)
    {
        success = betas_in_fp ?
             fq_nmod_mpoly_from_zipp(B + i, Z + i, H + i, Bdegs[i], m, ctx, M) :
             fq_nmod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m,
                                                       ctx, M, St->poly_stack);

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

    for (i = 0; i < 3*ctx->minfo->nvars; i++)
        n_poly_clear(caches + i);
    flint_free(caches);

    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

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
    FLINT_ASSERT(_n_fq_is_one(A->coeffs + d*0, d));
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
        _n_fq_poly_normalise(alphabetas + i, d);
    }

    _fq_nmod_mpoly_eval_rest_to_n_fq_bpoly(Ab, A, alphabetas, ctx);
    success = n_fq_bpoly_factor_smprime(Abfc, Abfp, Ab, 0, ctx->fqctx);
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
        fq_nmod_t q, qt;
        fq_nmod_init(q, ctx->fqctx);
        fq_nmod_init(qt, ctx->fqctx);
        FLINT_ASSERT(fq_nmod_mpoly_is_fq_nmod(new_lcs->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fq_nmod_mpoly_length(new_lcs->coeffs + 0*r + i, ctx) == 1);
        _fq_nmod_mpoly_set_n_fq_bpoly_gen1_zero(fac->coeffs + i, newA->bits, Abfp->coeffs + i, 0, ctx);
        FLINT_ASSERT(fac->coeffs[i].length > 0);
        n_fq_get_fq_nmod(qt, fac->coeffs[i].coeffs + d*0, ctx->fqctx);
        fq_nmod_inv(q, qt, ctx->fqctx);
        n_fq_get_fq_nmod(qt, new_lcs->coeffs[0*r + i].coeffs + 0, ctx->fqctx);
        fq_nmod_mul(q, q, qt, ctx->fqctx);
        fq_nmod_mpoly_scalar_mul_fq_nmod(fac->coeffs + i, fac->coeffs + i, q, ctx);
        fq_nmod_clear(q, ctx->fqctx);
        fq_nmod_clear(qt, ctx->fqctx);
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
        for (i = 0; i < r; i++)
        {
            /* hlift should not have returned any large bits */
            FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);

            if (!fq_nmod_mpolyl_content(t, fac->coeffs + i, 1, ctx))
            {
                success = -1;
                goto cleanup;
            }

            success = fq_nmod_mpoly_divides(fac->coeffs + i, fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
        }
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

#ifdef FLINT_WANT_ASSERT
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
