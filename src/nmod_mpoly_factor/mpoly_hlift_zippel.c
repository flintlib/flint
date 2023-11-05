/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


static void _delete_duplicates(
    nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(A->coeffs[j] == A->coeffs[i]);
            continue;
        }
        j++;
        A->coeffs[j] = A->coeffs[i];
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}


static int nmod_mpoly_from_zip(
    nmod_mpoly_t B,
    const n_polyun_t Z,
    nmod_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const nmod_mpoly_ctx_t ctx,
    n_polyun_t M)   /* temp */
{
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
    nmod_mpoly_struct * Hc;
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
        nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        /* fill in root product if it is missing */
        if (M->coeffs[Hi].length < 1)
            n_poly_mod_product_roots_nmod_vec(M->coeffs + Hi,
                                             Hc->coeffs, Hc->length, ctx->mod);

        n_poly_fit_length(M->coeffs + Hlen, Hc->length);

        success = _nmod_zip_vand_solve(Bcoeffs + Bi, Hc->coeffs, Hc->length,
                                Z->coeffs[Zi].coeffs, Z->coeffs[Zi].length,
                                M->coeffs[Hi].coeffs, M->coeffs[Hlen].coeffs,
                                                                     ctx->mod);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (Bcoeffs[j] == 0)
                continue;

            FLINT_ASSERT(Bi < B->coeffs_alloc);
            FLINT_ASSERT(N*Bi < B->exps_alloc);

            Bcoeffs[Bi] = Bcoeffs[j];
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    nmod_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(nmod_mpoly_is_canonical(B, ctx));

    return 1;
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

static void nmod_mpoly_set_eval_helper3(
    n_polyun_t EH,
    const nmod_mpoly_t A,
    slong yvar,
    n_poly_struct * caches,
    const nmod_mpoly_ctx_t ctx)
{
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

    T->exps   = FLINT_ARRAY_ALLOC(W->length, ulong);
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
        n_poly_fit_length(EH->coeffs + i, 3*n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Ai = ind[j];
            mp_limb_t meval = 1;

            for (k = 2; k < yvar; k++)
            {
                ulong ei = (Aexps[N*Ai + off[k]] >> shift[k]) & mask;
                meval = nmod_pow_cache_mulpow_ui(meval, ei, caches + 3*k + 0,
                                 caches + 3*k + 1, caches + 3*k + 2, ctx->mod);
            }

            /* set cur = monomial eval */
            p[j] = meval;

            /* copy cur to inc */
            p[j + n] = meval;

            /* copy coeff */
            p[j + 2*n] = Acoeffs[Ai];
        }
    }

    n_polyun_clear(T);

    TMP_END;
}

/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets coming from different y's (shouldn't happen)
*/
static slong nmod_mpoly_set_eval_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    nmod_mpolyu_t H,
    const nmod_mpoly_t B,
    n_poly_struct * caches,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const nmod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    slong * off, * shift;
    ulong y, x, z;
    mp_limb_t * p;
    nmod_mpoly_struct * Hc;
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

        T->exps   = FLINT_ARRAY_ALLOC(W->length, ulong);
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
        n_poly_fit_length(EH->coeffs + i, 3*n);
        EH->coeffs[i].length = n;
        p = EH->coeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Bi = ind[j];
            mp_limb_t meval = 1;

            for (k = 2; k < yvar; k++)
            {
                ulong ei = (Bexps[N*Bi + off[k]] >> shift[k]) & mask;
                meval = nmod_pow_cache_mulpow_ui(meval, ei, caches + 3*k + 0,
                                 caches + 3*k + 1, caches + 3*k + 2, ctx->mod);
            }

            /* set cur = monomial eval */
            p[j] = meval;

            /* copy cur to inc */
            p[j + n] = meval;

            /* copy coeff */
            p[j + 2*n] = Bcoeffs[Bi];
        }

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            flint_mpn_copyi(Hc->coeffs + old_len, p, n);
            for (j = 0; j < n; j++)
                mpoly_monomial_set(Hc->exps + N*(old_len + j), Bexps + N*ind[j], N);

            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                nmod_mpoly_sort_terms(Hc, ctx);
                _delete_duplicates(Hc, ctx);
            }
        }
    }

    n_polyun_clear(T);

    TMP_END;

    *deg_ = deg;

    return zip_length;
}


static void n_polyu_mod_eval_step(n_polyu_t E, n_polyun_t A, nmod_t ctx)
{
    slong Ai, Ei, n;
    mp_limb_t * p;

    n_polyu_fit_length(E, A->length);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = A->exps[Ai];

        n = A->coeffs[Ai].length;
        p = A->coeffs[Ai].coeffs;
        FLINT_ASSERT(A->coeffs[Ai].alloc >= 3*n);
        E->coeffs[Ei] = _nmod_zip_eval_step(p + 0*n, p + 1*n, p + 2*n, n, ctx);

        Ei += (E->coeffs[Ei] != 0);
    }
    E->length = Ei;
}


static void n_polyu3_add_zip_limit1(
    n_polyun_t Z,
    const n_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length)
{
    const n_poly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    n_poly_struct * Zcoeffs = Z->coeffs;
    ulong * Zexps = Z->exps;
    slong Ai, ai, Zi, j;

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
            n_poly_fit_length(Zcoeffs + Zi, fit_length);
            Zcoeffs[Zi].length = cur_length;
            flint_mpn_zero(Zcoeffs[Zi].coeffs, cur_length);
            goto in_both;
        }
        else if (Aexps[Ai] + ai < Zexps[Zi])
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
            Zcoeffs[Zi].coeffs[cur_length] = 0;
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
        }
        else
        {
    in_both:
            FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
            FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
            Zcoeffs[Zi].coeffs[cur_length] = Acoeffs[Ai].coeffs[ai];
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
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
        n_poly_fit_length(Zcoeffs + Zi, fit_length);
        Zcoeffs[Zi].length = cur_length;
        flint_mpn_zero(Zcoeffs[Zi].coeffs, cur_length);
        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
        Zcoeffs[Zi].coeffs[cur_length] = Acoeffs[Ai].coeffs[ai];
        Zcoeffs[Zi].length = cur_length + 1;

        Z->length = ++Zi;

        do {
            ai--;
        } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
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
        FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
        Zcoeffs[Zi].coeffs[cur_length] = 0;
        Zcoeffs[Zi].length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->coeffs[Zi].length == cur_length + 1);
    }
}


int nmod_mpoly_hlift_zippel(
    slong m,
    nmod_mpoly_struct * B,
    slong r,
    const mp_limb_t * alpha,
    const nmod_mpoly_t A,
    const slong * degs,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    flint_bitcnt_t bits = A->bits;
    int success;
    slong i, zip_fails_remaining, req_zip_images, cur_zip_image;
    nmod_mpolyu_struct * H;
    n_polyun_struct M[1], Aeh[1], * Beh, * BBeval, * Z;
    n_polyu_struct Aeval[1], * Beval;
    mp_limb_t * beta;
    n_poly_struct * caches;
    nmod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

#ifdef FLINT_WANT_ASSERT
    {
        nmod_mpoly_t T;
        slong j, * check_degs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);

        nmod_mpoly_init(T, ctx);

        nmod_mpoly_degrees_si(check_degs, A, ctx);
        for (j = 0; j < ctx->minfo->nvars; j++)
            FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);

        nmod_mpoly_one(T, ctx);
        for (i = 0; i < r; i++)
        {
            nmod_mpoly_degrees_si(check_degs, B + i, ctx);
            for (j = 0; j < ctx->minfo->nvars; j++)
                FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);
            nmod_mpoly_mul(T, T, B + i, ctx);
        }
        nmod_mpoly_sub(T, A, T, ctx);

        nmod_mpoly_evaluate_one_ui(T, T, m, alpha[m - 1], ctx);
        FLINT_ASSERT(nmod_mpoly_is_zero(T, ctx));

        nmod_mpoly_clear(T, ctx);
        flint_free(check_degs);
    }
#endif

    beta = FLINT_ARRAY_ALLOC(ctx->minfo->nvars,mp_limb_t);

    /* caches for powers of the betas */
    caches = FLINT_ARRAY_ALLOC(3*ctx->minfo->nvars, n_poly_struct);
    for (i = 0; i < 3*ctx->minfo->nvars; i++)
        n_poly_init(caches + i);

    Bdegs = FLINT_ARRAY_ALLOC(r, ulong);
    H = FLINT_ARRAY_ALLOC(r, nmod_mpolyu_struct);
    Beh = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    Beval = FLINT_ARRAY_ALLOC(r, n_polyu_struct);
    BBeval = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    Z = FLINT_ARRAY_ALLOC(r, n_polyun_struct);

    n_polyun_init(Aeh);
    n_polyu_init(Aeval);
    n_polyun_init(M);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyu_init(H + i, bits, ctx);
        n_polyun_init(Beh + i);
        n_polyu_init(Beval + i);
        n_polyun_init(BBeval + i);
        n_polyun_init(Z + i);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = nmod_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    FLINT_ASSERT(ctx->mod.n > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        beta[i] = n_urandint(state, ctx->mod.n - 3) + 2;
        nmod_pow_cache_start(beta[i], caches + 3*i + 0,
                                      caches + 3*i + 1, caches + 3*i + 2);
    }

    nmod_mpoly_set_eval_helper3(Aeh, A, m, caches, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_images;
        this_images = nmod_mpoly_set_eval_helper_and_zip_form3(Bdegs + i,
                                        Beh + i, H + i, B + i, caches, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_images);
        FLINT_ASSERT(Bdegs[i] > 0);
        Z[i].length = 0;
    }

    cur_zip_image = 0;

next_zip_image:

    n_polyu_mod_eval_step(Aeval, Aeh, ctx->mod);
    for (i = 0; i < r; i++)
        n_polyu_mod_eval_step(Beval + i, Beh + i, ctx->mod);

    success = n_polyu3_mod_hlift(r, BBeval, Aeval, Beval,
                                             alpha[m - 1], degs0, ctx->mod);
    if (success < 1)
    {
        if (--zip_fails_remaining >= 0)
            goto choose_betas;
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        n_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                                cur_zip_image, req_zip_images);
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = nmod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx, M);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    nmod_mpoly_init3(T1, A->length, bits, ctx);
    nmod_mpoly_init3(T2, A->length, bits, ctx);
    nmod_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        nmod_mpoly_mul(T2, T1, B + i, ctx);
        nmod_mpoly_swap(T1, T2, ctx);
    }

    success = nmod_mpoly_equal(T1, A, ctx);
    nmod_mpoly_clear(T1, ctx);
    nmod_mpoly_clear(T2, ctx);

cleanup:

    n_polyun_clear(Aeh);
    n_polyu_clear(Aeval);
    n_polyun_clear(M);
    for (i = 0; i < r; i++)
    {
        nmod_mpolyu_clear(H + i, ctx);
        n_polyun_clear(Beh + i);
        n_polyu_clear(Beval + i);
        n_polyun_clear(BBeval + i);
        n_polyun_clear(Z + i);
    }

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

    return success;
}

