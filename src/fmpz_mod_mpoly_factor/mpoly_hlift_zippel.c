/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod_mpoly_factor.h"
#include "n_poly.h"

void fmpz_mod_polyu_init(fmpz_mod_polyu_t A)
{
    A->length = 0;
    A->alloc = 0;
    A->exps = NULL;
    A->coeffs = NULL;
}

void fmpz_mod_polyu_clear(fmpz_mod_polyu_t A)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_clear(A->coeffs + i);
    flint_free(A->exps);
    flint_free(A->coeffs);
}

void fmpz_mod_polyu_realloc(fmpz_mod_polyu_t A, slong len)
{
    slong i, old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, old_alloc + 1 + old_alloc/2);

    FLINT_ASSERT(A->alloc >= 0);
    if (len <= A->alloc)
        return;

    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, ulong);
    A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, fmpz);

    for (i = old_alloc; i < new_alloc; i++)
        fmpz_init(A->coeffs + i);

    A->alloc = new_alloc;
}

void fmpz_mod_polyu3_degrees(
    slong * deg0,
    slong * deg1,
    slong * deg2,
    const fmpz_mod_polyu_t A)
{
    slong i;
    ulong m;
    ulong mask = mpoly_overflow_mask_sp(FLINT_BITS/3);

    if (A->length <= 0)
    {
        *deg0 = *deg1 = *deg2 = -1;
        return;
    }

    m = A->exps[0];
    for (i = 1; i < A->length; i++)
        m = mpoly_monomial_max1(m, A->exps[i], FLINT_BITS/3, mask);

    *deg0 = extract_exp(m, 2, 3);
    *deg1 = extract_exp(m, 1, 3);
    *deg2 = extract_exp(m, 0, 3);
}



/***************************************************************************/

void fmpz_mod_mpolyu_init(
    fmpz_mod_mpolyu_t A,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}


void fmpz_mod_mpolyu_clear(
    fmpz_mod_mpolyu_t A,
    const fmpz_mod_mpoly_ctx_t uctx)
{
    slong i;
    for (i = 0; i < A->alloc; i++)
        fmpz_mod_mpoly_clear(A->coeffs + i, uctx);
    flint_free(A->coeffs);
    flint_free(A->exps);
}


void fmpz_mod_mpolyu_swap(
    fmpz_mod_mpolyu_t A,
    fmpz_mod_mpolyu_t B,
    const fmpz_mod_mpoly_ctx_t uctx)
{
   fmpz_mod_mpolyu_struct t = *A;
   *A = *B;
   *B = t;
}

void fmpz_mod_mpolyu_zero(
    fmpz_mod_mpolyu_t A,
    const fmpz_mod_mpoly_ctx_t uctx)
{
    A->length = 0;
}

int fmpz_mod_mpolyu_is_one(
    fmpz_mod_mpolyu_t A,
    const fmpz_mod_mpoly_ctx_t uctx)
{
    if (A->length != 1 || A->exps[0] != UWORD(0))
        return 0;

    return fmpz_mod_mpoly_is_one(A->coeffs + 0, uctx);
}

void fmpz_mod_mpolyu_fit_length(
    fmpz_mod_mpolyu_t A,
    slong length,
    const fmpz_mod_mpoly_ctx_t uctx)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(length, 2*A->alloc);

    if (length > old_alloc)
    {
        A->exps = (ulong *) flint_realloc(A->exps, new_alloc*sizeof(ulong));
        A->coeffs = (fmpz_mod_mpoly_struct *) flint_realloc(A->coeffs,
                                      new_alloc*sizeof(fmpz_mod_mpoly_struct));

        for (i = old_alloc; i < new_alloc; i++)
            fmpz_mod_mpoly_init3(A->coeffs + i, 0, A->bits, uctx);

        A->alloc = new_alloc;
    }
}

void fmpz_mod_mpolyu_one(fmpz_mod_mpolyu_t A, const fmpz_mod_mpoly_ctx_t uctx)
{
    fmpz_mod_mpolyu_fit_length(A, WORD(1), uctx);
    A->exps[0] = UWORD(0);
    fmpz_mod_mpoly_one(A->coeffs + 0, uctx);
    A->length = WORD(1);
}

void fmpz_mod_mpolyu_repack_bits_inplace(
    fmpz_mod_mpolyu_t A,
    flint_bitcnt_t bits,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;

    if (bits == A->bits)
        return;

    A->bits = bits;

    for (i = 0; i < A->alloc; i++)
        fmpz_mod_mpoly_repack_bits_inplace(A->coeffs + i, bits, ctx);
}


/* if the coefficient doesn't exist, a new one is created (and set to zero) */
fmpz_mod_mpoly_struct * _fmpz_mod_mpolyu_get_coeff(
    fmpz_mod_mpolyu_t A,
    ulong pow,
    const fmpz_mod_mpoly_ctx_t uctx)
{
    slong i, j;
    fmpz_mod_mpoly_struct * xk;

    for (i = 0; i < A->length && A->exps[i] >= pow; i++)
    {
        if (A->exps[i] == pow)
        {
            return A->coeffs + i;
        }
    }

    fmpz_mod_mpolyu_fit_length(A, A->length + 1, uctx);

    for (j = A->length; j > i; j--)
    {
        A->exps[j] = A->exps[j - 1];
        fmpz_mod_mpoly_swap(A->coeffs + j, A->coeffs + j - 1, uctx);
    }

    A->length++;
    A->exps[i] = pow;
    xk = A->coeffs + i;
    xk->length = 0;
    FLINT_ASSERT(xk->bits == A->bits);

    return xk;
}


/****************************************************************/


static void _delete_duplicates(
    fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    j = -1;
    for (i = 0; i < A->length; i++)
    {
        if (j >= 0 && mpoly_monomial_equal(A->exps + N*j, A->exps + N*i, N))
        {
            FLINT_ASSERT(fmpz_equal(A->coeffs + j, A->coeffs + i));
            continue;
        }
        j++;
        fmpz_set(A->coeffs + j, A->coeffs + i);
        mpoly_monomial_set(A->exps + N*j, A->exps + N*i, N);
    }
    j++;
    A->length = j;
}


static int fmpz_mod_mpoly_from_zip(
    fmpz_mod_mpoly_t B,
    const fmpz_mod_polyun_t Z,
    fmpz_mod_mpolyu_t H,
    ulong deg,
    slong yvar,     /* Y = gen(yvar) */
    const fmpz_mod_mpoly_ctx_t ctx,
    fmpz_mod_polyun_t M)   /* temp */
{
    int success;
    slong Hi, Zi, Bi, i, j;
    slong xvar = 0;
    slong zvar = 1;
    ulong x, y, z;
    flint_bitcnt_t bits = B->bits;
    fmpz * Bcoeffs;
    ulong * Bexps;
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong xoff, xshift, yoff, yshift, zoff, zshift;
    fmpz_mod_poly_struct * Zcoeffs = Z->coeffs;
    ulong * Zexps = Z->exps;
    fmpz_mod_mpoly_struct * Hc;
    slong Hlen = H->length;

    FLINT_ASSERT(bits == H->bits);

    fmpz_mod_polyun_fit_length(M, Hlen + 1, ctx->ffinfo);
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
        y = extract_exp(Zexps[Zi], 2, 3);
        x = extract_exp(Zexps[Zi], 1, 3);
        z = extract_exp(Zexps[Zi], 0, 3);
        FLINT_ASSERT(x < deg);
        Hi = mpoly_monomial_index1_nomask(H->exps, H->length, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        fmpz_mod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->coeffs[Hi].length < 1)
            fmpz_mod_poly_product_roots_fmpz_vec(M->coeffs + Hi,
                                          Hc->coeffs, Hc->length, ctx->ffinfo);

        fmpz_mod_poly_fit_length(M->coeffs + Hlen, Hc->length, ctx->ffinfo);

        success = _fmpz_mod_zip_vand_solve(Bcoeffs + Bi, Hc->coeffs,
                    Hc->length, Zcoeffs[Zi].coeffs, Zcoeffs[Zi].length,
                    M->coeffs[Hi].coeffs, M->coeffs[Hlen].coeffs, ctx->ffinfo);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (fmpz_is_zero(Bcoeffs + j))
                continue;

            FLINT_ASSERT(Bi < B->coeffs_alloc);
            FLINT_ASSERT(N*Bi < B->exps_alloc);

            fmpz_set(Bcoeffs + Bi, Bcoeffs + j);
            mpoly_monomial_set(Bexps + N*Bi, Hc->exps + N*i, N);
            (Bexps + N*Bi)[yoff] += y << yshift;
            Bi++;
        }
    }
    B->length = Bi;
    fmpz_mod_mpoly_sort_terms(B, ctx);
    FLINT_ASSERT(fmpz_mod_mpoly_is_canonical(B, ctx));

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

static void fmpz_mod_mpoly_set_eval_helper3(
    fmpz_mod_polyun_t EH,
    const fmpz_mod_mpoly_t A,
    slong yvar,
    const fmpz * alphas,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    ulong y, x, z;
    slong yoff, xoff, zoff, * off;
    slong yshift, xshift, zshift, * shift;
    fmpz_mod_poly_struct * EHcoeffs;
    ulong * EHexps;
    fmpz * p;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    const fmpz * Acoeffs = A->coeffs;
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

    fmpz_mod_polyun_fit_length(EH, T->length, ctx->ffinfo);
    EH->length = T->length;
    EHcoeffs = EH->coeffs;
    EHexps = EH->exps;

    for (i = 0; i < T->length; i++)
    {
        EHexps[i] = T->exps[i];
        n = T->coeffs[i].length;
        fmpz_mod_poly_fit_length(EHcoeffs + i, 3*n, ctx->ffinfo);
        EHcoeffs[i].length = n;
        p = EHcoeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Ai = ind[j];

            fmpz_one(p + j);
            for (k = 2; k < yvar; k++)
            {
                fmpz_t tt;
                ulong ei = (Aexps[N*Ai + off[k]] >> shift[k]) & mask;
                fmpz_init(tt);
                fmpz_mod_pow_ui(tt, alphas + k, ei, ctx->ffinfo);
                fmpz_mod_mul(p + j, p + j, tt, ctx->ffinfo);
                fmpz_clear(tt);
            }

            /* set cur = monomial eval */

            /* copy cur to inc */
            fmpz_set(p + j + n, p + j);

            /* copy coeff */
            fmpz_set(p + j + 2*n, Acoeffs + Ai);
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
static slong fmpz_mod_mpoly_set_eval_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    fmpz_mod_polyun_t EH,
    fmpz_mod_mpolyu_t H,
    const fmpz_mod_mpoly_t B,
    fmpz * alphas,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, k, n;
    slong * off, * shift;
    ulong y, x, z;
    fmpz_mod_poly_struct * EHcoeffs;
    ulong * EHexps;
    fmpz * p;
    fmpz_mod_mpoly_struct * Hc;
    slong old_len, zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    const fmpz * Bcoeffs = B->coeffs;
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

    fmpz_mod_polyun_fit_length(EH, T->length, ctx->ffinfo);
    EH->length = T->length;
    EHcoeffs = EH->coeffs;
    EHexps = EH->exps;

    H->length = 0;

    for (i = 0; i < T->length; i++)
    {
        EHexps[i] = T->exps[i];
        y = extract_exp(EHexps[i], 2, 3);
        x = extract_exp(EHexps[i], 1, 3);
        z = extract_exp(EHexps[i], 0, 3);
        n = T->coeffs[i].length;
        fmpz_mod_poly_fit_length(EHcoeffs + i, 3*n, ctx->ffinfo);
        EHcoeffs[i].length = n;
        p = EHcoeffs[i].coeffs;
        ind = T->coeffs[i].coeffs;

        for (j = 0; j < n; j++)
        {
            slong Bi = ind[j];

            fmpz_one(p + j);

            for (k = 2; k < yvar; k++)
            {
                fmpz_t tt;
                ulong ei = (Bexps[N*Bi + off[k]] >> shift[k]) & mask;
                fmpz_init(tt);
                fmpz_init(tt);
                fmpz_mod_pow_ui(tt, alphas + k, ei, ctx->ffinfo);
                fmpz_mod_mul(p + j, p + j, tt, ctx->ffinfo);
                fmpz_clear(tt);
            }

            /* set cur = monomial eval */

            /* copy cur to inc */
            fmpz_set(p + j + n, p + j);

            /* copy coeff */
            fmpz_set(p + j + 2*n, Bcoeffs + Bi);
        }

        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _fmpz_mod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            fmpz_mod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            _fmpz_vec_set(Hc->coeffs + old_len, p, n);
            for (j = 0; j < n; j++)
                mpoly_monomial_set(Hc->exps + N*(old_len + j), Bexps + N*ind[j], N);

            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                fmpz_mod_mpoly_sort_terms(Hc, ctx);
                _delete_duplicates(Hc, ctx);
            }
        }
    }

    n_polyun_clear(T);

    TMP_END;

    *deg_ = deg;

    return zip_length;
}


static void fmpz_mod_polyu_eval_step(
    fmpz_mod_polyu_t E,
    fmpz_mod_polyun_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong Ai, Ei, n;
    fmpz * p;

    fmpz_mod_polyu_fit_length(E, A->length, ctx);

    Ei = 0;
    for (Ai = 0; Ai < A->length; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        E->exps[Ei] = A->exps[Ai];

        n = A->coeffs[Ai].length;
        p = A->coeffs[Ai].coeffs;
        FLINT_ASSERT(A->coeffs[Ai].alloc >= 3*n);
        _fmpz_mod_zip_eval_step(E->coeffs + Ei, p + 0*n, p + 1*n, p + 2*n, n, ctx);

        Ei += !fmpz_is_zero(E->coeffs + Ei);
    }
    E->length = Ei;
}



static void fmpz_mod_polyu3_add_zip_limit1(
    fmpz_mod_polyun_t Z,
    const fmpz_mod_polyun_t A,
    const ulong deg1,
    slong cur_length,
    slong fit_length,
    const fmpz_mod_ctx_t ctx)
{
    const fmpz_mod_poly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    fmpz_mod_poly_struct * Zcoeffs = Z->coeffs;
    ulong * Zexps = Z->exps;
    slong Ai, ai, Zi, j;

    Ai = -1;
    ai = -1;
    do {
        Ai++;
    } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
    if (Ai < A->length)
        ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx);

    Zi = 0;

    while (Ai < A->length && Zi < Z->length)
    {
        if (Aexps[Ai] + ai > Zexps[Zi])
        {
            /* missing from Z */
            fmpz_mod_polyun_fit_length(Z, Z->length + 1, ctx);
            Zcoeffs = Z->coeffs;
            Zexps = Z->exps;

            for (j = Z->length; j > Zi; j--)
            {
                fmpz_mod_poly_swap(Zcoeffs + j, Zcoeffs + j - 1, ctx);
                FLINT_SWAP(ulong, Zexps[j], Zexps[j - 1]);
            }

            Z->length++;

            Zexps[Zi] = Aexps[Ai] + ai;
            fmpz_mod_poly_fit_length(Zcoeffs + Zi, fit_length, ctx);
            Zcoeffs[Zi].length = cur_length;
            _fmpz_vec_zero(Zcoeffs[Zi].coeffs, cur_length);
            goto in_both;
        }
        else if (Aexps[Ai] + ai < Zexps[Zi])
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
            fmpz_zero(Zcoeffs[Zi].coeffs + cur_length);
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
        }
        else
        {
    in_both:
            FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
            FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
            fmpz_set(Zcoeffs[Zi].coeffs + cur_length, Acoeffs[Ai].coeffs + ai);
            Zcoeffs[Zi].length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
            if (ai < 0)
            {
                do {
                    Ai++;
                } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
                if (Ai < A->length)
                    ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx);
            }
        }
    }

    /* everything in A must be put on the end of Z */
    while (Ai < A->length)
    {
        Zi = Z->length;

        fmpz_mod_polyun_fit_length(Z, Zi + A->length - Ai, ctx);
        Zcoeffs = Z->coeffs;
        Zexps = Z->exps;

        Zexps[Zi] = Aexps[Ai] + ai;
        fmpz_mod_poly_fit_length(Zcoeffs + Zi, fit_length, ctx);
        Zcoeffs[Zi].length = cur_length;
        _fmpz_vec_zero(Zcoeffs[Zi].coeffs, cur_length);
        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
        fmpz_set(Zcoeffs[Zi].coeffs + cur_length, Acoeffs[Ai].coeffs + ai);
        Zcoeffs[Zi].length = cur_length + 1;

        Z->length = ++Zi;

        do {
            ai--;
        } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
        if (ai < 0)
        {
            do {
                Ai++;
            } while (Ai < A->length && extract_exp(Aexps[Ai], 1, 3) >= deg1);
            if (Ai < A->length)
                ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx);
        }
    }

    /* everything in Z must have a zero appended */
    while (Zi < Z->length)
    {
        FLINT_ASSERT(cur_length == Zcoeffs[Zi].length);
        FLINT_ASSERT(cur_length + 1 <= Zcoeffs[Zi].alloc);
        fmpz_zero(Zcoeffs[Zi].coeffs + cur_length);
        Zcoeffs[Zi].length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->coeffs[Zi].length == cur_length + 1);
    }
}


int fmpz_mod_mpoly_hlift_zippel(
    slong m,
    fmpz_mod_mpoly_struct * B,
    slong r,
    const fmpz * alpha,
    const fmpz_mod_mpoly_t A,
    const slong * degs,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    flint_bitcnt_t bits = A->bits;
    int success;
    slong i, zip_fails_remaining, req_zip_images, cur_zip_image;
    fmpz_mod_mpolyu_struct * H;
    fmpz_mod_polyun_struct M[1], Aeh[1], * Beh, * BBeval, * Z;
    fmpz_mod_polyu_struct Aeval[1], * Beval;
    fmpz * beta;
    fmpz_mod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

#ifdef FLINT_WANT_ASSERT
    {
        fmpz_mod_mpoly_t T;
        slong j, * check_degs = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, slong);

        fmpz_mod_mpoly_init(T, ctx);

        fmpz_mod_mpoly_degrees_si(check_degs, A, ctx);
        for (j = 0; j < ctx->minfo->nvars; j++)
            FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);

        fmpz_mod_mpoly_one(T, ctx);
        for (i = 0; i < r; i++)
        {
            fmpz_mod_mpoly_degrees_si(check_degs, B + i, ctx);
            for (j = 0; j < ctx->minfo->nvars; j++)
                FLINT_ASSERT(FLINT_BIT_COUNT(check_degs[j]) < FLINT_BITS/3);
            fmpz_mod_mpoly_mul(T, T, B + i, ctx);
        }
        fmpz_mod_mpoly_sub(T, A, T, ctx);

        fmpz_mod_mpoly_evaluate_one_fmpz(T, T, m, alpha + m - 1, ctx);
        FLINT_ASSERT(fmpz_mod_mpoly_is_zero(T, ctx));

        fmpz_mod_mpoly_clear(T, ctx);
        flint_free(check_degs);
    }
#endif

    beta = _fmpz_vec_init(ctx->minfo->nvars);

    Bdegs = FLINT_ARRAY_ALLOC(r, ulong);
    H = FLINT_ARRAY_ALLOC(r, fmpz_mod_mpolyu_struct);
    Beh = FLINT_ARRAY_ALLOC(r, fmpz_mod_polyun_struct);
    Beval = FLINT_ARRAY_ALLOC(r, fmpz_mod_polyu_struct);
    BBeval = FLINT_ARRAY_ALLOC(r, fmpz_mod_polyun_struct);
    Z = FLINT_ARRAY_ALLOC(r, fmpz_mod_polyun_struct);

    fmpz_mod_polyun_init(Aeh, ctx->ffinfo);
    fmpz_mod_polyu_init(Aeval);
    fmpz_mod_polyun_init(M, ctx->ffinfo);
    for (i = 0; i < r; i++)
    {
        fmpz_mod_mpolyu_init(H + i, bits, ctx);
        fmpz_mod_polyun_init(Beh + i, ctx->ffinfo);
        fmpz_mod_polyu_init(Beval + i);
        fmpz_mod_polyun_init(BBeval + i, ctx->ffinfo);
        fmpz_mod_polyun_init(Z + i, ctx->ffinfo);
    }

    /* init done */

    for (i = 0; i < r; i++)
    {
        success = fmpz_mod_mpoly_repack_bits_inplace(B + i, bits, ctx);
        if (!success)
            goto cleanup;
    }

    zip_fails_remaining = 3;

choose_betas:

    /* only beta[2], beta[3], ..., beta[m - 1] will be used */
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_mod_rand_not_zero(beta + i, state, ctx->ffinfo);

    fmpz_mod_mpoly_set_eval_helper3(Aeh, A, m, beta, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_images;
        this_images = fmpz_mod_mpoly_set_eval_helper_and_zip_form3(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_images);
        FLINT_ASSERT(Bdegs[i] > 0);
        Z[i].length = 0;
    }

    cur_zip_image = 0;

next_zip_image:

    fmpz_mod_polyu_eval_step(Aeval, Aeh, ctx->ffinfo);
    for (i = 0; i < r; i++)
        fmpz_mod_polyu_eval_step(Beval + i, Beh + i, ctx->ffinfo);

    success = fmpz_mod_polyu3_hlift(r, BBeval, Aeval, Beval,
                                            alpha + m - 1, degs0, ctx->ffinfo);
    if (success < 1)
    {
        if (--zip_fails_remaining >= 0)
            goto choose_betas;
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        fmpz_mod_polyu3_add_zip_limit1(Z + i, BBeval + i, Bdegs[i],
                                   cur_zip_image, req_zip_images, ctx->ffinfo);
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = fmpz_mod_mpoly_from_zip(B + i, Z + i, H + i, Bdegs[i], m, ctx, M);
        if (success < 1)
        {
FLINT_ASSERT(0);
            success = 0;
            goto cleanup;
        }
    }

    fmpz_mod_mpoly_init3(T1, A->length, bits, ctx);
    fmpz_mod_mpoly_init3(T2, A->length, bits, ctx);
    fmpz_mod_mpoly_mul(T1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fmpz_mod_mpoly_mul(T2, T1, B + i, ctx);
        fmpz_mod_mpoly_swap(T1, T2, ctx);
    }

    success = fmpz_mod_mpoly_equal(T1, A, ctx);
    fmpz_mod_mpoly_clear(T1, ctx);
    fmpz_mod_mpoly_clear(T2, ctx);

cleanup:

    fmpz_mod_polyun_clear(Aeh, ctx->ffinfo);
    fmpz_mod_polyu_clear(Aeval);
    fmpz_mod_polyun_clear(M, ctx->ffinfo);
    for (i = 0; i < r; i++)
    {
        fmpz_mod_mpolyu_clear(H + i, ctx);
        fmpz_mod_polyun_clear(Beh + i, ctx->ffinfo);
        fmpz_mod_polyu_clear(Beval + i);
        fmpz_mod_polyun_clear(BBeval + i, ctx->ffinfo);
        fmpz_mod_polyun_clear(Z + i, ctx->ffinfo);
    }

    _fmpz_vec_clear(beta, ctx->minfo->nvars);

    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

FLINT_ASSERT(success);

    return success;
}
