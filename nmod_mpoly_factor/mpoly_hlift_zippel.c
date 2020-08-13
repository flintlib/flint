/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


static void nmod_mpoly_delete_duplicate_terms(
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


static slong nmod_mpolyu_find_term(const nmod_mpolyu_t A, ulong e)
{
    slong start, i, stop;
    ulong * Aexps = A->exps;

    start = 0;
    stop = A->length;

again:

    if (stop - start < 8)
    {
        for (i = start; i < stop; i++)
        {
            if (Aexps[i] == e)
                return i;
        }
        return -1;
    }

    i = start + (stop - start)/2;

    FLINT_ASSERT(Aexps[start] > Aexps[i]);
    FLINT_ASSERT(stop >= A->length || Aexps[stop] < Aexps[i]);

    if (Aexps[i] < e)
    {
        stop = i;
        goto again;
    }
    else if (Aexps[i] > e)
    {
        start = i;
        goto again;
    }

    return i;
}


void _nmod_mpoly_monomial_evals(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong nvars = ctx->minfo->nvars;
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    TMP_INIT;

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
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

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if (((Aexps + N*i)[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}

static void _nmod_mpoly_monomial_evals_indirect(
    mp_limb_t * E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    ulong * Aind,
    slong Alen,
    const mp_limb_t * alpha,
    slong vstart,
    slong vstop,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    slong offset, shift;
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * LUToffset;
    ulong * LUTmask;
    mp_limb_t * LUTvalue;
    slong LUTlen;
    mp_limb_t xpoweval;
    ulong * inputexpmask;
    const ulong * thisAexp;
    TMP_INIT;

    FLINT_ASSERT(0 <= vstart);
    FLINT_ASSERT(vstart < vstop);
    FLINT_ASSERT(vstop <= ctx->minfo->nvars);

    TMP_START;

    LUToffset = (slong *) TMP_ALLOC(N*FLINT_BITS*sizeof(slong));
    LUTmask   = (ulong *) TMP_ALLOC(N*FLINT_BITS*sizeof(ulong));
    LUTvalue  = (mp_limb_t *) TMP_ALLOC(N*FLINT_BITS*sizeof(mp_limb_t));

    inputexpmask = (ulong *) TMP_ALLOC(N*sizeof(ulong));
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

        xpoweval = alpha[j]; /* xpoweval = alpha[i]^(2^i) */
        for (i = 0; i < Abits; i++)
        {
            LUToffset[LUTlen] = offset;
            LUTmask[LUTlen] = (UWORD(1) << (shift + i));
            LUTvalue[LUTlen] = xpoweval;
            if ((inputexpmask[offset] & LUTmask[LUTlen]) != 0)
                LUTlen++;
            xpoweval = nmod_mul(xpoweval, xpoweval, ctx->ffinfo->mod);
        }
    }
    FLINT_ASSERT(LUTlen < N*FLINT_BITS);

    for (i = 0; i < Alen; i++)
    {
        thisAexp = Aexps + N*Aind[i];
        xpoweval = 1;
        for (j = 0; j < LUTlen; j++)
        {
            if ((thisAexp[LUToffset[j]] & LUTmask[j]) != 0)
            {
                xpoweval = nmod_mul(xpoweval, LUTvalue[j], ctx->ffinfo->mod);
            }
        }
        E[i] = xpoweval;
    }

    TMP_END;
}


/*
    return 
        -1: singular
        0:  inconsistent
        1:  success
*/
int nmod_zip_find_coeffs_new(
    mp_limb_t * coeffs,             /* length mlength */
    const mp_limb_t * monomials,    /* length mlength */
    slong mlength,
    const mp_limb_t * evals,        /* length elength */
    slong elength,
    const mp_limb_t * master,       /* length mlength + 1 */
    mp_limb_t * scratch,            /* length mlength */
    nmod_t ctx)
{
    slong i, j;
    mp_limb_t V, V0, V1, V2, T, S, r, p0, p1;

    FLINT_ASSERT(elength >= mlength);

    for (i = 0; i < mlength; i++)
    {
        /* coeffs[i] is (coeffs(P).values)/P(roots[i]) =: V/S
            where P(x) = master(x)/(x-roots[i])     */
        V0 = V1 = V2 = T = S = 0;
        r = monomials[i];
        for (j = mlength; j > 0; j--)
        {
            T = nmod_add(nmod_mul(r, T, ctx), master[j], ctx);
            S = nmod_add(nmod_mul(r, S, ctx), T, ctx);
            umul_ppmm(p1, p0, evals[j - 1], T);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        /* roots[i] should be a root of master */
        FLINT_ASSERT(nmod_add(nmod_mul(r, T, ctx), master[0], ctx) == 0);
        NMOD_RED3(V, V2, V1, V0, ctx);
        S = nmod_mul(S, r, ctx); /* shift is one */
        if (S == 0)
            return -1;
        coeffs[i] = nmod_mul(V, nmod_inv(S, ctx), ctx);
    }

    /* check that the remaining points match */
    for (j = 0; j < mlength; j++)
        scratch[j] = nmod_pow_ui(monomials[j], mlength, ctx);

    for (i = mlength; i < elength; i++)
    {
        V0 = V1 = V2 = S = 0;
        for (j = 0; j < mlength; j++)
        {
            scratch[j] = nmod_mul(scratch[j], monomials[j], ctx);
            umul_ppmm(p1, p0, coeffs[j], scratch[j]);
            add_sssaaaaaa(V2, V1, V0, V2, V1, V0, 0, p1, p0);
        }
        NMOD_RED3(V, V2, V1, V0, ctx);
        if (V != evals[i])
            return 0;
    }
    return 1;
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
    n_polyun_term_struct * Zt = Z->terms;
    nmod_mpoly_struct * Hc;
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
        Hi = nmod_mpolyu_find_term(H, pack_exp3(0, x, z));
        if (Hi < 0)
            return 0;

        FLINT_ASSERT(Hi < Hlen);
        FLINT_ASSERT(H->exps[Hi] == pack_exp3(0, x, z));

        Hc = H->coeffs + Hi;
        FLINT_ASSERT(bits == Hc->bits);
        FLINT_ASSERT(Hc->length > 0);
        nmod_mpoly_fit_length(B, Bi + Hc->length, ctx);
        Bcoeffs = B->coeffs;

        if (M->terms[Hi].coeff->length < 1)
        {
            n_poly_mod_product_roots_nmod_vec(M->terms[Hi].coeff,
                                     Hc->coeffs, Hc->length, ctx->ffinfo->mod);
        }

        n_poly_fit_length(M->terms[Hlen].coeff, Hc->length);

        success = nmod_zip_find_coeffs_new(Bcoeffs + Bi, Hc->coeffs,
                    Hc->length, Zt[Zi].coeff->coeffs, Zt[Zi].coeff->length,
                    M->terms[Hi].coeff->coeffs, M->terms[Hlen].coeff->coeffs,
                                                             ctx->ffinfo->mod);
        if (success < 1)
            return success;

        Bexps = B->exps;
        for (j = Bi, i = 0; i < Hc->length; j++, i++)
        {
            if (Bcoeffs[j] == 0)
                continue;
            Bcoeffs[Bi] = Bcoeffs[j];
            FLINT_ASSERT(Bi < B->alloc);
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
    W->terms[W->length].exp = nodes[idx].key;
    W->terms[W->length].coeff[0] = ((n_poly_struct *) T->data)[idx];
    W->length++;

    if (nodes[idx].left >= 0)
        _clearit(W, T, nodes[idx].left);
}

static void nmod_mpoly_set_eval_helper3(
    n_polyun_t EH,
    const nmod_mpoly_t A,
    slong yvar,
    const mp_limb_t * alpha,
    const nmod_mpoly_ctx_t ctx)
{
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
    const mp_limb_t * Acoeffs = A->coeffs;
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
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _nmod_mpoly_monomial_evals_indirect(p, Aexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = Acoeffs[ind[j]];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    n_polyun_clear(T);
}

/*
    for each term Y^y*X^x*Z^z * pol(x1,...) in B with j < deg
    set Y^0*X^x*Z^z in H as the monomials with the monomial evals as coeffs
        merge monomial sets comming from different y's (shouldn't happen)
*/
static slong nmod_mpoly_set_eval_helper_and_zip_form3(
    ulong * deg_,       /* deg_X(B), output */
    n_polyun_t EH,
    nmod_mpolyu_t H,
    const nmod_mpoly_t B,
    const mp_limb_t * alpha,
    slong yvar,         /* Y = gen(yvar) (X = gen(0), Z = gen(1))*/
    const nmod_mpoly_ctx_t ctx)
{
    slong xvar = 0;
    slong zvar = 1;
    slong i, j, n;
    ulong y, x, z;
    n_polyun_term_struct * EHterms;
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
        n_poly_fit_length(EHterms[i].coeff, 3*n);
        EHterms[i].coeff->length = n;
        p = EHterms[i].coeff->coeffs;
        ind = T->terms[i].coeff->coeffs;
        _nmod_mpoly_monomial_evals_indirect(p, Bexps, bits, ind, n, alpha,
                                                                2, yvar, ctx);
        if (x < deg)
        {
            FLINT_ASSERT(y == 0 && "strange but ok");
            Hc = _nmod_mpolyu_get_coeff(H, pack_exp3(0, x, z), ctx);
            nmod_mpoly_fit_length(Hc, n, ctx);
            old_len = Hc->length;
            flint_mpn_copyi(Hc->coeffs + old_len, p, n);
            for (j = 0; j < n; j++)
            {
                mpoly_monomial_set(Hc->exps + N*(old_len + j),
                                   Bexps + N*ind[j], N);
            }
            Hc->length += n;
            zip_length = FLINT_MAX(zip_length, Hc->length);
            if (old_len > 0)
            {
                FLINT_ASSERT(0 && "strange but ok");
                nmod_mpoly_sort_terms(Hc, ctx);
                nmod_mpoly_delete_duplicate_terms(Hc, ctx);
            }
        }

        for (j = n - 1; j >= 0; j--)
        {
            mp_limb_t t1 = p[j];
            mp_limb_t t2 = Bcoeffs[ind[j]];
            p[3*j + 0] = t1;
            p[3*j + 1] = t2;
            p[3*j + 2] = t1;
        }
    }

    n_polyun_clear(T);

    *deg_ = deg;

    return zip_length;
}


mp_limb_t n_poly_mod_eval_step(n_poly_t A, nmod_t mod)
{
    slong i, Alen = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    ulong t0, t1, t2, p0, p1;

    FLINT_ASSERT(3*Alen <= A->alloc);

    t2 = t1 = t0 = 0;
    for (i = 0; i < Alen; i++)
    {
        umul_ppmm(p1, p0, Acoeffs[3*i + 0], Acoeffs[3*i + 1]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        Acoeffs[3*i + 0] = nmod_mul(Acoeffs[3*i + 0], Acoeffs[3*i + 2], mod);
    }
    NMOD_RED3(t0, t2, t1, t0, mod);
    return t0;
}


static void n_polyu_mod_eval_step(n_polyu_t E, n_polyun_t A, nmod_t mod)
{
    slong Ai, Ei;
    ulong * Eexps;
    mp_limb_t * Ecoeffs;
    n_polyun_term_struct * Aterms;
    slong Alen = A->length;

    n_polyu_fit_length(E, Alen);

    Eexps = E->exps;
    Ecoeffs = E->coeffs;
    Aterms = A->terms;
    Ei = 0;
    for (Ai = 0; Ai < Alen; Ai++)
    {
        FLINT_ASSERT(Ei < E->alloc);
        Eexps[Ei] = Aterms[Ai].exp;
        Ecoeffs[Ei] = n_poly_mod_eval_step(Aterms[Ai].coeff, mod);
        Ei += (Ecoeffs[Ei] != 0);
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
    const n_polyun_term_struct * At = A->terms;
    const n_polyun_term_struct * Ait;
    n_polyun_term_struct * Zit;
    slong Ai, ai, Zi, j;

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
            n_poly_fit_length(Zit->coeff, fit_length);
            Zit->coeff->length = cur_length;
            mpn_zero(Zit->coeff->coeffs, cur_length);
            goto in_both;            
        }
        else if (Ait->exp + ai < Zit->exp)
        {
            /* missing from A */
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            Zit->coeff->coeffs[cur_length] = 0;
            Zit->coeff->length = cur_length + 1;
            Zi++;
        }
        else
        {
in_both:
            FLINT_ASSERT(cur_length == Zit->coeff->length);
            FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
            Zit->coeff->coeffs[cur_length] = Ait->coeff->coeffs[ai];
            Zit->coeff->length = cur_length + 1;
            Zi++;
            do {
                ai--;
            } while (ai >= 0 && Ait->coeff->coeffs[ai] == 0);
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
        n_poly_fit_length(Zit->coeff, fit_length);
        Zit->coeff->length = cur_length;
        mpn_zero(Zit->coeff->coeffs, cur_length);
        Z->length = ++Zi;
        FLINT_ASSERT(cur_length == Zit->coeff->length);
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        Zit->coeff->coeffs[cur_length] = Ait->coeff->coeffs[ai];
        Zit->coeff->length = cur_length + 1;
        do {
            ai--;
        } while (ai >= 0 && Ait->coeff->coeffs[ai] == 0);
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
        FLINT_ASSERT(cur_length + 1 <= Zit->coeff->alloc);
        Zit->coeff->coeffs[cur_length] = 0;
        Zit->coeff->length = cur_length + 1;
        Zi++;
    }

    for (Zi = 0; Zi < Z->length; Zi++)
    {
        FLINT_ASSERT(Z->terms[Zi].coeff->length == cur_length + 1);
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
    nmod_mpoly_t T1, T2;
    ulong * Bdegs;
    const slong degs0 = degs[0];

    FLINT_ASSERT(m > 2);
    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(bits <= FLINT_BITS);

#if WANT_ASSERT
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

    beta = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
    Bdegs = (ulong *) flint_malloc(r*sizeof(ulong));
    H = (nmod_mpolyu_struct *) flint_malloc(r*sizeof(nmod_mpolyu_struct));
    Beh = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Beval = (n_polyu_struct *) flint_malloc(r*sizeof(n_polyu_struct));
    BBeval = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));
    Z = (n_polyun_struct *) flint_malloc(r*sizeof(n_polyun_struct));

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
    FLINT_ASSERT(ctx->ffinfo->mod.n > 3);
    for (i = 0; i < ctx->minfo->nvars; i++)
        beta[i] = n_urandint(state, ctx->ffinfo->mod.n - 3) + 2;

    nmod_mpoly_set_eval_helper3(Aeh, A, m, beta, ctx);

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        slong this_zip_images;
        this_zip_images = nmod_mpoly_set_eval_helper_and_zip_form3(Bdegs + i,
                                          Beh + i, H + i, B + i, beta, m, ctx);
        req_zip_images = FLINT_MAX(req_zip_images, this_zip_images);
        FLINT_ASSERT(Bdegs[i] > 0);
        Z[i].length = 0;
    }

    cur_zip_image = 0;

next_zip_image:

    n_polyu_mod_eval_step(Aeval, Aeh, ctx->ffinfo->mod);
    for (i = 0; i < r; i++)
        n_polyu_mod_eval_step(Beval + i, Beh + i, ctx->ffinfo->mod);

    success = n_polyu3_mod_hlift(r, BBeval, Aeval, Beval,
                                             alpha[m - 1], degs0, ctx->ffinfo);
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
    flint_free(Bdegs);
    flint_free(H);
    flint_free(Beh);
    flint_free(Beval);
    flint_free(BBeval);
    flint_free(Z);

    return success;
}

