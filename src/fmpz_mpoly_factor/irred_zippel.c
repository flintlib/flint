/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"
#include "nmod_mpoly_factor.h"


static void nmod_mpoly_get_eval_helper2(
    n_polyun_t EH,
    const nmod_mpoly_t A,
    n_poly_struct * caches,
    const nmod_mpoly_ctx_t ctx)
{
    slong start, Ai, j, k, n;
    slong e0, e1, EHi;
    mp_limb_t * p;
    flint_bitcnt_t bits = A->bits;
    slong Alen = A->length;
    const ulong * Aexps = A->exps;
    slong nvars = ctx->minfo->nvars;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong off0, off1, shift0, shift1;
    slong * off, * shift;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, bits, ctx->minfo);

    off = (slong *) TMP_ALLOC(2*nvars*sizeof(slong));
    shift = off + nvars;
    for (k = 2; k < nvars; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, bits, ctx->minfo);

    Ai = 0;
    EHi = 0;

    while (Ai < Alen)
    {
        start = Ai;
        e0 = (Aexps[N*Ai + off0] >> shift0) & mask;
        e1 = (Aexps[N*Ai + off1] >> shift1) & mask;
        while (1)
        {
            Ai++;
            if (Ai >= Alen)
                break;
            if (((Aexps[N*Ai + off0] >> shift0) & mask) != e0)
                break;
            if (((Aexps[N*Ai + off1] >> shift1) & mask) != e1)
                break;
        }

        n = Ai - start;

        n_polyun_fit_length(EH, EHi + 1);
        EH->exps[EHi] = pack_exp2(e0, e1);
        n_poly_fit_length(EH->coeffs + EHi, 2*n);
        EH->coeffs[EHi].length = n;
        p = EH->coeffs[EHi].coeffs;
        EHi++;

        for (j = 0; j < n; j++)
        {
            mp_limb_t meval = 1;

            for (k = 2; k < nvars; k++)
            {
                ulong ei = (Aexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                meval = nmod_pow_cache_mulpow_ui(meval, ei, caches + 3*k + 0,
                                 caches + 3*k + 1, caches + 3*k + 2, ctx->mod);
            }

            /* set cur = monomial eval */
            p[j] = meval;

            /* copy cur to inc */
            p[j + n] = meval;
        }
    }

    EH->length = EHi;

    TMP_END;
}

static slong nmod_mpoly_set_eval_helper_and_zip_form2(
    slong * deg1_, /* degree of B wrt main var 1 */
    n_polyun_t EH,
    n_polyun_t H,
    n_polyun_t M,
    const nmod_mpoly_t B,
    n_poly_struct * caches,
    const nmod_mpoly_ctx_t ctx)
{
    slong start, Bi, j, k, n;
    slong e0, e1, Hi, EHi;
    mp_limb_t * p;
    slong zip_length = 0;
    flint_bitcnt_t bits = B->bits;
    slong Blen = B->length;
    const ulong * Bexps = B->exps;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    slong N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong off0, off1, shift0, shift1;
    slong deg0, deg1 = -1;
    slong * off, * shift;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(Blen > 0);

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, bits, ctx->minfo);

    off = (slong *) TMP_ALLOC(2*ctx->minfo->nvars*sizeof(slong));
    shift = off + ctx->minfo->nvars;
    for (k = 2; k < ctx->minfo->nvars; k++)
        mpoly_gen_offset_shift_sp(&off[k], &shift[k], k, bits, ctx->minfo);

    Bi = 0;
    deg0 = (Bexps[N*Bi + off0] >> shift0) & mask;

    EHi = 0;
    Hi = 0;

    while (Bi < Blen)
    {
        start = Bi;
        e0 = (Bexps[N*Bi + off0] >> shift0) & mask;
        e1 = (Bexps[N*Bi + off1] >> shift1) & mask;
        deg1 = FLINT_MAX(deg1, e1);
        while (1)
        {
            Bi++;
            if (Bi >= Blen)
                break;
            if (((Bexps[N*Bi + off0] >> shift0) & mask) != e0)
                break;
            if (((Bexps[N*Bi + off1] >> shift1) & mask) != e1)
                break;
        }

        n = Bi - start;

        n_polyun_fit_length(EH, EHi + 1);
        EH->exps[EHi] = pack_exp2(e0, e1);
        n_poly_fit_length(EH->coeffs + EHi, 2*n);
        EH->coeffs[EHi].length = n;
        p = EH->coeffs[EHi].coeffs;
        EHi++;

        for (j = 0; j < n; j++)
        {
            mp_limb_t meval = 1;

            for (k = 2; k < ctx->minfo->nvars; k++)
            {
                ulong ei = (Bexps[N*(start + j) + off[k]] >> shift[k]) & mask;
                meval = nmod_pow_cache_mulpow_ui(meval, ei, caches + 3*k + 0,
                               caches + 3*k + 1, caches + 3*k + 2, ctx->mod);
            }

            /* set cur = monomial eval */
            p[j] = meval;

            /* copy cur to inc */
            p[j + n] = meval;
        }

        if (e0 < deg0)
        {
            n_polyun_fit_length(H, Hi + 1);
            n_polyun_fit_length(M, Hi + 1);
            H->exps[Hi] = pack_exp2(e0, e1);
            M->exps[Hi] = pack_exp2(e0, e1);
            n_poly_fit_length(H->coeffs + Hi, n);
            zip_length = FLINT_MAX(zip_length, n);
            H->coeffs[Hi].length = n;
            flint_mpn_copyi(H->coeffs[Hi].coeffs, p, n);
            n_poly_mod_product_roots_nmod_vec(M->coeffs + Hi, p, n, ctx->mod);
            Hi++;
        }
    }

    TMP_END;

    EH->length = EHi;
    H->length = Hi;
    M->length = Hi;

    *deg1_ = deg1;
    return zip_length;
}


static int _fmpz_mpoly_modpk_update_zip(
    fmpz_t pk,
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    n_polyun_t Z,
    const n_polyun_t H,
    const n_polyun_t M,
    const nmod_mpoly_ctx_t ctxp)
{
    slong i, j, Ai, n;
    int success;
    slong off, shift;
    slong N = mpoly_words_per_exp_sp(A->bits, ctx->minfo);
    ulong start, mask = (-UWORD(1)) >> (FLINT_BITS - A->bits);
    n_poly_t c, t;
    mp_limb_t * ccoeffs;

    mpoly_gen_offset_shift_sp(&off, &shift, 0, A->bits, ctx->minfo);

    mask = mask << shift;

    Ai = 1;
    start = (A->exps + N*0)[off] & mask;

    while (Ai < A->length && start == ((A->exps + N*Ai)[off] & mask))
    {
        Ai++;
    }

    FLINT_ASSERT(Ai < A->length);

    if (Ai >= A->length)
        return 1;

    n_poly_init(c);
    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    for (i = 0; i < Z->length; i++)
    {
        n = H->coeffs[i].length;
        FLINT_ASSERT(M->coeffs[i].length == n + 1);
        FLINT_ASSERT(Z->coeffs[i].length >= n);

        n_poly_fit_length(c, n);
        n_poly_fit_length(t, n);

        ccoeffs = c->coeffs;

        success = _nmod_zip_vand_solve(c->coeffs,
                              H->coeffs[i].coeffs, n,
                              Z->coeffs[i].coeffs, Z->coeffs[i].length,
                              M->coeffs[i].coeffs, t->coeffs, ctxp->mod);
        if (success < 1)
        {
            n_poly_clear(t);
            n_poly_clear(c);
            return success;
        }

        FLINT_ASSERT(Ai + n <= A->length);

        for (j = 0; j < n; j++)
        {
            if (ctxp->mod.n - ccoeffs[j] < ccoeffs[j])
                fmpz_submul_ui(A->coeffs + Ai + j, pk, ctxp->mod.n - ccoeffs[j]);
            else
                fmpz_addmul_ui(A->coeffs + Ai + j, pk, ccoeffs[j]);
        }

        Ai += n;
    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    n_poly_clear(c);

    return 1;
}


static void _nmod_mpoly_set_fmpz_mpoly(
    nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
    slong Ap_len, i;
    FLINT_ASSERT(ctx->minfo->nvars == ctx->minfo->nvars);
    FLINT_ASSERT(ctx->minfo->ord == ctx->minfo->ord);
    nmod_mpoly_fit_length_reset_bits(Ap, A->length, A->bits, ctxp);
    Ap_len = 0;
    for (i = 0; i < A->length; i++)
    {
        Ap->coeffs[Ap_len] = fmpz_fdiv_ui(A->coeffs + i, ctxp->mod.n);
        if (Ap->coeffs[Ap_len] == 0)
            continue;
        mpoly_monomial_set(Ap->exps + N*Ap_len, A->exps + N*i, N);
        Ap_len++;
    }
    Ap->length = Ap_len;
}


static void _fmpz_mpoly_modpk_taylor_coeff(
    const fmpz_t pk,
    nmod_mpoly_t T,
    const nmod_mpoly_ctx_t ctxp,
    const fmpz_mpoly_t E,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i, Tlen;
    slong N = mpoly_words_per_exp(E->bits, ctx->minfo);
    fmpz_t t;

    fmpz_init(t);

    nmod_mpoly_fit_length_reset_bits(T, E->length, E->bits, ctxp);
    Tlen = 0;
    for (i = 0; i < E->length; i++)
    {
        FLINT_ASSERT(fmpz_divisible(E->coeffs + i, pk)); /* TODO !!! */
        fmpz_divexact(t, E->coeffs + i, pk);
        T->coeffs[Tlen] = fmpz_fdiv_ui(t, ctxp->mod.n);
        if (T->coeffs[Tlen] == 0)
            continue;
        mpoly_monomial_set(T->exps + N*Tlen, E->exps + N*i, N);
        Tlen++;
    }
    T->length = Tlen;

    fmpz_clear(t);
}


static void _fmpz_mpoly_set_nmod_mpoly_smod(
    fmpz_mpoly_t A,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_t Ap,
    const nmod_mpoly_ctx_t ctxp)
{
    slong N = mpoly_words_per_exp(Ap->bits, ctxp->minfo);

    FLINT_ASSERT(ctx->minfo->ord == ctxp->minfo->ord);
    FLINT_ASSERT(ctx->minfo->nvars == ctxp->minfo->nvars);

    A->length = 0;
    fmpz_mpoly_fit_bits(A, Ap->bits, ctx);
    A->bits = Ap->bits;

    fmpz_mpoly_fit_length(A, Ap->length, ctx);
    A->length = Ap->length;

    mpoly_copy_monomials(A->exps, Ap->exps, Ap->length, N);
    _fmpz_vec_set_nmod_vec(A->coeffs, Ap->coeffs, Ap->length, ctxp->mod);
}


static void n_poly_eval_reset(n_poly_t A)
{
    _nmod_vec_set(A->coeffs, A->coeffs + A->length, A->length);
}


static void n_polyun_eval_reset(n_polyun_t A)
{
    slong Ai;
    for (Ai = 0; Ai < A->length; Ai++)
        n_poly_eval_reset(A->coeffs + Ai);
}


static void n_bpoly_mod_eval_step(
    n_bpoly_t E,
    n_polyun_t EH,
    const nmod_mpoly_t A,
    nmod_t ctx)
{
    slong i, n, Ai;
    mp_limb_t * p;
    mp_limb_t c;
    ulong e0, e1;
    slong EHlen = EH->length;

    Ai = 0;

    n_bpoly_zero(E);
    for (i = 0; i < EHlen; i++)
    {
        n = EH->coeffs[i].length;
        p = EH->coeffs[i].coeffs;
        FLINT_ASSERT(EH->coeffs[i].alloc >= 2*n);
        FLINT_ASSERT(Ai + n <= A->length);

        c = _nmod_zip_eval_step(p + 0*n, p + 1*n, A->coeffs + Ai, n, ctx);
        Ai += n;
        e0 = extract_exp(EH->exps[i], 1, 2);
        e1 = extract_exp(EH->exps[i], 0, 2);
        if (c == 0)
            continue;
        n_bpoly_set_coeff_nonzero(E, e0, e1, c);
    }

    FLINT_ASSERT(Ai == A->length);
}


/*
    n_polyun_t Beh has all x0^i*x1^j*poly(x2, ...) with
        coeffs triples suitable for sequential eval at x2,... = beta^i

    n_polyun_t H has x0^i*x1^j*poly(x2, ...) with i < deg_x0(B[i]) with
        the monomial evals at beta in the coeffs

    n_polyun_t M has x0^i*x1^j*poly(x2, ...) with i < deg_x0(B[i])
        the master poly of H[x0^i*x1^j] in the coeff

    n_polyun_t Z has x0^i*x1^j with i < deg_x0(B[i])
        ready to collect images
*/
static int fmpz_mfactor_lift_prime_power_zippel(
    slong r,
    fmpz_mpoly_struct * B,
    flint_rand_t state,
    const nmod_mpoly_struct * Bp,
    const fmpz_mpoly_t A,
    const mp_limb_t * alphap,
    const fmpz_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t ctxp,
    slong L)
{
    slong req_zip_images, cur_zip_image, this_images;
    slong n = ctxp->minfo->nvars;
    int success;
    slong i, j, k;
    n_polyun_struct * H, * M, * Z, * Beh;
    n_bpoly_struct * Ceval, * Beval, Teval[1];
    slong * Cdegs1;
    nmod_mpoly_t Tp;
    n_polyun_t Teh;
    fmpz_t pk;
    n_poly_struct * beta_caches;
    fmpz_mpoly_t e, t1, t2;

    FLINT_ASSERT(r > 1);
    FLINT_ASSERT(ctxp->mod.n > 3);

    fmpz_init(pk);

    fmpz_mpoly_init(e, ctx);
    fmpz_mpoly_init(t1, ctx);
    fmpz_mpoly_init(t2, ctx);

    nmod_mpoly_init(Tp, ctxp);
    n_polyun_init(Teh);

    Cdegs1 = FLINT_ARRAY_ALLOC(r, slong);

    beta_caches = FLINT_ARRAY_ALLOC(3*n, n_poly_struct);
    for (i = 0; i < 3*n; i++)
        n_poly_init(beta_caches + i);

    Beh = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    H = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    M = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    Z = FLINT_ARRAY_ALLOC(r, n_polyun_struct);
    Ceval = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    Beval = FLINT_ARRAY_ALLOC(r, n_bpoly_struct);
    for (i = 0; i < r; i++)
    {
        n_polyun_init(Beh + i);
        n_polyun_init(H + i);
        n_polyun_init(M + i);
        n_polyun_init(Z + i);
        n_bpoly_init(Ceval + i);
        n_bpoly_init(Beval + i);
    }

    n_bpoly_init(Teval);

    /* init done */

    /* choose betas */
    for (i = 2; i < n; i++)
    {
        mp_limb_t bb = n_urandint(state, ctxp->mod.n - 3) + 2;
        nmod_pow_cache_start(bb, beta_caches + 3*i + 0,
                                 beta_caches + 3*i + 1, beta_caches + 3*i + 2);
    }

    req_zip_images = 1;
    for (i = 0; i < r; i++)
    {
        if (Bp[i].bits > FLINT_BITS)
        {
            success = 0;
            goto cleanup;
        }

        this_images = nmod_mpoly_set_eval_helper_and_zip_form2(Cdegs1 + i,
                             Beh + i, H + i, M + i, Bp + i, beta_caches, ctxp);
        req_zip_images = FLINT_MAX(req_zip_images, this_images);
    }

    for (i = 0; i < r; i++)
    {
        n_polyun_fit_length(Z + i, H[i].length);
        Z[i].length = H[i].length;
        for (j = 0; j < H[i].length; j++)
        {
            Z[i].exps[j] = H[i].exps[j];
            n_poly_fit_length(Z[i].coeffs + j, req_zip_images);
            Z[i].coeffs[j].length = 0;
        }
    }

    fmpz_one(pk);

    k = 1;

next_power:

    fmpz_mul_ui(pk, pk, ctxp->mod.n);

    fmpz_mpoly_mul(t1, B + 0, B + 1, ctx);
    for (i = 2; i < r; i++)
    {
        fmpz_mpoly_mul(t2, t1, B + i, ctx);
        fmpz_mpoly_swap(t1, t2, ctx);
    }
    fmpz_mpoly_sub(e, A, t1, ctx);

    if (fmpz_mpoly_is_zero(e, ctx))
    {
        success = 1;
        goto cleanup;
    }

    if (k > L || e->bits > FLINT_BITS)
    {
        success = 0;
        goto cleanup;
    }

    _fmpz_mpoly_modpk_taylor_coeff(pk, Tp, ctxp, e, ctx);
    nmod_mpoly_get_eval_helper2(Teh, Tp, beta_caches, ctxp);

    if (fmpz_cmp_ui(pk, ctxp->mod.n) > 0)
    {
        for (i = 0; i < r; i++)
            n_polyun_eval_reset(Beh + i);
    }

    cur_zip_image = 0;

next_zip_image:

    n_bpoly_mod_eval_step(Teval, Teh, Tp, ctxp->mod);
    for (i = 0; i < r; i++)
        n_bpoly_mod_eval_step(Beval + i, Beh + i, Bp + i, ctxp->mod);

    success = n_bpoly_mod_pfrac(r, Ceval, Cdegs1, Teval, Beval, ctxp->mod);
    if (success < 1)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 0; i < r; i++)
    {
        success = n_polyu2n_add_zip_must_match(Z + i, Ceval + i, cur_zip_image);
        if (!success)
        {
            success = 0;
            goto cleanup;
        }
    }

    cur_zip_image++;
    if (cur_zip_image < req_zip_images)
        goto next_zip_image;

    for (i = 0; i < r; i++)
    {
        success = _fmpz_mpoly_modpk_update_zip(pk, B + i, ctx,
                                                    Z + i, H + i, M + i, ctxp);
        if (success < 1)
        {
            success = 0;
            goto cleanup;
        }
    }

    goto next_power;

cleanup:

    fmpz_clear(pk);

    fmpz_mpoly_clear(e, ctx);
    fmpz_mpoly_clear(t1, ctx);
    fmpz_mpoly_clear(t2, ctx);

    nmod_mpoly_clear(Tp, ctxp);
    n_polyun_clear(Teh);

    flint_free(Cdegs1);

    for (i = 0; i < 3*n; i++)
        n_poly_clear(beta_caches + i);
    flint_free(beta_caches);

    for (i = 0; i < r; i++)
    {
        n_polyun_clear(Beh + i);
        n_polyun_clear(H + i);
        n_polyun_clear(M + i);
        n_polyun_clear(Z + i);
        n_bpoly_clear(Ceval + i);
        n_bpoly_clear(Beval + i);
    }
    flint_free(Beh);
    flint_free(H);
    flint_free(M);
    flint_free(Z);
    flint_free(Ceval);
    flint_free(Beval);

    n_bpoly_clear(Teval);

    FLINT_ASSERT(success == 0 || success == 1);
    return success;
}


void nmod_poly_set_fmpz_poly(nmod_poly_t a, const fmpz_poly_t b)
{
    slong i;
    nmod_poly_fit_length(a, b->length);
    for (i = 0; i < b->length; i++)
        a->coeffs[i] = fmpz_fdiv_ui(b->coeffs + i, a->mod.n);
    a->length = b->length;
    _nmod_poly_normalise(a);
}

/*
    return 1: success
           0: failed
          -1: exception large exponents
*/
int fmpz_mpoly_factor_irred_zippel(
    fmpz_mpolyv_t fac,
    const fmpz_mpoly_t A,
    const fmpz_mpoly_factor_t lcAfac,
    int lcAfac_irred,
    const fmpz_mpoly_t lcA,
    const fmpz_mpoly_ctx_t ctx,
    flint_rand_t state,
    zassenhaus_prune_t zas)
{
    int success, kfails_left = 4;
    const slong n = ctx->minfo->nvars - 1;
    slong i, j, k;
    fmpz * alpha;
    fmpz_mpoly_struct * Aevals;
    slong * degs, * tdegs;
    fmpz_mpolyv_t tfac;
    fmpz_mpoly_t t, Acopy;
    fmpz_mpoly_struct * newA;
    fmpz_poly_t Au;
    fmpz_poly_factor_t Aufac;
    slong alpha_bits, alpha_count;
    fmpz_mpoly_t m, mpow;
    fmpz_mpolyv_t Alc, lc_divs;
    fmpz_t q, facBound;
    mp_limb_t p;
    nmod_mpoly_ctx_t ctxp;
    nmod_mpolyv_t facp, tfacp;
    nmod_mpolyv_t Aevalp, Alcp;
    nmod_poly_t Aup;
    mp_limb_t * alphap;
    slong r, L;

    FLINT_ASSERT(n > 1);
    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(fmpz_mpoly_factor_matches(lcA, lcAfac, ctx));

    fmpz_init(facBound);
    fmpz_init(q);

    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(m, ctx);
    fmpz_mpoly_init(mpow, ctx);

    fmpz_mpolyv_init(Alc, ctx);
    fmpz_mpolyv_init(lc_divs, ctx);

    fmpz_poly_factor_init(Aufac);
    fmpz_poly_init(Au);

    alpha = _fmpz_vec_init(n);
    alphap = (mp_limb_t *) flint_malloc(n*sizeof(mp_limb_t));

    degs  = (slong *) flint_malloc(2*(n + 1)*sizeof(slong));
    tdegs = degs + (n + 1);

    Aevals = (fmpz_mpoly_struct *) flint_malloc(n*sizeof(fmpz_mpoly_struct));
    for (i = 0; i < n; i++)
        fmpz_mpoly_init(Aevals + i, ctx);

    fmpz_mpolyv_init(tfac, ctx);
    fmpz_mpoly_init(t, ctx);

    nmod_mpoly_ctx_init(ctxp, n + 1, ORD_LEX, 2); /* modulus no care */
    nmod_mpolyv_init(facp, ctxp);
    nmod_mpolyv_init(tfacp, ctxp);
    nmod_mpolyv_init(Aevalp, ctxp);
    nmod_mpolyv_init(Alcp, ctxp);
    nmod_poly_init_mod(Aup, ctxp->mod);

    /* init done */

    fmpz_mpoly_degrees_si(degs, A, ctx);

    alpha_count = 0;
    alpha_bits = 10;
    p = UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX);

next_alpha:

    alpha_count++;
    if (alpha_count >= alpha_bits)
    {
        alpha_count = 0;
        alpha_bits++;
        if (alpha_bits > FLINT_BITS/2)
        {
            success = 0;
            goto cleanup;
        }
    }

    for (i = 0; i < n; i++)
    {
        ulong l = n_randlimb(state);
        ulong mask = UWORD(1) << alpha_bits;
        if (l & mask)
            fmpz_neg_ui(alpha + i, 1 + (l & (mask - 1)));
        else
            fmpz_set_ui(alpha + i, 1 + (l & (mask - 1)));
    }

#if FLINT_WANT_ASSERT
    fmpz_mpoly_degrees_si(tdegs, A, ctx);
    for (i = 0; i < n + 1; i++)
        FLINT_ASSERT(degs[i] == tdegs[i]);
#endif

    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(Aevals + i,
                       i == n - 1 ? A : Aevals + i + 1, i + 1, alpha + i, ctx);
        fmpz_mpoly_degrees_si(tdegs, Aevals + i, ctx);
        for (j = 0; j <= i; j++)
        {
            if (tdegs[j] != degs[j])
                goto next_alpha;
        }
    }

    success = fmpz_mpoly_get_fmpz_poly(Au, Aevals + 0, 0, ctx);
    FLINT_ASSERT(success);
    fmpz_poly_factor(Aufac, Au);
    r = Aufac->num;

    zassenhaus_prune_start_add_factors(zas);
    for (j = 0; j < r; j++)
        zassenhaus_prune_add_factor(zas, fmpz_poly_degree(Aufac->p + j),
                                                                Aufac->exp[j]);
    zassenhaus_prune_end_add_factors(zas);

    if ((r < 2 && Aufac->exp[0] == 1) ||
        zassenhaus_prune_must_be_irreducible(zas))
    {
        fmpz_mpolyv_fit_length(fac, 1, ctx);
        fac->length = 1;
        fmpz_mpoly_set(fac->coeffs + 0, A, ctx);
        success = 1;
        goto cleanup;
    }

    for (j = 0; j < r; j++)
    {
        if (Aufac->exp[j] != 1)
            goto next_alpha;
    }

    fmpz_mpolyv_fit_length(lc_divs, r, ctx);
    lc_divs->length = r;
    if (lcAfac->num > 0)
    {
        success = 0;
        if (lcAfac_irred)
            success = fmpz_mpoly_factor_lcc_wang(lc_divs->coeffs, lcAfac,
                                           &Aufac->c, Aufac->p, r, alpha, ctx);
        if (!success)
        {
            success = fmpz_mpoly_factor_lcc_kaltofen(lc_divs->coeffs, lcAfac,
                                                A, r, alpha, degs, Aufac, ctx);
            if (success < 0 || (success == 0 && --kfails_left >= 0))
                goto next_alpha;
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            FLINT_ASSERT(Aufac->p[i].length > 0);
            fmpz_mpoly_set_fmpz(lc_divs->coeffs + i,
                             Aufac->p[i].coeffs + Aufac->p[i].length - 1, ctx);
        }
    }

    FLINT_ASSERT(r > 1);

    if (!fmpz_mpoly_divides(m, lcA, lc_divs->coeffs + 0, ctx))
        goto next_alpha;
    for (i = 1; i < r; i++)
    {
        if (!fmpz_mpoly_divides(m, m, lc_divs->coeffs + i, ctx))
            goto next_alpha;
    }

    fmpz_mpoly_pow_ui(mpow, m, r - 1, ctx);
    if (fmpz_mpoly_is_one(mpow, ctx))
    {
        newA = (fmpz_mpoly_struct *) A;
    }
    else
    {
        newA = Acopy;
        fmpz_mpoly_mul(newA, A, mpow, ctx);
    }

    if (newA->bits > FLINT_BITS)
    {
        success = -1;
        goto cleanup;
    }

    fmpz_mpoly_degrees_si(tdegs, newA, ctx);

    for (i = 0; i < n + 1; i++)
    {
        if (FLINT_BIT_COUNT(tdegs[i]) >= FLINT_BITS/3)
        {
            success = -1;
            goto cleanup;
        }
    }

    _fmpz_vec_height(facBound, newA->coeffs, newA->length);
    if (!fmpz_mpoly_factor_bound_si(facBound, facBound, tdegs, n + 1))
    {
        success = -1;
        goto cleanup;
    }
    fmpz_mul_2exp(facBound, facBound, 1);
    fmpz_add_ui(facBound, facBound, 1);

    fmpz_mpoly_set(t, mpow, ctx);
    for (i = n - 1; i >= 0; i--)
    {
        fmpz_mpoly_evaluate_one_fmpz(t, mpow, i + 1, alpha + i, ctx);
        fmpz_mpoly_swap(t, mpow, ctx);
        fmpz_mpoly_mul(Aevals + i, Aevals + i, mpow, ctx);
    }

    fmpz_mpolyv_fit_length(Alc, (n + 1)*r, ctx);

    i = n;
    for (j = 0; j < r; j++)
    {
        fmpz_mpoly_mul(Alc->coeffs + i*r + j, lc_divs->coeffs + j, m, ctx);
    }
    for (i = n - 1; i >= 0; i--)
    {
        for (j = 0; j < r; j++)
        {
            fmpz_mpoly_evaluate_one_fmpz(Alc->coeffs + i*r + j,
                 Alc->coeffs + (i + 1)*r + j, i + 1, alpha + i, ctx);
        }
    }

    fmpz_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        FLINT_ASSERT(fmpz_mpoly_is_fmpz(Alc->coeffs + 0*r + i, ctx));
        FLINT_ASSERT(fmpz_mpoly_length(Alc->coeffs + 0*r + i, ctx) == 1);
        FLINT_ASSERT(fmpz_divisible(Alc->coeffs[i].coeffs + 0,
                                 Aufac->p[i].coeffs + Aufac->p[i].length - 1));
        fmpz_divexact(q, Alc->coeffs[i].coeffs + 0,
                                  Aufac->p[i].coeffs + Aufac->p[i].length - 1);
        _fmpz_mpoly_set_fmpz_poly(fac->coeffs + i, newA->bits,
                               Aufac->p[i].coeffs, Aufac->p[i].length, 0, ctx);
        fmpz_mpoly_scalar_mul_fmpz(fac->coeffs + i, fac->coeffs + i, q, ctx);
    }

next_zip_prime:

    if (p >= UWORD_MAX_PRIME)
    {
        success = 0;
        goto cleanup;
    }

    p = n_nextprime(p, 1);
    if ((p % 4) != 1 && p < UWORD_MAX_PRIME)
    {
        p = n_nextprime(p, 1);
        if ((p % 4) != 1 && p < UWORD_MAX_PRIME)
        {
            p = n_nextprime(p, 1);
        }
    }

    nmod_mpoly_ctx_set_modulus(ctxp, p);
    nmod_mpolyv_fit_length(facp, r, ctxp);
    nmod_mpolyv_fit_length(tfacp, r, ctxp);
    facp->length = r;
    tfacp->length = r;

    nmod_mpolyv_fit_length(Aevalp, n + 1, ctxp);
    nmod_mpolyv_fit_length(Alcp, (n + 1)*r, ctxp);

    for (i = 0; i < r; i++)
    {
        _nmod_mpoly_set_fmpz_mpoly(facp->coeffs + i, ctxp, fac->coeffs + i, ctx);
        if (facp->coeffs[i].length != fac->coeffs[i].length)
            goto next_zip_prime;
    }

    Aup->mod = ctxp->mod;
    nmod_poly_set_fmpz_poly(Aup, Au);
    if (Aup->length != Au->length || !nmod_poly_is_squarefree(Aup))
        goto next_zip_prime;

    for (k = 0; k <= n; k++)
    {
        _nmod_mpoly_set_fmpz_mpoly(Aevalp->coeffs + k, ctxp,
                                               k < n ? Aevals + k : newA, ctx);
        if (Aevalp->coeffs[k].length != (k < n ? Aevals + k : newA)->length)
            goto next_zip_prime;

        for (i = 0; i < r; i++)
        {
            _nmod_mpoly_set_fmpz_mpoly(Alcp->coeffs + k*r + i, ctxp,
                                                   Alc->coeffs + k*r + i, ctx);
            if (Alcp->coeffs[k*r + i].length != Alc->coeffs[k*r + i].length)
                goto next_zip_prime;
        }
    }

    for (i = 0; i < n; i++)
    {
        alphap[i] = fmpz_fdiv_ui(alpha + i, ctxp->mod.n);
        if (alphap[i] == 0)
            goto next_zip_prime;
    }

    for (k = 1; k <= n; k++)
    {
        for (i = 0; i < r; i++)
        {
            _nmod_mpoly_set_lead0(tfacp->coeffs + i, facp->coeffs + i,
                                                 Alcp->coeffs + k*r + i, ctxp);
        }

        if (k > 2)
        {
            success = nmod_mpoly_hlift_zippel(k, tfacp->coeffs, r, alphap,
                                       Aevalp->coeffs + k, tdegs, ctxp, state);
        }
        else
        {
            success = nmod_mpoly_hlift(k, tfacp->coeffs, r, alphap,
                                              Aevalp->coeffs + k, tdegs, ctxp);
        }

        if (!success)
            goto next_alpha;

        nmod_mpolyv_swap(tfacp, facp, ctxp);
    }

    fmpz_mpolyv_fit_length(fac, r, ctx);
    fac->length = r;
    for (i = 0; i < r; i++)
    {
        _fmpz_mpoly_set_nmod_mpoly_smod(fac->coeffs + i, ctx,
                                                       facp->coeffs + i, ctxp);
        _fmpz_mpoly_set_lead0(fac->coeffs + i, fac->coeffs + i,
                                                   Alc->coeffs + n*r + i, ctx);
    }

    L = fmpz_clog_ui(facBound, ctxp->mod.n);
    success = fmpz_mfactor_lift_prime_power_zippel(r, fac->coeffs, state,
                                     facp->coeffs, newA, alphap, ctx, ctxp, L);

    if (!success)
        goto next_alpha;

    if (fmpz_mpoly_is_fmpz(m, ctx))
    {
        for (i = 0; i < r; i++)
        {
            _fmpz_vec_content(q, fac->coeffs[i].coeffs, fac->coeffs[i].length);
            if (fmpz_sgn(fac->coeffs[i].coeffs + 0) < 0)
                fmpz_neg(q, q);
            fmpz_mpoly_scalar_divexact_fmpz(fac->coeffs + i,
                                            fac->coeffs + i, q, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            /* hlift should not have returned any large bits */
            FLINT_ASSERT(fac->coeffs[i].bits <= FLINT_BITS);

            if (!fmpz_mpolyl_content(t, fac->coeffs + i, 1, ctx))
            {
                success = -1;
                goto cleanup;
            }

            success = fmpz_mpoly_divides(fac->coeffs + i, fac->coeffs + i, t, ctx);
            FLINT_ASSERT(success);
            fmpz_mpoly_unit_normalize(fac->coeffs, ctx);
        }
    }

    success = 1;

cleanup:

    fmpz_clear(facBound);
    fmpz_clear(q);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(m, ctx);
    fmpz_mpoly_clear(mpow, ctx);

    fmpz_mpolyv_clear(Alc, ctx);
    fmpz_mpolyv_clear(lc_divs, ctx);

    fmpz_poly_factor_clear(Aufac);
    fmpz_poly_clear(Au);

    _fmpz_vec_clear(alpha, n);
    flint_free(alphap);
    flint_free(degs); /* and tdegs */

    for (i = 0; i < n; i++)
        fmpz_mpoly_clear(Aevals + i, ctx);
    flint_free(Aevals);

    fmpz_mpolyv_clear(tfac, ctx);
    fmpz_mpoly_clear(t, ctx);

    nmod_mpolyv_clear(facp, ctxp);
    nmod_mpolyv_clear(tfacp, ctxp);
    nmod_mpolyv_clear(Aevalp, ctxp);
    nmod_mpolyv_clear(Alcp, ctxp);
    nmod_poly_clear(Aup);
    nmod_mpoly_ctx_clear(ctxp);

    return success;
}

