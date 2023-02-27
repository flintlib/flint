/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"


void _fq_nmod_mpoly_monomial_evals_cache(
    n_fq_poly_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, Ai;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * off, * shift;
    n_poly_struct * caches;
    mp_limb_t * c;
    slong num = stop - start;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(num > 0);

    caches = FLINT_ARRAY_ALLOC(3*num, n_fq_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (i = 0; i < num; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i + start, Abits, ctx->minfo);
        n_fq_poly_init(caches + 3*i + 0);
        n_fq_poly_init(caches + 3*i + 1);
        n_fq_poly_init(caches + 3*i + 2);
        n_fq_pow_cache_start_fq_nmod(betas + i, caches + 3*i + 0,
                                                caches + 3*i + 1,
                                                caches + 3*i + 2, ctx->fqctx);
    }

    n_poly_fit_length(E, d*Alen);
    E->length = Alen;

    for (Ai = 0; Ai < Alen; Ai++)
    {
        c = E->coeffs + d*Ai;
        _n_fq_one(c, d);
        for (i = 0; i < num; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            n_fq_pow_cache_mulpow_ui(c, c, ei, caches + 3*i + 0,
                                               caches + 3*i + 1,
                                               caches + 3*i + 2, ctx->fqctx);
        }
    }

    for (i = 0; i < num; i++)
    {
        n_fq_poly_clear(caches + 3*i + 0);
        n_fq_poly_clear(caches + 3*i + 1);
        n_fq_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);
}



void _fq_nmod_mpoly_monomial_evals2_cache(
    n_fq_polyun_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const fq_nmod_struct * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i, Ai, Ei;
    ulong e0, e1, e01;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, ctx->minfo);
    slong * off, * shift;
    n_fq_poly_struct * caches;
    mp_limb_t * c;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(m > 2);

    caches = FLINT_ARRAY_ALLOC(3*(m - 2), n_fq_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (i = 0; i < m; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, Abits, ctx->minfo);
        if (i >= 2)
        {
            n_fq_poly_init(caches + 3*(i - 2) + 0);
            n_fq_poly_init(caches + 3*(i - 2) + 1);
            n_fq_poly_init(caches + 3*(i - 2) + 2);
            n_fq_pow_cache_start_fq_nmod(betas + i - 2, caches + 3*(i - 2) + 0,
                                           caches + 3*(i - 2) + 1,
                                           caches + 3*(i - 2) + 2, ctx->fqctx);
        }
    }

    Ai = 0;
    Ei = 0;

    e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
    e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
    e01 = pack_exp2(e0, e1);
    n_polyun_fit_length(E, Ei + 1);
    E->exps[Ei] = e01;
    n_poly_fit_length(E->coeffs + Ei, d*1);
    c = E->coeffs[Ei].coeffs + d*0;
    E->coeffs[Ei].length = 1;
    Ei++;
    _n_fq_one(c, d);
    for (i = 2; i < m; i++)
    {
        ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
        n_fq_pow_cache_mulpow_ui(c, c, ei, caches + 3*(i - 2) + 0,
                                           caches + 3*(i - 2) + 1,
                                           caches + 3*(i - 2) + 2, ctx->fqctx);
    }

    for (Ai++; Ai < Alen; Ai++)
    {
        e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
        e01 = pack_exp2(e0, e1);
        if (e01 == E->exps[Ei - 1])
        {
            slong len = E->coeffs[Ei - 1].length;
            n_poly_fit_length(E->coeffs + Ei - 1, d*(len + 1));
            c = E->coeffs[Ei - 1].coeffs + d*len;
            E->coeffs[Ei - 1].length = len + 1;
        }
        else
        {
            n_polyun_fit_length(E, Ei + 1);
            E->exps[Ei] = e01;
            n_poly_fit_length(E->coeffs + Ei, d*1);
            c = E->coeffs[Ei].coeffs + d*0;
            E->coeffs[Ei].length = 1;
            Ei++;
        }

        _n_fq_one(c, d);
        for (i = 2; i < m; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            n_fq_pow_cache_mulpow_ui(c, c, ei, caches + 3*(i - 2) + 0,
                                               caches + 3*(i - 2) + 1,
                                               caches + 3*(i - 2) + 2, ctx->fqctx);
        }
    }

    E->length = Ei;

    for (i = 0; i < m - 2; i++)
    {
        n_fq_poly_clear(caches + 3*i + 0);
        n_fq_poly_clear(caches + 3*i + 1);
        n_fq_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);

#if FLINT_WANT_ASSERT
    Ai = 0;
    for (i = 0; i < E->length; i++)
        Ai += E->coeffs[i].length;
    FLINT_ASSERT(Ai == Alen);
#endif
}



void n_fq_bpoly_eval_step_sep(
    n_fq_bpoly_t E,
    n_fq_polyun_t cur,
    const n_fq_polyun_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Ai;
    slong e0, e1;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_bpoly_zero(E);

    Ai = 0;
    for (i = 0; i < cur->length; i++)
    {
        slong this_len = cur->coeffs[i].length;

        _n_fq_zip_eval_step(c, cur->coeffs[i].coeffs, inc->coeffs[i].coeffs,
                                  A->coeffs + d*Ai, this_len, ctx);

        Ai += this_len;

        e0 = extract_exp(cur->exps[i], 1, 2);
        e1 = extract_exp(cur->exps[i], 0, 2);
        if (_n_fq_is_zero(c, d))
            continue;

        n_fq_bpoly_set_coeff_n_fq(E, e0, e1, c, ctx);
    }

    FLINT_ASSERT(Ai == A->length);

    flint_free(c);
}

static void n_fq_poly_eval_step_sep(
    mp_limb_t * res,
    n_fq_poly_t cur,
    const n_fq_poly_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    _n_fq_zip_eval_step(res, cur->coeffs, inc->coeffs, A->coeffs,
                                                               A->length, ctx);
}

static void n_fq_bpoly_evalp_step_sep(
    n_fq_bpoly_t E,
    n_polyun_t cur,
    const n_polyun_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Ai;
    slong e0, e1;
    mp_limb_t * c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_bpoly_zero(E);

    Ai = 0;
    for (i = 0; i < cur->length; i++)
    {
        slong this_len = cur->coeffs[i].length;

        _n_fqp_zip_eval_step(c, cur->coeffs[i].coeffs, inc->coeffs[i].coeffs,
                                  A->coeffs + d*Ai, this_len, d, ctx->mod);

        Ai += this_len;

        e0 = extract_exp(cur->exps[i], 1, 2);
        e1 = extract_exp(cur->exps[i], 0, 2);
        if (_n_fq_is_zero(c, d))
            continue;

        n_fq_bpoly_set_coeff_n_fq(E, e0, e1, c, ctx);
    }

    FLINT_ASSERT(Ai == A->length);

    flint_free(c);
}


static void n_fq_poly_evalp_step_sep(
    mp_limb_t * res,
    n_poly_t cur,
    const n_poly_t inc,
    const fq_nmod_mpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    _n_fqp_zip_eval_step(res, cur->coeffs, inc->coeffs, A->coeffs,
                                 A->length, fq_nmod_ctx_degree(ctx), ctx->mod);
}

static ulong _fq_nmod_mpoly_bidegree(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}

static void fq_nmod_mpoly_monomial_evals2(
    n_fq_polyun_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    _fq_nmod_mpoly_monomial_evals2_cache(E, A->exps, A->bits, A->length, betas, m, ctx);
}

static void fq_nmod_mpoly_monomial_evalsp2(
    n_polyun_t E,
    const fq_nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong m,
    const fq_nmod_mpoly_ctx_t ctx)
{
    _nmod_mpoly_monomial_evals2_cache(E, A->exps, A->bits, A->length, betas, m,
                                                  ctx->minfo, ctx->fqctx->mod);
}

static void fq_nmod_mpoly_monomial_evalsp(
    n_poly_t E,
    const fq_nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    _nmod_mpoly_monomial_evals_cache(E, A->exps, A->bits, A->length, betas,
                                     start, stop, ctx->minfo, ctx->fqctx->mod);
}

static void fq_nmod_mpoly_monomial_evals(
    n_fq_poly_t E,
    const fq_nmod_mpoly_t A,
    const fq_nmod_struct * betas,
    slong start,
    slong stop,
    const fq_nmod_mpoly_ctx_t ctx)
{
    _fq_nmod_mpoly_monomial_evals_cache(E, A->exps, A->bits, A->length,
                                                      betas, start, stop, ctx);
}


int n_fq_polyun_zip_solve(
    fq_nmod_mpoly_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong Ai, i, n;
    n_poly_t t;

    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    /* A has the right length, but possibly from a smaller ctx */
    if (A->length*d > A->coeffs_alloc)
    {
        slong new_alloc = FLINT_MAX(A->coeffs_alloc + A->coeffs_alloc/2, A->length*d);
        A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, new_alloc*sizeof(mp_limb_t));
        A->coeffs_alloc = new_alloc;
    }

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->coeffs[i].length;
        FLINT_ASSERT(M->coeffs[i].length == n + 1);
        FLINT_ASSERT(Z->coeffs[i].length >= n);
        FLINT_ASSERT(Ai + n <= A->length);

        n_poly_fit_length(t, d*n);

        success = _n_fq_zip_vand_solve(A->coeffs + d*Ai,
                         H->coeffs[i].coeffs, n,
                         Z->coeffs[i].coeffs, Z->coeffs[i].length,
                         M->coeffs[i].coeffs, t->coeffs, ctx->fqctx);
        if (success < 1)
        {
            n_poly_clear(t);
            return success;
        }

        Ai += n;
        FLINT_ASSERT(Ai <= A->length);
    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    return 1;
}



static int n_fq_polyun_zip_solvep(
    fq_nmod_mpoly_t A,
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    n_fq_polyun_t M,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success;
    slong Ai, i, n;
    n_poly_t t;

    n_poly_init(t);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    /* A has the right length, but possibly from a smaller ctx */
    if (A->length*d > A->coeffs_alloc)
    {
        slong new_alloc = FLINT_MAX(A->coeffs_alloc + A->coeffs_alloc/2, A->length*d);
        A->coeffs = FLINT_ARRAY_REALLOC(A->coeffs, new_alloc, mp_limb_t);
        A->coeffs_alloc = new_alloc;
    }

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->coeffs[i].length;
        FLINT_ASSERT(M->coeffs[i].length == n + 1);
        FLINT_ASSERT(Z->coeffs[i].length >= n);
        FLINT_ASSERT(Ai + n <= A->length);

        n_poly_fit_length(t, n);

        success = _n_fqp_zip_vand_solve(A->coeffs + d*Ai,
                         H->coeffs[i].coeffs, n,
                         Z->coeffs[i].coeffs, Z->coeffs[i].length,
                         M->coeffs[i].coeffs, t->coeffs, ctx->fqctx);
        if (success < 1)
        {
            n_poly_clear(t);
            return success;
        }

        Ai += n;
        FLINT_ASSERT(Ai <= A->length);
    }

    FLINT_ASSERT(Ai == A->length);

    n_poly_clear(t);
    return 1;
}


void n_fq_polyun_zip_start(
    n_fq_polyun_t Z,
    n_fq_polyun_t H,
    slong req_images,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong j;
    n_polyun_fit_length(Z, H->length);
    Z->length = H->length;
    for (j = 0; j < H->length; j++)
    {
        Z->exps[j] = H->exps[j];
        n_poly_fit_length(Z->coeffs + j, d*req_images);
        Z->coeffs[j].length = 0;
    }
}


int n_fq_polyu2n_add_zip_must_match(
    n_fq_polyun_t Z,
    const n_fq_bpoly_t A,
    slong cur_length,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Ai, ai;
    n_poly_struct * Zcoeffs = Z->coeffs;
    ulong * Zexps = Z->exps;
    const n_poly_struct * Acoeffs = A->coeffs;

    Ai = A->length - 1;
    ai = (Ai < 0) ? 0 : n_fq_poly_degree(A->coeffs + Ai);

    for (i = 0; i < Z->length; i++)
    {
        if (Ai >= 0 && Zexps[i] == pack_exp2(Ai, ai))
        {
            /* Z present, A present */
            _n_fq_set(Zcoeffs[i].coeffs + d*cur_length, Acoeffs[Ai].coeffs + d*ai, d);
            Zcoeffs[i].length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && _n_fq_is_zero(Acoeffs[Ai].coeffs + d*ai, d));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_fq_poly_degree(Acoeffs + Ai);
            }
        }
        else if (Ai < 0 || Zexps[i] > pack_exp2(Ai, ai))
        {
            /* Z present, A missing */
            _n_fq_zero(Zcoeffs[i].coeffs + d*cur_length, d);
            Zcoeffs[i].length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return 1;
}



#define USE_G    1
#define USE_ABAR 2
#define USE_BBAR 4

/*
    gamma = gcd(lc(A), lc(B))

    fake answers G, Abar, Bbar with

        G*Abar = gamma*A
        G*Bbar = gamma*B
        lc(G) = gamma
        lc(Abar) = lc(A)
        lc(Bbar) = lc(B)

    real answers

        rG = pp(G)
        rAbar = Abar/lc(rG)
        rBbar = Bbar/lc(rG)
*/
int fq_nmod_mpolyl_gcd_zippel_smprime(
    fq_nmod_mpoly_t rG, const slong * rGdegs,
    fq_nmod_mpoly_t rAbar,
    fq_nmod_mpoly_t rBbar,
    const fq_nmod_mpoly_t A, const slong * Adegs,
    const fq_nmod_mpoly_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    int success, use, betas_in_fp, main_tries_left;
    slong i, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    mp_limb_t * betasp;
    flint_rand_t state;
    fq_nmod_mpoly_t cont;
    fq_nmod_mpoly_t T, G, Abar, Bbar;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    fq_nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_fq_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_fq_poly_t gammaeh_cur, gammaeh_inc;
    n_fq_poly_t modulus, alphapow;
    fq_nmod_mpoly_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_poly_bpoly_stack_t St;
    mp_limb_t * c;
    fq_nmod_t start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    fq_nmod_mpoly_fit_length_reset_bits(rG, 1, bits, ctx);
    fq_nmod_mpoly_fit_length_reset_bits(rAbar, 1, bits, ctx);
    fq_nmod_mpoly_fit_length_reset_bits(rBbar, 1, bits, ctx);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fq_nmod_mpoly_degrees_si(tmp_degs, A, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, B, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    FLINT_ASSERT(gammadegs[0] == 0);
    FLINT_ASSERT(gammadegs[1] == 0);

    fq_nmod_init(start_alpha, ctx->fqctx);
    c = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_polyun_init(HG);
    n_fq_polyun_init(HAbar);
    n_fq_polyun_init(HBbar);
    n_fq_polyun_init(MG);
    n_fq_polyun_init(MAbar);
    n_fq_polyun_init(MBbar);
    n_fq_polyun_init(ZG);
    n_fq_polyun_init(ZAbar);
    n_fq_polyun_init(ZBbar);
    n_fq_bpoly_init(Aev);
    n_fq_bpoly_init(Bev);
    n_fq_bpoly_init(Gev);
    n_fq_bpoly_init(Abarev);
    n_fq_bpoly_init(Bbarev);
    n_fq_poly_init2(alphapow, 4, ctx->fqctx);
    fq_nmod_mpoly_init3(cont, 1, bits, ctx);
    fq_nmod_mpoly_init3(T, 0, bits, ctx);
    fq_nmod_mpoly_init3(G, 0, bits, ctx);
    fq_nmod_mpoly_init3(Abar, 0, bits, ctx);
    fq_nmod_mpoly_init3(Bbar, 0, bits, ctx);
    fq_nmod_mpolyn_init(Tn, bits, ctx);
    fq_nmod_mpolyn_init(Gn, bits, ctx);
    fq_nmod_mpolyn_init(Abarn, bits, ctx);
    fq_nmod_mpolyn_init(Bbarn, bits, ctx);
    n_fq_polyun_init(Aeh_cur);
    n_fq_polyun_init(Aeh_inc);
    n_fq_polyun_init(Beh_cur);
    n_fq_polyun_init(Beh_inc);
    n_fq_poly_init(gammaeh_cur);
    n_fq_poly_init(gammaeh_inc);
    n_fq_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betasp = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, ctx->fqctx);
        fq_nmod_init(alphas + i, ctx->fqctx);
    }

    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_init3(Aevals + i, 0, bits, ctx);
        fq_nmod_mpoly_init3(Bevals + i, 0, bits, ctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    Abideg = _fq_nmod_mpoly_bidegree(A, ctx);
    Bbideg = _fq_nmod_mpoly_bidegree(B, ctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    main_tries_left = 3;

choose_main:

    if (--main_tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 2; i < nvars; i++)
        fq_nmod_rand_not_zero(alphas + i, state, ctx->fqctx);

    for (i = nvars - 1; i >= 2; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + i, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + i, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, ctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_main;
        if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, ctx) != Abideg)
            goto choose_main;
        if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, ctx) != Bbideg)
            goto choose_main;
    }

    m = 2;
    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, ctx);
    fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, ctx);

    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                     Aev, Bev, ctx->fqctx, St);
    if (!success)
        goto cleanup;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_main;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, ctx);
    n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, ctx->fqctx);

    fq_nmod_mpolyn_interp_lift_sm_bpoly(Gn, Gev, ctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Abarn, Abarev, ctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Bbarn, Bbarev, ctx);

    n_fq_poly_one(modulus, ctx->fqctx);
    n_fq_set_fq_nmod(c, alphas + m, ctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, c, ctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, ctx->fqctx);
    while (1)
    {
    choose_alpha_2:

        fq_nmod_next_not_zero(alphas + m, ctx->fqctx);

        if (fq_nmod_equal(alphas + m, start_alpha, ctx->fqctx))
            goto choose_main;

        FLINT_ASSERT(alphapow->alloc >= d*2);
        alphapow->length = 2;
        _n_fq_one(alphapow->coeffs + d*0, d);
        n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
            goto choose_alpha_2;

        fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, ctx);
        fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, ctx);

        success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                     Aev, Bev, ctx->fqctx, St);
        if (!success)
            goto cleanup;

        newdegXY = n_bpoly_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_2;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_main;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, ctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, ctx->fqctx);

        n_fq_poly_eval_pow(c, modulus, alphapow, ctx->fqctx);
        n_fq_inv(c, c, ctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, c, ctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                success = fq_nmod_mpolyl_content(cont, rG, 2, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rG, rG, cont, ctx);
                fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                if (fq_nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(G, Gn, m, ctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpoly_divides(Abar, T, G, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                success = fq_nmod_mpolyl_content(cont, rAbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                if (fq_nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, ctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpoly_divides(G, T, Abar, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                success = fq_nmod_mpolyl_content(cont, rBbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                if (fq_nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                    fq_nmod_mpoly_divides(rAbar, A, rG, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, ctx);
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (fq_nmod_mpoly_divides(G, T, Bbar, ctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(Abar, T, G, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_main;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + d*1, ctx->fqctx);
    }

    for (m = 3; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Gn, G, ctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, ctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_fq_poly_one(modulus, ctx->fqctx);
        n_fq_set_fq_nmod(c, alphas + m, ctx->fqctx);
        n_fq_poly_shift_left_scalar_submul(modulus, 1, c, ctx->fqctx);

        fq_nmod_set(start_alpha, alphas + m, ctx->fqctx);

    choose_betas:

        /* only beta[0], beta[1], ..., beta[m - 1] will be used */
        betas_in_fp = (ctx->fqctx->mod.norm < FLINT_BITS/4);
        if (betas_in_fp)
        {
            for (i = 2; i < ctx->minfo->nvars; i++)
                betasp[i] = n_urandint(state, ctx->fqctx->mod.n - 3) + 2;

            fq_nmod_mpoly_monomial_evalsp2(HG, G, betasp + 2, m, ctx);
            fq_nmod_mpoly_monomial_evalsp2(HAbar, Abar, betasp + 2, m, ctx);
            fq_nmod_mpoly_monomial_evalsp2(HBbar, Bbar, betasp + 2, m, ctx);
        }
        else
        {
            for (i = 2; i < ctx->minfo->nvars; i++)
                fq_nmod_rand_not_zero(betas + i, state, ctx->fqctx);

            FLINT_ASSERT(G->bits == bits);
            FLINT_ASSERT(Abar->bits == bits);
            FLINT_ASSERT(Bbar->bits == bits);

            fq_nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, ctx);
            fq_nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, ctx);
            fq_nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, ctx);
        }

        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                  gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        req_zip_images = 1;
        if (use & USE_G)
        {
            this_length = betas_in_fp ?
                         n_polyun_product_roots(MG, HG, ctx->fqctx->mod) :
                         n_fq_polyun_product_roots(MG, HG, ctx->fqctx, St->poly_stack);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = betas_in_fp ?
                         n_polyun_product_roots(MAbar, HAbar, ctx->fqctx->mod) :
                         n_fq_polyun_product_roots(MAbar, HAbar, ctx->fqctx, St->poly_stack);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = betas_in_fp ?
                         n_polyun_product_roots(MBbar, HBbar, ctx->fqctx->mod) :
                         n_fq_polyun_product_roots(MBbar, HBbar, ctx->fqctx, St->poly_stack);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, ctx->fqctx);

            if (fq_nmod_equal(alphas + m, start_alpha, ctx->fqctx))
                goto choose_main;

            FLINT_ASSERT(alphapow->alloc >= d*2);
            alphapow->length = 2;
            _n_fq_one(alphapow->coeffs + d*0, d);
            n_fq_set_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx);

            fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
            fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
                goto choose_alpha_m;

            if (betas_in_fp)
            {
                fq_nmod_mpoly_monomial_evalsp2(Aeh_inc, Aevals + m, betasp + 2, m, ctx);
                fq_nmod_mpoly_monomial_evalsp2(Beh_inc, Bevals + m, betasp + 2, m, ctx);
                fq_nmod_mpoly_monomial_evalsp(gammaeh_inc, gammaevals + m, betasp + 2, 2, m, ctx);
                n_polyun_set(Aeh_cur, Aeh_inc);
                n_polyun_set(Beh_cur, Beh_inc);
                n_poly_set(gammaeh_cur, gammaeh_inc);
            }
            else
            {
                fq_nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, ctx);
                fq_nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, ctx);
                fq_nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, ctx);
                n_fq_polyun_set(Aeh_cur, Aeh_inc, ctx->fqctx);
                n_fq_polyun_set(Beh_cur, Beh_inc, ctx->fqctx);
                n_fq_poly_set(gammaeh_cur, gammaeh_inc, ctx->fqctx);
            }

            n_fq_polyun_zip_start(ZG, HG, req_zip_images, ctx->fqctx);
            n_fq_polyun_zip_start(ZAbar, HAbar, req_zip_images, ctx->fqctx);
            n_fq_polyun_zip_start(ZBbar, HBbar, req_zip_images, ctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                if (betas_in_fp)
                {
                    n_fq_bpoly_evalp_step_sep(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx->fqctx);
                    n_fq_bpoly_evalp_step_sep(Bev, Beh_cur, Beh_inc, Bevals + m, ctx->fqctx);
                    n_fq_poly_evalp_step_sep(c, gammaeh_cur, gammaeh_inc, gammaevals + m, ctx->fqctx);
                }
                else
                {
                    n_fq_bpoly_eval_step_sep(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx->fqctx);
                    n_fq_bpoly_eval_step_sep(Bev, Beh_cur, Beh_inc, Bevals + m, ctx->fqctx);
                    n_fq_poly_eval_step_sep(c, gammaeh_cur, gammaeh_inc, gammaevals + m, ctx->fqctx);
                }

                if (_n_fq_is_zero(c, d))
                    goto choose_betas;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                     Aev, Bev, ctx->fqctx, St);
                if (!success)
                    goto cleanup;

                newdegXY = n_bpoly_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_main;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, c, ctx->fqctx);
                if ((use & USE_G) && !n_fq_polyu2n_add_zip_must_match(ZG, Gev, cur_zip_image, ctx->fqctx))
                    goto choose_main;
                if ((use & USE_ABAR) && !n_fq_polyu2n_add_zip_must_match(ZAbar, Abarev, cur_zip_image, ctx->fqctx))
                    goto choose_main;
                if ((use & USE_BBAR) && !n_fq_polyu2n_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, ctx->fqctx))
                    goto choose_main;
            }

            if (betas_in_fp)
            {
                if ((use & USE_G) && n_fq_polyun_zip_solvep(G, ZG, HG, MG, ctx) < 1)
                    goto choose_main;
                if ((use & USE_ABAR) && n_fq_polyun_zip_solvep(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                    goto choose_main;
                if ((use & USE_BBAR) && n_fq_polyun_zip_solvep(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                    goto choose_main;
            }
            else
            {
                if ((use & USE_G) && n_fq_polyun_zip_solve(G, ZG, HG, MG, ctx) < 1)
                    goto choose_main;
                if ((use & USE_ABAR) && n_fq_polyun_zip_solve(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                    goto choose_main;
                if ((use & USE_BBAR) && n_fq_polyun_zip_solve(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                    goto choose_main;
            }

            n_fq_poly_eval_pow(c, modulus, alphapow, ctx->fqctx);
            n_fq_inv(c, c, ctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, c, ctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyl_content(cont, rG, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpoly_divides(rG, rG, cont, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    if (fq_nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                        fq_nmod_mpoly_divides(rBbar, B, rG,  ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(rAbar, T, rG, ctx))
                    {
                        fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fq_nmod_mpoly_swap(G, rG, ctx);
                            fq_nmod_mpoly_swap(Abar, rAbar, ctx);
                            fq_nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyl_content(cont, rAbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    success = fq_nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                    FLINT_ASSERT(success);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    if (fq_nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                        fq_nmod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(rG, T, rAbar, ctx))
                    {
                        fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fq_nmod_mpoly_swap(G, rG, ctx);
                            fq_nmod_mpoly_swap(Abar, rAbar, ctx);
                            fq_nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fq_nmod_mpolyl_content(cont, rBbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fq_nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    if (fq_nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                        fq_nmod_mpoly_divides(rAbar, A, rG, ctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fq_nmod_mpoly_divides(rG, T, rBbar, ctx))
                    {
                        fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (fq_nmod_mpoly_divides(rAbar, T, rG, ctx))
                        {
                            fq_nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fq_nmod_mpoly_swap(G, rG, ctx);
                            fq_nmod_mpoly_swap(Abar, rAbar, ctx);
                            fq_nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + d*1, alphas + m, ctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + d*1, ctx->fqctx);
        }
    }

    success = 1;

cleanup:

    n_fq_polyun_clear(HG);
    n_fq_polyun_clear(HAbar);
    n_fq_polyun_clear(HBbar);
    n_fq_polyun_clear(MG);
    n_fq_polyun_clear(MAbar);
    n_fq_polyun_clear(MBbar);
    n_fq_polyun_clear(ZG);
    n_fq_polyun_clear(ZAbar);
    n_fq_polyun_clear(ZBbar);
    n_fq_bpoly_clear(Aev);
    n_fq_bpoly_clear(Bev);
    n_fq_bpoly_clear(Gev);
    n_fq_bpoly_clear(Abarev);
    n_fq_bpoly_clear(Bbarev);
    n_fq_poly_clear(alphapow);
    fq_nmod_mpoly_clear(cont, ctx);
    fq_nmod_mpoly_clear(T, ctx);
    fq_nmod_mpoly_clear(G, ctx);
    fq_nmod_mpoly_clear(Abar, ctx);
    fq_nmod_mpoly_clear(Bbar, ctx);
    fq_nmod_mpolyn_clear(Tn, ctx);
    fq_nmod_mpolyn_clear(Gn, ctx);
    fq_nmod_mpolyn_clear(Abarn, ctx);
    fq_nmod_mpolyn_clear(Bbarn, ctx);
    n_fq_polyun_clear(Aeh_cur);
    n_fq_polyun_clear(Aeh_inc);
    n_fq_polyun_clear(Beh_cur);
    n_fq_polyun_clear(Beh_inc);
    n_fq_poly_clear(gammaeh_cur);
    n_fq_poly_clear(gammaeh_inc);
    n_fq_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    fq_nmod_clear(start_alpha, ctx->fqctx);
    flint_free(c);

    flint_free(betasp);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_clear(betas + i, ctx->fqctx);
        fq_nmod_clear(alphas + i, ctx->fqctx);
    }
    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_clear(Aevals + i, ctx);
        fq_nmod_mpoly_clear(Bevals + i, ctx);
        fq_nmod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

    return success;

gcd_is_trivial:

    fq_nmod_mpoly_one(rG, ctx);
    fq_nmod_mpoly_set(rAbar, A, ctx);
    fq_nmod_mpoly_set(rBbar, B, ctx);

    success = 1;

    goto cleanup;
}


/* TODO: make this the real fq_nmod_mpolyn_interp_mcrt_lg_mpoly */
static int newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
    slong * lastdeg,
    fq_nmod_mpolyn_t H,
    const fq_nmod_mpoly_ctx_t ctx,
    const n_fq_poly_t m,
    const mp_limb_t * inv_m_eval,
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ectx,
    const bad_fq_nmod_embed_t emb)
{
    slong lgd = fq_nmod_ctx_degree(ectx->fqctx);
    slong i;
#if FLINT_WANT_ASSERT
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);
#endif
    int changed = 0;
    mp_limb_t * u, * v, * tmp;
    n_fq_poly_struct * w, * u_sm;
    n_poly_stack_t St;

    FLINT_ASSERT(H->length == A->length);
    FLINT_ASSERT(H->bits == A->bits);

    n_poly_stack_init(St);
    n_poly_stack_fit_request(St, 3);

    w = n_poly_stack_take_top(St);
    u_sm = n_poly_stack_take_top(St);

    tmp = n_poly_stack_vec_init(St, lgd*(2 + N_FQ_MUL_ITCH));
    u = tmp + lgd*N_FQ_MUL_ITCH;
    v = u + lgd;

    for (i = 0; i < A->length; i++)
    {
        FLINT_ASSERT(mpoly_monomial_equal(H->exps + N*i, A->exps + N*i, N));
        bad_n_fq_embed_sm_to_lg(u, H->coeffs + i, emb);
        _n_fq_sub(v, A->coeffs + lgd*i, u, lgd, ectx->fqctx->mod);
        if (!_n_fq_is_zero(v, lgd))
        {
            changed = 1;
            _n_fq_mul(u, v, inv_m_eval, ectx->fqctx, tmp);
            bad_n_fq_embed_lg_to_sm(u_sm, u, emb);
            n_fq_poly_mul_(w, u_sm, m, ctx->fqctx, St);
            n_fq_poly_add(H->coeffs + i, H->coeffs + i, w, ctx->fqctx);
        }

        *lastdeg = FLINT_MAX(*lastdeg, n_fq_poly_degree(H->coeffs + i));

        FLINT_ASSERT(n_fq_poly_degree(H->coeffs + i) < n_fq_poly_degree(m) +
                                      fq_nmod_poly_degree(emb->h, ctx->fqctx));
    }

    n_poly_stack_vec_clear(St);

    n_poly_stack_give_back(St, 2);

    n_poly_stack_clear(St);

    return changed;
}



int fq_nmod_mpolyl_gcd_zippel_lgprime(
    fq_nmod_mpoly_t rG, const slong * rGdegs,
    fq_nmod_mpoly_t rAbar,
    fq_nmod_mpoly_t rBbar,
    const fq_nmod_mpoly_t A, const slong * Adegs,
    const fq_nmod_mpoly_t B, const slong * Bdegs,
    const fq_nmod_mpoly_t gamma, const slong * gammadegs,
    const fq_nmod_mpoly_ctx_t smctx)
{
    slong lgd;
    int success, use;
    slong i, m;
    slong nvars = smctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    flint_rand_t state;
    fq_nmod_mpoly_t cont;
    fq_nmod_mpoly_t T, G, Abar, Bbar;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    fq_nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_fq_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_fq_poly_t gammaeh_cur, gammaeh_inc;
    fq_nmod_mpoly_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_fq_poly_t modulus, alphapow;
    n_poly_bpoly_stack_t St;
    fq_nmod_t start_alpha;
    n_poly_t tmp;  /* tmp arithmetic space */
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;
    bad_fq_nmod_mpoly_embed_chooser_t embc;
    bad_fq_nmod_embed_struct * cur_emb;
    fq_nmod_mpoly_ctx_t lgctx;
    fq_nmod_mpolyn_t gamman;
    fq_nmod_mpolyn_t An, Bn;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    fq_nmod_mpoly_fit_length_reset_bits(rG, 1, bits, smctx);
    fq_nmod_mpoly_fit_length_reset_bits(rAbar, 1, bits, smctx);
    fq_nmod_mpoly_fit_length_reset_bits(rBbar, 1, bits, smctx);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fq_nmod_mpoly_degrees_si(tmp_degs, A, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, B, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        fq_nmod_mpoly_degrees_si(tmp_degs, gamma, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    FLINT_ASSERT(gammadegs[0] == 0);
    FLINT_ASSERT(gammadegs[1] == 0);

    n_fq_polyun_init(HG);
    n_fq_polyun_init(HAbar);
    n_fq_polyun_init(HBbar);
    n_fq_polyun_init(MG);
    n_fq_polyun_init(MAbar);
    n_fq_polyun_init(MBbar);
    n_fq_polyun_init(ZG);
    n_fq_polyun_init(ZAbar);
    n_fq_polyun_init(ZBbar);
    n_fq_bpoly_init(Aev);
    n_fq_bpoly_init(Bev);
    n_fq_bpoly_init(Gev);
    n_fq_bpoly_init(Abarev);
    n_fq_bpoly_init(Bbarev);
    fq_nmod_mpoly_init3(cont, 1, bits, smctx);
    fq_nmod_mpoly_init3(T, 0, bits, smctx);
    fq_nmod_mpoly_init3(G, 0, bits, smctx);
    fq_nmod_mpoly_init3(Abar, 0, bits, smctx);
    fq_nmod_mpoly_init3(Bbar, 0, bits, smctx);
    fq_nmod_mpolyn_init(Tn, bits, smctx);
    fq_nmod_mpolyn_init(Gn, bits, smctx);
    fq_nmod_mpolyn_init(Abarn, bits, smctx);
    fq_nmod_mpolyn_init(Bbarn, bits, smctx);
    n_fq_polyun_init(Aeh_cur);
    n_fq_polyun_init(Aeh_inc);
    n_fq_polyun_init(Beh_cur);
    n_fq_polyun_init(Beh_inc);
    n_fq_poly_init(gammaeh_cur);
    n_fq_poly_init(gammaeh_inc);
    n_fq_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);
    fq_nmod_mpolyn_init(An, bits, smctx);
    fq_nmod_mpolyn_init(Bn, bits, smctx);
    fq_nmod_mpolyn_init(gamman, bits, smctx);

    /* alphas[nvars - 1] not used - it is replaced by cur_emb */
    betas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    alphas = FLINT_ARRAY_ALLOC(nvars, fq_nmod_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, smctx->fqctx);
        fq_nmod_init(alphas + i, smctx->fqctx);
    }

    flint_randinit(state);
    cur_emb = bad_fq_nmod_mpoly_embed_chooser_init(embc, lgctx, smctx, state);
    lgd = fq_nmod_ctx_degree(lgctx->fqctx);
    n_poly_init2(tmp, lgd);
    n_poly_init2(alphapow, 2*lgd);
    fq_nmod_init(start_alpha, lgctx->fqctx);

    /* Aevals[nvars] does not exist - it is replaced by An */
    Aevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_init3(Aevals + i, 0, bits, smctx);
        fq_nmod_mpoly_init3(Bevals + i, 0, bits, smctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, smctx);
    }
    fq_nmod_mpoly_cvtto_mpolyn(An, A, nvars - 1, smctx);
    fq_nmod_mpoly_cvtto_mpolyn(Bn, B, nvars - 1, smctx);
    fq_nmod_mpoly_cvtto_mpolyn(gamman, gamma, nvars - 1, smctx);

    Abideg = _fq_nmod_mpoly_bidegree(A, smctx);
    Bbideg = _fq_nmod_mpoly_bidegree(B, smctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    goto got_alpha_m;

increase_degree:

    do {
        cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
        if (cur_emb == NULL)
        {
            /* ran out of primes */
            success = 0;
            goto cleanup;
        }
    } while (fq_nmod_ctx_degree(lgctx->fqctx) <= lgd);

    lgd = fq_nmod_ctx_degree(lgctx->fqctx);

    goto got_alpha_m;

choose_main:

    cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
    if (cur_emb == NULL)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    lgd = fq_nmod_ctx_degree(lgctx->fqctx);

got_alpha_m:

    n_poly_fit_length(tmp, lgd);
    n_poly_fit_length(alphapow, 2*lgd);

    for (i = 2; i < nvars - 1; i++)
        fq_nmod_rand_not_zero(alphas + i, state, lgctx->fqctx);

    i = nvars - 1;
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + i, An, lgctx, smctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + i, Bn, lgctx, smctx, cur_emb);
    fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + i, gamman, lgctx, smctx, cur_emb);
    if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
        goto choose_main;
    if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
        goto choose_main;
    if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
        goto choose_main;
    for (i--; i >= 2; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
            goto choose_main;
        if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
            goto choose_main;
        if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
            goto choose_main;
    }

    m = 2;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, lgctx);
    fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, lgctx);

    success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
    if (!success)
        goto increase_degree;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_main;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
    n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

    if (nvars == 3)
    {
        fq_nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Gn, smctx, Gev, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Abarn, smctx, Abarev, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Bbarn, smctx, Bbarev, lgctx, cur_emb);

        n_fq_poly_set(modulus, cur_emb->h_as_n_fq_poly, smctx->fqctx);

        while (1)
        {
        choose_alpha_2_last:

            cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
            if (cur_emb == NULL)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            lgd = fq_nmod_ctx_degree(lgctx->fqctx);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx, cur_emb);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_2_last;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_2_last;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_2_last;

            fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, lgctx);
            fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, lgctx);

            success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
            if (!success)
                goto increase_degree;

            newdegXY = n_bpoly_bidegree(Gev);
            if (newdegXY > GdegboundXY)
                goto choose_alpha_2_last;
            if (newdegXY < GdegboundXY)
            {
                GdegboundXY = newdegXY;
                if (GdegboundXY == 0)
                    goto gcd_is_trivial;
                goto choose_main;
            }

            gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
            n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_lg_bpoly(
                        &lastdeg, Gn, Tn, modulus, smctx, Gev, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rG, rG, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (fq_nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_lg_bpoly(
                  &lastdeg, Abarn, Tn, modulus, smctx, Abarev, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_lg_bpoly(
                  &lastdeg, Bbarn, Tn, modulus, smctx, Bbarev, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    fq_nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            n_fq_poly_mul_(modulus, modulus, cur_emb->h_as_n_fq_poly,
                                                 smctx->fqctx, St->poly_stack);
        }

        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyn_interp_lift_sm_bpoly(Gn, Gev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Abarn, Abarev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(Bbarn, Bbarev, lgctx);

    FLINT_ASSERT(tmp->alloc >= lgd);
    n_fq_poly_one(modulus, lgctx->fqctx);
    n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

    fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

    while (1)
    {
    choose_alpha_2:

        fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
        if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
            goto increase_degree;

        FLINT_ASSERT(alphapow->alloc >= lgd*2);
        alphapow->length = 2;
        _n_fq_one(alphapow->coeffs + lgd*0, lgd);
        n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);

        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
            goto choose_alpha_2;

        fq_nmod_mpoly_get_n_fq_bpoly(Aev, Aevals + m, 0, 1, lgctx);
        fq_nmod_mpoly_get_n_fq_bpoly(Bev, Bevals + m, 0, 1, lgctx);

        success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
        if (!success)
            goto increase_degree;

        newdegXY = n_bpoly_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_2;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_main;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
        n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                              &lastdeg, Gn, Tn, Gev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(G, Gn, m, lgctx);
            fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(Abar, T, G, lgctx))
            {
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(Bbar, T, G, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(Abar, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                        &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, lgctx);
            fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(G, T, Abar, lgctx))
            {
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(Bbar, T, G, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(Bbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                        &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, lgctx);
            fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(G, T, Bbar, lgctx))
            {
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(Abar, T, G, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(G, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(Abar, bits, lgctx);
                    break;
                }
            }
        }

        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_main;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
    }

    for (m = 3; m < nvars - 1; m++)
    {
        /* G, Abar, Bbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Gn, G, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, lgctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_one(modulus, lgctx->fqctx);
        n_fq_set_fq_nmod(tmp->coeffs, alphas + m, lgctx->fqctx);
        n_fq_poly_shift_left_scalar_submul(modulus, 1, tmp->coeffs, lgctx->fqctx);

        fq_nmod_set(start_alpha, alphas + m, lgctx->fqctx);

    choose_betas_m:

        /* only beta[2], beta[1], ..., beta[m - 1] will be used */
        for (i = 2; i < nvars; i++)
            fq_nmod_rand_not_zero(betas + i, state, lgctx->fqctx);

        FLINT_ASSERT(G->bits == bits);
        FLINT_ASSERT(Abar->bits == bits);
        FLINT_ASSERT(Bbar->bits == bits);

        fq_nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, lgctx);
        fq_nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, lgctx);
        fq_nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, lgctx);
        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                   gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        req_zip_images = 1;
        if (use & USE_G)
        {
            this_length = n_fq_polyun_product_roots(MG, HG, lgctx->fqctx, St->poly_stack);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = n_fq_polyun_product_roots(MAbar, HAbar, lgctx->fqctx, St->poly_stack);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = n_fq_polyun_product_roots(MBbar, HBbar, lgctx->fqctx, St->poly_stack);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            fq_nmod_next_not_zero(alphas + m, lgctx->fqctx);
            if (fq_nmod_equal(alphas + m, start_alpha, lgctx->fqctx))
                goto choose_main;

            FLINT_ASSERT(alphapow->alloc >= lgd*2);
            alphapow->length = 2;
            _n_fq_one(alphapow->coeffs + lgd*0, lgd);
            n_fq_set_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx);

            fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + m, Aevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_repack_bits_inplace(Aevals + m, bits, lgctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + m, Bevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_repack_bits_inplace(Bevals + m, bits, lgctx);
            fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + m, gammaevals + m + 1, m, alphas + m, lgctx);
            fq_nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, lgctx);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_m;

            fq_nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, lgctx);
            n_fq_polyun_set(Aeh_cur, Aeh_inc, lgctx->fqctx);
            n_fq_polyun_set(Beh_cur, Beh_inc, lgctx->fqctx);
            n_fq_poly_set(gammaeh_cur, gammaeh_inc, lgctx->fqctx);

            n_fq_polyun_zip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            n_fq_polyun_zip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            n_fq_polyun_zip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step_sep(Aev, Aeh_cur, Aeh_inc, Aevals + m, lgctx->fqctx);
                n_fq_bpoly_eval_step_sep(Bev, Beh_cur, Beh_inc, Bevals + m, lgctx->fqctx);
                FLINT_ASSERT(tmp->alloc >= lgd);
                n_fq_poly_eval_step_sep(tmp->coeffs, gammaeh_cur, gammaeh_inc, gammaevals + m, lgctx->fqctx);
                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_m;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas_m;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas_m;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
                if (!success)
                    goto increase_degree;

                newdegXY = n_bpoly_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas_m;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_main;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2n_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_ABAR) && !n_fq_polyu2n_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_BBAR) && !n_fq_polyu2n_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
            }

            if ((use & USE_G) && n_fq_polyun_zip_solve(G, ZG, HG, MG, lgctx) < 1)
                goto choose_main;
            if ((use & USE_ABAR) && n_fq_polyun_zip_solve(Abar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_main;
            if ((use & USE_BBAR) && n_fq_polyun_zip_solve(Bbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_main;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                    &lastdeg, Gn, G, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, lgctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(rAbar, T, rG, lgctx))
                {
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(rBbar, T, rG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(G, rG, lgctx);
                        fq_nmod_mpoly_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpoly_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                              &lastdeg, Abarn, Abar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, lgctx);
                fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(rG, T, rAbar, lgctx))
                {
                    fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(rBbar, T, rG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(G, rG, lgctx);
                        fq_nmod_mpoly_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpoly_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                              &lastdeg, Bbarn, Bbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, lgctx);
                fq_nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(rG, T, rBbar, lgctx))
                {
                    fq_nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(rAbar, T, rG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(rG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, lgctx);
                        fq_nmod_mpoly_swap(G, rG, lgctx);
                        fq_nmod_mpoly_swap(Abar, rAbar, lgctx);
                        fq_nmod_mpoly_swap(Bbar, rBbar, lgctx);
                        break;
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
        }
    }

    m = nvars - 1;
    {
        /* G, Abar, Bbar are in Fq/alpha(gen(m-1))[gen(0), ..., gen(m - 1)] */

        fq_nmod_mpolyn_interp_lift_lg_mpoly(Gn, smctx, G, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_mpoly(Abarn, smctx, Abar, lgctx, cur_emb);
        fq_nmod_mpolyn_interp_lift_lg_mpoly(Bbarn, smctx, Bbar, lgctx, cur_emb);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_fq_poly_set(modulus, cur_emb->h_as_n_fq_poly, smctx->fqctx);

        while (1)
        {
        choose_alpha_last:

            cur_emb = bad_fq_nmod_mpoly_embed_chooser_next(embc, lgctx, smctx, state);
            if (cur_emb == NULL)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            lgd = fq_nmod_ctx_degree(lgctx->fqctx);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx, cur_emb);
            fq_nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx, cur_emb);
            if (fq_nmod_mpoly_is_zero(gammaevals + m, lgctx))
                goto choose_alpha_last;
            if (Aevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + m, lgctx) != Abideg)
                goto choose_alpha_last;
            if (Bevals[m].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + m, lgctx) != Bbideg)
                goto choose_alpha_last;

        choose_betas_last:

            /* only beta[2], ..., beta[m - 1] will be used */
            for (i = 2; i < nvars; i++)
                fq_nmod_rand_not_zero(betas + i, state, lgctx->fqctx);

            fq_nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, lgctx);

            if (use == 0)
            {
                this_length = gamma->length + A->length + B->length;
                use = nmod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                       gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
            }

            req_zip_images = 1;
            if (use & USE_G)
            {
                this_length = n_fq_polyun_product_roots(MG, HG, lgctx->fqctx, St->poly_stack);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);
            }
            if (use & USE_ABAR)
            {
                this_length = n_fq_polyun_product_roots(MAbar, HAbar, lgctx->fqctx, St->poly_stack);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);
            }
            if (use & USE_BBAR)
            {
                this_length = n_fq_polyun_product_roots(MBbar, HBbar, lgctx->fqctx, St->poly_stack);
                req_zip_images = FLINT_MAX(req_zip_images, this_length);
            }

            fq_nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, lgctx);
            n_fq_polyun_set(Aeh_cur, Aeh_inc, lgctx->fqctx);
            n_fq_polyun_set(Beh_cur, Beh_inc, lgctx->fqctx);
            n_fq_poly_set(gammaeh_cur, gammaeh_inc, lgctx->fqctx);

            n_fq_polyun_zip_start(ZG, HG, req_zip_images, lgctx->fqctx);
            n_fq_polyun_zip_start(ZAbar, HAbar, req_zip_images, lgctx->fqctx);
            n_fq_polyun_zip_start(ZBbar, HBbar, req_zip_images, lgctx->fqctx);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_fq_bpoly_eval_step_sep(Aev, Aeh_cur, Aeh_inc, Aevals + m, lgctx->fqctx);
                n_fq_bpoly_eval_step_sep(Bev, Beh_cur, Beh_inc, Bevals + m, lgctx->fqctx);
                FLINT_ASSERT(tmp->alloc >= lgd);
                n_fq_poly_eval_step_sep(tmp->coeffs, gammaeh_cur, gammaeh_inc, gammaevals + m, lgctx->fqctx);
                if (_n_fq_is_zero(tmp->coeffs, lgd))
                    goto choose_betas_last;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas_last;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas_last;

                success = n_fq_bpoly_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                   Aev, Bev, lgctx->fqctx, St);
                if (!success)
                    goto cleanup;

                newdegXY = n_bpoly_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas_last;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_main;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2n_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_ABAR) && !n_fq_polyu2n_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
                if ((use & USE_BBAR) && !n_fq_polyu2n_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_main;
            }
            if ((use & USE_G) && n_fq_polyun_zip_solve(G, ZG, HG, MG, lgctx) < 1)
                goto choose_main;
            if ((use & USE_ABAR) && n_fq_polyun_zip_solve(Abar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_main;
            if ((use & USE_BBAR) && n_fq_polyun_zip_solve(Bbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_main;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            bad_n_fq_embed_sm_to_lg(tmp->coeffs, modulus, cur_emb);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
                                    &lastdeg, Gn, smctx, modulus, tmp->coeffs,
                                                            G, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rG, rG, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (fq_nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
                                 &lastdeg, Abarn, smctx, modulus, tmp->coeffs,
                                                         Abar, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    fq_nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !newfq_nmod_mpolyn_interp_mcrt_lg_mpoly(
                                 &lastdeg, Bbarn, smctx, modulus, tmp->coeffs,
                                                         Bbar, lgctx, cur_emb))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = fq_nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                fq_nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                fq_nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (fq_nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    fq_nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    fq_nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_main;
            }

            n_fq_poly_mul_(modulus, modulus, cur_emb->h_as_n_fq_poly,
                                                 smctx->fqctx, St->poly_stack);
        }
    }

    success = 1;

cleanup:

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_clear(betas + i, smctx->fqctx);
        fq_nmod_clear(alphas + i, smctx->fqctx);
    }
    flint_free(betas);
    flint_free(alphas);

    fq_nmod_mpolyn_clear(An, smctx);
    fq_nmod_mpolyn_clear(Bn, smctx);
    fq_nmod_mpolyn_clear(gamman, smctx);
    bad_fq_nmod_mpoly_embed_chooser_clear(embc, lgctx, smctx, state);

    n_fq_polyun_clear(HG);
    n_fq_polyun_clear(HAbar);
    n_fq_polyun_clear(HBbar);
    n_fq_polyun_clear(MG);
    n_fq_polyun_clear(MAbar);
    n_fq_polyun_clear(MBbar);
    n_fq_polyun_clear(ZG);
    n_fq_polyun_clear(ZAbar);
    n_fq_polyun_clear(ZBbar);
    n_fq_bpoly_clear(Aev);
    n_fq_bpoly_clear(Bev);
    n_fq_bpoly_clear(Gev);
    n_fq_bpoly_clear(Abarev);
    n_fq_bpoly_clear(Bbarev);
    fq_nmod_mpoly_clear(cont, smctx);
    fq_nmod_mpoly_clear(T, smctx);
    fq_nmod_mpoly_clear(G, smctx);
    fq_nmod_mpoly_clear(Abar, smctx);
    fq_nmod_mpoly_clear(Bbar, smctx);
    fq_nmod_mpolyn_clear(Tn, smctx);
    fq_nmod_mpolyn_clear(Gn, smctx);
    fq_nmod_mpolyn_clear(Abarn, smctx);
    fq_nmod_mpolyn_clear(Bbarn, smctx);
    n_fq_polyun_clear(Aeh_cur);
    n_fq_polyun_clear(Aeh_inc);
    n_fq_polyun_clear(Beh_cur);
    n_fq_polyun_clear(Beh_inc);
    n_fq_poly_clear(gammaeh_cur);
    n_fq_poly_clear(gammaeh_inc);
    n_fq_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    n_poly_clear(tmp);
    n_fq_poly_clear(alphapow);
    fq_nmod_clear(start_alpha, smctx->fqctx);

    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_clear(Aevals + i, smctx);
        fq_nmod_mpoly_clear(Bevals + i, smctx);
        fq_nmod_mpoly_clear(gammaevals + i, smctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

    return success;

gcd_is_trivial:

    fq_nmod_mpoly_one(rG, smctx);
    fq_nmod_mpoly_set(rAbar, A, smctx);
    fq_nmod_mpoly_set(rBbar, B, smctx);

    success = 1;

    goto cleanup;
}


/* should find its way back here in interesting cases */
int fq_nmod_mpoly_gcd_zippel2(
    fq_nmod_mpoly_t G,
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    if (fq_nmod_mpoly_is_zero(A, ctx) || fq_nmod_mpoly_is_zero(B, ctx))
        return fq_nmod_mpoly_gcd(G, A, B, ctx);

    return _fq_nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL2);
}

