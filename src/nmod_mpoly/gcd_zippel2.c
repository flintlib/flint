/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "nmod_mpoly_factor.h"
#include "fq_nmod_mpoly.h"
#include "fq_nmod_mpoly_factor.h"

/*
    evaluation at
        gen(start) -> betas[0]
        gen(start+1) -> betas[1]
        ...
        gen(stop-1) -> betas[stop-start-1]

    the other gen are assumed to not appear in A
*/
void _nmod_mpoly_monomial_evals_cache(
    n_poly_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * betas,
    slong start,
    slong stop,
    const mpoly_ctx_t mctx,
    nmod_t mod)
{
    slong i, Ai;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    n_poly_struct * caches;
    mp_limb_t * c;
    slong num = stop - start;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(num > 0);

    caches = FLINT_ARRAY_ALLOC(3*num, n_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*num, slong);
    shift = off + num;
    for (i = 0; i < num; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i + start, Abits, mctx);
        n_poly_init(caches + 3*i + 0);
        n_poly_init(caches + 3*i + 1);
        n_poly_init(caches + 3*i + 2);
        nmod_pow_cache_start(betas[i], caches + 3*i + 0,
                                       caches + 3*i + 1,
                                       caches + 3*i + 2);
    }

    n_poly_fit_length(E, Alen);
    E->length = Alen;

    for (Ai = 0; Ai < Alen; Ai++)
    {
        c = E->coeffs + Ai;
        *c = 1;
        for (i = 0; i < num; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            *c = nmod_pow_cache_mulpow_ui(*c, ei, caches + 3*i + 0,
                                                  caches + 3*i + 1,
                                                  caches + 3*i + 2, mod);
        }
    }

    for (i = 0; i < num; i++)
    {
        n_poly_clear(caches + 3*i + 0);
        n_poly_clear(caches + 3*i + 1);
        n_poly_clear(caches + 3*i + 2);
    }
    flint_free(caches);
    flint_free(off);
}


/*
    evaluation at
        gen(0) -> x
        gen(1) -> y
        gen(2) -> betas[0]
        gen(3) -> betas[1]
        ...
        gen(m-1) -> betas[m-3]
*/
void _nmod_mpoly_monomial_evals2_cache(
    n_polyun_t E,
    const ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mp_limb_t * betas,
    slong m,
    const mpoly_ctx_t mctx,
    nmod_t mod)
{
    slong i, Ai, Ei;
    ulong e0, e1, e01;
    ulong mask = (-UWORD(1)) >> (FLINT_BITS - Abits);
    slong N = mpoly_words_per_exp_sp(Abits, mctx);
    slong * off, * shift;
    n_poly_struct * caches;
    mp_limb_t * c;

    FLINT_ASSERT(Abits <= FLINT_BITS);
    FLINT_ASSERT(Alen > 0);
    FLINT_ASSERT(m > 2);

    caches = FLINT_ARRAY_ALLOC(3*(m - 2), n_poly_struct);
    off = FLINT_ARRAY_ALLOC(2*m, slong);
    shift = off + m;
    for (i = 0; i < m; i++)
    {
        mpoly_gen_offset_shift_sp(&off[i], &shift[i], i, Abits, mctx);
        if (i >= 2)
        {
            n_poly_init(caches + 3*(i - 2) + 0);
            n_poly_init(caches + 3*(i - 2) + 1);
            n_poly_init(caches + 3*(i - 2) + 2);
            nmod_pow_cache_start(betas[i - 2], caches + 3*(i - 2) + 0,
                                           caches + 3*(i - 2) + 1,
                                           caches + 3*(i - 2) + 2);
        }
    }

    Ai = 0;
    Ei = 0;

    e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
    e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
    e01 = pack_exp2(e0, e1);
    n_polyun_fit_length(E, Ei + 1);
    E->exps[Ei] = e01;
    n_poly_fit_length(E->coeffs + Ei, 1);
    c = E->coeffs[Ei].coeffs + 0;
    E->coeffs[Ei].length = 1;
    Ei++;
    *c = 1;
    for (i = 2; i < m; i++)
    {
        ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
        *c = nmod_pow_cache_mulpow_ui(*c, ei, caches + 3*(i - 2) + 0,
                                              caches + 3*(i - 2) + 1,
                                              caches + 3*(i - 2) + 2, mod);
    }

    for (Ai++; Ai < Alen; Ai++)
    {
        e0 = (Aexps[N*Ai + off[0]] >> shift[0]) & mask;
        e1 = (Aexps[N*Ai + off[1]] >> shift[1]) & mask;
        e01 = pack_exp2(e0, e1);
        if (e01 == E->exps[Ei - 1])
        {
            slong len = E->coeffs[Ei - 1].length;
            n_poly_fit_length(E->coeffs + Ei - 1, len + 1);
            c = E->coeffs[Ei - 1].coeffs + len;
            E->coeffs[Ei - 1].length = len + 1;
        }
        else
        {
            n_polyun_fit_length(E, Ei + 1);
            E->exps[Ei] = e01;
            n_poly_fit_length(E->coeffs + Ei, 1);
            c = E->coeffs[Ei].coeffs + 0;
            E->coeffs[Ei].length = 1;
            Ei++;
        }

        *c = 1;
        for (i = 2; i < m; i++)
        {
            ulong ei = (Aexps[N*Ai + off[i]] >> shift[i]) & mask;
            *c = nmod_pow_cache_mulpow_ui(*c, ei, caches + 3*(i - 2) + 0,
                                                  caches + 3*(i - 2) + 1,
                                                  caches + 3*(i - 2) + 2, mod);
        }
    }

    E->length = Ei;

    for (i = 0; i < m - 2; i++)
    {
        n_poly_clear(caches + 3*i + 0);
        n_poly_clear(caches + 3*i + 1);
        n_poly_clear(caches + 3*i + 2);
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

static ulong _nmod_mpoly_bidegree(
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}

static ulong _fq_nmod_mpoly_bidegree(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}


#define USE_G    1
#define USE_ABAR 2
#define USE_BBAR 4

/*
    The cost of interpolating gen(m) into
        G mod (gen(m) - alpha) = sum_i x^e_i*y^f_i c_i(gen(0), ..., gen(m-1))
    is (# = length, deg = deg_(gen(m))):
        (deg G)*(#A + #B + #gamma)                      eval setup
      + (deg G)*(max_i #c_i)*(#A + #B + #gamma)         zip eval
      + (deg G)*(max_i #c_i)*(deg_x AB)^2 (deg_y AB)^2  base gcd
      + (deg G)*(max_i #c_i)*(sum_i # c_i)              zip solve
      + (deg(G))^2*(sum_i #c_i)                         final interp
*/
static double interp_cost(
    double degG,
    double numABgamma,
    double maxnumci,
    double totnumci,
    double degxAB,
    double degyAB)
{
    return degG*(degG*totnumci + numABgamma + 0.01*maxnumci*(
                     numABgamma + totnumci + (degxAB*degyAB)*(degxAB*degyAB)));
}

int mpoly_gcd_get_use_first(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg)
{
    int use = 0;
    slong lower = FLINT_MAX(gammadeg, rGdeg);
    slong upper = gammadeg + FLINT_MIN(FLINT_MIN(Adeg, Bdeg), rGdeg);
    if (lower <= upper)
    {
        slong Gdeg = ((ulong)upper + (ulong)lower)/2;
        slong Abardeg = gammadeg + Adeg - Gdeg;
        slong Bbardeg = gammadeg + Bdeg - Gdeg;

        if (Gdeg <= Abardeg && Gdeg <= Bbardeg)
            use |= USE_G;

        if (Abardeg <= Gdeg && Abardeg <= Bbardeg)
            use |= USE_ABAR;

        if (Bbardeg <= Gdeg && Bbardeg <= Abardeg)
            use |= USE_BBAR;
    }

    if (use == 0)
        use = USE_G | USE_ABAR | USE_BBAR;

    return use;
}

int nmod_mpoly_gcd_get_use_new(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const n_polyun_t G,
    const n_polyun_t Abar,
    const n_polyun_t Bbar)
{
    int use = 0;
    slong i, lower = FLINT_MAX(gammadeg, rGdeg);
    slong upper = gammadeg + FLINT_MIN(FLINT_MIN(Adeg, Bdeg), rGdeg);
    if (lower <= upper)
    {
        slong Gdeg = ((ulong)upper + (ulong)lower)/2;
        slong maxnumci, totnumci;
        double Gcost, Abarcost, Bbarcost;

        maxnumci = totnumci = 0;
        for (i = 0; i < G->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, G->coeffs[i].length);
            totnumci += G->coeffs[i].length;
        }
        FLINT_ASSERT(Gdeg >= 0);
        Gcost = interp_cost(Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        maxnumci = totnumci = 0;
        for (i = 0; i < Abar->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, Abar->coeffs[i].length);
            totnumci += Abar->coeffs[i].length;
        }
        FLINT_ASSERT(gammadeg + Adeg - Gdeg >= 0);
        Abarcost = interp_cost(gammadeg + Adeg - Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        maxnumci = totnumci = 0;
        for (i = 0; i < Bbar->length; i++)
        {
            maxnumci = FLINT_MAX(maxnumci, Bbar->coeffs[i].length);
            totnumci += Bbar->coeffs[i].length;
        }
        FLINT_ASSERT(gammadeg + Bdeg - Gdeg >= 0);
        Bbarcost = interp_cost(gammadeg + Bdeg - Gdeg,
                       numABgamma, maxnumci, totnumci, degxAB, degyAB);

        if (Gcost <= FLINT_MIN(Abarcost, Bbarcost)*1.125)
            use |= USE_G;

        if (Abarcost <= FLINT_MIN(Gcost, Bbarcost)*1.125)
            use |= USE_ABAR;

        if (Bbarcost <= FLINT_MIN(Gcost, Abarcost)*1.125)
            use |= USE_BBAR;
    }

    if (use == 0)
        use = USE_G | USE_ABAR | USE_BBAR;

    return use;
}

mp_limb_t n_poly_mod_eval_step_sep(
    n_poly_t cur,
    const n_poly_t inc,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    return _nmod_zip_eval_step(cur->coeffs, inc->coeffs, A->coeffs, A->length, ctx->mod);
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
    _n_fq_zip_eval_step(res, cur->coeffs, inc->coeffs, A->coeffs, A->length, ctx);
}

void static n_bpoly_mod_eval_step_sep(
    n_bpoly_t E,
    n_polyun_t cur,
    const n_polyun_t inc,
    const nmod_mpoly_t A,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, Ai;
    slong e0, e1;
    mp_limb_t c;

    n_bpoly_zero(E);

    Ai = 0;
    for (i = 0; i < cur->length; i++)
    {
        slong this_len = cur->coeffs[i].length;

        c = _nmod_zip_eval_step(cur->coeffs[i].coeffs, inc->coeffs[i].coeffs,
                                  A->coeffs + Ai, this_len, ctx->mod);
        Ai += this_len;

        e0 = extract_exp(cur->exps[i], 1, 2);
        e1 = extract_exp(cur->exps[i], 0, 2);
        if (c == 0)
            continue;

        n_bpoly_set_coeff_nonzero(E, e0, e1, c);
    }

    FLINT_ASSERT(Ai == A->length);
}



static void nmod_mpoly_monomial_evals(
    n_poly_t E,
    const nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong start,
    slong stop,
    const nmod_mpoly_ctx_t ctx)
{
    _nmod_mpoly_monomial_evals_cache(E, A->exps, A->bits, A->length,
                                     betas, start, stop, ctx->minfo, ctx->mod);
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

static void nmod_mpoly_monomial_evals2(
    n_polyun_t E,
    const nmod_mpoly_t A,
    const mp_limb_t * betas,
    slong m,
    const nmod_mpoly_ctx_t ctx)
{
    _nmod_mpoly_monomial_evals2_cache(E, A->exps, A->bits, A->length, betas, m,
                                                         ctx->minfo, ctx->mod);
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

    The degrees of A, B, and gamma wrt the minor vars must be passed in.
    A guess of the degrees of rG wrt the minor vars can be passed in.


    deg(G) = deg(gamma) - deg(lc(rG)) +  deg(rG)
    deg(G) <= deg(gamma) + deg(rG)
    deg(G) <= deg(gamma) + deg(A)
    deg(G) <= deg(gamma) + deg(B)
    deg(G) >= deg(gamma)
    deg(G) >= deg(rG)

    deg(A) = deg(gamma) + deg(A) - deg(G)
    deg(B) = deg(gamma) + deg(B) - deg(G)
*/
int nmod_mpolyl_gcd_zippel_smprime(
    nmod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpoly_t rAbar,
    nmod_mpoly_t rBbar,
    const nmod_mpoly_t A, const slong * Adegs,
    const nmod_mpoly_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t ctx)
{
    int success, use, alpha_tries_left;
    slong i, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    mp_limb_t * alphas, * betas;
    flint_rand_t state;
    nmod_mpoly_t cont;
    nmod_mpoly_t T, G, Abar, Bbar;
    n_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    mp_limb_t gammaev;
    nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_poly_t gammaeh_cur, gammaeh_inc;
    n_poly_t modulus, alphapow;
    nmod_mpoly_struct * Aevals, * Bevals;
    nmod_mpoly_struct * gammaevals;
    n_poly_bpoly_stack_t St;
    mp_limb_t c, start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    nmod_mpoly_fit_length_reset_bits(rG, 1, bits, ctx);
    nmod_mpoly_fit_length_reset_bits(rAbar, 1, bits, ctx);
    nmod_mpoly_fit_length_reset_bits(rBbar, 1, bits, ctx);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        nmod_mpoly_degrees_si(tmp_degs, A, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        nmod_mpoly_degrees_si(tmp_degs, B, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        nmod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    FLINT_ASSERT(gammadegs[0] == 0);
    FLINT_ASSERT(gammadegs[1] == 0);

    if (ctx->mod.n < 7)
        return 0;

    n_polyun_init(HG);
    n_polyun_init(HAbar);
    n_polyun_init(HBbar);
    n_polyun_init(MG);
    n_polyun_init(MAbar);
    n_polyun_init(MBbar);
    n_polyun_init(ZG);
    n_polyun_init(ZAbar);
    n_polyun_init(ZBbar);
    n_bpoly_init(Aev);
    n_bpoly_init(Bev);
    n_bpoly_init(Gev);
    n_bpoly_init(Abarev);
    n_bpoly_init(Bbarev);
    n_poly_init2(alphapow, 4);
    nmod_mpoly_init3(cont, 1, bits, ctx);
    nmod_mpoly_init3(T, 1, bits, ctx);
    nmod_mpoly_init3(G, 1, bits, ctx);
    nmod_mpoly_init3(Abar, 1, bits, ctx);
    nmod_mpoly_init3(Bbar, 1, bits, ctx);
    nmod_mpolyn_init(Tn, bits, ctx);
    nmod_mpolyn_init(Gn, bits, ctx);
    nmod_mpolyn_init(Abarn, bits, ctx);
    nmod_mpolyn_init(Bbarn, bits, ctx);
    n_polyun_init(Aeh_cur);
    n_polyun_init(Aeh_inc);
    n_polyun_init(Beh_cur);
    n_polyun_init(Beh_inc);
    n_poly_init(gammaeh_cur);
    n_poly_init(gammaeh_inc);
    n_poly_init(modulus);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    betas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    alphas = FLINT_ARRAY_ALLOC(nvars, mp_limb_t);
    flint_randinit(state);

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        nmod_mpoly_init3(Aevals + i, 0, bits, ctx);
        nmod_mpoly_init3(Bevals + i, 0, bits, ctx);
        nmod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    Abideg = _nmod_mpoly_bidegree(A, ctx);
    Bbideg = _nmod_mpoly_bidegree(B, ctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    alpha_tries_left = 20;

choose_alphas:

    if (--alpha_tries_left < 0)
    {
        success = 0;
        goto cleanup;
    }

    for (i = 2; i < nvars; i++)
        alphas[i] = n_urandint(state, ctx->mod.n - 2) + 1;

    for (i = nvars - 1; i >= 2; i--)
    {
        nmod_mpoly_evaluate_one_ui(Aevals + i, Aevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + i, bits, ctx);
        nmod_mpoly_evaluate_one_ui(Bevals + i, Bevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(Bevals + i, bits, ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + i, gammaevals + i + 1, i, alphas[i], ctx);
        nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, ctx);
        if (nmod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || _nmod_mpoly_bidegree(Aevals + i, ctx) != Abideg)
            goto choose_alphas;
        if (Bevals[i].length < 1 || _nmod_mpoly_bidegree(Bevals + i, ctx) != Bbideg)
            goto choose_alphas;
    }

    m = 2;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    nmod_mpoly_get_bpoly(Aev, Aevals + m, 0, 1, ctx);
    nmod_mpoly_get_bpoly(Bev, Bevals + m, 0, 1, ctx);

    success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->mod, St);
    if (!success)
        goto cleanup;

    newdegXY = n_bpoly_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_alphas;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    gammaev = nmod_mpoly_get_ui(gammaevals + m, ctx);
    n_bpoly_scalar_mul_nmod(Gev, gammaev, ctx->mod);

    nmod_mpolyn_interp_lift_sm_bpoly(Gn, Gev, ctx);
    nmod_mpolyn_interp_lift_sm_bpoly(Abarn, Abarev, ctx);
    nmod_mpolyn_interp_lift_sm_bpoly(Bbarn, Bbarev, ctx);

    n_poly_one(modulus);
    c = nmod_neg(alphas[m], ctx->mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, ctx->mod);

    start_alpha = alphas[m];
    while (1)
    {
    choose_alpha_2:

        alphas[m] = (alphas[m] < 2) ? ctx->mod.n - 1 : alphas[m] - 1;
        if (alphas[m] == start_alpha)
            goto choose_alphas;

        FLINT_ASSERT(alphapow->alloc >= 2);
        alphapow->coeffs[0] = 1;
        alphapow->coeffs[1] = alphas[m];
        alphapow->length = 2;

        nmod_mpoly_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
        nmod_mpoly_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
        nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);
        nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);
        if (nmod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
            goto choose_alpha_2;

        nmod_mpoly_get_bpoly(Aev, Aevals + m, 0, 1, ctx);
        nmod_mpoly_get_bpoly(Bev, Bevals + m, 0, 1, ctx);

        success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev, ctx->mod, St);
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
            goto choose_alphas;
        }

        gammaev = nmod_mpoly_get_ui(gammaevals + m, ctx);
        n_bpoly_scalar_mul_nmod(Gev, gammaev, ctx->mod);

        c = n_poly_mod_eval_pow(modulus, alphapow, ctx->mod);
        c = nmod_inv(c, ctx->mod);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, ctx->mod);

        if ((use & USE_G) && !nmod_mpolyn_interp_crt_sm_bpoly(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                success = nmod_mpolyl_content(cont, rG, 2, ctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rG, rG, cont, ctx);
                nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                if (nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }
            }
            else
            {
                nmod_mpoly_cvtfrom_mpolyn(G, Gn, m, ctx);
                nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpoly_divides(Abar, T, G, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                    nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                success = nmod_mpolyl_content(cont, rAbar, 2, ctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                if (nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                nmod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, ctx);
                nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpoly_divides(G, T, Abar, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_BBAR) && !nmod_mpolyn_interp_crt_sm_bpoly(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                success = nmod_mpolyl_content(cont, rBbar, 2, ctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                if (nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                    nmod_mpoly_divides(rAbar, A, rG, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                nmod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, ctx);
                nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (nmod_mpoly_divides(G, T, Bbar, ctx))
                {
                    nmod_mpoly_repack_bits_inplace(G, bits, ctx);
                    nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(Abar, T, G, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if (n_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        c = nmod_neg(alphas[m], ctx->mod);
        n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, ctx->mod);
    }

    for (m = 3; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fp[gen(0), ..., gen(m - 1)] */
        nmod_mpolyn_interp_lift_sm_mpoly(Gn, G, ctx);
        nmod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, ctx);
        nmod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_poly_one(modulus);
        c = nmod_neg(alphas[m], ctx->mod);
        n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, ctx->mod);

        start_alpha = alphas[m];

    choose_betas:

        /* only beta[2], beta[1], ..., beta[m - 1] will be used */
        for (i = 2; i < ctx->minfo->nvars; i++)
            betas[i] = n_urandint(state, ctx->mod.n - 3) + 2;

        nmod_mpoly_monomial_evals2(HG, G, betas + 2, m, ctx);
        nmod_mpoly_monomial_evals2(HAbar, Abar, betas + 2, m, ctx);
        nmod_mpoly_monomial_evals2(HBbar, Bbar, betas + 2, m, ctx);

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
            this_length = n_polyun_product_roots(MG, HG, ctx->mod);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = n_polyun_product_roots(MAbar, HAbar, ctx->mod);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = n_polyun_product_roots(MBbar, HBbar, ctx->mod);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            alphas[m] = (alphas[m] < 2) ? ctx->mod.n - 1 : alphas[m] - 1;
            if (alphas[m] == start_alpha)
                goto choose_alphas;

            FLINT_ASSERT(alphapow->alloc >= 2);
            alphapow->coeffs[0] = 1;
            alphapow->coeffs[1] = alphas[m];
            alphapow->length = 2;

            nmod_mpoly_evaluate_one_ui(Aevals + m, Aevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_repack_bits_inplace(Aevals + m, bits, ctx);
            nmod_mpoly_evaluate_one_ui(Bevals + m, Bevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_repack_bits_inplace(Bevals + m, bits, ctx);
            nmod_mpoly_evaluate_one_ui(gammaevals + m, gammaevals + m + 1, m, alphas[m], ctx);
            nmod_mpoly_repack_bits_inplace(gammaevals + m, bits, ctx);
            if (nmod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _nmod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _nmod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
                goto choose_alpha_m;

            nmod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, betas + 2, m, ctx);
            nmod_mpoly_monomial_evals2(Beh_inc, Bevals + m, betas + 2, m, ctx);
            nmod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, betas + 2, 2, m, ctx);
            n_polyun_set(Aeh_cur, Aeh_inc);
            n_polyun_set(Beh_cur, Beh_inc);
            n_poly_set(gammaeh_cur, gammaeh_inc);

            n_polyun_zip_start(ZG, HG, req_zip_images);
            n_polyun_zip_start(ZAbar, HAbar, req_zip_images);
            n_polyun_zip_start(ZBbar, HBbar, req_zip_images);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                n_bpoly_mod_eval_step_sep(Aev, Aeh_cur, Aeh_inc, Aevals + m, ctx);
                n_bpoly_mod_eval_step_sep(Bev, Beh_cur, Beh_inc, Bevals + m, ctx);
                gammaev = n_poly_mod_eval_step_sep(gammaeh_cur, gammaeh_inc, gammaevals + m, ctx);
                if (gammaev == 0)
                    goto choose_betas;
                if (Aev->length < 1 || n_bpoly_bidegree(Aev) != Abideg)
                    goto choose_betas;
                if (Bev->length < 1 || n_bpoly_bidegree(Bev) != Bbideg)
                    goto choose_betas;

                success = n_bpoly_mod_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                            Aev, Bev, ctx->mod, St);        
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
                    goto choose_alphas;
                }

                n_bpoly_scalar_mul_nmod(Gev, gammaev, ctx->mod);
                if ((use & USE_G) && !n_polyu2n_add_zip_must_match(ZG, Gev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_polyu2n_add_zip_must_match(ZAbar, Abarev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_polyu2n_add_zip_must_match(ZBbar, Bbarev, cur_zip_image))
                    goto choose_alphas;
            }

            if ((use & USE_G) && n_polyun_zip_solve(G, ZG, HG, MG, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && n_polyun_zip_solve(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && n_polyun_zip_solve(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_poly_degree(modulus) > 0);
            c = n_poly_mod_eval_pow(modulus, alphapow, ctx->mod);
            c = nmod_inv(c, ctx->mod);
            _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, ctx->mod);

            if ((use & USE_G) && !nmod_mpolyn_interp_mcrt_sm_mpoly(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyl_content(cont, rG, 2, ctx);
                    if (!success)
                        goto cleanup;
                    nmod_mpoly_divides(rG, rG, cont, ctx);
                    nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    if (nmod_mpoly_divides(rAbar, A, rG, ctx) &&
                        nmod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(rAbar, T, rG, ctx))
                    {
                        nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            nmod_mpoly_swap(G, rG, ctx);
                            nmod_mpoly_swap(Abar, rAbar, ctx);
                            nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyl_content(cont, rAbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    nmod_mpoly_divides(rAbar, rAbar, cont, ctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    if (nmod_mpoly_divides(rG, A, rAbar, ctx) &&
                        nmod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(rG, T, rAbar, ctx))
                    {
                        nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            nmod_mpoly_swap(G, rG, ctx);
                            nmod_mpoly_swap(Abar, rAbar, ctx);
                            nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !nmod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = nmod_mpolyl_content(cont, rBbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    nmod_mpoly_divides(rBbar, rBbar, cont, ctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    if (nmod_mpoly_divides(rG, B, rBbar, ctx) &&
                        nmod_mpoly_divides(rAbar, A, rG, ctx))
                    {
                        nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    nmod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (nmod_mpoly_divides(rG, T, rBbar, ctx))
                    {
                        nmod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (nmod_mpoly_divides(rAbar, T, rG, ctx))
                        {
                            nmod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            nmod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            nmod_mpoly_swap(G, rG, ctx);
                            nmod_mpoly_swap(Abar, rAbar, ctx);
                            nmod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            if (n_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            c = nmod_neg(alphas[m], ctx->mod);
            n_poly_mod_shift_left_scalar_addmul(modulus, 1, c, ctx->mod);
        }
    }

    success = 1;

cleanup:

    n_polyun_clear(HG);
    n_polyun_clear(HAbar);
    n_polyun_clear(HBbar);
    n_polyun_clear(MG);
    n_polyun_clear(MAbar);
    n_polyun_clear(MBbar);
    n_polyun_clear(ZG);
    n_polyun_clear(ZAbar);
    n_polyun_clear(ZBbar);
    n_bpoly_clear(Aev);
    n_bpoly_clear(Bev);
    n_bpoly_clear(Gev);
    n_bpoly_clear(Abarev);
    n_bpoly_clear(Bbarev);
    n_poly_clear(alphapow);
    nmod_mpoly_clear(cont, ctx);
    nmod_mpoly_clear(T, ctx);
    nmod_mpoly_clear(G, ctx);
    nmod_mpoly_clear(Abar, ctx);
    nmod_mpoly_clear(Bbar, ctx);
    nmod_mpolyn_clear(Tn, ctx);
    nmod_mpolyn_clear(Gn, ctx);
    nmod_mpolyn_clear(Abarn, ctx);
    nmod_mpolyn_clear(Bbarn, ctx);
    n_polyun_clear(Aeh_cur);
    n_polyun_clear(Aeh_inc);
    n_polyun_clear(Beh_cur);
    n_polyun_clear(Beh_inc);
    n_poly_clear(gammaeh_cur);
    n_poly_clear(gammaeh_inc);
    n_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    flint_free(betas);
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        nmod_mpoly_clear(Aevals + i, ctx);
        nmod_mpoly_clear(Bevals + i, ctx);
        nmod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

    return success;

gcd_is_trivial:

    nmod_mpoly_one(rG, ctx);
    nmod_mpoly_set(rAbar, A, ctx);
    nmod_mpoly_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}


int nmod_mpolyl_gcd_zippel_lgprime(
    nmod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    nmod_mpoly_t rAbar,
    nmod_mpoly_t rBbar,
    const nmod_mpoly_t A, const slong * Adegs,
    const nmod_mpoly_t B, const slong * Bdegs,
    const nmod_mpoly_t gamma, const slong * gammadegs,
    const nmod_mpoly_ctx_t smctx)
{
    slong lgd;
    int success, use;
    slong i, m;
    slong nvars = smctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fq_nmod_struct * alphas, * betas;
    flint_rand_t state;
    nmod_mpoly_t cont;
    nmod_mpoly_t T, G, Abar, Bbar;
    nmod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    fq_nmod_mpoly_t qT, qG, qAbar, qBbar;
    fq_nmod_mpoly_t qrG, qrAbar, qrBbar;
    fq_nmod_mpolyn_t qTn, qGn, qAbarn, qBbarn;
    n_fq_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    n_fq_bpoly_t Aev, Bev, Gev, Abarev, Bbarev;
    const mp_limb_t * gammaev;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    n_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    n_fq_poly_t gammaeh_cur, gammaeh_inc;
    n_fq_poly_t alphapow;
    fq_nmod_mpoly_struct * Aevals, * Bevals;
    fq_nmod_mpoly_struct * gammaevals;
    n_poly_t modulus;
    n_poly_bpoly_stack_t St;
    n_poly_t tmp;  /* tmp arithmetic space */
    fq_nmod_t start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;
    fq_nmod_mpoly_ctx_t lgctx;
    nmod_mpolyn_t gamman;
    nmod_mpolyn_t An, Bn;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    nmod_mpoly_fit_length_reset_bits(rG, 1, bits, smctx);
    nmod_mpoly_fit_length_reset_bits(rAbar, 1, bits, smctx);
    nmod_mpoly_fit_length_reset_bits(rBbar, 1, bits, smctx);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        nmod_mpoly_degrees_si(tmp_degs, A, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        nmod_mpoly_degrees_si(tmp_degs, B, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        nmod_mpoly_degrees_si(tmp_degs, gamma, smctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    FLINT_ASSERT(gammadegs[0] == 0);
    FLINT_ASSERT(gammadegs[1] == 0);

    flint_randinit(state);

    lgd = WORD(20)/(FLINT_BIT_COUNT(smctx->mod.n));
    lgd = FLINT_MAX(WORD(2), lgd);
    fq_nmod_mpoly_ctx_init_deg(lgctx, nvars, ORD_LEX, smctx->mod.n, lgd);
    n_poly_init2(tmp, lgd);
    n_poly_init2(alphapow, 2*lgd);

    fq_nmod_init(start_alpha, lgctx->fqctx);
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
    nmod_mpoly_init3(cont, 1, bits, smctx);

    nmod_mpoly_init3(T, 1, bits, smctx);
    nmod_mpoly_init3(G, 1, bits, smctx);
    nmod_mpoly_init3(Abar, 1, bits, smctx);
    nmod_mpoly_init3(Bbar, 1, bits, smctx);
    nmod_mpolyn_init(Tn, bits, smctx);
    nmod_mpolyn_init(Gn, bits, smctx);
    nmod_mpolyn_init(Abarn, bits, smctx);
    nmod_mpolyn_init(Bbarn, bits, smctx);

    fq_nmod_mpoly_init3(qT, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qG, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qAbar, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qBbar, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qrG, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qrAbar, 1, bits, lgctx);
    fq_nmod_mpoly_init3(qrBbar, 1, bits, lgctx);
    fq_nmod_mpolyn_init(qTn, bits, lgctx);
    fq_nmod_mpolyn_init(qGn, bits, lgctx);
    fq_nmod_mpolyn_init(qAbarn, bits, lgctx);
    fq_nmod_mpolyn_init(qBbarn, bits, lgctx);

    n_fq_polyun_init(Aeh_cur);
    n_fq_polyun_init(Aeh_inc);
    n_fq_polyun_init(Beh_cur);
    n_fq_polyun_init(Beh_inc);
    n_fq_poly_init(gammaeh_cur);
    n_fq_poly_init(gammaeh_inc);

    n_poly_init(modulus);

    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);

    nmod_mpolyn_init(An, bits, smctx);
    nmod_mpolyn_init(Bn, bits, smctx);
    nmod_mpolyn_init(gamman, bits, smctx);

    /* alphas[nvars - 1] not used - it is replaced lgctx->fqctx->modulus */
    alphas = FLINT_ARRAY_ALLOC(2*nvars, fq_nmod_struct);
    betas = alphas + nvars;
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_init(betas + i, lgctx->fqctx);
        fq_nmod_init(alphas + i, lgctx->fqctx);
    }

    /* Aevals[nvars] does not exist - it is replaced by An */
    Aevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars, fq_nmod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_init3(Aevals + i, 0, bits, lgctx);
        fq_nmod_mpoly_init3(Bevals + i, 0, bits, lgctx);
        fq_nmod_mpoly_init3(gammaevals + i, 0, bits, lgctx);
    }
    nmod_mpoly_cvtto_mpolyn(An, A, nvars - 1, smctx);
    nmod_mpoly_cvtto_mpolyn(Bn, B, nvars - 1, smctx);
    nmod_mpoly_cvtto_mpolyn(gamman, gamma, nvars - 1, smctx);

    Abideg = _nmod_mpoly_bidegree(A, smctx);
    Bbideg = _nmod_mpoly_bidegree(B, smctx);

    degxAB = FLINT_MAX(Adegs[0], Bdegs[0]);
    degyAB = FLINT_MAX(Adegs[1], Bdegs[1]);

    GdegboundXY = pack_exp2(FLINT_MIN(Adegs[0], Bdegs[0]),
                            FLINT_MIN(Adegs[1], Bdegs[1]));
    if (GdegboundXY == 0)
        goto gcd_is_trivial;

    goto got_alpha_last;

increase_degree:

choose_alphas:

    /* TODO: don't necessarily increase degree here */
    lgd++;
    if (lgd > 10000)
    {
        /* ran out of primes */
        success = 0;
        goto cleanup;
    }

    fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
    n_poly_fit_length(tmp, lgd);
    n_poly_fit_length(alphapow, 2*lgd);

got_alpha_last:

    for (i = 2; i < nvars - 1; i++)
        fq_nmod_rand_not_zero(alphas + i, state, lgctx->fqctx);

    i = nvars - 1;
    nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + i, An, lgctx, smctx);
    nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + i, Bn, lgctx, smctx);
    nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + i, gamman, lgctx, smctx);
    if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
        goto choose_alphas;
    if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
        goto choose_alphas;
    if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
        goto choose_alphas;
    for (i--; i >= 2; i--)
    {
        fq_nmod_mpoly_evaluate_one_fq_nmod(Aevals + i, Aevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Aevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(Bevals + i, Bevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(Bevals + i, bits, lgctx);
        fq_nmod_mpoly_evaluate_one_fq_nmod(gammaevals + i, gammaevals + i + 1, i, alphas + i, lgctx);
        fq_nmod_mpoly_repack_bits_inplace(gammaevals + i, bits, lgctx);
        if (fq_nmod_mpoly_is_zero(gammaevals + i, lgctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Aevals + i, lgctx) != Abideg)
            goto choose_alphas;
        if (Bevals[i].length < 1 || _fq_nmod_mpoly_bidegree(Bevals + i, lgctx) != Bbideg)
            goto choose_alphas;
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
        goto choose_alphas;
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
        nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Gn, smctx, Gev, lgctx);
        nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Abarn, smctx, Abarev, lgctx);
        nmod_mpolyn_interp_lift_lg_bpoly(&lastdeg, Bbarn, smctx, Bbarev, lgctx);

        n_poly_set_nmod_poly(modulus, lgctx->fqctx->modulus);

        while (1)
        {
        choose_alpha_2_last:
            lgd++;
            if (lgd > 10000)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx);
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
                goto choose_alphas;
            }

            gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
            n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

            if ((use & USE_G) && !nmod_mpolyn_interp_crt_lg_bpoly(
                                 &lastdeg, Gn, Tn, modulus, smctx, Gev, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rG, rG, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !nmod_mpolyn_interp_crt_lg_bpoly(
                           &lastdeg, Abarn, Tn, modulus, smctx, Abarev, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !nmod_mpolyn_interp_crt_lg_bpoly(
                          &lastdeg, Bbarn, Tn, modulus, smctx, Bbarev, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            {
                n_poly_t h_mock;
                n_poly_mock(h_mock, lgctx->fqctx->modulus);
                n_poly_mod_mul(modulus, modulus, h_mock, smctx->mod);
            }
        }

        success = 1;
        goto cleanup;
    }

    fq_nmod_mpolyn_interp_lift_sm_bpoly(qGn, Gev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(qAbarn, Abarev, lgctx);
    fq_nmod_mpolyn_interp_lift_sm_bpoly(qBbarn, Bbarev, lgctx);

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
            goto choose_alphas;
        }

        gammaev = fq_nmod_mpoly_get_nonzero_n_fq(gammaevals + m, lgctx);
        n_fq_bpoly_scalar_mul_n_fq(Gev, gammaev, lgctx->fqctx);

        FLINT_ASSERT(tmp->alloc >= lgd);
        n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
        n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

        if ((use & USE_G) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                              &lastdeg, qGn, qTn, Gev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(qG, qGn, m, lgctx);
            fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(qAbar, qT, qG, lgctx))
            {
                fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qBbar, qT, qG, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(qAbar, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(qBbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                      &lastdeg, qAbarn, qTn, Abarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(qAbar, qAbarn, m, lgctx);
            fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(qG, qT, qAbar, lgctx))
            {
                fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qBbar, qT, qG, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(qG, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(qBbar, bits, lgctx);
                    break;
                }
            }
        }

        if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_crt_sm_bpoly(
                      &lastdeg, qBbarn, qTn, Bbarev, modulus, alphapow, lgctx))
        {
            fq_nmod_mpoly_cvtfrom_mpolyn(qBbar, qBbarn, m, lgctx);
            fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
            if (fq_nmod_mpoly_divides(qG, qT, qBbar, lgctx))
            {
                fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qAbar, qT, qG, lgctx))
                {
                    fq_nmod_mpoly_repack_bits_inplace(qG, bits, lgctx);
                    fq_nmod_mpoly_repack_bits_inplace(qAbar, bits, lgctx);
                    break;
                }
            }
        }

        if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
            n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
        n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
    }

    for (m = 3; m < nvars - 1; m++)
    {
        /* qG, qAbar, qBbar are in Fq[gen(0), ..., gen(m - 1)] */
        fq_nmod_mpolyn_interp_lift_sm_mpoly(qGn, qG, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(qAbarn, qAbar, lgctx);
        fq_nmod_mpolyn_interp_lift_sm_mpoly(qBbarn, qBbar, lgctx);

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

        FLINT_ASSERT(qG->bits == bits);
        FLINT_ASSERT(qAbar->bits == bits);
        FLINT_ASSERT(qBbar->bits == bits);

        fq_nmod_mpoly_monomial_evals2(HG, qG, betas + 2, m, lgctx);
        fq_nmod_mpoly_monomial_evals2(HAbar, qAbar, betas + 2, m, lgctx);
        fq_nmod_mpoly_monomial_evals2(HBbar, qBbar, betas + 2, m, lgctx);
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
                    goto choose_alphas;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2n_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_fq_polyu2n_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_fq_polyu2n_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
            }

            if ((use & USE_G) && n_fq_polyun_zip_solve(qG, ZG, HG, MG, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && n_fq_polyun_zip_solve(qAbar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && n_fq_polyun_zip_solve(qBbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            FLINT_ASSERT(tmp->alloc >= lgd);
            n_fq_poly_eval_pow(tmp->coeffs, modulus, alphapow, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);
            n_fq_poly_scalar_mul_n_fq(modulus, modulus, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                                  &lastdeg, qGn, qG, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(qrG, qGn, m, lgctx);
                fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qrAbar, qT, qrG, lgctx))
                {
                    fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(qrBbar, qT, qrG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(qrAbar, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(qrBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(qG, qrG, lgctx);
                        fq_nmod_mpoly_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpoly_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_ABAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                              &lastdeg, qAbarn, qAbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(qrAbar, qAbarn, m, lgctx);
                fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qrG, qT, qrAbar, lgctx))
                {
                    fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(qrBbar, qT, qrG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(qrG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(qrBbar, bits, lgctx);
                        fq_nmod_mpoly_swap(qG, qrG, lgctx);
                        fq_nmod_mpoly_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpoly_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if ((use & USE_BBAR) && !fq_nmod_mpolyn_interp_mcrt_sm_mpoly(
                            &lastdeg, qBbarn, qBbar, modulus, alphapow, lgctx))
            {
                fq_nmod_mpoly_cvtfrom_mpolyn(qrBbar, qBbarn, m, lgctx);
                fq_nmod_mpoly_mul(qT, Bevals + m + 1, gammaevals + m + 1, lgctx);
                if (fq_nmod_mpoly_divides(qrG, qT, qrBbar, lgctx))
                {
                    fq_nmod_mpoly_mul(qT, Aevals + m + 1, gammaevals + m + 1, lgctx);
                    if (fq_nmod_mpoly_divides(qrAbar, qT, qrG, lgctx))
                    {
                        fq_nmod_mpoly_repack_bits_inplace(qrG, bits, lgctx);
                        fq_nmod_mpoly_repack_bits_inplace(qrAbar, bits, lgctx);
                        fq_nmod_mpoly_swap(qG, qrG, lgctx);
                        fq_nmod_mpoly_swap(qAbar, qrAbar, lgctx);
                        fq_nmod_mpoly_swap(qBbar, qrBbar, lgctx);
                        break;
                    }
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            FLINT_ASSERT(n_fq_equal_fq_nmod(alphapow->coeffs + lgd*1, alphas + m, lgctx->fqctx));
            n_fq_poly_shift_left_scalar_submul(modulus, 1, alphapow->coeffs + lgd*1, lgctx->fqctx);
        }
    }

    m = nvars - 1;
    {
        /* G, Abar, Bbar are in Fq/alpha(gen(m-1))[gen(0), ..., gen(m - 1)] */
        nmod_mpolyn_interp_lift_lg_mpoly(Gn, smctx, qG, lgctx);
        nmod_mpolyn_interp_lift_lg_mpoly(Abarn, smctx, qAbar, lgctx);
        nmod_mpolyn_interp_lift_lg_mpoly(Bbarn, smctx, qBbar, lgctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        n_poly_set_nmod_poly(modulus, lgctx->fqctx->modulus);

        while (1)
        {
        choose_alpha_last:

            lgd++;
            if (lgd > 10000)
            {
                /* ran out of primes */
                success = 0;
                goto cleanup;
            }

            fq_nmod_mpoly_ctx_change_modulus(lgctx, lgd);
            n_poly_fit_length(tmp, lgd);
            n_poly_fit_length(alphapow, 2*lgd);

            nmod_mpolyn_interp_reduce_lg_mpoly(Aevals + m, An, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(Bevals + m, Bn, lgctx, smctx);
            nmod_mpolyn_interp_reduce_lg_mpoly(gammaevals + m, gamman, lgctx, smctx);
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

            fq_nmod_mpoly_monomial_evals2(HG, qG, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(HAbar, qAbar, betas + 2, m, lgctx);
            fq_nmod_mpoly_monomial_evals2(HBbar, qBbar, betas + 2, m, lgctx);

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
                    goto choose_alphas;
                }

                n_fq_bpoly_scalar_mul_n_fq(Gev, tmp->coeffs, lgctx->fqctx);

                if ((use & USE_G) && !n_fq_polyu2n_add_zip_must_match(ZG, Gev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !n_fq_polyu2n_add_zip_must_match(ZAbar, Abarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !n_fq_polyu2n_add_zip_must_match(ZBbar, Bbarev, cur_zip_image, lgctx->fqctx))
                    goto choose_alphas;
            }

            if ((use & USE_G) && n_fq_polyun_zip_solve(qG, ZG, HG, MG, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && n_fq_polyun_zip_solve(qAbar, ZAbar, HAbar, MAbar, lgctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && n_fq_polyun_zip_solve(qBbar, ZBbar, HBbar, MBbar, lgctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(n_fq_poly_degree(modulus) > 0);

            n_poly_fit_length(tmp, lgd);
            _n_fq_set_n_poly(tmp->coeffs, modulus->coeffs, modulus->length, lgctx->fqctx);
            n_fq_inv(tmp->coeffs, tmp->coeffs, lgctx->fqctx);

            if ((use & USE_G) && !nmod_mpolyn_interp_mcrt_lg_mpoly(
                         &lastdeg, Gn, smctx, modulus, tmp->coeffs, qG, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rG, Gn, m, smctx);
                success = nmod_mpolyl_content(cont, rG, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rG, rG, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                if (nmod_mpoly_divides(rAbar, A, rG, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_ABAR) && !nmod_mpolyn_interp_mcrt_lg_mpoly(
                   &lastdeg, Abarn, smctx, modulus, tmp->coeffs, qAbar, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, smctx);
                success = nmod_mpolyl_content(cont, rAbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rAbar, rAbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                if (nmod_mpoly_divides(rG, A, rAbar, smctx) &&
                    nmod_mpoly_divides(rBbar, B, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                    break;
                }
            }

            if ((use & USE_BBAR) && !nmod_mpolyn_interp_mcrt_lg_mpoly(
                   &lastdeg, Bbarn, smctx, modulus, tmp->coeffs, qBbar, lgctx))
            {
                nmod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, smctx);
                success = nmod_mpolyl_content(cont, rBbar, 2, smctx);
                if (!success)
                    goto cleanup;
                nmod_mpoly_divides(rBbar, rBbar, cont, smctx);
                nmod_mpoly_repack_bits_inplace(rBbar, bits, smctx);
                if (nmod_mpoly_divides(rG, B, rBbar, smctx) &&
                    nmod_mpoly_divides(rAbar, A, rG, smctx))
                {
                    nmod_mpoly_repack_bits_inplace(rG, bits, smctx);
                    nmod_mpoly_repack_bits_inplace(rAbar, bits, smctx);
                    break;
                }
            }

            if (n_fq_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                n_fq_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            {
                n_poly_t h_mock;
                n_poly_mock(h_mock, lgctx->fqctx->modulus);
                n_poly_mod_mul(modulus, modulus, h_mock, smctx->mod);
            }
        }
    }

    success = 1;

cleanup:

    nmod_mpolyn_clear(An, smctx);
    nmod_mpolyn_clear(Bn, smctx);
    nmod_mpolyn_clear(gamman, smctx);

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
    nmod_mpoly_clear(cont, smctx);
    nmod_mpoly_clear(T, smctx);
    nmod_mpoly_clear(G, smctx);
    nmod_mpoly_clear(Abar, smctx);
    nmod_mpoly_clear(Bbar, smctx);
    nmod_mpolyn_clear(Tn, smctx);
    nmod_mpolyn_clear(Gn, smctx);
    nmod_mpolyn_clear(Abarn, smctx);
    nmod_mpolyn_clear(Bbarn, smctx);
    fq_nmod_mpoly_clear(qT, lgctx);
    fq_nmod_mpoly_clear(qG, lgctx);
    fq_nmod_mpoly_clear(qAbar, lgctx);
    fq_nmod_mpoly_clear(qBbar, lgctx);
    fq_nmod_mpoly_clear(qrG, lgctx);
    fq_nmod_mpoly_clear(qrAbar, lgctx);
    fq_nmod_mpoly_clear(qrBbar, lgctx);
    fq_nmod_mpolyn_clear(qTn, lgctx);
    fq_nmod_mpolyn_clear(qGn, lgctx);
    fq_nmod_mpolyn_clear(qAbarn, lgctx);
    fq_nmod_mpolyn_clear(qBbarn, lgctx);
    n_fq_polyun_clear(Aeh_cur);
    n_fq_polyun_clear(Aeh_inc);
    n_fq_polyun_clear(Beh_cur);
    n_fq_polyun_clear(Beh_inc);
    n_fq_poly_clear(gammaeh_cur);
    n_fq_poly_clear(gammaeh_inc);
    n_poly_clear(modulus);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);

    fq_nmod_clear(start_alpha, lgctx->fqctx);

    n_poly_clear(tmp);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_clear(alphas + i, lgctx->fqctx);
        fq_nmod_clear(betas + i, lgctx->fqctx);
    }
    flint_free(alphas);
    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fq_nmod_mpoly_clear(Aevals + i, lgctx);
        fq_nmod_mpoly_clear(Bevals + i, lgctx);
        fq_nmod_mpoly_clear(gammaevals + i, lgctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    fq_nmod_mpoly_ctx_clear(lgctx);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

    return success;

gcd_is_trivial:

    nmod_mpoly_one(rG, smctx);
    nmod_mpoly_set(rAbar, A, smctx);
    nmod_mpoly_set(rBbar, B, smctx);

    success = 1;
    
    goto cleanup;
}


/* should find its way back here in interesting cases */
int nmod_mpoly_gcd_zippel2(
    nmod_mpoly_t G,
    const nmod_mpoly_t A,
    const nmod_mpoly_t B,
    const nmod_mpoly_ctx_t ctx)
{
    if (nmod_mpoly_is_zero(A, ctx) || nmod_mpoly_is_zero(B, ctx))
        return nmod_mpoly_gcd(G, A, B, ctx);

    return _nmod_mpoly_gcd_algo(G, NULL, NULL, A, B, ctx, MPOLY_GCD_USE_ZIPPEL2);
}

