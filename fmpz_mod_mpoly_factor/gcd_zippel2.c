/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mpoly.h"
#include "fmpz_mod_mpoly_factor.h"


void fmpz_mod_mpolyn_interp_lift_sm_polyu1n(
    fmpz_mod_mpolyn_t F,
    fmpz_mod_polyun_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp_sp(F->bits, ctx->minfo);
    slong i, j, Fi;
    slong off0, shift0, off1, shift1;

    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    Fi = 0;
    for (i = 0; i < A->length; i++)
    {
        fmpz * Aicoeffs = A->coeffs[i].coeffs;
        ulong e0 = A->exps[i] << shift0;

        for (j = A->coeffs[i].length - 1; j >= 0; j--)
        {
            if (fmpz_is_zero(Aicoeffs + j))
                continue;

            fmpz_mod_mpolyn_fit_length(F, Fi + 1, ctx);
            mpoly_monomial_zero(F->exps + N*Fi, N);
            (F->exps + N*Fi)[off0] = e0;
            (F->exps + N*Fi)[off1] += (j << shift1);
            fmpz_mod_poly_set_fmpz(F->coeffs + Fi, Aicoeffs + j, ctx->ffinfo);
            Fi++;
        }
    }

    F->length = Fi;
}

int fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
    slong * lastdeg,
    fmpz_mod_mpolyn_t F,
    fmpz_mod_mpolyn_t T,
    fmpz_mod_polyun_t A,
    fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int changed = 0;
    slong N = mpoly_words_per_exp(F->bits, ctx->minfo);
    slong off0, shift0, off1, shift1;
    fmpz_mod_poly_struct * Acoeffs = A->coeffs;
    ulong * Aexps = A->exps;
    slong Fi, Ti, Ai, ai;
    slong Alen = A->length;
    slong Flen = F->length;
    ulong * Fexps = F->exps;
    fmpz_mod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Texps = T->exps;
    fmpz_mod_poly_struct * Tcoeffs = T->coeffs;
    fmpz_t v;
    ulong Fexpi, mask;

    fmpz_init(v);

    mask = (-UWORD(1)) >> (FLINT_BITS - F->bits);
    mpoly_gen_offset_shift_sp(&off0, &shift0, 0, F->bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off1, &shift1, 1, F->bits, ctx->minfo);

    FLINT_ASSERT(T->bits == F->bits);

    *lastdeg = -1;

    Ti = Fi = 0;
    Ai = 0;
    ai = 0;
    if (Ai < Alen)
        ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx->ffinfo);

    while (Fi < Flen || Ai < Alen)
    {
        if (Ti >= T->alloc)
        {
            slong extra = FLINT_MAX(Flen - Fi, Alen - Ai);
            fmpz_mod_mpolyn_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        if (Fi < Flen)
            Fexpi = pack_exp2(((Fexps + N*Fi)[off0]>>shift0)&mask,
                              ((Fexps + N*Fi)[off1]>>shift1)&mask);
        else
            Fexpi = 0;

        if (Fi < Flen && Ai < Alen && Fexpi == pack_exp2(Aexps[Ai], ai))
        {
            /* F term ok, A term ok */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            fmpz_mod_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->ffinfo);
            fmpz_mod_sub(v, Acoeffs[Ai].coeffs + ai, v, ctx->ffinfo);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeffs + Ti,
                                        Fcoeffs + Fi, modulus, v, ctx->ffinfo);
            Fi++;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx->ffinfo);
            }
        }
        else if (Ai < Alen && (Fi >= Flen || Fexpi < pack_exp2(Aexps[Ai], ai)))
        {
            /* F term missing, A term ok */
            mpoly_monomial_zero(Texps + N*Ti, N);
            (Texps + N*Ti)[off0] += (Ai << shift0);
            (Texps + N*Ti)[off1] += (ai << shift1);

            changed = 1;
            fmpz_mod_poly_scalar_mul_fmpz(Tcoeffs + Ti, modulus,
                                         Acoeffs[Ai].coeffs + ai, ctx->ffinfo);

            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx->ffinfo);
            }
        }
        else
        {
            FLINT_ASSERT(Fi < Flen && (Ai >= Alen || Fexpi > pack_exp2(Aexps[Ai], ai)));
            /* F term ok, Aterm missing */
            mpoly_monomial_set(Texps + N*Ti, Fexps + N*Fi, N);

            fmpz_mod_poly_eval_pow(v, Fcoeffs + Fi, alphapow, ctx->ffinfo);
            fmpz_mod_neg(v, v, ctx->ffinfo);
            changed |= !fmpz_is_zero(v);
            fmpz_mod_poly_scalar_addmul_fmpz_mod(Tcoeffs + Ti,
                                        Fcoeffs + Fi, modulus, v, ctx->ffinfo);
            Fi++;
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeffs + Ti, ctx->ffinfo));
        *lastdeg = FLINT_MAX(*lastdeg, fmpz_mod_poly_degree(Tcoeffs + Ti, ctx->ffinfo));
        Ti++;
    }

    T->length = Ti;

    if (changed)
        fmpz_mod_mpolyn_swap(T, F, ctx);

    fmpz_clear(v);

    return changed;
}


static ulong _fmpz_mod_mpoly_bidegree(
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length > 0);
    return _mpoly_bidegree(A->exps, A->bits, ctx->minfo);
}

int fmpz_mod_polyun_zip_solve(
    fmpz_mod_mpoly_t A,
    fmpz_mod_polyun_t Z,
    fmpz_mod_polyun_t H,
    fmpz_mod_polyun_t M,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success;
    slong Ai, i, n;
    fmpz * Acoeffs = A->coeffs;
    fmpz_mod_poly_t t;

    fmpz_mod_poly_init(t, ctx->ffinfo);

    FLINT_ASSERT(Z->length == H->length);
    FLINT_ASSERT(Z->length == M->length);

    Ai = 0;
    for (i = 0; i < H->length; i++)
    {
        n = H->coeffs[i].length;
        FLINT_ASSERT(M->coeffs[i].length == n + 1);
        FLINT_ASSERT(Z->coeffs[i].length >= n);
        FLINT_ASSERT(Ai + n <= A->length);

        fmpz_mod_poly_fit_length(t, n, ctx->ffinfo);

        success = _fmpz_mod_zip_vand_solve(Acoeffs + Ai,
                         H->coeffs[i].coeffs, n,
                         Z->coeffs[i].coeffs, Z->coeffs[i].length,
                         M->coeffs[i].coeffs, t->coeffs, ctx->ffinfo);
        if (success < 1)
        {
            fmpz_mod_poly_clear(t, ctx->ffinfo);
            return success;
        }

        Ai += n;
        FLINT_ASSERT(Ai <= A->length);

    }

    FLINT_ASSERT(Ai == A->length);

    fmpz_mod_poly_clear(t, ctx->ffinfo);
    return 1;
}


#define USE_G    1
#define USE_ABAR 2
#define USE_BBAR 4

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

int fmpz_mod_mpoly_gcd_get_use_new(
    slong rGdeg,
    slong Adeg,
    slong Bdeg,
    slong gammadeg,
    slong degxAB,
    slong degyAB,
    slong numABgamma,
    const fmpz_mod_polyun_t G,
    const fmpz_mod_polyun_t Abar,
    const fmpz_mod_polyun_t Bbar)
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


static void fmpz_mod_poly_zip_eval_cur_inc_coeff(
    fmpz_t eval,
    fmpz_mod_poly_t cur,
    const fmpz_mod_poly_t inc,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(A->length == cur->length);
    FLINT_ASSERT(A->length == inc->length);
    _fmpz_mod_zip_eval_step(eval, cur->coeffs, inc->coeffs, A->coeffs, A->length, ctx->ffinfo);
}


void fmpz_mod_mpoly_mock_eval_coeff(
    fmpz_mod_polyun_t mock,
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_polyun_t Aeh_inc,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, k;

    if (mock->alloc < Aeh_inc->length)
    {
        mock->alloc = FLINT_MAX(Aeh_inc->length, mock->alloc + mock->alloc/2);
        mock->coeffs = FLINT_ARRAY_REALLOC(mock->coeffs, mock->alloc,
                                                         fmpz_mod_poly_struct);        
    }

    mock->length = Aeh_inc->length;

    k = 0;
    for (i = 0; i < Aeh_inc->length; i++)
    {
        slong l = Aeh_inc->coeffs[i].length;
        mock->coeffs[i].coeffs = A->coeffs + k;
        mock->coeffs[i].alloc = l;
        mock->coeffs[i].length = l;
        k += l;
    }

    FLINT_ASSERT(k == A->length);
}


static void fmpz_mod_mpoly_monomial_evals(
    fmpz_mod_poly_t E,
    const fmpz_mod_mpoly_t A,
    fmpz_mod_poly_struct * beta_caches,
    slong start,
    slong stop,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    mpoly_monomial_evals_fmpz_mod(E, A->exps, A->length, A->bits,
                            beta_caches, start, stop, ctx->minfo, ctx->ffinfo);
}

static void fmpz_mod_mpoly_monomial_evals2(
    fmpz_mod_polyun_t E,
    const fmpz_mod_mpoly_t A,
    fmpz_mod_poly_struct * betas,
    slong m,
    const fmpz_mod_mpoly_ctx_t ctx,
    n_poly_t tmark) /* temp space */
{
    mpoly2_fill_marks(&tmark->coeffs, &tmark->length, &tmark->alloc,
                                      A->exps, A->length, A->bits, ctx->minfo);

    mpoly2_monomial_evals_fmpz_mod(E, A->exps, A->bits, tmark->coeffs,
                            tmark->length, betas, m, ctx->minfo, ctx->ffinfo);
}


static int fmpz_mod_polyun_add_zip_must_match(
    fmpz_mod_polyun_t Z,
    const fmpz_mod_polyun_t A,
    slong cur_length)
{
    slong i, Ai, ai;
    slong Alen = A->length;
    ulong * Zexps = Z->exps;
    fmpz_mod_poly_struct * Zcoeffs = Z->coeffs;
    ulong * Aexps = A->exps;
    fmpz_mod_poly_struct * Acoeffs = A->coeffs;

    Ai = 0;
    ai = 0;
    if (Ai < Alen)
        ai = Acoeffs[Ai].length - 1;

    for (i = 0; i < Z->length; i++)
    {
        if (Ai < Alen && Zexps[i] == pack_exp2(Aexps[Ai], ai))
        {
            /* Z present, A present */
            fmpz_set(Zcoeffs[i].coeffs + cur_length,
                     Acoeffs[Ai].coeffs + ai);
            Zcoeffs[i].length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
            if (ai < 0)
            {
                Ai++;
                if (Ai < Alen)
                    ai = Acoeffs[Ai].length - 1;
            }
        }
        else if (Ai >= Alen || Zexps[i] > pack_exp2(Aexps[Ai], ai))
        {
            /* Z present, A missing */
            fmpz_zero(Zcoeffs[i].coeffs + cur_length);
            Zcoeffs[i].length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return Ai >= Alen;
}


void fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(
    fmpz_mod_polyun_t E,
    fmpz_mod_polyun_t Acur,
    const fmpz_mod_polyun_t Ainc,
    const fmpz_mod_polyun_t Acoeff,
    const fmpz_mod_ctx_t ctx)
{
    slong i, Ei;
    slong e0, e1;
    fmpz_t c;

    FLINT_ASSERT(Acur->length > 0);
    FLINT_ASSERT(Acur->length == Ainc->length);
    FLINT_ASSERT(Acur->length == Acoeff->length);

    fmpz_init(c);

    e0 = extract_exp(Acur->exps[0], 1, 2);
    e1 = extract_exp(Acur->exps[0], 0, 2);

    fmpz_mod_polyun_fit_length(E, 4, ctx);
    Ei = 0;
    E->exps[Ei] = e1;
    fmpz_mod_poly_zero(E->coeffs + Ei, ctx);

    for (i = 0; i < Acur->length; i++)
    {
        slong this_len = Acur->coeffs[i].length;
        FLINT_ASSERT(this_len == Ainc->coeffs[i].length);
        FLINT_ASSERT(this_len == Acoeff->coeffs[i].length);

        _fmpz_mod_zip_eval_step(c, Acur->coeffs[i].coeffs,
              Ainc->coeffs[i].coeffs, Acoeff->coeffs[i].coeffs, this_len, ctx);

        e0 = extract_exp(Acur->exps[i], 1, 2);
        e1 = extract_exp(Acur->exps[i], 0, 2);

        if (E->exps[Ei] != e0)
        {
            fmpz_mod_polyun_fit_length(E, Ei + 2, ctx);
            Ei += !fmpz_mod_poly_is_zero(E->coeffs + Ei, ctx);
            E->exps[Ei] = e0;
            fmpz_mod_poly_zero(E->coeffs + Ei, ctx);
        }

        fmpz_mod_poly_set_coeff_fmpz(E->coeffs + Ei, e1, c, ctx);
    }

    Ei += !fmpz_mod_poly_is_zero(E->coeffs + Ei, ctx);
    E->length = Ei;

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(E, ctx));

    fmpz_clear(c);
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
int fmpz_mod_mpolyl_gcd_zippel2_smprime(
    fmpz_mod_mpoly_t rG, const slong * rGdegs, /* guess at rG degrees, could be NULL */
    fmpz_mod_mpoly_t rAbar,
    fmpz_mod_mpoly_t rBbar,
    const fmpz_mod_mpoly_t A, const slong * Adegs,
    const fmpz_mod_mpoly_t B, const slong * Bdegs,
    const fmpz_mod_mpoly_t gamma, const slong * gammadegs,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int success, use, alpha_tries_left;
    slong i, m;
    slong nvars = ctx->minfo->nvars;
    flint_bitcnt_t bits = A->bits;
    fmpz * alphas, * betas;
    fmpz_mod_poly_struct * alpha_caches, * beta_caches;
    flint_rand_t state;
    n_poly_t tmark;
    fmpz_mod_mpoly_t cont;
    fmpz_mod_mpoly_t T, G, Abar, Bbar;
    fmpz_mod_polyun_t HG, HAbar, HBbar, MG, MAbar, MBbar, ZG, ZAbar, ZBbar;
    fmpz_mod_polyun_t Aev, Bev, Gev, Abarev, Bbarev;
    fmpz_t gammaev;
    fmpz_mod_mpolyn_t Tn, Gn, Abarn, Bbarn;
    slong lastdeg;
    slong cur_zip_image, req_zip_images, this_length;
    fmpz_mod_polyun_t Aeh_coeff_mock, Beh_coeff_mock;
    fmpz_mod_polyun_t Aeh_cur, Aeh_inc, Beh_cur, Beh_inc;
    fmpz_mod_poly_t gammaeh_cur, gammaeh_inc;
    fmpz_mod_poly_t modulus, alphapow;
    fmpz_mod_mpoly_struct * Aevals, * Bevals;
    fmpz_mod_mpoly_struct * gammaevals;
    fmpz_mod_poly_polyun_stack_t St;
    fmpz_t c, start_alpha;
    ulong GdegboundXY, newdegXY, Abideg, Bbideg;
    slong degxAB, degyAB;

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(bits == A->bits);
    FLINT_ASSERT(bits == B->bits);
    FLINT_ASSERT(bits == gamma->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(gamma->length > 0);

    fmpz_mod_mpoly_fit_length_reset_bits(rG, 1, bits, ctx);
    fmpz_mod_mpoly_fit_length_reset_bits(rAbar, 1, bits, ctx);
    fmpz_mod_mpoly_fit_length_reset_bits(rBbar, 1, bits, ctx);

    fmpz_init(gammaev);
    fmpz_init(c);
    fmpz_init(start_alpha);

#if FLINT_WANT_ASSERT
    {
        slong * tmp_degs = FLINT_ARRAY_ALLOC(nvars, slong);

        fmpz_mod_mpoly_degrees_si(tmp_degs, A, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Adegs[i]);

        fmpz_mod_mpoly_degrees_si(tmp_degs, B, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == Bdegs[i]);

        fmpz_mod_mpoly_degrees_si(tmp_degs, gamma, ctx);
        for (i = 0; i < nvars; i++)
            FLINT_ASSERT(tmp_degs[i] == gammadegs[i]);

        flint_free(tmp_degs);
    }
#endif

    FLINT_ASSERT(gammadegs[0] == 0);
    FLINT_ASSERT(gammadegs[1] == 0);

    n_poly_init(tmark);

    fmpz_mod_polyun_init(HG, ctx->ffinfo);
    fmpz_mod_polyun_init(HAbar, ctx->ffinfo);
    fmpz_mod_polyun_init(HBbar, ctx->ffinfo);
    fmpz_mod_polyun_init(MG, ctx->ffinfo);
    fmpz_mod_polyun_init(MAbar, ctx->ffinfo);
    fmpz_mod_polyun_init(MBbar, ctx->ffinfo);
    fmpz_mod_polyun_init(ZG, ctx->ffinfo);
    fmpz_mod_polyun_init(ZAbar, ctx->ffinfo);
    fmpz_mod_polyun_init(ZBbar, ctx->ffinfo);
    fmpz_mod_polyun_init(Aev, ctx->ffinfo);
    fmpz_mod_polyun_init(Bev, ctx->ffinfo);
    fmpz_mod_polyun_init(Gev, ctx->ffinfo);
    fmpz_mod_polyun_init(Abarev, ctx->ffinfo);
    fmpz_mod_polyun_init(Bbarev, ctx->ffinfo);
    fmpz_mod_poly_init2(alphapow, 4, ctx->ffinfo);
    fmpz_mod_mpoly_init3(cont, 1, bits, ctx);
    fmpz_mod_mpoly_init3(T, 1, bits, ctx);
    fmpz_mod_mpoly_init3(G, 1, bits, ctx);
    fmpz_mod_mpoly_init3(Abar, 1, bits, ctx);
    fmpz_mod_mpoly_init3(Bbar, 1, bits, ctx);
    fmpz_mod_mpolyn_init(Tn, bits, ctx);
    fmpz_mod_mpolyn_init(Gn, bits, ctx);
    fmpz_mod_mpolyn_init(Abarn, bits, ctx);
    fmpz_mod_mpolyn_init(Bbarn, bits, ctx);
    fmpz_mod_polyun_init(Aeh_cur, ctx->ffinfo);
    fmpz_mod_polyun_init(Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_init(Beh_cur, ctx->ffinfo);
    fmpz_mod_polyun_init(Beh_inc, ctx->ffinfo);
    fmpz_mod_poly_init(gammaeh_cur, ctx->ffinfo);
    fmpz_mod_poly_init(gammaeh_inc, ctx->ffinfo);
    fmpz_mod_poly_init(modulus, ctx->ffinfo);
    fmpz_mod_poly_stack_init(St->poly_stack);
    fmpz_mod_polyun_stack_init(St->polyun_stack);

    Aeh_coeff_mock->exps = NULL;
    Aeh_coeff_mock->coeffs = NULL;
    Aeh_coeff_mock->length = 0;
    Aeh_coeff_mock->alloc = 0;

    Beh_coeff_mock->exps = NULL;
    Beh_coeff_mock->coeffs = NULL;
    Beh_coeff_mock->length = 0;
    Beh_coeff_mock->alloc = 0;

    betas = _fmpz_vec_init(nvars);
    alphas = _fmpz_vec_init(nvars);
    flint_randinit(state);

    beta_caches = FLINT_ARRAY_ALLOC(nvars, fmpz_mod_poly_struct);
    alpha_caches = FLINT_ARRAY_ALLOC(nvars, fmpz_mod_poly_struct);
    for (i = 0; i < nvars; i++)
    {
        fmpz_mod_poly_init(beta_caches + i, ctx->ffinfo);
        fmpz_mod_poly_init(alpha_caches + i, ctx->ffinfo);
    }

    Aevals = FLINT_ARRAY_ALLOC(nvars + 1, fmpz_mod_mpoly_struct);
    Bevals = FLINT_ARRAY_ALLOC(nvars + 1, fmpz_mod_mpoly_struct);
    gammaevals = FLINT_ARRAY_ALLOC(nvars + 1, fmpz_mod_mpoly_struct);
    for (i = 0; i < nvars; i++)
    {
        fmpz_mod_mpoly_init3(Aevals + i, 0, bits, ctx);
        fmpz_mod_mpoly_init3(Bevals + i, 0, bits, ctx);
        fmpz_mod_mpoly_init3(gammaevals + i, 0, bits, ctx);
    }
    Aevals[nvars] = *A;
    Bevals[nvars] = *B;
    gammaevals[nvars] = *gamma;

    Abideg = _fmpz_mod_mpoly_bidegree(A, ctx);
    Bbideg = _fmpz_mod_mpoly_bidegree(B, ctx);

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
        fmpz_mod_rand_not_zero(alphas + i, state, ctx->ffinfo);

    for (i = nvars - 1; i >= 2; i--)
    {
        fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + i, Aevals + i + 1, i, alphas + i, ctx);
        fmpz_mod_mpoly_evaluate_one_fmpz(Bevals + i, Bevals + i + 1, i, alphas + i, ctx);
        fmpz_mod_mpoly_evaluate_one_fmpz(gammaevals + i, gammaevals + i + 1, i, alphas + i, ctx);
        if (fmpz_mod_mpoly_is_zero(gammaevals + i, ctx))
            goto choose_alphas;
        if (Aevals[i].length < 1 || _fmpz_mod_mpoly_bidegree(Aevals + i, ctx) != Abideg)
            goto choose_alphas;
        if (Bevals[i].length < 1 || _fmpz_mod_mpoly_bidegree(Bevals + i, ctx) != Bbideg)
            goto choose_alphas;
    }

    m = 2;

    if (rGdegs == NULL)
        use = USE_G | USE_ABAR | USE_BBAR;
    else
        use = mpoly_gcd_get_use_first(rGdegs[m], Adegs[m], Bdegs[m], gammadegs[m]);

    fmpz_mod_mpoly_get_polyu1n(Aev, Aevals + m, 0, 1, ctx);
    fmpz_mod_mpoly_get_polyu1n(Bev, Bevals + m, 0, 1, ctx);

    success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev, Aev, Bev,
                                 ctx->ffinfo, St->poly_stack, St->polyun_stack);
    if (!success)
        goto cleanup;

    newdegXY = fmpz_mod_polyu1n_bidegree(Gev);
    if (newdegXY > GdegboundXY)
        goto choose_alphas;
    if (newdegXY < GdegboundXY)
    {
        GdegboundXY = newdegXY;
        if (GdegboundXY == 0)
            goto gcd_is_trivial;
    }

    fmpz_mod_mpoly_get_fmpz(gammaev, gammaevals + m, ctx);
    _fmpz_mod_poly_vec_mul_fmpz_mod(Gev->coeffs, Gev->length, gammaev, ctx->ffinfo);

    fmpz_mod_mpolyn_interp_lift_sm_polyu1n(Gn, Gev, ctx);
    fmpz_mod_mpolyn_interp_lift_sm_polyu1n(Abarn, Abarev, ctx);
    fmpz_mod_mpolyn_interp_lift_sm_polyu1n(Bbarn, Bbarev, ctx);

    fmpz_mod_poly_one(modulus, ctx->ffinfo);
    fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);

    fmpz_set(start_alpha, alphas + m);
    while (1)
    {
    choose_alpha_2:

        if (fmpz_cmp_ui(alphas + m, 2) < 0)
            fmpz_set(alphas + m, fmpz_mod_ctx_modulus(ctx->ffinfo));
        fmpz_sub_ui(alphas + m, alphas + m, 1);
        if (fmpz_equal(alphas + m, start_alpha))
            goto choose_alphas;

        FLINT_ASSERT(alphapow->alloc >= 2);
        fmpz_one(alphapow->coeffs + 0);
        fmpz_set(alphapow->coeffs + 1, alphas + m);
        alphapow->length = 2;

        fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
        fmpz_mod_mpoly_evaluate_one_fmpz(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);

        fmpz_mod_mpoly_evaluate_one_fmpz(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
        if (fmpz_mod_mpoly_is_zero(gammaevals + m, ctx))
            goto choose_alpha_2;
        if (Aevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
            goto choose_alpha_2;
        if (Bevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
            goto choose_alpha_2;

        fmpz_mod_mpoly_get_polyu1n(Aev, Aevals + m, 0, 1, ctx);
        fmpz_mod_mpoly_get_polyu1n(Bev, Bevals + m, 0, 1, ctx);

        success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev,
                      Aev, Bev, ctx->ffinfo, St->poly_stack, St->polyun_stack);
        if (!success)
            goto cleanup;

        newdegXY = fmpz_mod_polyu1n_bidegree(Gev);
        if (newdegXY > GdegboundXY)
            goto choose_alpha_2;
        if (newdegXY < GdegboundXY)
        {
            GdegboundXY = newdegXY;
            if (GdegboundXY == 0)
                goto gcd_is_trivial;
            goto choose_alphas;
        }

        fmpz_mod_mpoly_get_fmpz(gammaev, gammaevals + m, ctx);
        _fmpz_mod_poly_vec_mul_fmpz_mod(Gev->coeffs, Gev->length, gammaev, ctx->ffinfo);

        fmpz_mod_poly_eval_pow(c, modulus, alphapow, ctx->ffinfo);
        fmpz_mod_inv(c, c, ctx->ffinfo);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, c, ctx->ffinfo);

        if ((use & USE_G) && !fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
                                &lastdeg, Gn, Tn, Gev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                success = fmpz_mod_mpolyl_content(cont, rG, 2, ctx);
                if (!success)
                    goto cleanup;
                fmpz_mod_mpoly_divides(rG, rG, cont, ctx);
                fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                if (fmpz_mod_mpoly_divides(rAbar, A, rG, ctx) &&
                    fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }
            }
            else
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(G, Gn, m, ctx);
                fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fmpz_mod_mpoly_divides(Abar, T, G, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                    fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_ABAR) && !fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
                          &lastdeg, Abarn, Tn, Abarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                success = fmpz_mod_mpolyl_content(cont, rAbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fmpz_mod_mpoly_divides(rAbar, rAbar, cont, ctx);
                fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                if (fmpz_mod_mpoly_divides(rG, A, rAbar, ctx) &&
                    fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(Abar, Abarn, m, ctx);
                fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                if (fmpz_mod_mpoly_divides(G, T, Abar, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(Bbar, T, G, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(Bbar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if ((use & USE_BBAR) && !fmpz_mod_mpolyn_interp_crt_sm_polyu1n(
                          &lastdeg, Bbarn, Tn, Bbarev, modulus, alphapow, ctx))
        {
            if (m == nvars - 1)
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                success = fmpz_mod_mpolyl_content(cont, rBbar, 2, ctx);
                if (!success)
                    goto cleanup;
                fmpz_mod_mpoly_divides(rBbar, rBbar, cont, ctx);
                fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                if (fmpz_mod_mpoly_divides(rG, B, rBbar, ctx) &&
                    fmpz_mod_mpoly_divides(rAbar, A, rG, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    break;
                }                
            }
            else
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(Bbar, Bbarn, m, ctx);
                fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                if (fmpz_mod_mpoly_divides(G, T, Bbar, ctx))
                {
                    fmpz_mod_mpoly_repack_bits_inplace(G, bits, ctx);
                    fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(Abar, T, G, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(Abar, bits, ctx);
                        break;
                    }
                }
            }
        }

        if (fmpz_mod_poly_degree(modulus, ctx->ffinfo) > gammadegs[m] + Adegs[m] &&
            fmpz_mod_poly_degree(modulus, ctx->ffinfo) > gammadegs[m] + Bdegs[m])
        {
            goto choose_alphas;
        }

        fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
        fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);
    }

    for (m = 3; m < nvars; m++)
    {
        /* G, Abar, Bbar are in Fp[gen(0), ..., gen(m - 1)] */
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(Gn, G, ctx);
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(Abarn, Abar, ctx);
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(Bbarn, Bbar, ctx);

        if (rGdegs == NULL)
            use = USE_G | USE_ABAR | USE_BBAR;
        else
            use = 0;

        fmpz_mod_poly_one(modulus, ctx->ffinfo);
        fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
        fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);

        fmpz_set(start_alpha, alphas + m);

    choose_betas:

        /* only beta[2], beta[1], ..., beta[m - 1] will be used */
        for (i = 2; i < nvars; i++)
        {
            fmpz_mod_rand_not_zero(betas + i, state, ctx->ffinfo);
            fmpz_mod_pow_cache_start(betas + i, beta_caches + i, ctx->ffinfo);
        }

        fmpz_mod_mpoly_monomial_evals2(HG, G, beta_caches + 2, m, ctx, tmark);
        fmpz_mod_mpoly_monomial_evals2(HAbar, Abar, beta_caches + 2, m, ctx, tmark);
        fmpz_mod_mpoly_monomial_evals2(HBbar, Bbar, beta_caches + 2, m, ctx, tmark);

        if (use == 0)
        {
            this_length = gammaevals[m + 1].length + Aevals[m + 1].length +
                                                     Bevals[m + 1].length;

            use = fmpz_mod_mpoly_gcd_get_use_new(rGdegs[m], Adegs[m], Bdegs[m],
                  gammadegs[m], degxAB, degyAB, this_length, HG, HAbar, HBbar);
        }

        req_zip_images = 1;
        if (use & USE_G)
        {
            this_length = fmpz_mod_polyun_product_roots(MG, HG, ctx->ffinfo);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_ABAR)
        {
            this_length = fmpz_mod_polyun_product_roots(MAbar, HAbar, ctx->ffinfo);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }
        if (use & USE_BBAR)
        {
            this_length = fmpz_mod_polyun_product_roots(MBbar, HBbar, ctx->ffinfo);
            req_zip_images = FLINT_MAX(req_zip_images, this_length);
        }

        while (1)
        {
        choose_alpha_m:

            if (fmpz_cmp_ui(alphas + m, 2) < 0)
                fmpz_set(alphas + m, fmpz_mod_ctx_modulus(ctx->ffinfo));
            fmpz_sub_ui(alphas + m, alphas + m, 1);
            if (fmpz_equal(alphas + m, start_alpha))
                goto choose_alphas;

            FLINT_ASSERT(alphapow->alloc >= 2);
            fmpz_one(alphapow->coeffs + 0);
            fmpz_set(alphapow->coeffs + 1, alphas + m);
            alphapow->length = 2;

            fmpz_mod_mpoly_evaluate_one_fmpz(Aevals + m, Aevals + m + 1, m, alphas + m, ctx);
            fmpz_mod_mpoly_evaluate_one_fmpz(Bevals + m, Bevals + m + 1, m, alphas + m, ctx);
            fmpz_mod_mpoly_evaluate_one_fmpz(gammaevals + m, gammaevals + m + 1, m, alphas + m, ctx);
            if (fmpz_mod_mpoly_is_zero(gammaevals + m, ctx))
                goto choose_alpha_m;
            if (Aevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Aevals + m, ctx) != Abideg)
                goto choose_alpha_m;
            if (Bevals[m].length < 1 || _fmpz_mod_mpoly_bidegree(Bevals + m, ctx) != Bbideg)
                goto choose_alpha_m;

            fmpz_mod_mpoly_monomial_evals2(Aeh_inc, Aevals + m, beta_caches + 2, m, ctx, tmark);
            fmpz_mod_mpoly_monomial_evals2(Beh_inc, Bevals + m, beta_caches + 2, m, ctx, tmark);
            fmpz_mod_mpoly_monomial_evals(gammaeh_inc, gammaevals + m, beta_caches + 2, 2, m, ctx);
            fmpz_mod_polyun_set(Aeh_cur, Aeh_inc, ctx->ffinfo);
            fmpz_mod_polyun_set(Beh_cur, Beh_inc, ctx->ffinfo);
            fmpz_mod_poly_set(gammaeh_cur, gammaeh_inc, ctx->ffinfo);
            fmpz_mod_mpoly_mock_eval_coeff(Aeh_coeff_mock, Aevals + m, Aeh_inc, ctx);
            fmpz_mod_mpoly_mock_eval_coeff(Beh_coeff_mock, Bevals + m, Beh_inc, ctx);

            fmpz_mod_polyun_zip_start(ZG, HG, req_zip_images, ctx->ffinfo);
            fmpz_mod_polyun_zip_start(ZAbar, HAbar, req_zip_images, ctx->ffinfo);
            fmpz_mod_polyun_zip_start(ZBbar, HBbar, req_zip_images, ctx->ffinfo);
            for (cur_zip_image = 0; cur_zip_image < req_zip_images; cur_zip_image++)
            {
                fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(Aev, Aeh_cur, Aeh_inc, Aeh_coeff_mock, ctx->ffinfo);
                fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(Bev, Beh_cur, Beh_inc, Beh_coeff_mock, ctx->ffinfo);
                fmpz_mod_poly_zip_eval_cur_inc_coeff(gammaev, gammaeh_cur, gammaeh_inc, gammaevals + m, ctx);
                if (fmpz_is_zero(gammaev))
                    goto choose_betas;
                if (Aev->length < 1 || fmpz_mod_polyu1n_bidegree(Aev) != Abideg)
                    goto choose_betas;
                if (Bev->length < 1 || fmpz_mod_polyu1n_bidegree(Bev) != Bbideg)
                    goto choose_betas;

                success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev,
                      Aev, Bev, ctx->ffinfo, St->poly_stack, St->polyun_stack);        
                if (!success)
                    goto cleanup;

                newdegXY = fmpz_mod_polyu1n_bidegree(Gev);
                if (newdegXY > GdegboundXY)
                    goto choose_betas;
                if (newdegXY < GdegboundXY)
                {
                    GdegboundXY = newdegXY;
                    if (GdegboundXY == 0)
                        goto gcd_is_trivial;
                    goto choose_alphas;
                }

                _fmpz_mod_poly_vec_mul_fmpz_mod(Gev->coeffs, Gev->length, gammaev, ctx->ffinfo);
                if ((use & USE_G) && !fmpz_mod_polyun_add_zip_must_match(ZG, Gev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_ABAR) && !fmpz_mod_polyun_add_zip_must_match(ZAbar, Abarev, cur_zip_image))
                    goto choose_alphas;
                if ((use & USE_BBAR) && !fmpz_mod_polyun_add_zip_must_match(ZBbar, Bbarev, cur_zip_image))
                    goto choose_alphas;
            }

            if ((use & USE_G) && fmpz_mod_polyun_zip_solve(G, ZG, HG, MG, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_ABAR) && fmpz_mod_polyun_zip_solve(Abar, ZAbar, HAbar, MAbar, ctx) < 1)
                goto choose_alphas;
            if ((use & USE_BBAR) && fmpz_mod_polyun_zip_solve(Bbar, ZBbar, HBbar, MBbar, ctx) < 1)
                goto choose_alphas;

            FLINT_ASSERT(fmpz_mod_poly_degree(modulus, ctx->ffinfo) > 0);
            fmpz_mod_poly_eval_pow(c, modulus, alphapow, ctx->ffinfo);
            fmpz_mod_inv(c, c, ctx->ffinfo);
            fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, c, ctx->ffinfo);

            if ((use & USE_G) && !fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
                                      &lastdeg, Gn, G, modulus, alphapow, ctx))
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rG, Gn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fmpz_mod_mpolyl_content(cont, rG, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fmpz_mod_mpoly_divides(rG, rG, cont, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                    if (fmpz_mod_mpoly_divides(rAbar, A, rG, ctx) &&
                        fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(rAbar, T, rG, ctx))
                    {
                        fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fmpz_mod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fmpz_mod_mpoly_swap(G, rG, ctx);
                            fmpz_mod_mpoly_swap(Abar, rAbar, ctx);
                            fmpz_mod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_ABAR) && !fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Abarn, Abar, modulus, alphapow, ctx))
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rAbar, Abarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fmpz_mod_mpolyl_content(cont, rAbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fmpz_mod_mpoly_divides(rAbar, rAbar, cont, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                    if (fmpz_mod_mpoly_divides(rG, A, rAbar, ctx) &&
                        fmpz_mod_mpoly_divides(rBbar, B, rG, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(rG, T, rAbar, ctx))
                    {
                        fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                        if (fmpz_mod_mpoly_divides(rBbar, T, rG, ctx))
                        {
                            fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                            fmpz_mod_mpoly_swap(G, rG, ctx);
                            fmpz_mod_mpoly_swap(Abar, rAbar, ctx);
                            fmpz_mod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }
            if ((use & USE_BBAR) && !fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(
                                &lastdeg, Bbarn, Bbar, modulus, alphapow, ctx))
            {
                fmpz_mod_mpoly_cvtfrom_mpolyn(rBbar, Bbarn, m, ctx);
                if (m == nvars - 1)
                {
                    success = fmpz_mod_mpolyl_content(cont, rBbar, 2, ctx);
                    if (!success)
                        goto cleanup;
                    fmpz_mod_mpoly_divides(rBbar, rBbar, cont, ctx);
                    fmpz_mod_mpoly_repack_bits_inplace(rBbar, bits, ctx);
                    if (fmpz_mod_mpoly_divides(rG, B, rBbar, ctx) &&
                        fmpz_mod_mpoly_divides(rAbar, A, rG, ctx))
                    {
                        fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                        fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                        break;
                    }
                }
                else
                {
                    fmpz_mod_mpoly_mul(T, Bevals + m + 1, gammaevals + m + 1, ctx);
                    if (fmpz_mod_mpoly_divides(rG, T, rBbar, ctx))
                    {
                        fmpz_mod_mpoly_mul(T, Aevals + m + 1, gammaevals + m + 1, ctx);
                        if (fmpz_mod_mpoly_divides(rAbar, T, rG, ctx))
                        {
                            fmpz_mod_mpoly_repack_bits_inplace(rG, bits, ctx);
                            fmpz_mod_mpoly_repack_bits_inplace(rAbar, bits, ctx);
                            fmpz_mod_mpoly_swap(G, rG, ctx);
                            fmpz_mod_mpoly_swap(Abar, rAbar, ctx);
                            fmpz_mod_mpoly_swap(Bbar, rBbar, ctx);
                            break;
                        }
                    }
                }
            }

            if (_fmpz_mod_poly_degree(modulus) > gammadegs[m] + Adegs[m] &&
                _fmpz_mod_poly_degree(modulus) > gammadegs[m] + Bdegs[m])
            {
                goto choose_alphas;
            }

            fmpz_mod_neg(c, alphas + m, ctx->ffinfo);
            fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, c, ctx->ffinfo);
        }
    }

    success = 1;

cleanup:

    n_poly_clear(tmark);

    fmpz_mod_polyun_clear(HG, ctx->ffinfo);
    fmpz_mod_polyun_clear(HAbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(HBbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(MG, ctx->ffinfo);
    fmpz_mod_polyun_clear(MAbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(MBbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZG, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZAbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZBbar, ctx->ffinfo);
    fmpz_mod_polyun_clear(Aev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Bev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Gev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Abarev, ctx->ffinfo);
    fmpz_mod_polyun_clear(Bbarev, ctx->ffinfo);
    fmpz_mod_poly_clear(alphapow, ctx->ffinfo);
    fmpz_mod_mpoly_clear(cont, ctx);
    fmpz_mod_mpoly_clear(T, ctx);
    fmpz_mod_mpoly_clear(G, ctx);
    fmpz_mod_mpoly_clear(Abar, ctx);
    fmpz_mod_mpoly_clear(Bbar, ctx);
    fmpz_mod_mpolyn_clear(Tn, ctx);
    fmpz_mod_mpolyn_clear(Gn, ctx);
    fmpz_mod_mpolyn_clear(Abarn, ctx);
    fmpz_mod_mpolyn_clear(Bbarn, ctx);
    fmpz_mod_polyun_clear(Aeh_cur, ctx->ffinfo);
    fmpz_mod_polyun_clear(Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_clear(Beh_cur, ctx->ffinfo);
    fmpz_mod_polyun_clear(Beh_inc, ctx->ffinfo);
    fmpz_mod_poly_clear(gammaeh_cur, ctx->ffinfo);
    fmpz_mod_poly_clear(gammaeh_inc, ctx->ffinfo);
    fmpz_mod_poly_clear(modulus, ctx->ffinfo);
    fmpz_mod_poly_stack_clear(St->poly_stack);
    fmpz_mod_polyun_stack_clear(St->polyun_stack);

    flint_free(Aeh_coeff_mock->coeffs);
    flint_free(Beh_coeff_mock->coeffs);

    for (i = 0; i < nvars; i++)
    {
        fmpz_mod_poly_clear(beta_caches + i, ctx->ffinfo);
        fmpz_mod_poly_clear(alpha_caches + i, ctx->ffinfo);
    }
    flint_free(beta_caches);
    flint_free(alpha_caches);

    _fmpz_vec_clear(betas, nvars);
    _fmpz_vec_clear(alphas, nvars);

    flint_randclear(state);

    for (i = 0; i < nvars; i++)
    {
        fmpz_mod_mpoly_clear(Aevals + i, ctx);
        fmpz_mod_mpoly_clear(Bevals + i, ctx);
        fmpz_mod_mpoly_clear(gammaevals + i, ctx);
    }
    flint_free(Aevals);
    flint_free(Bevals);
    flint_free(gammaevals);

    fmpz_clear(gammaev);
    fmpz_clear(c);
    fmpz_clear(start_alpha);

    FLINT_ASSERT(!success || rG->bits == bits);
    FLINT_ASSERT(!success || rAbar->bits == bits);
    FLINT_ASSERT(!success || rBbar->bits == bits);

    FLINT_ASSERT(success || fmpz_abs_fits_ui(fmpz_mod_mpoly_ctx_modulus(ctx)));

    return success;

gcd_is_trivial:

    fmpz_mod_mpoly_one(rG, ctx);
    fmpz_mod_mpoly_set(rAbar, A, ctx);
    fmpz_mod_mpoly_set(rBbar, B, ctx);

    success = 1;
    
    goto cleanup;
}

