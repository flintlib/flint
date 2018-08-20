/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"
#include "fmpz_mpoly.h"


void nmod_mpoly_ctx_change_modulus(nmod_mpoly_ctx_t ctx, mp_limb_t modulus)
{
    nmodf_ctx_clear(ctx->ffinfo);
    nmodf_ctx_init(ctx->ffinfo, modulus);
}

/*
    Find a bound on the bits of the coefficients of gcd(A,B).
    If this overflows a mp_bitcnt_t, the max mp_bitcnt_t is returned.
    We will apply a Kronecker substitution and use the Landau-Mignotte bound
        for univariates.

         min(deg A, deg B)
        2                  * gcd(A[0], B[0])

                 * min( L2Norm(A) / A[0] , L2Norm(B) / B[0] )

    This number is almost certainly not going to be used,
        so it just needs to be correct and not nececessary as tight as possible.
*/
mp_bitcnt_t fmpz_mpolyu_gcd_bitbound(const fmpz_t gcdlc,
                            const fmpz_mpolyu_t A, const fmpz_mpolyu_t B,
                       const fmpz_mpoly_ctx_t ctx, const mpoly_zipinfo_t zinfo)
{
    slong i;
    fmpz_t n, an, bn;
    ulong max, len;
    mp_bitcnt_t r;

    /* find the degree of the kronecker substitution in the lesser variables */
    fmpz_init_set_ui(n, UWORD(1));
    for (i = 0; i < zinfo->nvars - 1; i++)
    {
        fmpz_mul_ui(n, n, FLINT_MAX((ulong)(zinfo->Adegs[i]),
                                    (ulong)(zinfo->Bdegs[i])) + UWORD(1));
    }

    /* n = min(deg A, deg B) after the KS */
    i = zinfo->nvars - 1;
    fmpz_addmul_ui(n, n, FLINT_MIN((ulong)(zinfo->Adegs[i]),
                                   (ulong)(zinfo->Bdegs[i])));

    /* n += log2(gcd(A[0], B[0])) */
    fmpz_add_ui(n, n, fmpz_bits(gcdlc));

    len = max = UWORD(0);
    for (i = 0; i < A->length; i++)
    {
        len += (A->coeffs + i)->length;
        max = FLINT_MAX(max, FLINT_ABS(
                    _fmpz_vec_max_bits((A->coeffs + i)->coeffs,
                                       (A->coeffs + i)->length)
              ));
    }
    fmpz_init_set_ui(an, n_clog(len, UWORD(2))/UWORD(2));
    fmpz_add_ui(an, an, max);
    fmpz_sub_ui(an, an, fmpz_bits(fmpz_mpolyu_leadcoeff_ref(A)));
    FLINT_ASSERT(fmpz_sgn(an) >= 0);

    len = max = UWORD(0);
    for (i = 0; i < B->length; i++)
    {
        len += (B->coeffs + i)->length;
        max = FLINT_MAX(max, FLINT_ABS(
                    _fmpz_vec_max_bits((B->coeffs + i)->coeffs,
                                       (B->coeffs + i)->length)
              ));
    }
    fmpz_init_set_ui(bn, n_clog(len, UWORD(2))/UWORD(2));
    fmpz_add_ui(bn, bn, max);
    fmpz_sub_ui(bn, bn, fmpz_bits(fmpz_mpolyu_leadcoeff_ref(B)));
    FLINT_ASSERT(fmpz_sgn(bn) >= 0);

    /* n += log2( min( L2Norm(A) / A[0] , L2Norm(B) / B[0] ) ) */
    fmpz_add(n, n, fmpz_cmp(an, bn) < 0 ? an : bn);

    FLINT_ASSERT(fmpz_sgn(n) > 0);

    r = fmpz_abs_fits_ui(n) ? fmpz_get_ui(n) : -UWORD(1);

    fmpz_clear(n);
    fmpz_clear(an);
    fmpz_clear(bn);

    return r;
}

int fmpz_mpolyu_gcdm_zippel(fmpz_mpolyu_t G, fmpz_mpolyu_t A, fmpz_mpolyu_t B,
                            const fmpz_mpoly_ctx_t ctx, mpoly_zipinfo_t zinfo,
                                                        flint_rand_t randstate)
{
    mp_bitcnt_t coeffbitbound;
    mp_bitcnt_t coeffbits;
    slong degbound;
    int success, changed;
    mp_limb_t p = UWORD(1) << (FLINT_BITS - 1), old_p, t, gammap;
    fmpz_t gamma, pp, gammapp, modulus;
    nmod_mpolyu_t Ap, Bp, Gp, Gform;
    fmpz_mpolyu_t H;
    nmod_mpoly_ctx_t ctxp;

    fmpz_init(pp);
    fmpz_init(gammapp);
    fmpz_init_set_si(modulus, 1);
    fmpz_init(gamma);
    fmpz_gcd(gamma, fmpz_mpolyu_leadcoeff_ref(A), fmpz_mpolyu_leadcoeff_ref(B));

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->bits == G->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->exps[A->length - 1] == 0);
    FLINT_ASSERT(B->exps[B->length - 1] == 0);

    degbound = FLINT_MIN(A->exps[0], B->exps[0]);

    coeffbitbound = fmpz_mpolyu_gcd_bitbound(gamma, A, B, ctx, zinfo);

    nmod_mpoly_ctx_init(ctxp, ctx->minfo->nvars, ORD_LEX, 2);

    nmod_mpolyu_init(Ap, A->bits, ctxp);
    nmod_mpolyu_init(Bp, A->bits, ctxp);
    nmod_mpolyu_init(Gp, A->bits, ctxp);
    nmod_mpolyu_init(Gform, A->bits, ctxp);

    fmpz_mpolyu_init(H, A->bits, ctx);

choose_prime_outer:
    old_p = p;
    p = n_nextprime(p, 1);
    if (p <= old_p) {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    /* make sure mod p reduction does not kill both lc(A) and lc(B) */
    fmpz_set_ui(pp, p);
    fmpz_mod(gammapp, gamma, pp);
    gammap = fmpz_get_ui(gammapp);
    if (gammap == UWORD(0))
        goto choose_prime_outer;

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    /* make sure mod p reduction does not kill either A or B */
    fmpz_mpolyu_to_nmod_mpolyu(Ap, ctxp, A, ctx);
    fmpz_mpolyu_to_nmod_mpolyu(Bp, ctxp, B, ctx);
    if (Ap->length == 0 || Bp->length == 0)
        goto choose_prime_outer;

    success = nmod_mpolyu_gcdp_zippel(Gp, Ap, Bp, ctx->minfo->nvars - 1,
                                                       ctxp, zinfo, randstate);
    if (!success || Gp->exps[0] > degbound)
        goto choose_prime_outer;
    degbound = Gp->exps[0];

    t = nmod_mpolyu_leadcoeff(Gp, ctxp);
    t = nmod_inv(t, ctxp->ffinfo->mod);
    t = nmod_mul(t, gammap, ctxp->ffinfo->mod);
    nmod_mpolyu_scalar_mul_nmod(Gp, t, ctxp);

    if (Gp->length == 1 && Gp->exps[0] == 0)
    {
        FLINT_ASSERT(nmod_mpoly_is_nmod(Gp->coeffs + 0, ctxp));
        FLINT_ASSERT(!nmod_mpoly_is_zero(Gp->coeffs + 0, ctxp));
        fmpz_mpolyu_one(G, ctx);
        
        success = 1;
        goto finished;
    }

    nmod_mpolyu_setform(Gform, Gp, ctxp);
    fmpz_mpolyu_set_nmod_mpolyu(H, ctx, Gp, ctxp);
    fmpz_set_ui(modulus, p);

choose_prime_inner:
    old_p = p;
    p = n_nextprime(p, 1);
    if (p <= old_p) {
        /* ran out of primes */
        success = 0;
        goto finished;
    }

    /* make sure mod p reduction does not kill both lc(A) and lc(B) */
    fmpz_set_ui(pp, p);
    fmpz_mod(gammapp, gamma, pp);
    gammap = fmpz_get_ui(gammapp);
    if (gammap == UWORD(0))
        goto choose_prime_inner;

    nmod_mpoly_ctx_change_modulus(ctxp, p);

    /* make sure mod p reduction does not kill either A or B */
    fmpz_mpolyu_to_nmod_mpolyu(Ap, ctxp, A, ctx);
    fmpz_mpolyu_to_nmod_mpolyu(Bp, ctxp, B, ctx);
    if (Ap->length == 0 || Bp->length == 0)
        goto choose_prime_inner;

    switch (nmod_mpolyu_gcds_zippel(Gp, Ap, Bp, Gform, ctx->minfo->nvars,
                                                   ctxp, randstate, &degbound))
    {
        default:
            FLINT_ASSERT(0);
        case nmod_gcds_form_main_degree_too_high:
        case nmod_gcds_form_wrong:
        case nmod_gcds_no_solution:
            goto choose_prime_outer;
        case nmod_gcds_scales_not_found:
        case nmod_gcds_eval_point_not_found:
        case nmod_gcds_eval_gcd_deg_too_high:
            goto choose_prime_inner;
        case nmod_gcds_success:
            NULL;
    }

    if (nmod_mpolyu_leadcoeff(Gp, ctxp) == UWORD(0))
        goto choose_prime_inner;


    t = nmod_mpolyu_leadcoeff(Gp, ctxp);
    t = nmod_inv(t, ctxp->ffinfo->mod);
    t = nmod_mul(t, gammap, ctxp->ffinfo->mod);
    nmod_mpolyu_scalar_mul_nmod(Gp, t, ctxp);

    changed = fmpz_mpolyu_CRT_nmod_mpolyu(&coeffbits, H, ctx, modulus, Gp, ctxp);
    fmpz_mul_ui(modulus, modulus, p);

    if (changed)
    {
        if (coeffbits > coeffbitbound)
        {
            goto choose_prime_outer;
        }
        goto choose_prime_inner;
    }

    fmpz_mpolyu_fmpz_content(pp, H, ctx);
    fmpz_mpolyu_scalar_divexact_fmpz(G, H, pp, ctx);

    if (!fmpz_mpolyu_divides(A, G, ctx) || !fmpz_mpolyu_divides(B, G, ctx))
        goto choose_prime_inner;

    success = 1;

finished:

    nmod_mpolyu_clear(Ap, ctxp);
    nmod_mpolyu_clear(Bp, ctxp);
    nmod_mpolyu_clear(Gp, ctxp);
    nmod_mpolyu_clear(Gform, ctxp);

    fmpz_mpolyu_clear(H, ctx);

    nmod_mpoly_ctx_clear(ctxp);

    fmpz_clear(gammapp);
    fmpz_clear(gamma);
    fmpz_clear(modulus);
    fmpz_clear(pp);

    return success;
}



int fmpz_mpolyu_gcd_zippel(
    fmpz_mpolyu_t G,
    fmpz_mpolyu_t A,
    fmpz_mpolyu_t B,
    fmpz_mpoly_ctx_t ctx,
    mpoly_zipinfo_t zinfo,
    flint_rand_t randstate)
{
    int success = 0;
    slong i;
    slong ABminshift;
    fmpz_mpoly_t content;
    fmpz_mpolyu_t Abar, Bbar, Gbar;

    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    fmpz_mpoly_init(content, ctx);
    fmpz_mpolyu_init(Abar, A->bits, ctx);
    fmpz_mpolyu_init(Bbar, A->bits, ctx);
    fmpz_mpolyu_init(Gbar, A->bits, ctx);

    fmpz_mpoly_set(content, A->coeffs + 0, ctx);
    for (i = 1; i < A->length; i++)
    {
        if (fmpz_mpoly_is_one(content, ctx))
            break;
        success = _fmpz_mpoly_gcd_zippel(content, content, A->coeffs + i, ctx);
        if (!success)
            goto finished;
    }
    for (i = 0; i < B->length; i++)
    {
        if (fmpz_mpoly_is_one(content, ctx))
            break;
        success = _fmpz_mpoly_gcd_zippel(content, content, B->coeffs + i, ctx);
        if (!success)
            goto finished;
    }

    fmpz_mpolyu_divexact_mpoly(Abar, A, content, ctx);
    fmpz_mpolyu_divexact_mpoly(Bbar, B, content, ctx);

    ABminshift = FLINT_MIN(Abar->exps[Abar->length - 1], Bbar->exps[Bbar->length - 1]);
    fmpz_mpolyu_shift_right(Abar, Abar->exps[Abar->length - 1]);
    fmpz_mpolyu_shift_right(Bbar, Bbar->exps[Bbar->length - 1]);

    success = fmpz_mpolyu_gcdm_zippel(Gbar, Abar, Bbar, ctx, zinfo, randstate);
    if (!success)
        goto finished;

    fmpz_mpolyu_shift_left(Gbar, ABminshift);
    fmpz_mpolyu_mul_mpoly(G, Gbar, content, ctx);

    success = 1;

finished:

    fmpz_mpolyu_clear(Abar, ctx);
    fmpz_mpolyu_clear(Bbar, ctx);
    fmpz_mpolyu_clear(Gbar, ctx);
    fmpz_mpoly_clear(content, ctx);

    return 1;
}










void fmpz_mpoly_to_fmpz_poly_keepbits(fmpz_poly_t A, slong * Ashift,
               const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, shift, off, N;
    slong _Ashift = 0, len = B->length;
    fmpz * coeff = B->coeffs;
    ulong * exp = B->exps;
    mp_bitcnt_t bits = B->bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift(&off, &shift, var, N, bits, ctx->minfo);

    fmpz_poly_zero(A);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _Ashift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            ulong k = ((exp[N*i + off] >> shift) & mask) - _Ashift;
            FLINT_ASSERT(((slong)k) >= 0);
            fmpz_poly_set_coeff_fmpz(A, k, coeff + i);
        }
    }

    *Ashift = _Ashift;
}

void fmpz_mpoly_from_fmpz_poly_keepbits(fmpz_mpoly_t A, const fmpz_poly_t B,
                           slong Bshift, slong var, mp_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)
{
    slong shift, off, N;
    slong k;
    slong Alen;
    fmpz * Acoeff;
    ulong * Aexp;
    slong Aalloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    FLINT_ASSERT(!fmpz_poly_is_zero(B));
    FLINT_ASSERT(Bshift >= 0);
    FLINT_ASSERT(Bshift + fmpz_poly_degree(B) >= 0);
    FLINT_ASSERT(1 + FLINT_BIT_COUNT(Bshift + fmpz_poly_degree(B)) <= bits);
    
    N = mpoly_words_per_exp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_oneexp_offset_shift(one, &off, &shift, var, N, bits, ctx->minfo);

    fmpz_mpoly_fit_bits(A, bits, ctx);
    A->bits = bits;

    Acoeff = A->coeffs;
    Aexp = A->exps;
    Aalloc = A->alloc;
    Alen = 0;
    for (k = fmpz_poly_degree(B); k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&Acoeff, &Aexp, &Aalloc, Alen + 1, N);
        mpoly_monomial_mul_si(Aexp + N*Alen, one, N, k + Bshift);
        fmpz_poly_get_coeff_fmpz(Acoeff + Alen, B, k);
        Alen += !fmpz_is_zero(Acoeff + Alen);
    }

    A->coeffs = Acoeff;
    A->exps = Aexp;
    A->alloc = Aalloc;
    _fmpz_mpoly_set_length(A, Alen, ctx);

    TMP_END;
}


/*
    like fmpz_mpoly_gcd_zippel, but G and A and B all have the same bits
*/
int _fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G, fmpz_mpoly_t A,
                                    fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    flint_rand_t randstate;
    int ret, success = 0;
    mp_bitcnt_t new_bits;
    mpoly_zipinfo_t zinfo;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu;

    fmpz_mpoly_t Acopy, Bcopy;
    fmpz_mpoly_init(Acopy, ctx);
    fmpz_mpoly_init(Bcopy, ctx);
    fmpz_mpoly_set(Acopy, A, ctx);
    fmpz_mpoly_set(Bcopy, B, ctx);

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);
    FLINT_ASSERT(A->bits == B->bits);
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));
    new_bits = A->bits;

    if (ctx->minfo->nvars == 1) {
        slong shiftA, shiftB;
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_mpoly_to_fmpz_poly_keepbits(a, &shiftA, A, 0, ctx);
        fmpz_mpoly_to_fmpz_poly_keepbits(b, &shiftB, B, 0, ctx);
        fmpz_poly_gcd(g, a, b);
        fmpz_mpoly_from_fmpz_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
        fmpz_mpoly_clear(Acopy, ctx);
        fmpz_mpoly_clear(Bcopy, ctx);
        return 1;
    }


    flint_randinit(randstate);

    /* TODO: choose nvars and perm more cleverly */
    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fmpz_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fmpz_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = i;
    }

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);


    fmpz_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

    ret = fmpz_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret) {
        fmpz_mpoly_from_mpolyu_perm(G, Gu, A->bits, zinfo->perm, uctx, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, G, ctx);
        success = 1;
    }

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    fmpz_mpoly_clear(Acopy, ctx);
    fmpz_mpoly_clear(Bcopy, ctx);

    return success;
}




int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G,
        const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    mp_bitcnt_t new_bits;
    flint_rand_t randstate;
    int ret, success = 0;
    mpoly_zipinfo_t zinfo;
    fmpz_mpoly_ctx_t uctx;
    fmpz_mpolyu_t Au, Bu, Gu;

    if (fmpz_mpoly_is_zero(A, ctx))
    {
        if (fmpz_mpoly_is_zero(B, ctx))
        {
            fmpz_mpoly_zero(G, ctx);
            return 1;
        }
        if (fmpz_sgn(B->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, B, ctx);
            return 1;
        } else
        {
            fmpz_mpoly_set(G, B, ctx);
            return 1;
        }
    }

    if (fmpz_mpoly_is_zero(B, ctx))
    {
        if (fmpz_sgn(A->coeffs + 0) < 0)
        {
            fmpz_mpoly_neg(G, A, ctx);
            return 1;
        } else
        {
            fmpz_mpoly_set(G, A, ctx);
            return 1;
        }
    }

    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        return 0;

    if (ctx->minfo->nvars == 1) {
        slong shiftA, shiftB;
        fmpz_poly_t a, b, g;
        fmpz_poly_init(a);
        fmpz_poly_init(b);
        fmpz_poly_init(g);
        fmpz_mpoly_to_fmpz_poly_keepbits(a, &shiftA, A, 0, ctx);
        fmpz_mpoly_to_fmpz_poly_keepbits(b, &shiftB, B, 0, ctx);
        fmpz_poly_gcd(g, a, b);
        fmpz_mpoly_from_fmpz_poly_keepbits(G, g, FLINT_MIN(shiftA, shiftB), 0, A->bits, ctx);
        fmpz_poly_clear(a);
        fmpz_poly_clear(b);
        fmpz_poly_clear(g);
        return 1;
    }

    FLINT_ASSERT(A->bits <= FLINT_BITS);
    FLINT_ASSERT(B->bits <= FLINT_BITS);

    FLINT_ASSERT(ctx->minfo->nvars > WORD(1));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(A, ctx));
    FLINT_ASSERT(!fmpz_mpoly_is_zero(B, ctx));

    flint_randinit(randstate);

    /* TODO: choose nvars and perm more cleverly */
    mpoly_zipinfo_init(zinfo, ctx->minfo->nvars);
    fmpz_mpoly_degrees_si(zinfo->Adegs, A, ctx);
    fmpz_mpoly_degrees_si(zinfo->Bdegs, B, ctx);
    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        zinfo->perm[i] = ctx->minfo->nvars - 1 - i;
        zinfo->perm[i] = i;
    }

    new_bits = FLINT_MAX(A->bits, B->bits);

    fmpz_mpoly_ctx_init(uctx, ctx->minfo->nvars - 1, ORD_LEX);
    fmpz_mpolyu_init(Au, new_bits, uctx);
    fmpz_mpolyu_init(Bu, new_bits, uctx);
    fmpz_mpolyu_init(Gu, new_bits, uctx);

    fmpz_mpoly_to_mpolyu_perm(Au, A, zinfo->perm, uctx, ctx);
    fmpz_mpoly_to_mpolyu_perm(Bu, B, zinfo->perm, uctx, ctx);

    ret = fmpz_mpolyu_gcd_zippel(Gu, Au, Bu, uctx, zinfo, randstate);
    if (ret) {
        fmpz_mpoly_from_mpolyu_perm(G, Gu, 0, zinfo->perm, uctx, ctx);
        if (fmpz_sgn(G->coeffs + 0) < 0)
            fmpz_mpoly_neg(G, G, ctx);
        success = 1;
    }

    fmpz_mpolyu_clear(Au, uctx);
    fmpz_mpolyu_clear(Bu, uctx);
    fmpz_mpolyu_clear(Gu, uctx);
    fmpz_mpoly_ctx_clear(uctx);

    mpoly_zipinfo_clear(zinfo);

    flint_randclear(randstate);

    return success;
}
