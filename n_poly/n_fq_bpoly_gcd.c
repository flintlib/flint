/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "fq_nmod_mpoly_factor.h"


/* g is in Fq[gen(1)] */
static void n_fq_bpoly_content_var0(
    n_fq_poly_t g,
    const n_fq_bpoly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong i;

    n_fq_poly_zero(g);
    for (i = 0; i < A->length; i++)
    {
        n_fq_poly_gcd(g, g, A->coeffs + i, ctx);
        if (n_fq_poly_degree(g) == 0)
            break;
    }
}

static void n_fq_bpoly_divexact_poly_var1(
    n_fq_bpoly_t A,
    const n_fq_poly_t b, /* in Fq[gen(1)] */
    const fq_nmod_ctx_t ctx)
{
    slong i;
    n_fq_poly_t t, r;

    n_fq_poly_init(t);
    n_fq_poly_init(r);

    for (i = 0; i < A->length; i++)
    {
        if (n_fq_poly_is_zero(A->coeffs + i))
            continue;

        n_fq_poly_divrem(t, r, A->coeffs + i, b, ctx);
        n_fq_poly_swap(A->coeffs + i, t);
    }

    n_fq_poly_clear(t);
    n_fq_poly_clear(r);
}


void n_fq_bpoly_mul_last(n_bpoly_t A, const n_poly_t b, const fq_nmod_ctx_t ctx)
{
    slong i;
    n_fq_poly_t t;

    n_fq_poly_init(t);

    for (i = 0; i < A->length; i++)
    {
        if (n_fq_poly_is_zero(A->coeffs + i))
            continue;

        n_fq_poly_mul(t, A->coeffs + i, b, ctx);
        n_fq_poly_set(A->coeffs + i, t, ctx);
    }

    n_fq_poly_clear(t);
}

/*****************************************************************************/

void n_fq_poly_eval2p_pow(
    mp_limb_t * vp,
    mp_limb_t * vm,
    const n_fq_poly_t P, 
    n_poly_t alphapow,
    slong d,
    nmod_t ctx)
{
    const mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong i, k;

    FLINT_ASSERT(P->alloc >= d*Plen);

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                               alphapow->coeffs[1], ctx);
        }
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
    }

    for (i = 0; i < d; i++)
    {
        a0 = a1 = a2 = 0;
        b0 = b1 = b2 = 0;

        for (k = 0; k + 2 <= Plen; k += 2)
        {
            umul_ppmm(p1, p0, Pcoeffs[d*(k + 0) + i], alpha_powers[k + 0]);
            umul_ppmm(q1, q0, Pcoeffs[d*(k + 1) + i], alpha_powers[k + 1]);
            add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
            add_sssaaaaaa(b2, b1, b0, b2, b1, b0, 0, q1, q0);
        }

        if (k < Plen)
        {
            umul_ppmm(p1, p0, Pcoeffs[d*(k + 0) + i], alpha_powers[k + 0]);
            add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
            k++;
        }

        FLINT_ASSERT(k == Plen);

        NMOD_RED3(p0, a2, a1, a0, ctx);
        NMOD_RED3(q0, b2, b1, b0, ctx);

        vp[i] = nmod_add(p0, q0, ctx);
        vm[i] = nmod_sub(p0, q0, ctx);
    }
}

void n_fq_bpoly_interp_reduce_2psm_poly(
    n_fq_poly_t Ap,
    n_fq_poly_t Am,
    const n_fq_bpoly_t A,
    n_poly_t alphapow,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Alen = A->length;
    const n_poly_struct * Ac = A->coeffs;
    mp_limb_t * Apc, * Amc;

    n_poly_fit_length(Ap, d*Alen);
    n_poly_fit_length(Am, d*Alen);

    Apc = Ap->coeffs;
    Amc = Am->coeffs;

    for (i = 0; i < Alen; i++)
        n_fq_poly_eval2p_pow(Apc + d*i, Amc + d*i, Ac + i, alphapow,
                                                      d, fq_nmod_ctx_mod(ctx));

    Ap->length = Alen;
    _n_fq_poly_normalise(Ap, d);
    Am->length = Alen;
    _n_fq_poly_normalise(Am, d);

}

void n_fq_bpoly_interp_lift_2psm_poly(
    slong * deg1,
    n_fq_bpoly_t T,
    const n_fq_poly_t A,
    const n_fq_poly_t B,
    mp_limb_t alpha,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    slong i, j;
    slong lastlength = 0;
    const mp_limb_t * Acoeffs = A->coeffs;
    const mp_limb_t * Bcoeffs = B->coeffs;
    n_fq_poly_struct * Tcoeffs;
    slong Alen = A->length;
    slong Blen = B->length;
    slong Tlen = FLINT_MAX(Alen, Blen);
    mp_limb_t d0 = (1 + mod.n)/2;
    mp_limb_t d1 = nmod_inv(nmod_add(alpha, alpha, mod), mod);
    mp_limb_t * u, u1nonzero, u0nonzero;

    u = FLINT_ARRAY_ALLOC(2*d, mp_limb_t);

    n_bpoly_fit_length(T, Tlen);

    Tcoeffs = T->coeffs;

    for (i = 0; i < Tlen; i++)
    {
        _nmod_vec_zero(u, 2*d);

        if (i < Alen && i < Blen)
        {
            u0nonzero = u1nonzero = 0;
            for (j = 0; j < d; j++)
            {
                ulong t0 = nmod_add(Acoeffs[d*i + j], Bcoeffs[d*i + j], mod);
                ulong t1 = nmod_sub(Acoeffs[d*i + j], Bcoeffs[d*i + j], mod);
                u[d*0 + j] = t0;
                u[d*1 + j] = t1;
                u1nonzero |= t1;
                u0nonzero |= t0;
            }
        }
        else if (i < Alen)
        {
            u0nonzero = 0;
            for (j = 0; j < d; j++)
            {
                ulong t0 = Acoeffs[d*i + j];
                u0nonzero |= t0;
                u[d*0 + j] = t0;
                u[d*1 + j] = t0;
            }
            u1nonzero = u0nonzero;
        }
        else
        {
            u0nonzero = 0;
            for (j = 0; j < d; j++)
            {
                ulong t0 = Bcoeffs[d*i + j];
                u0nonzero |= t0;
                u[d*0 + j] = t0;
                u[d*1 + j] = nmod_neg(t0, mod);
            }
            u1nonzero = u0nonzero;
        }

        if (u1nonzero | u0nonzero)
        {
            n_poly_fit_length(Tcoeffs + i, d*2);
            _nmod_vec_scalar_mul_nmod(Tcoeffs[i].coeffs + d*0, u + d*0, d, d0, mod);
            if (u1nonzero)
            {
                _nmod_vec_scalar_mul_nmod(Tcoeffs[i].coeffs + d*1, u + d*1, d, d1, mod);
                Tcoeffs[i].length = 2;
            }
            else
            {
                Tcoeffs[i].length = 1;
            }
            lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
        }
        else
        {
            Tcoeffs[i].length = 0;
        }
    }

    *deg1 = lastlength - 1;

    flint_free(u);

    FLINT_ASSERT(Tlen < 1 || !n_fq_poly_is_zero(Tcoeffs + Tlen - 1));
    T->length = Tlen;

    FLINT_ASSERT(n_fq_bpoly_is_canonical(T, ctx));
}

void _n_fq_poly_addmul_plinear(
    n_fq_poly_t A,
    mp_limb_t * Bcoeffs, slong Blen,
    const n_poly_t C,
    mp_limb_t * s,
    slong d,
    nmod_t mod)
{
    slong i, j;
    mp_limb_t * Acoeffs;
    mp_limb_t * Ccoeffs = C->coeffs;
    slong Clen = C->length;
    slong Alen = FLINT_MAX(Blen, Clen + 1);

    n_poly_fit_length(A, d*Alen);
    Acoeffs = A->coeffs;

    for (i = 0; i < Alen; i++)
    {
        for (j = 0; j < d; j++)
        {
            ulong p1, p0, t0 = 0, t1 = 0, t2 = 0;

            if (i < Blen)
            {
                t0 = Bcoeffs[d*i + j];
            }
            if (i < Clen)
            {
                umul_ppmm(p1, p0, Ccoeffs[i], s[0*d + j]);
                add_ssaaaa(t1, t0, t1, t0, p1, p0);
            }
            if (0 < i && i - 1 < Clen)
            {
                umul_ppmm(p1, p0, Ccoeffs[i - 1], s[1*d + j]);
                add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
            }
            NMOD_RED3(Acoeffs[i*d + j], t2, t1, t0, mod);
        }
    }

    A->length = Alen;
    _n_fq_poly_normalise(A, d);
}

int n_fq_bpoly_interp_crt_2psm_poly(
    slong * deg1,
    n_fq_bpoly_t F,
    n_fq_bpoly_t T,
    n_fq_poly_t A,
    n_fq_poly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    int changed = 0;
    slong i, j, lastlength = 0;
    slong Alen = A->length;
    slong Blen = B->length;
    slong Flen = F->length;
    slong Tlen = FLINT_MAX(FLINT_MAX(Alen, Blen), Flen);
    n_fq_poly_struct * Tcoeffs, * Fcoeffs;
    mp_limb_t * Acoeffs, * Bcoeffs;
    mp_limb_t * u, unonzero;
    mp_limb_t malpha = mod.n - alphapow->coeffs[1];

    n_bpoly_fit_length(T, Tlen);
    Tcoeffs = T->coeffs;
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Fcoeffs = F->coeffs;

    u = FLINT_ARRAY_ALLOC(2*d, mp_limb_t);

    for (i = 0; i < Tlen; i++)
    {
        if (i < Flen)
            n_fq_poly_eval2p_pow(u + d*0, u + d*1, Fcoeffs + i, alphapow, d, mod);
        else
            _nmod_vec_zero(u, 2*d);

        if (i < Alen)
            _nmod_vec_sub(u + d*0, u + d*0, Acoeffs + d*i, d, mod);

        if (i < Blen)
            _nmod_vec_sub(u + d*1, u + d*1, Bcoeffs + d*i, d, mod);

        unonzero = 0;
        for (j = 0; j < d; j++)
        {
            mp_limb_t t1 = nmod_sub(u[d*1 + j], u[d*0 + j], mod);
            mp_limb_t t0 = nmod_add(u[d*1 + j], u[d*0 + j], mod);
            u[d*1 + j] = t1;
            unonzero |= u[d*1 + j];
            u[d*0 + j] = nmod_mul(malpha, t0, mod);
            unonzero |= u[d*0 + j];
        }                

        if (unonzero)
        {
            mp_limb_t * Ficoeffs = i < Flen ? Fcoeffs[i].coeffs : NULL;
            slong Filen = i < Flen ? Fcoeffs[i].length : 0;
            _n_fq_poly_addmul_plinear(Tcoeffs + i, Ficoeffs, Filen, modulus, u, d, mod);
            changed = 1;
        }
        else
        {
            if (i < Flen)
                n_fq_poly_set(Tcoeffs + i, Fcoeffs + i, ctx);
            else
                n_fq_poly_zero(Tcoeffs + i);
        }

        lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
    }

    FLINT_ASSERT(i < 1 || !n_fq_poly_is_zero(Tcoeffs + i - 1));
    T->length = i;

    if (changed)
        n_bpoly_swap(T, F);

    FLINT_ASSERT(n_fq_bpoly_is_canonical(F, ctx));

    *deg1 = lastlength - 1;

    flint_free(u);

    FLINT_ASSERT(n_fq_bpoly_is_canonical(T, ctx));

    return changed;
}


int n_fq_bpoly_gcd_brown_smprime2p(
    n_fq_bpoly_t G,
    n_fq_bpoly_t Abar,
    n_fq_bpoly_t Bbar,
    const n_fq_bpoly_t A,
    const n_fq_bpoly_t B,
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t Sp,
    n_fq_poly_t cA,
    n_fq_poly_t cB,
    n_fq_poly_t cG,
    n_fq_poly_t cAbar,
    n_fq_poly_t cBbar,
    n_fq_poly_t gamma,
    n_fq_poly_t r)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_mod(ctx);
    int success;
    slong bound;
    mp_limb_t alpha, temp, * gammaevalp, * gammaevalm;
    n_fq_poly_struct * Aevalp, * Bevalp, * Gevalp, * Abarevalp, * Bbarevalp;
    n_fq_poly_struct * Aevalm, * Bevalm, * Gevalm, * Abarevalm, * Bbarevalm;
    n_bpoly_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * modulus, * alphapow;
    int gstab, astab, bstab, use_stab;

    ldegA = n_fq_bpoly_degree1(A);
    ldegB = n_fq_bpoly_degree1(B);
    deggamma = n_fq_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    if (bound >= mod.n/2)
        return 0;

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    gammaevalp = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    gammaevalm = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_poly_stack_fit_request(Sp->poly_stack, 12);
    Aevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Bevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Gevalp      = n_poly_stack_take_top(Sp->poly_stack);
    Abarevalp   = n_poly_stack_take_top(Sp->poly_stack);
    Bbarevalp   = n_poly_stack_take_top(Sp->poly_stack);
    Aevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Bevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Gevalm      = n_poly_stack_take_top(Sp->poly_stack);
    Abarevalm   = n_poly_stack_take_top(Sp->poly_stack);
    Bbarevalm   = n_poly_stack_take_top(Sp->poly_stack);
    alphapow    = n_poly_stack_take_top(Sp->poly_stack);
    modulus     = n_poly_stack_take_top(Sp->poly_stack);

    n_bpoly_stack_fit_request(Sp->bpoly_stack, 1);
    T           = n_bpoly_stack_take_top(Sp->bpoly_stack);

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    use_stab = 1;
    gstab = bstab = astab = 0;

    if ((mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    alpha = (mod.n - UWORD(1))/UWORD(2);

choose_prime:

    if (alpha < 2)
    {
        success = 0;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_fq_poly_eval2p_pow(gammaevalp, gammaevalm, gamma, alphapow, d, mod);
    if (_n_fq_is_zero(gammaevalp, d) || _n_fq_is_zero(gammaevalm, d))
        goto choose_prime;

    n_fq_bpoly_interp_reduce_2psm_poly(Aevalp, Aevalm, A, alphapow, ctx);
    n_fq_bpoly_interp_reduce_2psm_poly(Bevalp, Bevalm, B, alphapow, ctx);
    FLINT_ASSERT(Aevalp->length > 0);
    FLINT_ASSERT(Aevalm->length > 0);
    FLINT_ASSERT(Bevalp->length > 0);
    FLINT_ASSERT(Bevalm->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        n_fq_bpoly_interp_reduce_2psm_poly(Gevalp, Gevalm, G, alphapow, ctx);
        Gdeg = n_fq_bpoly_degree0(G);
        success = 1;
        success = success && n_fq_poly_degree(Gevalp) == Gdeg;
        success = success && n_fq_poly_degree(Gevalm) == Gdeg;
        success = success && _n_fq_equal(Gevalp->coeffs + d*Gdeg, gammaevalp, d);
        success = success && _n_fq_equal(Gevalm->coeffs + d*Gdeg, gammaevalm, d);
        n_fq_poly_divrem_(Abarevalp, r, Aevalp, Gevalp, ctx, Sp->poly_stack);
        success = success && (r->length == 0);
        n_fq_poly_divrem_(Abarevalm, r, Aevalm, Gevalm, ctx, Sp->poly_stack);
        success = success && (r->length == 0);
        n_fq_poly_divrem_(Bbarevalp, r, Bevalp, Gevalp, ctx, Sp->poly_stack);
        success = success && (r->length == 0);
        n_fq_poly_divrem_(Bbarevalm, r, Bevalm, Gevalm, ctx, Sp->poly_stack);
        success = success && (r->length == 0);

        if (!success)
        {
            use_stab = 0;
            n_poly_one(modulus);
            alpha = (fq_nmod_ctx_mod(ctx).n - UWORD(1))/UWORD(2);
            goto choose_prime;
        }

        n_fq_poly_scalar_mul_n_fq(Abarevalp, Abarevalp, gammaevalp, ctx);
        n_fq_poly_scalar_mul_n_fq(Abarevalm, Abarevalm, gammaevalm, ctx);
        n_fq_poly_scalar_mul_n_fq(Bbarevalp, Bbarevalp, gammaevalp, ctx);
        n_fq_poly_scalar_mul_n_fq(Bbarevalm, Bbarevalm, gammaevalm, ctx);
    }
    else
    {
        n_fq_poly_gcd_(Gevalp, Aevalp, Bevalp, ctx, Sp->poly_stack);
        n_fq_poly_divrem_(Abarevalp, r, Aevalp, Gevalp, ctx, Sp->poly_stack);
        n_fq_poly_divrem_(Bbarevalp, r, Bevalp, Gevalp, ctx, Sp->poly_stack);
        n_fq_poly_gcd_(Gevalm, Aevalm, Bevalm, ctx, Sp->poly_stack);
        n_fq_poly_divrem_(Abarevalm, r, Aevalm, Gevalm, ctx, Sp->poly_stack);
        n_fq_poly_divrem_(Bbarevalm, r, Bevalm, Gevalm, ctx, Sp->poly_stack);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Gevalp->length > 0);
    FLINT_ASSERT(Abarevalp->length > 0);
    FLINT_ASSERT(Bbarevalp->length > 0);
    FLINT_ASSERT(Gevalm->length > 0);
    FLINT_ASSERT(Abarevalm->length > 0);
    FLINT_ASSERT(Bbarevalm->length > 0);

    if (n_fq_poly_degree(Gevalp) == 0 || n_fq_poly_degree(Gevalm) == 0)
    {
        n_fq_bpoly_one(G, ctx);
        n_fq_bpoly_set(Abar, A, ctx);
        n_fq_bpoly_set(Bbar, B, ctx);
        goto successful_put_content;    
    }

    if (n_fq_poly_degree(Gevalp) != n_fq_poly_degree(Gevalm))
    {
        goto choose_prime;
    }

    if (n_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);
        if (n_fq_poly_degree(Gevalp) > n_fq_bpoly_degree0(G))
        {
            goto choose_prime;
        }
        else if (n_fq_poly_degree(Gevalp) < n_fq_bpoly_degree0(G))
        {
            n_poly_one(modulus);
        }
    }

    n_fq_poly_scalar_mul_n_fq(Gevalp, Gevalp, gammaevalp, ctx);
    n_fq_poly_scalar_mul_n_fq(Gevalm, Gevalm, gammaevalm, ctx);
    if (n_poly_degree(modulus) > 0)
    {
        temp = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
        FLINT_ASSERT(temp == n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod));
        temp = nmod_mul(temp, alpha, mod);
        temp = nmod_add(temp, temp, mod);
        temp = nmod_inv(temp, mod);
        _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, mod);

        gstab = gstab || !n_fq_bpoly_interp_crt_2psm_poly(&ldegG, G, T,
                                       Gevalp, Gevalm, modulus, alphapow, ctx);
        n_fq_bpoly_interp_crt_2psm_poly(&ldegAbar, Abar, T,
                                 Abarevalp, Abarevalm, modulus, alphapow, ctx);
        n_fq_bpoly_interp_crt_2psm_poly(&ldegBbar, Bbar, T,
                                 Bbarevalp, Bbarevalm, modulus, alphapow, ctx);
    }
    else
    {
        n_fq_bpoly_interp_lift_2psm_poly(&ldegG, G, Gevalp, Gevalm, alpha, ctx);
        n_fq_bpoly_interp_lift_2psm_poly(&ldegAbar, Abar,
                                             Abarevalp, Abarevalm, alpha, ctx);
        n_fq_bpoly_interp_lift_2psm_poly(&ldegBbar, Bbar,
                                             Bbarevalp, Bbarevalm, alpha, ctx);
        gstab = astab = bstab = 0;
    }

    temp = mod.n - nmod_mul(alpha, alpha, mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, temp, mod);

    if (n_poly_degree(modulus) < bound)
        goto choose_prime;

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    n_poly_one(modulus);
    goto choose_prime;

successful:

    n_fq_bpoly_content_var0(r, G, ctx);
    n_fq_bpoly_divexact_poly_var1(G, r, ctx);
    n_fq_bpoly_divexact_poly_var1(Abar, G->coeffs + G->length - 1, ctx);
    n_fq_bpoly_divexact_poly_var1(Bbar, G->coeffs + G->length - 1, ctx);

successful_put_content:

    n_fq_bpoly_mul_last(G, cG, ctx);
    n_fq_bpoly_mul_last(Abar, cAbar, ctx);
    n_fq_bpoly_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

    FLINT_ASSERT(n_fq_bpoly_is_canonical(G, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(Abar, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(Bbar, ctx));

    flint_free(gammaevalp);
    flint_free(gammaevalm);

    n_poly_stack_give_back(Sp->poly_stack, 12);
    n_bpoly_stack_give_back(Sp->bpoly_stack, 1);

    return success;
}


/*****************************************************************************/

void n_fq_bpoly_interp_reduce_sm_poly(
    n_fq_poly_t E,
    const n_fq_bpoly_t A,
    n_fq_poly_t alphapow,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i, Alen = A->length;
    const n_fq_poly_struct * Ac = A->coeffs;
    mp_limb_t * Ec;

    n_poly_fit_length(E, d*Alen);
    Ec = E->coeffs;

    for (i = 0; i < Alen; i++)
        n_fq_poly_eval_pow(Ec + d*i, Ac + i, alphapow, ctx);

    E->length = Alen;
    _n_fq_poly_normalise(E, d);
}

void n_fq_bpoly_interp_lift_sm_poly(
    n_fq_bpoly_t T,
    const n_fq_poly_t A,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    slong i;
    const mp_limb_t * Acoeffs = A->coeffs;
    n_poly_struct * Tcoeffs;
    slong Alen = A->length;

    n_fq_bpoly_fit_length(T, Alen);

    Tcoeffs = T->coeffs;

    for (i = 0; i < Alen; i++)
    {
        n_fq_poly_set_n_fq(Tcoeffs + i, Acoeffs + d*i, ctx);
    }

    FLINT_ASSERT(i < 1 || !n_fq_poly_is_zero(Tcoeffs + i - 1));
    T->length = i;
}


int n_fq_bpoly_interp_crt_sm_poly(
    slong * deg1,
    n_fq_bpoly_t F,
    n_fq_bpoly_t T,
    n_fq_poly_t A,
    const n_fq_poly_t modulus,
    n_fq_poly_t alphapow,
    const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int changed = 0;
    slong i, lastlength = 0;
    slong Alen = A->length;
    slong Flen = F->length;
    n_fq_poly_struct * Tcoeffs, * Fcoeffs;
    mp_limb_t * Acoeffs;
    mp_limb_t * u, * v;

    FLINT_ASSERT(n_fq_bpoly_is_canonical(F, ctx));
    FLINT_ASSERT(n_fq_poly_is_canonical(A, ctx));

    u = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    v = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_fq_bpoly_fit_length(T, FLINT_MAX(Alen, Flen));
    Tcoeffs = T->coeffs;
    Acoeffs = A->coeffs;
    Fcoeffs = F->coeffs;

    for (i = 0; i < Flen; i++)
    {
        /* F term ok, A term ok/missing */
        n_fq_poly_eval_pow(u, Fcoeffs + i, alphapow, ctx);

        if (i < Alen)
            n_fq_sub(v, Acoeffs + d*i, u, ctx);
        else
            _n_fq_neg(v, u, d, ctx->mod);

        if (!_n_fq_is_zero(v, d))
        {
            changed = 1;
            n_fq_poly_scalar_addmul_n_fq(Tcoeffs + i, Fcoeffs + i, modulus, v, ctx);
        }
        else
        {
            n_fq_poly_set(Tcoeffs + i, Fcoeffs + i, ctx);
        }

        lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
    }

    for ( ; i < Alen; i++)
    {
        /* F term missing, A term ok */
        if (!_n_fq_is_zero(Acoeffs + d*i, d))
        {
            changed = 1;
            n_fq_poly_scalar_mul_n_fq(Tcoeffs + i, modulus, Acoeffs + d*i, ctx);
        }
        else
        {
            n_fq_poly_zero(Tcoeffs + i);
        }

        lastlength = FLINT_MAX(lastlength, Tcoeffs[i].length);
    }

    flint_free(u);
    flint_free(v);

    FLINT_ASSERT(i < 1 || !n_fq_poly_is_zero(Tcoeffs + i - 1));
    T->length = i;

    if (changed)
        n_bpoly_swap(T, F);

    FLINT_ASSERT(n_fq_bpoly_is_canonical(F, ctx));

    *deg1 = lastlength - 1;
    return changed;
}


int n_fq_bpoly_gcd_brown_smprime(
    n_fq_bpoly_t G,
    n_fq_bpoly_t Abar,
    n_fq_bpoly_t Bbar,
    n_fq_bpoly_t A,
    n_fq_bpoly_t B,
    const fq_nmod_ctx_t ctx,
    n_poly_bpoly_stack_t Sp)
{
    slong d = fq_nmod_ctx_degree(ctx);
    int success;
    slong bound;
    fq_nmod_t alpha;
    mp_limb_t * temp, * gammaeval;
    n_poly_struct * Aeval, * Beval, * Geval, * Abareval, * Bbareval;
    n_bpoly_struct * T;
    slong deggamma, ldegG, ldegAbar, ldegBbar, ldegA, ldegB;
    n_poly_struct * cA, * cB, * cG, * cAbar, * cBbar, * gamma;
    n_poly_struct * modulus, * alphapow, * r;
    int gstab, astab, bstab, use_stab;
#if FLINT_WANT_ASSERT
    n_fq_bpoly_t Asave, Bsave;
    const slong Sp_size_poly = n_poly_stack_size(Sp->poly_stack);
    const slong Sp_size_bpoly = n_bpoly_stack_size(Sp->bpoly_stack);

    n_fq_bpoly_init(Asave);
    n_fq_bpoly_init(Bsave);
    n_fq_bpoly_set(Asave, A, ctx);
    n_fq_bpoly_set(Bsave, B, ctx);
#endif

    FLINT_ASSERT(n_fq_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(n_fq_bpoly_is_canonical(B, ctx));

    n_poly_stack_fit_request(Sp->poly_stack, 7);
    cA      = n_poly_stack_take_top(Sp->poly_stack);
    cB      = n_poly_stack_take_top(Sp->poly_stack);
    cG      = n_poly_stack_take_top(Sp->poly_stack);
    cAbar   = n_poly_stack_take_top(Sp->poly_stack);
    cBbar   = n_poly_stack_take_top(Sp->poly_stack);
    gamma   = n_poly_stack_take_top(Sp->poly_stack);
    r       = n_poly_stack_take_top(Sp->poly_stack);

    n_fq_bpoly_content_var0(cA, A, ctx);
    n_fq_bpoly_content_var0(cB, B, ctx);
    n_fq_bpoly_divexact_poly_var1(A, cA, ctx);
    n_fq_bpoly_divexact_poly_var1(B, cB, ctx);

    n_fq_poly_gcd(cG, cA, cB, ctx);
    n_fq_poly_divrem(cAbar, r, cA, cG, ctx);
    n_fq_poly_divrem(cBbar, r, cB, cG, ctx);

    n_fq_poly_gcd(gamma, A->coeffs + A->length - 1, B->coeffs + B->length - 1, ctx);


    if (n_fq_bpoly_gcd_brown_smprime2p(G, Abar, Bbar, A, B, ctx, Sp,
                                           cA, cB, cG, cAbar, cBbar, gamma, r))
    {
#if FLINT_WANT_ASSERT
        {
            n_fq_poly_struct * Glead = G->coeffs + G->length - 1;
            n_fq_bpoly_t P;
            FLINT_ASSERT(_n_fq_is_one(Glead->coeffs + d*(Glead->length - 1), d));
            n_fq_bpoly_init(P);
            n_fq_bpoly_mul(P, G, Abar, ctx);
            FLINT_ASSERT(n_fq_bpoly_equal(P, Asave, ctx));
            n_fq_bpoly_mul(P, G, Bbar, ctx);
            FLINT_ASSERT(n_fq_bpoly_equal(P, Bsave, ctx));
            n_fq_bpoly_clear(P);
        }
        n_fq_bpoly_clear(Asave);
        n_fq_bpoly_clear(Bsave);
#endif

        n_poly_stack_give_back(Sp->poly_stack, 7);
        FLINT_ASSERT(Sp_size_poly == n_poly_stack_size(Sp->poly_stack));
        FLINT_ASSERT(Sp_size_bpoly == n_bpoly_stack_size(Sp->bpoly_stack));
        return 1;
    }

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    fq_nmod_init(alpha, ctx);
    temp = FLINT_ARRAY_ALLOC(d, mp_limb_t);
    gammaeval = FLINT_ARRAY_ALLOC(d, mp_limb_t);

    n_poly_stack_fit_request(Sp->poly_stack, 7);
    Aeval       = n_poly_stack_take_top(Sp->poly_stack);
    Beval       = n_poly_stack_take_top(Sp->poly_stack);
    Geval       = n_poly_stack_take_top(Sp->poly_stack);
    Abareval    = n_poly_stack_take_top(Sp->poly_stack);
    Bbareval    = n_poly_stack_take_top(Sp->poly_stack);
    alphapow    = n_poly_stack_take_top(Sp->poly_stack);
    modulus     = n_poly_stack_take_top(Sp->poly_stack);

    n_bpoly_stack_fit_request(Sp->bpoly_stack, 1);
    T           = n_bpoly_stack_take_top(Sp->bpoly_stack);

    ldegA = n_fq_bpoly_degree1(A);
    ldegB = n_fq_bpoly_degree1(B);
    deggamma = n_fq_poly_degree(gamma);
    bound = 1 + deggamma + FLINT_MAX(ldegA, ldegB);

    n_poly_fit_length(alphapow, d*FLINT_MAX(WORD(3), bound + 1));
    n_fq_poly_one(modulus, ctx);

    use_stab = 1;
    gstab = bstab = astab = 0;

    fq_nmod_zero(alpha, ctx);

choose_prime:

    if (fq_nmod_next(alpha, ctx) == 0)
    {
        success = 0;
        goto cleanup;
    }

    alphapow->length = 2;
    _n_fq_one(alphapow->coeffs + d*0, d);
    n_fq_set_fq_nmod(alphapow->coeffs + d*1, alpha, ctx);

    n_fq_poly_eval_pow(gammaeval, gamma, alphapow, ctx);
    if (_n_fq_is_zero(gammaeval, d))
        goto choose_prime;

    n_fq_bpoly_interp_reduce_sm_poly(Aeval, A, alphapow, ctx);
    n_fq_bpoly_interp_reduce_sm_poly(Beval, B, alphapow, ctx);
    FLINT_ASSERT(Aeval->length > 0);
    FLINT_ASSERT(Beval->length > 0);

    if (use_stab && gstab)
    {
        slong Gdeg;
        n_fq_bpoly_interp_reduce_sm_poly(Geval, G, alphapow, ctx);
        Gdeg = n_fq_bpoly_degree0(G);
        success = 1;
        success = success && n_fq_poly_degree(Geval) == Gdeg;
        success = success && _n_fq_equal(Geval->coeffs + d*Gdeg,  gammaeval, d);
        n_fq_poly_divrem(Abareval, r, Aeval, Geval, ctx);
        success = success && n_fq_poly_is_zero(r);
        n_fq_poly_divrem(Bbareval, r, Beval, Geval, ctx);
        success = success && n_fq_poly_is_zero(r);

        if (!success)
        {
            use_stab = 0;
            n_fq_poly_one(modulus, ctx);
            fq_nmod_zero(alpha, ctx);
            goto choose_prime;
        }

        n_fq_poly_scalar_mul_n_fq(Abareval, Abareval, gammaeval, ctx);
        n_fq_poly_scalar_mul_n_fq(Bbareval, Bbareval, gammaeval, ctx);
    }
    else
    {
        n_fq_poly_gcd(Geval, Aeval, Beval, ctx);
        n_fq_poly_divrem(Abareval, r, Aeval, Geval, ctx);
        n_fq_poly_divrem(Bbareval, r, Beval, Geval, ctx);
        gstab = astab = bstab = 0;
    }

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Abareval->length > 0);
    FLINT_ASSERT(Bbareval->length > 0);

    if (n_fq_poly_degree(Geval) == 0)
    {
        n_fq_bpoly_one(G, ctx);
        n_fq_bpoly_set(Abar, A, ctx);
        n_fq_bpoly_set(Bbar, B, ctx);
        goto successful_put_content;    
    }

    if (n_fq_poly_degree(modulus) > 0)
    {
        FLINT_ASSERT(G->length > 0);

        if (n_fq_poly_degree(Geval) > n_fq_bpoly_degree0(G))
            goto choose_prime;

        if (n_fq_poly_degree(Geval) < n_fq_bpoly_degree0(G))
            n_fq_poly_one(modulus, ctx);
    }

    n_fq_poly_scalar_mul_n_fq(Geval, Geval, gammaeval, ctx);

    if (n_fq_poly_degree(modulus) > 0)
    {
        n_fq_poly_eval_pow(temp, modulus, alphapow, ctx);
        n_fq_inv(temp, temp, ctx);
        n_fq_poly_scalar_mul_n_fq(modulus, modulus, temp, ctx);

        gstab = gstab || !n_fq_bpoly_interp_crt_sm_poly(&ldegG, G, T,
                                                Geval, modulus, alphapow, ctx);

        n_fq_bpoly_interp_crt_sm_poly(&ldegAbar, Abar, T,
                                             Abareval, modulus, alphapow, ctx);
        n_fq_bpoly_interp_crt_sm_poly(&ldegBbar, Bbar, T,
                                             Bbareval, modulus, alphapow, ctx);
    }
    else
    {
        n_fq_bpoly_interp_lift_sm_poly(G, Geval, ctx);
        n_fq_bpoly_interp_lift_sm_poly(Abar, Abareval, ctx);
        n_fq_bpoly_interp_lift_sm_poly(Bbar, Bbareval, ctx);
        ldegG = ldegAbar = ldegBbar = 0;
        gstab = astab = bstab = 0;
    }

    n_fq_set_fq_nmod(temp, alpha, ctx);
    n_fq_poly_shift_left_scalar_submul(modulus, 1, temp, ctx);

    if (n_fq_poly_degree(modulus) < bound)
        goto choose_prime;

    FLINT_ASSERT(ldegG >= 0);
    FLINT_ASSERT(ldegAbar >= 0);
    FLINT_ASSERT(ldegBbar >= 0);

    if (deggamma + ldegA == ldegG + ldegAbar &&
        deggamma + ldegB == ldegG + ldegBbar)
    {
        goto successful;
    }

    n_fq_poly_one(modulus, ctx);
    goto choose_prime;

successful:

    n_fq_bpoly_content_var0(modulus, G, ctx);
    n_fq_bpoly_divexact_poly_var1(G, modulus, ctx);
    n_fq_bpoly_divexact_poly_var1(Abar, G->coeffs + G->length - 1, ctx);
    n_fq_bpoly_divexact_poly_var1(Bbar, G->coeffs + G->length - 1, ctx);

successful_put_content:

    n_fq_bpoly_mul_last(G, cG, ctx);
    n_fq_bpoly_mul_last(Abar, cAbar, ctx);
    n_fq_bpoly_mul_last(Bbar, cBbar, ctx);

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success)
    {
        n_fq_poly_struct * Glead = G->coeffs + G->length - 1;
        n_fq_bpoly_t P;
        FLINT_ASSERT(_n_fq_is_one(Glead->coeffs + d*(Glead->length - 1), d));
        n_fq_bpoly_init(P);
        n_fq_bpoly_mul(P, G, Abar, ctx);
        FLINT_ASSERT(n_fq_bpoly_equal(P, Asave, ctx));
        n_fq_bpoly_mul(P, G, Bbar, ctx);
        FLINT_ASSERT(n_fq_bpoly_equal(P, Bsave, ctx));
        n_fq_bpoly_clear(P);
    }
    n_fq_bpoly_clear(Asave);
    n_fq_bpoly_clear(Bsave);
#endif

    fq_nmod_clear(alpha, ctx);
    flint_free(temp);
    flint_free(gammaeval);

    n_poly_stack_give_back(Sp->poly_stack, 14);
    n_bpoly_stack_give_back(Sp->bpoly_stack, 1);
    FLINT_ASSERT(Sp_size_poly == n_poly_stack_size(Sp->poly_stack));
    FLINT_ASSERT(Sp_size_bpoly == n_bpoly_stack_size(Sp->bpoly_stack));

    return success;
}
