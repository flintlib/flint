/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mod_vec.h"
#include "ulong_extras.h"

/* return n with ||a||_2 < 2^n */
static flint_bitcnt_t fmpz_poly_norm2_bits(const fmpz_poly_t a)
{
    fmpz_t t;
    flint_bitcnt_t bits;
    fmpz_init(t);
    _fmpz_vec_dot(t, a->coeffs, a->coeffs, a->length);
    bits = fmpz_bits(t);
    fmpz_clear(t);
    return (bits + 1)/2;
}

/* reduce A mod B, return new length of A */
static slong _reduce_inplace(
    fmpz * Acoeffs, slong Alen,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t Binv,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_t Q,  /* temps */
    fmpz_mod_poly_t R)
{
    slong Blen = B->length;
    fmpz * Bcoeffs = B->coeffs;
    fmpz * Qcoeffs, * Rcoeffs;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);

    if (Alen < Blen)
        return Alen;

    FLINT_ASSERT(Blen > 1);

    fmpz_mod_poly_fit_length(Q, Alen - Blen + 1, ctx);
    fmpz_mod_poly_fit_length(R, Blen - 1, ctx);
    Qcoeffs = Q->coeffs;
    Rcoeffs = R->coeffs;

    while (Alen >= Blen)
    {
        slong n = FLINT_MAX(0, (Alen-1) - 2*(Blen-1) + 1);
        slong Qlen = (Alen-1) - n - (Blen-1) + 1;
        FLINT_ASSERT(Qlen < Blen);

        _fmpz_mod_poly_div_newton_n_preinv(Qcoeffs + n, Acoeffs + n, Alen - n,
                                 Bcoeffs, Blen, Binv->coeffs, Binv->length, p);

        _fmpz_mod_poly_mullow(Rcoeffs, Bcoeffs, Blen - 1, Qcoeffs + n, Qlen,
                                                                  p, Blen - 1);

        _fmpz_mod_vec_sub(Acoeffs + n, Acoeffs + n, Rcoeffs, Blen - 1, ctx);

        Alen = n + Blen - 1;
        while (Alen > 0 && fmpz_is_zero(Acoeffs + Alen - 1))
            Alen--;
    }

    return Alen;
}

void fmpz_poly_pfrac_init(fmpz_poly_pfrac_t I)
{
    I->r = 0;
    I->bits = NULL;
    I->b = I->bprod = NULL;
    I->B = I->invBprod = I->inwBprod = I->B_inv = NULL;
    I->ctxs = NULL;
    I->halfpks = NULL;

    fmpz_poly_init(I->a);
    fmpz_poly_init(I->newa);
    fmpz_poly_init(I->t);

    fmpz_init(I->old_pk);
    fmpz_init(I->pk);
    fmpz_init(I->p);

    fmpz_mod_ctx_init_ui(I->ctxp, 2);

    fmpz_mod_poly_init(I->T, I->ctxp);
    fmpz_mod_poly_init(I->R, I->ctxp);
    fmpz_mod_poly_init(I->Q, I->ctxp);
}

static void _clear_arrays(fmpz_poly_pfrac_t I)
{
    slong i;

    for (i = 0; i < I->r; i++)
    {
        fmpz_poly_clear(I->b + i);
        fmpz_poly_clear(I->bprod + i);
        fmpz_mod_poly_clear(I->B + i, I->ctxs + i);
        fmpz_mod_poly_clear(I->invBprod + i, I->ctxs + i);
        fmpz_mod_poly_clear(I->inwBprod + i, I->ctxs + i);
        fmpz_mod_poly_clear(I->B_inv + i, I->ctxs + i);
        fmpz_clear(I->halfpks + i);
        fmpz_mod_ctx_clear(I->ctxs + i);
    }
    flint_free(I->halfpks);
    flint_free(I->ctxs);
    flint_free(I->bits);
    flint_free(I->b);
    flint_free(I->B);
    I->r = 0;
}

void fmpz_poly_pfrac_clear(fmpz_poly_pfrac_t I)
{
    _clear_arrays(I);

    fmpz_poly_clear(I->a);
    fmpz_poly_clear(I->newa);
    fmpz_poly_clear(I->t);

    fmpz_clear(I->old_pk);
    fmpz_clear(I->pk);
    fmpz_clear(I->p);

    fmpz_mod_poly_clear(I->T, I->ctxp);
    fmpz_mod_poly_clear(I->R, I->ctxp);
    fmpz_mod_poly_clear(I->Q, I->ctxp);

    fmpz_mod_ctx_clear(I->ctxp);
}

int fmpz_poly_pfrac_precompute(
    fmpz_poly_pfrac_t I,
    const fmpz_poly_struct * b,
    slong r)
{
    slong i;

    if (r < 2)
        return 0;

    for (i = 0; i < r; i++)
    {
        if (fmpz_poly_degree(b + i) < 1)
            return 0;
    }

    _clear_arrays(I);

    I->r = r;
    I->bits = FLINT_ARRAY_ALLOC(r, flint_bitcnt_t);

    I->ctxs = FLINT_ARRAY_ALLOC(r, fmpz_mod_ctx_struct);
    I->halfpks = FLINT_ARRAY_ALLOC(r, fmpz);
    for (i = 0; i < r; i++)
    {
        fmpz_init(I->halfpks + i);
        fmpz_mod_ctx_init_ui(I->ctxs + i, 2);
    }

    I->b = FLINT_ARRAY_ALLOC(2*r, fmpz_poly_struct);
    I->bprod = I->b + r;
    for (i = 0; i < r; i++)
    {
        fmpz_poly_init(I->bprod + i);
        fmpz_poly_init(I->b + i);
        fmpz_poly_set(I->b + i, b + i);
    }

    I->B = FLINT_ARRAY_ALLOC(4*r, fmpz_mod_poly_struct);
    I->invBprod = I->B + r;
    I->inwBprod = I->invBprod + r;
    I->B_inv = I->inwBprod + r;

    for (i = 0; i < r; i++)
    {
        fmpz_mod_poly_init(I->B + i, I->ctxs + i);
        fmpz_mod_poly_init(I->invBprod + i, I->ctxs + i);
        fmpz_mod_poly_init(I->inwBprod + i, I->ctxs + i);
        fmpz_mod_poly_init(I->B_inv + i, I->ctxs + i);
    }

    /* init done */

    fmpz_poly_one(I->bprod + r - 1);
    for (i = r - 2; i >= 0; i--)
    {
        fmpz_poly_mul(I->bprod + i, I->bprod + i + 1, I->b + i + 1);
        I->bits[i] = (fmpz_poly_degree(I->b + i) - 1)*
                                            fmpz_poly_norm2_bits(I->bprod + i);
        I->bits[i] += fmpz_poly_degree(I->bprod + i)*
                                                fmpz_poly_norm2_bits(I->b + i);

        fmpz_poly_resultant(I->pk, I->bprod + i, I->b + i);
        if (fmpz_is_zero(I->pk))
            return 0;

        if (n_sub_checked(&I->bits[i], I->bits[i] + 2, fmpz_bits(I->pk)))
            I->bits[i] = 1;
    }

    fmpz_set_ui(I->p, UWORD(1) << (SMALL_FMPZ_BITCOUNT_MAX));

next_p:

    fmpz_nextprime(I->p, I->p, 1);
    fmpz_mod_ctx_set_modulus(I->ctxp, I->p);
    fmpz_set(I->pk, I->p);

    for (i = 0; i < r; i++)
    {
        /* B[i] = make_monic(b[i] mod p^k) */
        fmpz_mod_poly_set_fmpz_poly(I->B + i, I->b + i, I->ctxp);
        if (I->B[i].length != I->b[i].length)
            goto next_p;

        fmpz_mod_poly_make_monic(I->B + i, I->B + i, I->ctxp);

        fmpz_mod_poly_reverse(I->B_inv + i, I->B + i,
                                                     I->B[i].length, I->ctxp);
        fmpz_mod_poly_inv_series_newton(I->B_inv + i, I->B_inv + i,
                                                     I->B[i].length, I->ctxp);
    }

    for (i = 0; i < r; i++)
    {
        fmpz_mod_poly_set_fmpz_poly(I->T, I->bprod + i, I->ctxp);
        fmpz_mod_poly_xgcd(I->R, I->invBprod + i, I->inwBprod + i,
                                                   I->T, I->B + i, I->ctxp);
        if (!fmpz_mod_poly_is_one(I->R, I->ctxp))
            goto next_p;

        /* now 1 = invBprod[i]*bprod[i] + inwBprod[i]*B[i] mod p^k */
    }

    /* set all ctx's to mod p*/
    for (i = 0; i < r; i++)
    {
        fmpz_mod_ctx_set_modulus(I->ctxs + i, I->p);
        fmpz_fdiv_q_2exp(I->halfpks + i, fmpz_mod_ctx_modulus(I->ctxs + i), 1);
    }

    return 1;
}

/*
set di = deg(Bi)
when solving for the Ci given A and Bi in
    A/(B1*B2) = C1/B1 + C2/B2
the coefficients of C2 and C1 come from a sylvester system
with B1 in the first d2 columns and B2 in the next d1 columns and
the coefficients of A on the rhs. Therefore cramer + hadamard
||C2||_infty <= ||B1||_2^(d2-1) * ||B2||_2^d1 * ||A||_2 / |res(B1,B2)|
||C1||_infty <= ||B2||_2^(d1-1) * ||B1||_2^d2 * ||A||_2 / |res(B1,B2)|

suppose ||B2||_2^(d1-1) * ||B1||_2^d2 / |res(B1,B2)| < 2^n
    and ||A||_2 < 2^a

Assume fmpz_bits(p^k) > n + 1 + a. Then, p^k/2 >= 2^(n + a).
If C1 (correct mod p^k with coeffs <= p^k/2 in abs) does not pass the
divisibility test, then C1 does not exist.
divisibility test is: C2 = (A - B2*C1)/B1
*/
int fmpz_poly_pfrac_precomp(
    fmpz_poly_struct * c,
    const fmpz_poly_t A,
    fmpz_poly_pfrac_t I)
{
    slong i, clen;
    const fmpz_poly_struct * a;

again:

    a = A;

    for (i = 0; i + 1 < I->r; i++)
    {
        /* T = a mod B[i] */
        fmpz_mod_poly_set_fmpz_poly(I->T, a, I->ctxs + i);
        I->T->length = _reduce_inplace(I->T->coeffs, I->T->length,
                              I->B + i, I->B_inv + i, I->ctxs + i, I->Q, I->R);

        /* c = T*invBprod[i] */
        if (I->T->length < 1)
        {
            clen = 0;
        }
        else
        {
            FLINT_ASSERT(I->invBprod[i].length > 0);

            clen = I->T->length + I->invBprod[i].length - 1;
            fmpz_poly_fit_length(c + i, clen);
            _fmpz_mod_poly_mul(c[i].coeffs, I->T->coeffs, I->T->length,
                               I->invBprod[i].coeffs, I->invBprod[i].length,
                                            fmpz_mod_ctx_modulus(I->ctxs + i));
            while (clen > 0 && fmpz_is_zero(c[i].coeffs + clen - 1))
                clen--;
        }

        /* c = c smod B[i] */
        clen = _reduce_inplace(c[i].coeffs, clen,
                              I->B + i, I->B_inv + i, I->ctxs + i, I->Q, I->R);
        c[i].length = clen;
        while (--clen >= 0)
        {
            if (fmpz_cmp(c[i].coeffs + clen, I->halfpks + i) > 0)
            {
                fmpz_sub(c[i].coeffs + clen, c[i].coeffs + clen,
                                            fmpz_mod_ctx_modulus(I->ctxs + i));
            }
        }

        /* now divisibility test */
        fmpz_poly_mul(I->t, c + i, I->bprod + i);
        fmpz_poly_sub(I->t, a, I->t);
        if (!fmpz_poly_divides(I->newa, I->t, I->b + i))
        {
            flint_bitcnt_t abits = fmpz_poly_norm2_bits(a);
            flint_bitcnt_t pkbits = fmpz_bits(fmpz_mod_ctx_modulus(I->ctxs + i));

            if (pkbits > abits && pkbits - abits > I->bits[i])
                return 0;

            goto more_prec;
        }

        a = I->a;
        fmpz_poly_swap(I->a, I->newa);
    }

    FLINT_ASSERT(a == I->a);
    fmpz_poly_swap(c + i, I->a);

    return 1;

more_prec:

    /* increase the precision of the i^th factor only */

    fmpz_set(I->old_pk, fmpz_mod_ctx_modulus(I->ctxs + i));
    fmpz_pow_ui(I->pk, I->p, 1 + fmpz_bits(I->old_pk)/512);
    fmpz_mul(I->halfpks + i, fmpz_mod_ctx_modulus(I->ctxs + i), I->pk);
    fmpz_mod_ctx_set_modulus(I->ctxs + i, I->halfpks + i);
    fmpz_fdiv_q_2exp(I->halfpks + i, fmpz_mod_ctx_modulus(I->ctxs + i), 1);

    fmpz_mod_poly_set_fmpz_poly(I->T, I->bprod + i, I->ctxs + i);
    /* lift_only_inverses wants monic bases */
    fmpz_mod_poly_scalar_div_fmpz(I->T, I->T,
                                    fmpz_poly_lead(I->bprod + i), I->ctxs + i);
    fmpz_mod_poly_scalar_mul_fmpz(I->invBprod + i, I->invBprod + i,
                                    fmpz_poly_lead(I->bprod + i), I->ctxs + i);

    fmpz_mod_poly_set_fmpz_poly(I->B + i, I->b + i, I->ctxs + i);
    fmpz_mod_poly_make_monic(I->B + i, I->B + i, I->ctxs + i);

    fmpz_mod_poly_fit_length(I->invBprod + i, I->B[i].length - 1, I->ctxs + i);
    fmpz_mod_poly_fit_length(I->inwBprod + i, I->T->length - 1, I->ctxs + i);

    _fmpz_poly_hensel_lift_only_inverse(
            I->invBprod[i].coeffs, I->inwBprod[i].coeffs,
            I->T->coeffs, I->T->length, I->B[i].coeffs, I->B[i].length,
            I->invBprod[i].coeffs, I->invBprod[i].length,
            I->inwBprod[i].coeffs, I->inwBprod[i].length,
            I->old_pk, I->pk);

    I->invBprod[i].length = I->B[i].length - 1;
    _fmpz_mod_poly_normalise(I->invBprod + i);

    I->inwBprod[i].length = I->T->length - 1;
    _fmpz_mod_poly_normalise(I->inwBprod + i);

    /* correct monic bases */
    fmpz_mod_poly_scalar_mul_fmpz(I->T, I->T,
                                    fmpz_poly_lead(I->bprod + i), I->ctxs + i);

    fmpz_mod_poly_scalar_div_fmpz(I->invBprod + i, I->invBprod + i,
                                    fmpz_poly_lead(I->bprod + i), I->ctxs + i);

    fmpz_mod_poly_reverse(I->B_inv + i, I->B + i, I->B[i].length, I->ctxs + i);
    fmpz_mod_poly_inv_series_newton(I->B_inv + i, I->B_inv + i,
                                                  I->B[i].length, I->ctxs + i);
    goto again;
}

