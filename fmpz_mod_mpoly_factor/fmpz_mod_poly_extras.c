/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_vec.h"
#include "fmpz_mod_mpoly.h"

/* multiply A by (x^k + c) */
void fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(
    fmpz_mod_poly_t A,
    slong k,
    const fmpz_t c,
    const fmpz_mod_ctx_t ctx)
{
    fmpz * Acoeffs;
    slong i;
    slong Alen = A->length;

    fmpz_mod_poly_fit_length(A, Alen + k, ctx);

    Acoeffs = A->coeffs;

    for (i = Alen - 1; i >= 0; i--)
        fmpz_set(Acoeffs + k + i, Acoeffs + i);
    for (i = 0; i < k; i++)
        fmpz_zero(Acoeffs + i);

    for (i = 0; i < Alen; i++)
        fmpz_mod_addmul(Acoeffs + i, Acoeffs + i, c, Acoeffs + i + k, ctx);

    A->length = Alen + k;
}


/* A = B + C*d0 */
void fmpz_mod_poly_scalar_addmul_fmpz_mod(
    fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t C,
    const fmpz_t d0,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = FLINT_MAX(B->length, C->length);

    fmpz_mod_poly_fit_length(A, Alen, ctx);
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Ccoeffs = C->coeffs;

    for (i = 0; i < FLINT_MIN(Blen, Clen); i++)
        fmpz_mod_addmul(Acoeffs + i, Bcoeffs + i, Ccoeffs + i, d0, ctx);

    for ( ; i < Clen; i++)
        fmpz_mod_mul(Acoeffs + i, Ccoeffs + i, d0, ctx);

    for ( ; i < Blen; i++)
        fmpz_set(Acoeffs + i, Bcoeffs + i);

    A->length = Alen;
    _fmpz_mod_poly_normalise(A);
}


/* A = B + C*(d1*x+d0) */
void fmpz_mod_poly_addmul_linear(
    fmpz_mod_poly_t A,
    const fmpz_mod_poly_t B,
    const fmpz_mod_poly_t C,
    const fmpz_t d1, const fmpz_t d0,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = FLINT_MAX(B->length, C->length + 1);

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_poly_fit_length(A, Alen, ctx);
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Ccoeffs = C->coeffs;

    for (i = 0; i < Alen; i++)
    {
        if (i < Blen)
            fmpz_set(Acoeffs + i, Bcoeffs + i);
        else
            fmpz_zero(Acoeffs + i);

        if (i < Clen)
            fmpz_addmul(Acoeffs + i, Ccoeffs + i, d0);

        if (0 < i && i - 1 < Clen)
            fmpz_addmul(Acoeffs + i, Ccoeffs + i - 1, d1);

        fmpz_mod_set_fmpz(Acoeffs + i, Acoeffs + i, ctx);
    }

    A->length = Alen;
    _fmpz_mod_poly_normalise(A);
}


/* evaluation at alphapow->coeffs[1] */
void fmpz_mod_poly_eval_pow(
    fmpz_t eval,
    const fmpz_mod_poly_t P,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx)
{
    slong Plen = P->length;
    if (Plen > alphapow->length)
    {
        slong i = alphapow->length;
        FLINT_ASSERT(2 <= i);
        fmpz_mod_poly_fit_length(alphapow, Plen, ctx);
        alphapow->length = Plen;
        for ( ; i < Plen; i++)
            fmpz_mod_mul(alphapow->coeffs + i, alphapow->coeffs + i - 1,
                                                    alphapow->coeffs + 1, ctx);
    }

    _fmpz_mod_vec_dot(eval, P->coeffs, alphapow->coeffs, Plen, ctx);
}

void fmpz_mod_poly_eval2_pow(
    fmpz_t vp,
    fmpz_t vm,
    const fmpz_mod_poly_t P, 
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx)
{
    const fmpz * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    fmpz * alpha_powers = alphapow->coeffs;
    fmpz_t a, b;
    slong k;

    fmpz_init(a);
    fmpz_init(b);

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        fmpz_mod_poly_fit_length(alphapow, Plen, ctx);
        for (k = oldlength; k < Plen; k++)
        {
            fmpz_mod_mul(alphapow->coeffs + k, alphapow->coeffs + k - 1,
                                               alphapow->coeffs + 1, ctx);
        }
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
    }

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        fmpz_addmul(a, Pcoeffs + k + 0, alpha_powers + k + 0);
        fmpz_addmul(b, Pcoeffs + k + 1, alpha_powers + k + 1);
    }

    if (k < Plen)
    {
        fmpz_addmul(a, Pcoeffs + k + 0, alpha_powers + k + 0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    fmpz_mod_set_fmpz(a, a, ctx);
    fmpz_mod_set_fmpz(b, b, ctx);

    fmpz_mod_add(vp, a, b, ctx);
    fmpz_mod_sub(vm, a, b, ctx);

    fmpz_clear(a);
    fmpz_clear(b);
}

