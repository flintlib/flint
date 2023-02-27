/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void 
_fmpz_mod_poly_div_series(fmpz * Q, const fmpz * A, slong Alen,
    const fmpz * B, slong Blen, const fmpz_t p, slong n)
{
    fmpz_t u, d;

    fmpz_init(d);
    fmpz_init(u);
     
    if (!fmpz_is_one(B + 0))
    {
       fmpz_gcdinv(d, u, B + 0, p);

       if (!fmpz_is_one(d)) /* check for invertibility */
       {
           printf("Exception (fmpz_mod_poly_div_series). Impossible inverse.");
               
           fmpz_clear(u);
           fmpz_clear(d);
               
           flint_abort();
       }
    } else
       fmpz_set_ui(u, 1);
      
    Alen = FLINT_MIN(Alen, n);
    Blen = FLINT_MIN(Blen, n);

    if (Blen == 1)
    {
        if (fmpz_is_one(B + 0))
            _fmpz_vec_set(Q, A, Alen);
        else
           _fmpz_mod_poly_scalar_mul_fmpz(Q, A, Alen, u, p);

        _fmpz_vec_zero(Q + Alen, n - Alen);
    }
    else if (n < 32 || Blen < 20)
    {
        slong i, j;

        if (fmpz_is_one(B + 0))
            fmpz_set(Q + 0, A + 0);
        else
        {
           fmpz_mul(Q + 0, u, A + 0);
           fmpz_mod(Q + 0, Q + 0, p);
        }

        for (i = 1; i < n; i++)
        {
            fmpz_mul(Q + i, B + 1, Q + i - 1);

            for (j = 2; j < FLINT_MIN(i + 1, Blen); j++)
                fmpz_addmul(Q + i, B + j, Q + i - j);

            if (i < Alen)
               fmpz_sub(Q + i, A + i, Q + i);
            else
               fmpz_neg(Q + i, Q + i);

            if (!fmpz_is_one(B + 0))
               fmpz_mul(Q + i, Q + i, u);
                           
            fmpz_mod(Q + i, Q + i, p);
        }
    }
    else
    {
        fmpz * B2, * Binv = _fmpz_vec_init(n);
        
        if (n > Blen)
        {
           B2 = _fmpz_vec_init(n);
           _fmpz_vec_set(B2, B, Blen);
        } else
           B2 = (fmpz *) B;

        _fmpz_mod_poly_inv_series(Binv, B2, n, u, p);
        _fmpz_mod_poly_mullow(Q, Binv, n, A, Alen, p, n);

        _fmpz_vec_clear(Binv, n);
        if (n > Blen)
           _fmpz_vec_clear(B2, n);
    }

    fmpz_clear(d);
    fmpz_clear(u);
}

void fmpz_mod_poly_div_series(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, 
                    const fmpz_mod_poly_t B, slong n, const fmpz_mod_ctx_t ctx)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_div_series). Division by zero.\n");
        flint_abort();
    }

    if (Alen == 0)
    {
        fmpz_mod_poly_zero(Q, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, n, ctx);
        _fmpz_mod_poly_div_series(t->coeffs, A->coeffs, Alen,
                                B->coeffs, Blen, fmpz_mod_ctx_modulus(ctx), n);
        fmpz_mod_poly_swap(Q, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, n, ctx);
        _fmpz_mod_poly_div_series(Q->coeffs, A->coeffs, Alen,
                                B->coeffs, Blen, fmpz_mod_ctx_modulus(ctx), n);
    }

    _fmpz_mod_poly_set_length(Q, n);
    _fmpz_mod_poly_normalise(Q);
}
