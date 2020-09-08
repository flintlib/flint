/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

void
unity_zpq_mul(unity_zpq f, const unity_zpq g, const unity_zpq h)
{
    slong i, j, k;
    ulong p, q;
    fmpz_mod_poly_t temp;

    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(h->ctx)));

    q = f->q;
    p = f->p;
    fmpz_mod_poly_init(temp, f->ctx);

    for (i = 0; i < p; i++)
    {
        fmpz_mod_poly_zero(f->polys[i], f->ctx);
    }

    for (i = 0; i < p; i++)
    {
        for (j = 0; j < p; j++)
        {
            ulong qpow;

            qpow = n_addmod(i, j, p);
            fmpz_mod_poly_mul(temp, g->polys[i], h->polys[j], f->ctx);

            if (temp->length == 0)
                continue;

            for (k = temp->length - 1; k >= q; k--)
            {
                fmpz_add(temp->coeffs + k - q,
                        temp->coeffs + k - q, temp->coeffs + k);
                fmpz_set_ui(temp->coeffs + k, 0);
                fmpz_mod(temp->coeffs + k - q, temp->coeffs + k - q,
                                                 fmpz_mod_ctx_modulus(f->ctx));
            }
            _fmpz_mod_poly_normalise(temp);

            fmpz_mod_poly_add(f->polys[qpow], f->polys[qpow], temp, f->ctx);
        }
    }

    fmpz_mod_poly_clear(temp, f->ctx);
}

