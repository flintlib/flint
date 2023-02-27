/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "fmpz_mat.h"

void
TEMPLATE(T, mat_mul_KS) (TEMPLATE(T, mat_t) C,
                         const TEMPLATE(T, mat_t) A,
                         const TEMPLATE(T, mat_t) B,
                         const TEMPLATE(T, ctx_t) ctx)
{
    slong bits;
    slong ar, bc, br;
    slong i, j;
    fmpz_mat_t fa, fb, fc;
    fmpz_t beta;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (br == 0)
    {
        TEMPLATE(T, mat_zero) (C, ctx);
        return;
    }

    /* Compute the number of bits needed */

    /* TODO: Make this generic based on say an TEMPLATE(T,
     * bits_needed) and TEMPLATE(T, bit_pack) */
    fmpz_init(beta);
    fmpz_set(beta, TEMPLATE(T, ctx_prime) (ctx));
    fmpz_sub_ui(beta, beta, 1);
    fmpz_mul(beta, beta, beta);
    fmpz_mul_si(beta, beta, A->c);
    fmpz_mul_si(beta, beta, TEMPLATE(T, ctx_degree) (ctx));
    bits = fmpz_bits(beta) + 1;

    fmpz_mat_init(fa, A->r, A->c);
    fmpz_mat_init(fb, B->r, B->c);
    fmpz_mat_init(fc, A->r, B->c);

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            TEMPLATE(T, bit_pack) (fmpz_mat_entry(fa, i, j),
                                   TEMPLATE(T, mat_entry) (A, i, j), bits,
                                   ctx);

    for (i = 0; i < B->r; i++)
        for (j = 0; j < B->c; j++)
            TEMPLATE(T, bit_pack) (fmpz_mat_entry(fb, i, j),
                                   TEMPLATE(T, mat_entry) (B, i, j), bits,
                                   ctx);

    fmpz_mat_mul(fc, fa, fb);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            TEMPLATE(T, bit_unpack) (TEMPLATE(T, mat_entry) (C, i, j),
                                     fmpz_mat_entry(fc, i, j), bits, ctx);
        }
    }

    fmpz_mat_clear(fa);
    fmpz_mat_clear(fb);
    fmpz_mat_clear(fc);

    fmpz_clear(beta);
}


#endif
