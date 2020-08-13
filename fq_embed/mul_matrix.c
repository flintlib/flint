/*
    Copyright (C) 2017 Luca De Feo

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_embed.h"

void fq_embed_mul_matrix(fmpz_mod_mat_t matrix,
                   const fq_t gen,
                   const fq_ctx_t ctx)
{
    slong i, j, len = fq_ctx_degree(ctx);
    fmpz_t lead;

    /* This is usually 1, unless the context is non-monic */
    fmpz_init(lead);
    fmpz_invmod(lead, ctx->modulus->coeffs + len, fq_ctx_prime(ctx));

    for (i = 0; i < gen->length; i++)
        fmpz_set(fmpz_mod_mat_entry(matrix, i, 0), gen->coeffs + i);
    for (i = gen->length; i < len; i++)
        fmpz_zero(fmpz_mod_mat_entry(matrix, i, 0));

    /* M[i, j] = M[i - 1, j - 1] - M[len - 1, j - 1] * lead * ctx->modulus->coeffs[i] */
    for (j = 1; j < len; j++)
    {
        fmpz_mul(fmpz_mod_mat_entry(matrix, len - 1, j),
                 fmpz_mod_mat_entry(matrix, len - 1, j - 1),
                 lead);
        for (i = 0; i < len; i++)
        {
            fmpz_mul(fmpz_mod_mat_entry(matrix, i, j),
                     fmpz_mod_mat_entry(matrix, len - 1, j),
                     ctx->modulus->coeffs + i);
            if (i > 0)
                fmpz_sub(fmpz_mod_mat_entry(matrix, i, j),
                         fmpz_mod_mat_entry(matrix, i, j),
                         fmpz_mod_mat_entry(matrix, i - 1, j - 1));
            fmpz_neg(fmpz_mod_mat_entry(matrix, i, j),
                     fmpz_mod_mat_entry(matrix, i, j));
        }
    }

    _fmpz_mod_mat_reduce(matrix);

    fmpz_clear(lead);
}
