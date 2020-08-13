/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_poly.h"
#include "padic_mat.h"

static void 
_padic_mat_canonicalise_fmpz(fmpz *vec, slong len, slong *val, const fmpz_t p)
{
    int nonzero = 0;
    slong i;

    for (i = 0; i < len; i++)
    {
        if (vec[i] != WORD(0))
        {
            nonzero = 1;
            if (!fmpz_divisible(vec + i, p))
            {
                return;
            }
        }
    }

    if (nonzero)
    {
        _fmpz_vec_scalar_divexact_fmpz(vec, vec, len, p);
        (*val)++;

        while (1)
        {
            for (i = 0; i < len; i++)
            {
                if (!fmpz_divisible(vec + i, p))
                {
                    return;
                }
            }

            _fmpz_vec_scalar_divexact_fmpz(vec, vec, len, p);
            (*val)++;
        }
    }
    else
    {
        *val = 0;
    }
}

static void 
_padic_mat_canonicalise_si(fmpz *vec, slong len, slong *val, slong p)
{
    int nonzero = 0;
    slong i;

    for (i = 0; i < len; i++)
    {
        if (vec[i] != WORD(0))
        {
            nonzero = 1;
            if (!fmpz_divisible_si(vec + i, p))
            {
                return;
            }
        }
    }

    if (nonzero)
    {
        _fmpz_vec_scalar_divexact_ui(vec, vec, len, p);
        (*val)++;

        while (1)
        {
            for (i = 0; i < len; i++)
            {
                if (!fmpz_divisible_si(vec + i, p))
                {
                    return;
                }
            }

            _fmpz_vec_scalar_divexact_ui(vec, vec, len, p);
            (*val)++;
        }
    }
    else
    {
        *val = 0;
    }
}

void _padic_mat_canonicalise(padic_mat_t A, const padic_ctx_t ctx)
{
    if (COEFF_IS_MPZ(*(ctx->p)))
    {
        _padic_mat_canonicalise_fmpz(padic_mat(A)->entries, 
                                     padic_mat(A)->r * padic_mat(A)->c, 
                                     &(A->val), ctx->p);
    }
    else
    {
        _padic_mat_canonicalise_si(padic_mat(A)->entries, 
                                   padic_mat(A)->r * padic_mat(A)->c, 
                                   &(A->val), *(ctx->p));
    }
}

