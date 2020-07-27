/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"


static void _fmpz_poly_product(
    fmpz_poly_t res,
    const fmpz_poly_struct * lifted_fac,
    const slong * subset,
    slong len,
    const fmpz_t P)
{
    slong i;
    fmpz_poly_one(res);
    for (i = 0; i < len; i++)
    {
        if (subset[i] >= 0)
        {
            fmpz_poly_mul(res, res, lifted_fac + subset[i]);
            fmpz_poly_scalar_smod_fmpz(res, res, P);
        }
    }
}


void fmpz_poly_factor_zassenhaus_recombination(
    fmpz_poly_factor_t final_fac,
	const fmpz_poly_factor_t lifted_fac, 
    const fmpz_poly_t F,
    const fmpz_t P,
    slong exp)
{
    slong * subset;
    slong k, len;
    fmpz_poly_t f, Q, tryme;

    len = lifted_fac->num;
    subset = (slong *) flint_malloc(len*sizeof(slong));
    for (k = 0; k < len; k++)
        subset[k] = k;

    fmpz_poly_init(f);
    fmpz_poly_init(Q);
    fmpz_poly_init(tryme);

    fmpz_poly_set(f, F);

    for (k = 1; k <= len/2; k++)
    {
        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            _fmpz_poly_product(tryme, lifted_fac->p, subset, len, P);
            fmpz_poly_scalar_mul_fmpz(tryme, tryme, fmpz_poly_lead(f));
            fmpz_poly_scalar_smod_fmpz(tryme, tryme, P);
            fmpz_poly_primitive_part(tryme, tryme);
            if (fmpz_poly_divides(Q, f, tryme))
            {
                fmpz_poly_factor_insert(final_fac, tryme, exp);
                fmpz_poly_swap(f, Q);
                len -= k;
                if (!zassenhaus_subset_next_disjoint(subset, len + k))
                    break;
            }
            else
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
            }
        }
    }

    if (fmpz_poly_degree(f) > 0)
    {
        fmpz_poly_factor_insert(final_fac, f, exp);
    }
    else
    {
        FLINT_ASSERT(fmpz_poly_is_one(f));
    }

    fmpz_poly_clear(f);
    fmpz_poly_clear(tryme);
    fmpz_poly_clear(Q);
    flint_free(subset);
}

