/*
    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "fmpz_poly.h"


static void _fmpz_poly_product(
    fmpz_poly_t res,
    const fmpz_poly_struct * lifted_fac,
    const slong * subset,
    slong len,
    const fmpz_t P,
    const fmpz_t leadf,
    fmpz_poly_struct ** stack,
    fmpz_poly_struct * tmp)
{
    slong i, j, k;
    fmpz_poly_struct * t;

    k = 0;
    for (i = 0; i < len; i++)
    {
        if (subset[i] < 0)
            continue;

        stack[k] = (fmpz_poly_struct *) lifted_fac + subset[i];
        k++;

        for (j = k - 1; j > 0 && stack[j - 1]->length < stack[j]->length; j--)
        {
            t = stack[j - 1];
            stack[j - 1] = stack[j];
            stack[j] = t;
        }
    }

    while (k > 1)
    {
        for (j = 1; j < k; j++)
            FLINT_ASSERT(stack[j - 1]->length >= stack[j]->length);

        fmpz_poly_mul(res, stack[k - 2], stack[k - 1]);
        fmpz_poly_scalar_smod_fmpz(res, res, P);

        k--;
        stack[k - 1] = tmp + k - 1; /* make sure stack[k - 1] is writeable */
        fmpz_poly_swap(res, stack[k - 1]);

        for (j = k - 1; j > 0 && stack[j - 1]->length < stack[j]->length; j--)
        {
            t = stack[j - 1];
            stack[j - 1] = stack[j];
            stack[j] = t;
        }
    }

    if (k == 1)
    {
        fmpz_poly_scalar_mul_fmpz(res, stack[0], leadf);
        fmpz_poly_scalar_smod_fmpz(res, res, P);
    }
    else
    {
        FLINT_ASSERT(0);
        fmpz_poly_one(res);
    }
}


void fmpz_poly_factor_zassenhaus_recombination(
    fmpz_poly_factor_t final_fac,
	const fmpz_poly_factor_t lifted_fac,
    const fmpz_poly_t F,
    const fmpz_t P,
    slong exp)
{
    const slong r = lifted_fac->num;
    slong * subset;
    slong k, len;
    fmpz_poly_t Fcopy, Q, tryme;
    fmpz_poly_struct * tmp;
    fmpz_poly_struct ** stack;
    fmpz_poly_struct * f;

    subset = (slong *) flint_malloc(r*sizeof(slong));
    for (k = 0; k < r; k++)
        subset[k] = k;

    stack = (fmpz_poly_struct **) flint_malloc(r*sizeof(fmpz_poly_struct *));

    tmp = (fmpz_poly_struct *) flint_malloc(r*sizeof(fmpz_poly_struct));
    for (k = 0; k < r; k++)
        fmpz_poly_init(tmp + k);

    fmpz_poly_init(Q);
    fmpz_poly_init(tryme);
    fmpz_poly_init(Fcopy);

    f = (fmpz_poly_struct *) F;

    len = r;
    for (k = 1; k <= len/2; k++)
    {
        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            _fmpz_poly_product(tryme, lifted_fac->p, subset, len, P,
                                                fmpz_poly_lead(f), stack, tmp);
            fmpz_poly_primitive_part(tryme, tryme);
            if (fmpz_poly_divides(Q, f, tryme))
            {
                fmpz_poly_factor_insert(final_fac, tryme, exp);
                f = Fcopy;  /* make sure f is writeable */
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

    fmpz_poly_clear(Fcopy);
    fmpz_poly_clear(tryme);
    fmpz_poly_clear(Q);

    flint_free(stack);

    for (k = 0; k < r; k++)
        fmpz_poly_clear(tmp + k);
    flint_free(tmp);

    flint_free(subset);
}


void fmpz_poly_factor_zassenhaus_recombination_with_prune(
    fmpz_poly_factor_t final_fac,
    const fmpz_poly_factor_t lifted_fac,
    const fmpz_poly_t F,
    const fmpz_t P,
    slong exp,
    const zassenhaus_prune_t Z)
{
    const slong r = lifted_fac->num;
    slong * subset;
    slong i, k, len, total;
    fmpz_poly_t Fcopy, Q, tryme;
    fmpz_poly_struct * tmp;
    fmpz_poly_struct ** stack;
    fmpz_poly_struct * f;

    subset = (slong *) flint_malloc(r*sizeof(slong));
    for (k = 0; k < r; k++)
        subset[k] = k;

    stack = (fmpz_poly_struct **) flint_malloc(r*sizeof(fmpz_poly_struct *));

    tmp = (fmpz_poly_struct *) flint_malloc(r*sizeof(fmpz_poly_struct));
    for (k = 0; k < r; k++)
        fmpz_poly_init(tmp + k);

    fmpz_poly_init(Q);
    fmpz_poly_init(tryme);
    fmpz_poly_init(Fcopy);

    f = (fmpz_poly_struct *) F;

    len = r;
    for (k = 1; k <= len/2; k++)
    {
        zassenhaus_subset_first(subset, len, k);
        while (1)
        {
            total = 0;
            for (i = 0; i < len; i++)
                if (subset[i] >= 0)
                    total += fmpz_poly_degree(lifted_fac->p + subset[i]);

            if (!zassenhaus_prune_degree_is_possible(Z, total))
            {
                if (!zassenhaus_subset_next(subset, len))
                    break;
                continue;
            }

            _fmpz_poly_product(tryme, lifted_fac->p, subset, len, P,
                                                fmpz_poly_lead(f), stack, tmp);
            fmpz_poly_primitive_part(tryme, tryme);
            if (fmpz_poly_divides(Q, f, tryme))
            {
                fmpz_poly_factor_insert(final_fac, tryme, exp);
                f = Fcopy;  /* make sure f is writeable */
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

    fmpz_poly_clear(Fcopy);
    fmpz_poly_clear(tryme);
    fmpz_poly_clear(Q);

    flint_free(stack);

    for (k = 0; k < r; k++)
        fmpz_poly_clear(tmp + k);
    flint_free(tmp);

    flint_free(subset);
}

