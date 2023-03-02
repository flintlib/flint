/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"


void zassenhaus_prune_clear(zassenhaus_prune_t Z)
{
    if (Z->alloc > 0)
    {
        flint_free(Z->pos_degs);
        flint_free(Z->new_degs);
    }
}

void zassenhaus_prune_set_degree(zassenhaus_prune_t Z, slong d)
{
    slong i;

    if (d < 1)
    {
        flint_throw(FLINT_ERROR, "zassenhaus_prune_set_degree");
        return;
    }

    if (Z->alloc > 0)
    {
        Z->pos_degs = (unsigned char *) flint_realloc(Z->pos_degs,
                                                (d + 1)*sizeof(unsigned char));
        Z->new_degs = (slong *) flint_realloc(Z->new_degs,
                                                        (d + 1)*sizeof(slong));
    }
    else
    {
        Z->pos_degs = (unsigned char *) flint_malloc(
                                                (d + 1)*sizeof(unsigned char));
        Z->new_degs = (slong *) flint_malloc((d + 1)*sizeof(slong));
    }
    Z->alloc = d + 1;
    Z->deg = d;

    for (i = 0; i <= d; i++)
        Z->pos_degs[i] = 1;

    Z->new_length = 0;
    Z->new_total = 0;
}

void zassenhaus_prune_add_factor(zassenhaus_prune_t Z, slong deg, slong exp)
{
    slong i;

    if (exp < 1 || deg < 1)
        return;

    for (i = 0; i < exp; i++)
    {
        if (Z->new_length >= Z->deg)
        {
            flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
            return;
        }
        Z->new_total += deg;
        Z->new_degs[Z->new_length] = deg;
        Z->new_length++;
    }
}

void zassenhaus_prune_end_add_factors(zassenhaus_prune_t Z)
{
    slong i, j;
    unsigned char * a = Z->pos_degs;
    unsigned char pos_mask = 1;
    unsigned char new_mask = 2;

    if (Z->new_total != Z->deg)
    {
        flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
    }

    a[0] |= new_mask;
    for (j = 1; j <= Z->deg; j++)
        a[j] &= ~new_mask;

    for (i = 0; i < Z->new_length; i++)
    {
        slong d = Z->new_degs[i];

        for (j = Z->deg; j >= 0; j--)
        {
            if ((a[j] & new_mask) != 0)
            {
                if (j + d > Z->deg)
                {
                    flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
                }
                a[j + d] |= new_mask;
            }
        }
    }

    /* merge new possibilities with old */
    for (j = 0; j <= Z->deg; j++)
        a[j] &= a[j] >> 1;

    /* 0 and deg should always be possible */
    if (a[0] != pos_mask || a[Z->deg] != pos_mask)
    {
        flint_throw(FLINT_ERROR, "zassenhaus_prune_add_factor");
    }
}

int zassenhaus_prune_must_be_irreducible(const zassenhaus_prune_t Z)
{
    slong i;

    for (i = 1; i < Z->deg; i++)
    {
        if (Z->pos_degs[i] != 0)
            return 0;
    }

    return 1;
}

