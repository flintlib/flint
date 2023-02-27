/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_poly.h"

void nmod_poly_factor_realloc(nmod_poly_factor_t fac, slong alloc)
{
    if (alloc == 0)             /* Clear up, reinitialise */
    {
        nmod_poly_factor_clear(fac);
        nmod_poly_factor_init(fac);
    }
    else if (fac->alloc)            /* Realloc */
    {
        if (fac->alloc > alloc)
        {
            slong i;

            for (i = alloc; i < fac->num; i++)
                nmod_poly_clear(fac->p + i);

            fac->p     = flint_realloc(fac->p, alloc * sizeof(nmod_poly_struct));
            fac->exp   = flint_realloc(fac->exp, alloc * sizeof(slong));
            fac->alloc = alloc;
        }
        else if (fac->alloc < alloc)
        {
            slong i;

            fac->p   = flint_realloc(fac->p, alloc * sizeof(nmod_poly_struct));
            fac->exp = flint_realloc(fac->exp, alloc * sizeof(slong));

            for (i = fac->alloc; i < alloc; i++)
            {
                nmod_poly_init_preinv(fac->p + i, 1, 0);
                fac->exp[i] = WORD(0);
            }
            fac->alloc = alloc;
        }
    }
    else                        /* Nothing allocated already so do it now */
    {
        slong i;

        fac->p   = flint_malloc(alloc * sizeof(nmod_poly_struct));
        fac->exp = flint_calloc(alloc, sizeof(slong));

        for (i = 0; i < alloc; i++)
            nmod_poly_init_preinv(fac->p + i, 1, 0);
        fac->num   = 0;
        fac->alloc = alloc;
    }
}

