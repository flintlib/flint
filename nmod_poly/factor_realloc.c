/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "nmod_poly.h"

void nmod_poly_factor_realloc(nmod_poly_factor_t fac, len_t alloc)
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
            len_t i;

            for (i = alloc; i < fac->num; i++)
                nmod_poly_clear(fac->p + i);

            fac->p     = flint_realloc(fac->p, alloc * sizeof(nmod_poly_struct));
            fac->exp   = flint_realloc(fac->exp, alloc * sizeof(len_t));
            fac->alloc = alloc;
        }
        else if (fac->alloc < alloc)
        {
            len_t i;

            fac->p   = flint_realloc(fac->p, alloc * sizeof(nmod_poly_struct));
            fac->exp = flint_realloc(fac->exp, alloc * sizeof(len_t));

            for (i = fac->alloc; i < alloc; i++)
            {
                nmod_poly_init_preinv(fac->p + i, 0, 0);
                fac->exp[i] = 0L;
            }
            fac->alloc = alloc;
        }
    }
    else                        /* Nothing allocated already so do it now */
    {
        len_t i;

        fac->p   = flint_malloc(alloc * sizeof(nmod_poly_struct));
        fac->exp = flint_calloc(alloc, sizeof(len_t));

        for (i = 0; i < alloc; i++)
            nmod_poly_init_preinv(fac->p + i, 0, 0);
        fac->num   = 0;
        fac->alloc = alloc;
    }
}

