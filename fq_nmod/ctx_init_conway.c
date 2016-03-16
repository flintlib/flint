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

    Copyright (C) 2012 Andres Goens
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "fq_nmod.h"

/* from qadic/ctx_init_conway.c */
extern int flint_conway_polynomials [];

int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    unsigned int position;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        return 0;
    }

    for (position = 0; flint_conway_polynomials[position] != 0; position += 3+flint_conway_polynomials[position+1])
    {
        /* Different prime? */
        if (fmpz_cmp_ui(p, flint_conway_polynomials[position]))
            continue;

        /* Same degree? */
        if (d == flint_conway_polynomials[position+1])
        {
            nmod_poly_t mod;
            slong i;

            nmod_poly_init(mod, fmpz_get_ui(p));
            
            /* Copy the polynomial */
            for (i = 0; i < d; i++)
            {
                int coeff = flint_conway_polynomials[position+2+i];                
                nmod_poly_set_coeff_ui(mod, i, coeff);
            }

            nmod_poly_set_coeff_ui(mod, d, 1);

            fq_nmod_ctx_init_modulus(ctx, mod, var);

            nmod_poly_clear(mod);
            return 1;
        }
    }

    return 0;
}


void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    int result;
    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        flint_printf("Exception (fq_nmod_ctx_init_conway).  Conway polynomials \n");
        flint_printf("are only available for primes up to 109987.\n");
        flint_abort();
    }

    result = _fq_nmod_ctx_init_conway(ctx, p, d, var);
    if (!result) {
        flint_printf("Exception (fq_nmod_ctx_init_conway).  The polynomial for \n(p,d) = (");
        fmpz_print(p), flint_printf(",%wd) is not present in the database.\n", d);
        flint_abort();
    }
}
