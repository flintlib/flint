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

int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, long d, const char *var)
{
    char *buf;
    FILE *file;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        return 0;
    }

    buf  = flint_malloc(832);
    file = fopen(FLINT_CPIMPORT, "r");

    if (!file)
    {
        file = fopen("../qadic/CPimport.txt", "r");

        if (!file)
        {
            printf("Exception (fq_nmod_ctx_init_conway).  File loading.\n");
            abort();
        }
    }

    while (fgets(buf, 832, file))
    {
        char *tmp = buf;

        /* Different prime? */
        if (fmpz_cmp_ui(p, atoi(tmp)))
            continue;

        while (*tmp++ != ' ') ;

        /* Same degree? */
        if (d == atoi(tmp))
        {
            nmod_poly_t mod;
            long i;
            char *ptr;

            nmod_poly_init(mod, fmpz_get_ui(p));
            
            /* Copy the polynomial */
            ptr = tmp;

            for (i = 0; i < d; i++)
            {
                int coeff;

                while (*ptr++ != ' ') ;

                coeff = atoi(ptr);
                
                nmod_poly_set_coeff_ui(mod, i, coeff);
            }

            nmod_poly_set_coeff_ui(mod, d, 1);

            fq_nmod_ctx_init_modulus(ctx, p, d, mod, var);

            fclose(file);
            flint_free(buf);
            return 1;
        }
    }

    fclose(file);
    flint_free(buf);

    return 0;
}


void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, const fmpz_t p, long d, const char *var)
{
    int result;
    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        printf("Exception (fq_nmod_ctx_init_conway).  Conway polynomials \n");
        printf("are only available for primes up to 109987.\n");
        abort();
    }

    result = _fq_nmod_ctx_init_conway(ctx, p, d, var);
    if (!result) {
        printf("Exception (fq_nmod_ctx_init_conway).  The polynomial for \n(p,d) = (");
        fmpz_print(p), printf(",%ld) is not present in the database.\n", d);
        abort();
    }
}
