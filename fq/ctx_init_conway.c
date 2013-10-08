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

    Copyright (C) 2012 Andres Goens, 2012 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "fq.h"

#ifndef FLINT_CPIMPORT
#define FLINT_CPIMPORT "../qadic/CPimport.txt"
#endif

void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, long d, const char *var)
{
    char *buf;
    FILE *file;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        printf("Exception (fq_ctx_init_conway).  Conway polynomials \n");
        printf("are only available for primes up to 109987.\n");
        abort();
    }

    buf  = flint_malloc(832);
    file = fopen(FLINT_CPIMPORT, "r");

    if (!file)
    {
        file = fopen("../qadic/CPimport.txt", "r");

        if (!file)
        {
            printf("Exception (fq_ctx_init_conway).  File loading.\n");
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
            fmpz_mod_poly_t mod;
            long i;
            char *ptr;

            fmpz_mod_poly_init(mod, p);

            /* Copy the polynomial */
            ptr = tmp;

            for (i = 0; i < d; i++)
            {
                int coeff;

                while (*ptr++ != ' ') ;

                coeff = atoi(ptr);
                fmpz_mod_poly_set_coeff_ui(mod, i, coeff);
            }
            fmpz_mod_poly_set_coeff_ui(mod, d, 1);

            fq_ctx_init_modulus(ctx, p, d, mod, var);

            fclose(file);
            flint_free(buf);
            return;
        }
    }

    fclose(file);
    flint_free(buf);

    printf("Exception (fq_ctx_init_conway).  The polynomial for \n(p,d) = (");
    fmpz_print(p), printf(",%ld) is not present in the database.\n", d);
    abort();
}

