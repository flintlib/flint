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

#include "fq.h"

#ifndef FLINT_CPIMPORT
#define FLINT_CPIMPORT "../qadic/CPimport.txt"
#endif

int
_fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    char *buf;
    FILE *file;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        return 0;
    }

    buf = flint_malloc(832);
    file = fopen(FLINT_CPIMPORT, "r");

    if (!file)
    {
        file = fopen("../qadic/CPimport.txt", "r");

        if (!file)
        {
            flint_printf("Exception (fq_ctx_init_conway).  File loading.\n");
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
            slong i;
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

            fq_ctx_init_modulus(ctx, mod, var);

            fmpz_mod_poly_clear(mod);
            fclose(file);
            flint_free(buf);
            return 1;
        }
    }
    fclose(file);
    flint_free(buf);
    return 0;
}

void
fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char *var)
{
    int result;
    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        flint_printf("Exception (fq_ctx_init_conway).  Conway polynomials \n");
        flint_printf("are only available for primes up to 109987.\n");
        abort();
    }

    result = _fq_ctx_init_conway(ctx, p, d, var);
    if (!result)
    {
        flint_printf
            ("Exception (fq_ctx_init_conway).  The polynomial for \n(p,d) = (");
        fmpz_print(p);
        flint_printf(",%wd) is not present in the database.\n", d);
        abort();
    }
}
