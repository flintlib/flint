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

#include <stdio.h>
#include <string.h>

#include "fmpz_vec.h"
#include "padic.h"
#include "qadic.h"

#ifndef FLINT_CPIMPORT
#define FLINT_CPIMPORT "/home/user/FLINT/flint-2/fq/CP.txt"
#endif

void qadic_ctx_init_conway(qadic_ctx_t ctx,
                           const fmpz_t p, long d, long N, const char *var, 
                           enum padic_print_mode mode)
{
    char *buf;
    FILE *file;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        printf("Exception (qadic_ctx_init_conway).  Conway polynomials \n");
        printf("are only available for primes up to 109987.\n");
        abort();
    }

    buf  = flint_malloc(832);
    file = fopen(FLINT_CPIMPORT, "r");

    if (!file)
    {
        printf("Exception (qadic_ctx_init_conway).  File loading.\n");
        abort();
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
            long i, j;
            char *ptr;

            /* Find number of non-zero coefficients */
            ctx->len = 1;
            ptr = tmp;

            for (i = 0; i < d; i++)
            {
                while (*ptr++ != ' ') ;

                if (atoi(ptr))
                    ctx->len ++;
            }

            ctx->a = _fmpz_vec_init(ctx->len);
            ctx->j = flint_malloc(ctx->len * sizeof(long));

            /* Copy the polynomial */
            j = 0;
            ptr = tmp;

            for (i = 0; i < d; i++)
            {
                int coeff;

                while (*ptr++ != ' ') ;

                coeff = atoi(ptr);

                if (coeff)
                {
                    fmpz_set_ui(ctx->a + j, coeff);
                    ctx->j[j] = i;
                    j++;
                }
            }

            fmpz_set_ui(ctx->a + j, 1);
            ctx->j[j] = d;

            /* Complete the initialisation of the context */
            padic_ctx_init(&ctx->pctx, p, N, mode);

            ctx->var = flint_malloc(strlen(var) + 1);
            strcpy(ctx->var, var);

            fclose(file);
            flint_free(buf);
            return;
        }
    }

    fclose(file);
    flint_free(buf);

    printf("Exception (qadic_ctx_init_conway).  The polynomial for \n");
    printf("(p,d) = (%ld,%ld) is not present in the database.\n", *p, d);
    abort();
}

