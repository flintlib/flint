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

#include <limits.h>

#include "padic.h"

int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        fprintf(file, "0 + O(");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);

        return 1;
    }
    
    if (ctx->mode == PADIC_TERSE)
    {
        if (padic_val(op) >= 0)
        {
            fmpz_t pow, x;

            fmpz_init(pow);
            fmpz_init(x);

            fmpz_pow_ui(pow, ctx->p, padic_val(op));
            fmpz_mul(x, pow, padic_unit(op));

            fmpz_fprint(file, x);

            fmpz_clear(pow);
            fmpz_clear(x);
        }
        else
        {
            fmpz_t pow;
            mpq_t y;

            fmpz_init(pow);
            mpq_init(y);

            fmpz_pow_ui(pow, ctx->p, - (ulong) padic_val(op));
            fmpz_get_mpz(mpq_numref(y), padic_unit(op));
            fmpz_get_mpz(mpq_denref(y), pow);

            gmp_fprintf(file, "%Qd", y);

            fmpz_clear(pow);
            mpq_clear(y);
        }

        fprintf(file, " + O(");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);
    }
    else if (ctx->mode == PADIC_SERIES)
    {
        fmpz_t x;
        fmpz_t d;
        long j;

        fmpz_init(d);
        fmpz_init(x);

        fmpz_set(x, padic_unit(op));
        
        for (j = 0; j < ctx->N - padic_val(op); j++)
        {
            fmpz_mod(d, x, ctx->p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);            /* x = x - d */
            fmpz_divexact(x, x, ctx->p);  /* x = x / p */

            if (!fmpz_is_zero(d))
            {
                if(j + padic_val(op) != 0)
                {
                    fmpz_fprint(file, d);
                    fprintf(file, "*");
                    fmpz_fprint(file, ctx->p);
                    fprintf(file, "^%ld + ", j + padic_val(op));
                }
                else
                {
                    fmpz_fprint(file, d);
                    fprintf(file, " + ");
                }
            }
        }

        fprintf(file, "O(");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);

        fmpz_clear(x);
        fmpz_clear(d);
    }
    else if (ctx->mode == PADIC_VAL_UNIT)
    {
        fmpz_t x, pow;
        int alloc;

        fmpz_init(x);

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(op), ctx);
        fmpz_mod(x, padic_unit(op), pow);
        if (alloc)
            fmpz_clear(pow);

        fmpz_fprint(file, x);
        fprintf(file, "*");
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld + O(", padic_val(op));
        fmpz_fprint(file, ctx->p);
        fprintf(file, "^%ld)", ctx->N);
    }
    else
    {
        printf("Exception (padic_fprint).  Unknown printing mode.\n");
        abort();
    }

    return 1;
}
