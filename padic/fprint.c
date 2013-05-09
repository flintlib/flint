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

    Copyright (C) 2011, 2012 Sebastian Pancratz
 
******************************************************************************/

#include <limits.h>

#include "padic.h"

int _padic_fprint(FILE * file, const fmpz_t u, len_t v, const padic_ctx_t ctx)
{
    const fmpz *p = ctx->p;

    if (fmpz_is_zero(u))
    {
        fputc('0', file);
        return 1;
    }

    if (ctx->mode == PADIC_TERSE)
    {
        if (v == 0)
        {
            fmpz_fprint(file, u);
        }
        else if (v > 0)
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, v);
            fmpz_mul(t, t, u);
            fmpz_fprint(file, t);
            fmpz_clear(t);
        }
        else  /* v < 0 */
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, -v);
            _fmpq_fprint(file, u, t);
            fmpz_clear(t);
        }
    }
    else if (ctx->mode == PADIC_SERIES)
    {
        fmpz_t x;
        fmpz_t d;
        len_t j;

        fmpz_init(d);
        fmpz_init(x);

        fmpz_set(x, u);

        /* Unroll first step */
        j = 0;
        {
            fmpz_mod(d, x, p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);       /* x = x - d */
            fmpz_divexact(x, x, p);  /* x = x / p */

            if (!fmpz_is_zero(d))
            {
                if (j + v != 0)
                {
                    fmpz_fprint(file, d);
                    fputc('*', file);
                    fmpz_fprint(file, p);
                    fprintf(file, "^%ld", j + v);
                }
                else
                {
                    fmpz_fprint(file, d);
                }
            }

            j++;
        }

        for ( ; !fmpz_is_zero(x); j++)
        {
            fmpz_mod(d, x, p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);       /* x = x - d */
            fmpz_divexact(x, x, p);  /* x = x / p */

            if (!fmpz_is_zero(d))
            {
                if (j + v != 0)
                {
                    fprintf(file, " + ");
                    fmpz_fprint(file, d);
                    fputc('*', file);
                    fmpz_fprint(file, p);
                    fprintf(file, "^%ld", j + v);
                }
                else
                {
                    fprintf(file, " + ");
                    fmpz_fprint(file, d);
                }
            }
        }
        
        fmpz_clear(x);
        fmpz_clear(d);
    }
    else if (ctx->mode == PADIC_VAL_UNIT)
    {
        if (v == 0)
        {
            fmpz_fprint(file, u);
        }
        else if (v == 1)
        {
            fmpz_fprint(file, u);
            fputc('*', file);
            fmpz_fprint(file, p);
        }
        else 
        {
            fmpz_fprint(file, u);
            fputc('*', file);
            fmpz_fprint(file, p);
            fprintf(file, "^%ld", v);
        }
    }
    else
    {
        printf("Exception (_padic_fprint).  Unknown print mode.\n");
        abort();
    }

    return 1;
}

int padic_fprint(FILE * file, const padic_t op, const padic_ctx_t ctx)
{
    return _padic_fprint(file, padic_unit(op), padic_val(op), ctx);
}

