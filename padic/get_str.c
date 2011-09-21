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
#include "long_extras.h"

char * _padic_get_str(char * str, const padic_t op, const padic_ctx_t ctx)
{
    const fmpz * u = padic_unit(op);
    const long v   = padic_val(op);

    if (fmpz_is_zero(u))
    {
        if (!str)
        {
            str = malloc(2);
        }
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    if (ctx->mode == PADIC_TERSE)
    {
        if (v == 0)
        {
            str = fmpz_get_str(str, 10, u);
        }
        else if (v > 0)
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, ctx->p, v);
            fmpz_mul(t, t, u);
            str = fmpz_get_str(str, 10, t);
            fmpz_clear(t);
        }
        else  /* v < 0 */
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, ctx->p, -v);
            str = _fmpq_get_str(str, 10, u, t);
            fmpz_clear(t);
        }
    }
    else if (ctx->mode == PADIC_SERIES)
    {
        char *s;
        fmpz_t x;
        fmpz_t d;
        long j, N;

        if (fmpz_sgn(u) < 0)
        {
            printf("ERROR (_padic_get_str).  u < 0 in SERIES mode.\n");
            abort();
        }

        N = fmpz_clog(u, ctx->p) + v;

        if (!str)
        {
            long b = (N - v) * (2 * fmpz_sizeinbase(ctx->p, 10) 
                     + z_sizeinbase(FLINT_MAX(FLINT_ABS(v), FLINT_ABS(N)), 10) 
                     + 5) + 1;

            str = malloc(b);
            if (!str)
            {
                printf("ERROR (_padic_get_str).  Memory allocation failed.\n");
                abort();
            }
        }

        s = str;

        fmpz_init(d);
        fmpz_init(x);

        fmpz_set(x, u);

        /* Unroll first step */
        j = 0;
        {
            fmpz_mod(d, x, ctx->p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);            /* x = x - d */
            fmpz_divexact(x, x, ctx->p);  /* x = x / p */

            if (!fmpz_is_zero(d))
            {
                if (j + v != 0)
                {
                    fmpz_get_str(s, 10, d);
                    while (*++s != '\0') ;
                    *s++ = '*';
                    fmpz_get_str(s, 10, ctx->p);
                    while (*++s != '\0') ;
                    *s++ = '^';
                    sprintf(s, "%ld", j + v);
                    while (*++s != '\0') ;
                }
                else
                {
                    fmpz_get_str(s, 10, d);
                    while (*++s != '\0') ;
                }
            }

            j++;
        }

        for ( ; !fmpz_is_zero(x); j++)
        {
            fmpz_mod(d, x, ctx->p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);            /* x = x - d */
            fmpz_divexact(x, x, ctx->p);  /* x = x / p */

            if (!fmpz_is_zero(d))
            {
                if (j + v != 0)
                {
                    *s++ = ' ';
                    *s++ = '+';
                    *s++ = ' ';
                    fmpz_get_str(s, 10, d);
                    while (*++s != '\0') ;
                    *s++ = '*';
                    fmpz_get_str(s, 10, ctx->p);
                    while (*++s != '\0') ;
                    *s++ = '^';
                    sprintf(s, "%ld", j + v);
                    while (*++s != '\0') ;
                }
                else
                {
                    *s++ = ' ';
                    *s++ = '+';
                    *s++ = ' ';
                    fmpz_get_str(s, 10, d);
                    while (*++s != '\0') ;
                }
            }
        }
        
        fmpz_clear(x);
        fmpz_clear(d);
    }
    else  /* ctx->mode == PADIC_VAL_UNIT */
    {
        if (!str)
        {
            long b = fmpz_sizeinbase(u, 10) + fmpz_sizeinbase(ctx->p, 10) 
                   + z_sizeinbase(v, 10) + 4;

            str = malloc(b);
            if (!str)
            {
                printf("ERROR (_padic_get_str).  Memory allocation failed.\n");
                abort();
            }
        }

        if (v == 0)
        {
            str = fmpz_get_str(str, 10, u);
        }
        else if (v == 1)
        {
            char *s = str;

            fmpz_get_str(s, 10, u);
            while (*++s != '\0') ;
            *s++ = '*';
            fmpz_get_str(s, 10, ctx->p);
        }
        else
        {
            char *s = str;

            fmpz_get_str(s, 10, u);
            while (*++s != '\0') ;
            *s++ = '*';
            fmpz_get_str(s, 10, ctx->p);
            while (*++s != '\0') ;
            *s++ = '^';
            sprintf(s, "%ld", v);
        }
    }

    return str;
}

char * padic_get_str(char * str, const padic_t op, const padic_ctx_t ctx)
{
    padic_t t;

    _padic_init(t);
    padic_set(t, op, ctx);
    str = _padic_get_str(str, t, ctx);
    _padic_clear(t);
    return str;
}

