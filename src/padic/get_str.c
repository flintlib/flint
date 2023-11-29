/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "fmpq.h"
#include "padic.h"

char * _padic_get_str(char *str, const padic_t op, const fmpz_t p, enum padic_print_mode mode)
{
    const fmpz *u = padic_unit(op);
    const slong v = padic_val(op);

    if (fmpz_is_zero(u))
    {
        if (!str)
        {
            str = flint_malloc(2);
        }
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    if (mode == PADIC_TERSE)
    {
        if (v == 0)
        {
            str = fmpz_get_str(str, 10, u);
        }
        else if (v > 0)
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, v);
            fmpz_mul(t, t, u);
            str = fmpz_get_str(str, 10, t);
            fmpz_clear(t);
        }
        else  /* v < 0 */
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, -v);
            str = _fmpq_get_str(str, 10, u, t);
            fmpz_clear(t);
        }
    }
    else if (mode == PADIC_SERIES)
    {
        char *s;
        fmpz_t x;
        fmpz_t d;
        slong j, N;

        N = fmpz_clog(u, p) + v + 1;

        if (!str)
        {
            slong b = (N - v) * (2 * fmpz_sizeinbase(p, 10)
                     + z_sizeinbase(FLINT_MAX(FLINT_ABS(v), FLINT_ABS(N)), 10)
                     + 5) + 1;

            str = flint_malloc(b);
            if (!str)
            {
                flint_throw(FLINT_ERROR, "Exception (padic_get_str).  Memory allocation failed.\n");
            }
        }

        s = str;

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
                    fmpz_get_str(s, 10, d);
                    while (*++s != '\0') ;
                    *s++ = '*';
                    fmpz_get_str(s, 10, p);
                    while (*++s != '\0') ;
                    *s++ = '^';
                    flint_sprintf(s, "%wd", j + v);
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
            fmpz_mod(d, x, p);       /* d = u mod p^{j+1} */
            fmpz_sub(x, x, d);       /* x = x - d */
            fmpz_divexact(x, x, p);  /* x = x / p */

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
                    fmpz_get_str(s, 10, p);
                    while (*++s != '\0') ;
                    *s++ = '^';
                    flint_sprintf(s, "%wd", j + v);
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
            slong b = fmpz_sizeinbase(u, 10) + fmpz_sizeinbase(p, 10)
                   + z_sizeinbase(v, 10) + 4;

            str = flint_malloc(b);
            if (!str)
            {
                flint_throw(FLINT_ERROR, "Exception (padic_get_str).  Memory allocation failed.\n");
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
            fmpz_get_str(s, 10, p);
        }
        else
        {
            char *s = str;

            fmpz_get_str(s, 10, u);
            while (*++s != '\0') ;
            *s++ = '*';
            fmpz_get_str(s, 10, p);
            while (*++s != '\0') ;
            *s++ = '^';
            flint_sprintf(s, "%wd", v);
        }
    }

    return str;
}

char * padic_get_str(char *str, const padic_t op, const padic_ctx_t ctx)
{
    return _padic_get_str(str, op, ctx->p, ctx->mode);
}

