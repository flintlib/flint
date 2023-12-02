/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "qadic.h"

/* printing *******************************************************************/

int _qadic_fprint_pretty(FILE * file, const fmpz * u, slong len, slong v, const qadic_ctx_t ctx)
{
    const fmpz *p = (&ctx->pctx)->p;

    if (_fmpz_vec_is_zero(u, len))
    {
        fputc('0', file);
        return 1;
    }

    if ((&ctx->pctx)->mode == PADIC_TERSE)
    {
        if (v == 0)
        {
            _fmpz_poly_fprint_pretty(file, u, len, ctx->var);
        }
        else if (v > 0)
        {
            fmpz *t = _fmpz_vec_init(len + 1);

            fmpz_pow_ui(t + len, p, v);
            _fmpz_vec_scalar_mul_fmpz(t, u, len, t + len);
            _fmpz_poly_fprint_pretty(file, t, len, ctx->var);
            _fmpz_vec_clear(t, len + 1);
        }
        else  /* v < 0 */
        {
            fmpz_t t;

            fmpz_init(t);
            fmpz_pow_ui(t, p, -v);
            _fmpq_poly_fprint_pretty(file, u, t, len, ctx->var);
            fmpz_clear(t);
        }
    }
    else if ((&ctx->pctx)->mode == PADIC_SERIES)
    {
        fmpz *x, *d;
        slong i, j;

        for (i = 0; i < len; i++)
            if (fmpz_sgn(u + i) < 0)
                break;

        if (i < len)
        {
            flint_throw(FLINT_ERROR, "ERROR (qadic_fprint_pretty).  u < 0 in SERIES mode.\n");
        }

        x = _fmpz_vec_init(len);
        d = _fmpz_vec_init(len);

        _fmpz_vec_set(x, u, len);

        /* Unroll first step */
        j = 0;
        {
            _fmpz_vec_scalar_mod_fmpz(d, x, len, p);       /* d = u mod p^{j+1} */
            _fmpz_vec_sub(x, x, d, len);                        /* x = x - d */
            _fmpz_vec_scalar_divexact_fmpz(x, x, len, p);  /* x = x / p */

            if (!_fmpz_vec_is_zero(d, len))
            {
                fputc('(', file);
                _fmpz_poly_fprint_pretty(file, d, len, ctx->var);
                fputc(')', file);
                if (j + v != 0)
                {
                    fputc('*', file);
                    fmpz_fprint(file, p);
                    if (j + v != 1)
                        flint_fprintf(file, "^%wd", j + v);
                }
            }

            j++;
        }

        for ( ; !_fmpz_vec_is_zero(x, len); j++)
        {
            _fmpz_vec_scalar_mod_fmpz(d, x, len, p);       /* d = u mod p^{j+1} */
            _fmpz_vec_sub(x, x, d, len);                        /* x = x - d */
            _fmpz_vec_scalar_divexact_fmpz(x, x, len, p);  /* x = x / p */

            if (!_fmpz_vec_is_zero(d, len))
            {
                flint_fprintf(file, " + ");
                fputc('(', file);
                _fmpz_poly_fprint_pretty(file, d, len, ctx->var);
                fputc(')', file);
                if (j + v != 0)
                {
                    fputc('*', file);
                    fmpz_fprint(file, p);
                    if (j + v != 1)
                        flint_fprintf(file, "^%wd", j + v);
                }
            }
        }
        _fmpz_vec_clear(x, len);
        _fmpz_vec_clear(d, len);
    }
    else if ((&ctx->pctx)->mode == PADIC_VAL_UNIT)
    {
        if (v == 0)
        {
            _fmpz_poly_fprint_pretty(file, u, len, ctx->var);
        }
        else if (v == 1)
        {
            fputc('(', file);
            _fmpz_poly_fprint_pretty(file, u, len, ctx->var);
            fputc(')', file);
            fputc('*', file);
            fmpz_fprint(file, p);
        }
        else
        {
            fputc('(', file);
            _fmpz_poly_fprint_pretty(file, u, len, ctx->var);
            fputc(')', file);
            fputc('*', file);
            fmpz_fprint(file, p);
            flint_fprintf(file, "^%wd", v);
        }
    }
    else
    {
        flint_throw(FLINT_ERROR, "Exception (qadic_fprint_pretty).  Unknown print mode.\n");
    }

    return 1;
}

int qadic_fprint_pretty(FILE *file, const qadic_t op, const qadic_ctx_t ctx) { return _qadic_fprint_pretty(file, op->coeffs, op->length, op->val, ctx); }
int qadic_print_pretty(const qadic_t op, const qadic_ctx_t ctx) { return qadic_fprint_pretty(stdout, op, ctx); }

/* debugging ******************************************************************/

int qadic_debug(const qadic_t op) { return padic_poly_debug(op); }
