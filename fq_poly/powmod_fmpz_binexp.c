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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fq_poly.h"

void
_fq_poly_powmod_fmpz_binexp(fq_struct * res, const fq_struct * poly,
                            const fmpz_t e, const fq_struct * f, slong lenf,
                            const fq_ctx_t ctx)
{
    fq_struct *T, *Q;
    fq_t invf;
    slong lenT, lenQ;
    slong i;

    if (lenf == 2)
    {
        fq_pow(res, poly, e, ctx);
        return;
    }

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _fq_vec_init(lenT + lenQ, ctx);
    Q = T + lenT;

    fq_init(invf, ctx);
    fq_inv(invf, f + lenf - 1, ctx);

    _fq_vec_set(res, poly, lenf - 1, ctx);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        _fq_poly_sqr(T, res, lenf - 1, ctx);
        _fq_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, invf, ctx);

        if (fmpz_tstbit(e, i))
        {
            _fq_poly_mul(T, res, lenf - 1, poly, lenf - 1, ctx);
            _fq_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, invf, ctx);
        }
    }

    fq_clear(invf, ctx);
    _fq_vec_clear(T, lenT + lenQ, ctx);
}


void
fq_poly_powmod_fmpz_binexp(fq_poly_t res, const fq_poly_t poly, const fmpz_t e,
                           const fq_poly_t f, const fq_ctx_t ctx)
{
    fq_struct *q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;

    if (lenf == 0)
    {
        flint_printf("Exception: fq_poly_powmod: divide by zero\n");
        abort();
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_printf
            ("Exception: fq_poly_powmod: negative exp not implemented\n");
        abort();
    }

    if (len >= lenf)
    {
        fq_poly_t t, r;
        fq_poly_init(t, ctx);
        fq_poly_init(r, ctx);
        fq_poly_divrem(t, r, poly, f, ctx);
        fq_poly_powmod_fmpz_binexp(res, r, e, f, ctx);
        fq_poly_clear(t, ctx);
        fq_poly_clear(r, ctx);
        return;
    }

    if (fmpz_abs_fits_ui(e))
    {
        ulong exp = fmpz_get_ui(e);

        if (exp <= 2)
        {
            if (exp == 0UL)
            {
                fq_poly_fit_length(res, 1, ctx);
                fq_one(res->coeffs, ctx);
                _fq_poly_set_length(res, 1, ctx);
            }
            else if (exp == 1UL)
            {
                fq_poly_set(res, poly, ctx);
            }
            else
                fq_poly_mulmod(res, poly, poly, f, ctx);
            return;
        }
    }

    if (lenf == 1 || len == 0)
    {
        fq_poly_zero(res, ctx);
        return;
    }

    if (poly->length < trunc)
    {
        q = _fq_vec_init(trunc, ctx);
        _fq_vec_set(q, poly->coeffs, len, ctx);
        _fq_vec_zero(q + len, trunc - len, ctx);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    if ((res == poly && !qcopy) || (res == f))
    {
        fq_poly_t t;
        fq_poly_init2(t, 2 * lenf - 3, ctx);
        _fq_poly_powmod_fmpz_binexp(t->coeffs, q, e, f->coeffs, lenf, ctx);
        fq_poly_swap(res, t, ctx);
        fq_poly_clear(t, ctx);
    }
    else
    {
        fq_poly_fit_length(res, 2 * lenf - 3, ctx);
        _fq_poly_powmod_fmpz_binexp(res->coeffs, q, e, f->coeffs, lenf, ctx);
    }

    if (qcopy)
        _fq_vec_clear(q, trunc, ctx);

    _fq_poly_set_length(res, trunc, ctx);
    _fq_poly_normalise(res, ctx);
}
