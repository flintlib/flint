/*============================================================================

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

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

void _padic_poly_add(fmpz *rop, long *val, 
                     const fmpz *op1, long val1, long len1, 
                     const fmpz *op2, long val2, long len2, 
                     const padic_ctx_t ctx)
{
    const long len = FLINT_MAX(len1, len2);
    fmpz_t pow;
    int alloc;

    *val = FLINT_MIN(val1, val2);

    alloc = _padic_ctx_pow_ui(pow, ctx->N - *val, ctx);

    if (val1 == val2)
    {
        _fmpz_mod_poly_add(rop, op1, len1, op2, len2, pow);
    }
    else  /* => (op1 != op2) */
    {
        long i;
        fmpz_t x;

        fmpz_init(x);
        if (val1 < val2)  /* F := p^g (G + p^{h-g} H) */
        {
            fmpz_pow_ui(x, ctx->p, val2 - val1);

            if (rop == op1)
            {
                _fmpz_vec_zero(rop + len1, len2 - len1);
                _fmpz_vec_scalar_addmul_fmpz(rop, op2, len2, x);
            }
            else
            {
                _fmpz_vec_scalar_mul_fmpz(rop, op2, len2, x);
                _fmpz_poly_add(rop, op1, len1, rop, len2);
            }
        }
        else  /* F := p^h (p^{g-h} G + H) */
        {
            fmpz_pow_ui(x, ctx->p, val1 - val2);

            if (rop == op2)
            {
                _fmpz_vec_zero(rop + len2, len1 - len2);
                _fmpz_vec_scalar_addmul_fmpz(rop, op1, len1, x);
            }
            else
            {
                _fmpz_vec_scalar_mul_fmpz(rop, op1, len1, x);
                _fmpz_poly_add(rop, rop, len1, op2, len2);
            }
        }
        fmpz_clear(x);

        for (i = 0; i < len; i++)
        {
            if (fmpz_cmpabs(rop + i, pow) >= 0)
                fmpz_sub(rop + i, rop + i, pow);
        }
    }

    if (alloc)
        fmpz_clear(pow);

    _padic_poly_canonicalise(rop, val, len, ctx->p);
}

void padic_poly_add(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)
{
    const long lenG = g->length;
    const long lenH = h->length;
    const long lenF = FLINT_MAX(lenG, lenH);

    if (lenG == 0 && lenH == 0)
    {
        padic_poly_zero(f);
        return;
    }

    padic_poly_fit_length(f, lenF);

    _padic_poly_add(f->coeffs, &(f->val), g->coeffs, g->val, lenG, 
                                          h->coeffs, h->val, lenH, ctx);

    _padic_poly_set_length(f, lenF);
    _padic_poly_normalise(f);
}

