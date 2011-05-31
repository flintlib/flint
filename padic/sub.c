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

#include "padic.h"

/*
    Computes the difference $u_1 p^{v_1} - u_2 p^{v_2}$ in canonical form.

    Supports aliasing.
 */
void _padic_sub(padic_t rop, const padic_t op1, const padic_t op2, 
                const padic_ctx_t ctx)
{
    if (_padic_is_zero(op1))
    {
        _padic_neg(rop, op2);
    }
    if (_padic_is_zero(op2))
    {
        _padic_set(rop, op1);
    }
    else
    {
        fmpz_t pow;
        int alloc = 0;

        if (padic_val(op1) < padic_val(op2))  /* u1 - p^{v2-v1} u2 */
        {
            _padic_ctx_pow_ui(pow, &alloc, padic_val(op2) - padic_val(op1), ctx);

            if (rop != op2)
            {
                fmpz_set(padic_unit(rop), padic_unit(op1));
                fmpz_submul(padic_unit(rop), pow, padic_unit(op2));
            }
            else
            {
                fmpz_mul(padic_unit(rop), pow, padic_unit(op2));
                fmpz_sub(padic_unit(rop), padic_unit(op1), padic_unit(rop));
            }

            padic_val(rop) = padic_val(op1);
        }
        else if (padic_val(op1) > padic_val(op2))  /* p^{v1-v2} u1 - u2 */
        {
            _padic_ctx_pow_ui(pow, &alloc, padic_val(op1) - padic_val(op2), ctx);

            if (rop != op1)
            {
                fmpz_neg(padic_unit(rop), padic_unit(op2));
                fmpz_addmul(padic_unit(rop), pow, padic_unit(op1));
            }
            else
            {
                fmpz_mul(padic_unit(rop), pow, padic_unit(op1));
                fmpz_sub(padic_unit(rop), padic_unit(rop), padic_unit(op2));
            }

            padic_val(rop) = padic_val(op2);
        }
        else  /* u1 - u2 */
        {
            fmpz_sub(padic_unit(rop), padic_unit(op1), padic_unit(op2));
            padic_val(rop) = padic_val(op1);
            _padic_canonicalise(rop, ctx);
        }

        if (alloc)
            fmpz_clear(pow);
    }
}

void padic_sub(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op1, ctx))
    {
        padic_neg(rop, op2, ctx);
        return;
    }
    if (padic_is_zero(op2, ctx))
    {
        padic_set(rop, op1, ctx);
        return;
    }

    _padic_sub(rop, op1, op2, ctx);
    _padic_reduce(rop, ctx);
}

