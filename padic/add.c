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

#include "padic.h"

void padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_prec(rop) <= FLINT_MIN(padic_val(op1), padic_val(op2)))
    {
        padic_zero(rop);
        return;
    }

    if (padic_is_zero(op1))
    {
        padic_set(rop, op2, ctx);
    }
    else if (padic_is_zero(op2))
    {
        padic_set(rop, op1, ctx);
    }
    else
    {
        if (padic_val(op1) == padic_val(op2))
        {

            fmpz_add(padic_unit(rop), padic_unit(op1), padic_unit(op2));
            padic_val(rop) = padic_val(op1);

            _padic_canonicalise(rop, ctx);

            if (padic_prec(rop) <= padic_val(rop))
            {
                padic_zero(rop);
                return;
            }
        }
        else if (padic_val(op1) < padic_val(op2))
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_pow_ui(f, ctx->p, padic_val(op2) - padic_val(op1));
            if (rop != op2)
            {
                fmpz_set(padic_unit(rop), padic_unit(op1));
                fmpz_addmul(padic_unit(rop), f, padic_unit(op2));
            }
            else
            {
                fmpz_mul(padic_unit(rop), f, padic_unit(op2));
                fmpz_add(padic_unit(rop), padic_unit(rop), padic_unit(op1));
            }
            fmpz_clear(f);

            padic_val(rop) = padic_val(op1);
        }
        else  /* padic_val(op1) > padic_val(op2) */
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_pow_ui(f, ctx->p, padic_val(op1) - padic_val(op2));
            if (rop != op1)
            {
                fmpz_set(padic_unit(rop), padic_unit(op2));
                fmpz_addmul(padic_unit(rop), f, padic_unit(op1));
            }
            else
            {
                fmpz_mul(padic_unit(rop), f, padic_unit(op1));
                fmpz_add(padic_unit(rop), padic_unit(rop), padic_unit(op2));
            }
            fmpz_clear(f);

            padic_val(rop) = padic_val(op2);
        }

        /* Reduce */
        {
            int alloc;
            fmpz_t pow;

            alloc = _padic_ctx_pow_ui(pow, padic_prec(rop) - padic_val(rop), ctx);

            if (padic_prec(rop) >= FLINT_MAX(padic_prec(op1), padic_prec(op2)))
            {
                if (fmpz_cmpabs(padic_unit(rop), pow) >= 0)
                    fmpz_sub(padic_unit(rop), padic_unit(rop), pow);
            }
            else
            {
                fmpz_mod(padic_unit(rop), padic_unit(rop), pow);
            }

            if (fmpz_is_zero(padic_unit(rop)))
                padic_val(rop) = 0;

            if (alloc)
                fmpz_clear(pow);
        }
    }
}

