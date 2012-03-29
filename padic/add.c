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

void _padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
                const padic_ctx_t ctx)
{
    if (_padic_is_zero(op1))
    {
        _padic_set(rop, op2);
    }
    else if (_padic_is_zero(op2))
    {
        _padic_set(rop, op1);
    }
    else if (padic_val(op1) == padic_val(op2))
    {
        fmpz_add(padic_unit(rop), padic_unit(op1), padic_unit(op2));
        padic_val(rop) = padic_val(op1);
        _padic_canonicalise(rop, ctx);
    }
    else
    {
        fmpz_t pow;

        fmpz_init(pow);
        if (padic_val(op1) < padic_val(op2))  /* u1 + p^{v2-v1} u2 */
        {
            fmpz_pow_ui(pow, ctx->p, padic_val(op2) - padic_val(op1));

            if (rop != op2)
            {
                fmpz_set(padic_unit(rop), padic_unit(op1));
                fmpz_addmul(padic_unit(rop), pow, padic_unit(op2));
            }
            else
            {
                fmpz_mul(padic_unit(rop), pow, padic_unit(op2));
                fmpz_add(padic_unit(rop), padic_unit(rop), padic_unit(op1));
            }

            padic_val(rop) = padic_val(op1);
        }
        else  /* p^{v1-v2} u1 + u2 */
        {
            fmpz_pow_ui(pow, ctx->p, padic_val(op1) - padic_val(op2));

            if (rop != op1)
            {
                fmpz_set(padic_unit(rop), padic_unit(op2));
                fmpz_addmul(padic_unit(rop), pow, padic_unit(op1));
            }
            else
            {
                fmpz_mul(padic_unit(rop), pow, padic_unit(op1));
                fmpz_add(padic_unit(rop), padic_unit(rop), padic_unit(op2));
            }

            padic_val(rop) = padic_val(op2);
        }
        fmpz_clear(pow);
    }
}

void padic_add(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (fmpz_is_zero(padic_unit(op1)))
    {
        _padic_set(rop, op2);
    }
    else if (fmpz_is_zero(padic_unit(op2)))
    {
        _padic_set(rop, op1);
    }
    else
    {
        if (padic_val(op1) == padic_val(op2))
        {
            int alloc;
            fmpz_t pow;

            fmpz_add(padic_unit(rop), padic_unit(op1), padic_unit(op2));
            padic_val(rop) = padic_val(op1);

            alloc = _padic_ctx_pow_ui(pow, ctx->N - padic_val(rop), ctx);

            if (fmpz_cmpabs(padic_unit(rop), pow) >= 0)
                fmpz_sub(padic_unit(rop), padic_unit(rop), pow);

            if (alloc)
                fmpz_clear(pow);

            _padic_canonicalise(rop, ctx);
        }
        else if (padic_val(op1) < padic_val(op2))
        {
            int alloc;
            fmpz_t f, pow;

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

            alloc = _padic_ctx_pow_ui(pow, ctx->N - padic_val(rop), ctx);

            if (fmpz_cmpabs(padic_unit(rop), pow) >= 0)
                fmpz_sub(padic_unit(rop), padic_unit(rop), pow);

            if (alloc)
                fmpz_clear(pow);
        }
        else  /* padic_val(op1) > padic_val(op2) */
        {
            int alloc;
            fmpz_t f, pow;

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

            alloc = _padic_ctx_pow_ui(pow, ctx->N - padic_val(rop), ctx);

            if (fmpz_cmpabs(padic_unit(rop), pow) >= 0)
                fmpz_sub(padic_unit(rop), padic_unit(rop), pow);

            if (alloc)
                fmpz_clear(pow);
        }
    }
}

