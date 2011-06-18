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

void _padic_reduce(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(padic_unit(rop)))
    {
        if (padic_val(rop) >= ctx->N)
        {
            fmpz_zero(padic_unit(rop));
            padic_val(rop) = 0;
        }
        else
        {
            int alloc;
            fmpz_t pow;

            _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(rop), ctx);
            fmpz_mod(padic_unit(rop), padic_unit(rop), pow);
            if (alloc)
                fmpz_clear(pow);
        }
    }
}

void padic_reduce(padic_t rop, const padic_ctx_t ctx)
{
    _padic_canonicalise(rop, ctx);
    _padic_reduce(rop, ctx);
}

