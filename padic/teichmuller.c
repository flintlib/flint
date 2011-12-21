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

void padic_teichmuller(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    int alloc;
    fmpz_t u, x, ppow;

    if (padic_val(op) < 0)
    {
        printf("ERROR (padic_teichmuller).  op is not a p-adic integer.\n");
        abort();
    }

    if (_padic_is_zero(op) || (padic_val(op) > 0))
    {
        padic_zero(rop);
        return;
    }

    fmpz_init(x);
    fmpz_init(u);

    alloc = _padic_ctx_pow_ui(ppow, ctx->N, ctx);

    /* Set x = op mod p^N */
    fmpz_mod(x, padic_unit(op), ppow);
    
    /* Let u be the inverse of 1-p mod p^N */
    fmpz_sub(u, ppow, ctx->p);
    fmpz_add_ui(u, u, 1);
    _padic_inv(u, u, ctx->p, ctx->N);

    /* Let rop = x + u * (x^p - x) mod p^N */
    fmpz_powm(padic_unit(rop), x, ctx->p, ppow);
    fmpz_sub(padic_unit(rop), padic_unit(rop), x);
    fmpz_mul(padic_unit(rop), u, padic_unit(rop));
    fmpz_add(padic_unit(rop), x, padic_unit(rop));
    fmpz_mod(padic_unit(rop), padic_unit(rop), ppow);
    
    /* Repeat this until rop == x mod p^N */
    while (!fmpz_equal(padic_unit(rop), x))
    {
        fmpz_swap(x, padic_unit(rop));
        fmpz_powm(padic_unit(rop), x, ctx->p, ppow);
        fmpz_sub(padic_unit(rop), padic_unit(rop), x);
        fmpz_mul(padic_unit(rop), u, padic_unit(rop));
        fmpz_add(padic_unit(rop), x, padic_unit(rop));
        fmpz_mod(padic_unit(rop), padic_unit(rop), ppow);
    }
    
    padic_val(rop) = 0;
    
    fmpz_clear(x);
    fmpz_clear(u);
    if (alloc)
        fmpz_clear(ppow);
}

