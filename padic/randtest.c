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

void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
{
    long min, max;
    fmpz_t pow;
    int alloc;

    if (ctx->N > 0)
    {
        min = - ((ctx->N + 9) / 10);
        max = ctx->N;
    }
    else if (ctx->N < 0)
    {
        min = ctx->N - ((-ctx->N + 9) / 10);
        max = ctx->N;
    }
    else  /* ctx->N == 0 */
    {
        min = -10;
        max = 0;
    }

    padic_val(rop) = n_randint(state, max - min) + min;

    _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(rop), ctx);
    fmpz_randm(padic_unit(rop), state, pow);
    _padic_canonicalise(rop, ctx);
    if (alloc)
        fmpz_clear(pow);
}

void padic_randtest_not_zero(padic_t rop, flint_rand_t state, 
                             const padic_ctx_t ctx)
{
    long i;

    padic_randtest(rop, state, ctx);

    for (i = 0; !_padic_is_zero(rop) && i < 10; i++)
        padic_randtest(rop, state, ctx);

    if (_padic_is_zero(rop))
    {
        fmpz_set_ui(padic_unit(rop), 1);
        padic_val(rop) = ctx->N - 1;
    }
}

