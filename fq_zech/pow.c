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

#include "fq_zech.h"

void
fq_zech_pow(fq_zech_t rop, const fq_zech_t op, const fmpz_t e,
            const fq_zech_ctx_t ctx)
{
    if (fmpz_sgn(e) < 0)
    {
        flint_printf("Exception (fq_zech_pow).  e < 0.\n");
        flint_abort();
    }
    if (fmpz_is_zero(e))
    {
        fq_zech_one(rop, ctx);
    }
    else if (fq_zech_is_zero(op, ctx))
    {
        fq_zech_zero(rop, ctx);
    }
    else if (fmpz_is_one(e))
    {
        fq_zech_set(rop, op, ctx);
    }
    else
    {
        fmpz_t new_e;
        fmpz_init(new_e);
        fmpz_mul_ui(new_e, e, op->value);
        fmpz_mod_ui(new_e, new_e, ctx->qm1);
        rop->value = fmpz_get_ui(new_e);
        fmpz_clear(new_e);
    }
}

/* TODO: Move into separate function and optimize */
void fq_zech_pow_ui(fq_zech_t rop, const fq_zech_t op,
                    const ulong e, const fq_zech_ctx_t ctx)
{
    fmpz_t exp;
    fmpz_init_set_ui(exp, e);
    fq_zech_pow(rop, op, exp, ctx);
    fmpz_clear(exp);
}
