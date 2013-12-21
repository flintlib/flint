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
fq_zech_add(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2,
            const fq_zech_ctx_t ctx)
{
    mp_limb_t index, c;
    if (op1->value == ctx->qm1)
    {
        rop->value = op2->value;
    }
    else if (op2->value == ctx->qm1)
    {
        rop->value = op1->value;
    }
    else
    {
        index = n_submod(op1->value, op2->value, ctx->qm1);

        c = ctx->zech_log_table[index];
        if (c != ctx->qm1)
        {
            c = n_addmod(c, op2->value, ctx->qm1);
        }
        rop->value = c;
    }
}
