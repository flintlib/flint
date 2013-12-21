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
#include <stdio.h>
#include <string.h>

#include "fq_zech.h"

void
fq_zech_inv(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    if (fq_zech_is_zero(op, ctx))
    {
        flint_printf("Exception (fq_inv).  Zero is not invertible.\n");
        abort();
    }
    if (fq_zech_is_one(op, ctx))
    {
        fq_zech_one(rop, ctx);
        return;
    }

    rop->value = ctx->qm1 - op->value;
}
