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
fq_zech_neg(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
{
    if (op1->value == ctx->qm1)
    {
        rop->value = ctx->qm1;
        return;
    }
    rop->value = op1->value + ctx->qm1o2;
    if (rop->value >= ctx->qm1)
    {
        rop->value -= ctx->qm1;
    }
}
