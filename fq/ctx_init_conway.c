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

#include <stdio.h>
#include <string.h>

#include "fmpz_vec.h"
#include "padic.h"
#include "qadic.h"
 /* I'm not sure if I should include all these above(I copied them from the qadic code), as I am including fq.h which should include them anyway */

#include "fq.h"

void
fq_ctx_init_conway(qadic_ctx_t ctx,
                   const fmpz_t p, long d, const char *var,
                   enum padic_print_mode mode)
{
    return qadic_ctx_init_conway(ctx, p, d, 1, *var, mode);
}
