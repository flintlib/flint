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

    Copyright (C) 2012 Andres Goens
 
******************************************************************************/

#include "padic.h"
#include "fmpz.h"

int
padic_ctx_equal(const padic_ctx_t ctx1,const padic_ctx_t ctx2)
{
    if(ctx1->p != ctx2->p) return 0;
    if(ctx1->N != ctx2->N) return 0;
    if(ctx1->p != ctx2->p) return 0;
    if(ctx1->min != ctx2->min) return 0;
    if(ctx1->max != ctx2->max) return 0;
    if(ctx1->mode != ctx2->mode) return 0;
    if(!fmpz_equal(ctx1->pow,ctx2->pow)) return 0;
    return 1;
}
