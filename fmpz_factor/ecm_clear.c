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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

void
fmpz_factor_ecm_clear(ecm_t ecm_inf)
{
    fmpz_clear(ecm_inf->t);
    fmpz_clear(ecm_inf->u);
    fmpz_clear(ecm_inf->v);
    fmpz_clear(ecm_inf->w);

    fmpz_clear(ecm_inf->x);
    fmpz_clear(ecm_inf->z);

    fmpz_clear(ecm_inf->a24);
    
} 
