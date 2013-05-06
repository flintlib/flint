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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_mod_poly_factor.h"

void
fmpz_mod_poly_factor_init(fmpz_mod_poly_factor_t fac)
{
    long i;
    fmpz_t p;

    fac->alloc = 5;
    fac->num = 0;
    fac->poly = flint_malloc(sizeof(fmpz_mod_poly_struct) * 5);
    fac->exp = flint_malloc(sizeof(long) * 5);

    fmpz_init_set_ui(p, 5);
    for (i = 0; i < 5; i++)
        fmpz_mod_poly_init(fac->poly + i, p);
    fmpz_clear(p);
}
