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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_sqrt(fmpz_t f, const fmpz_t g)
{
    if (fmpz_sgn(g) < 0)
    {
        printf("Exception (fmpz_sqrt). g is negative.\n");
        abort();
    }
    
    if (!COEFF_IS_MPZ(*g))
        fmpz_set_ui(f, n_sqrt(*g));
    else
    {
        __mpz_struct * mpz_ptr = _fmpz_promote(f);
        mpz_sqrt(mpz_ptr, COEFF_TO_PTR(*g));
        _fmpz_demote_val(f);
    }
}
