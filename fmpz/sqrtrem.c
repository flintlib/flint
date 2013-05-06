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

void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g)
{
    if (fmpz_sgn(g) < 0)
    {
        printf("Exception (fmpz_sqrtrem). g is negative.\n");
        abort();
    }
    
    if (!COEFF_IS_MPZ(*g))
    {
        if (COEFF_IS_MPZ(*r))
            _fmpz_clear_mpz(*r);
        fmpz_set_ui(f, n_sqrtrem((mp_limb_t *) r, *g));
    }
    else
    {
        __mpz_struct * r_mpz_ptr, * f_mpz_ptr;
        _fmpz_promote(f); /* must not hang on to pointer whilst promoting */
        r_mpz_ptr = _fmpz_promote(r);
		f_mpz_ptr = COEFF_TO_PTR(*f);

        mpz_sqrtrem(f_mpz_ptr, r_mpz_ptr, COEFF_TO_PTR(*g));
        _fmpz_demote_val(f);
        _fmpz_demote_val(r);
    }
}
