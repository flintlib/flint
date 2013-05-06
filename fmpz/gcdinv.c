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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)
{
    if (fmpz_is_zero(f))
    {
        fmpz_set(d, g);
        return;
    }

    if (!COEFF_IS_MPZ(*g))  /* g is small, hence f is small */
    {
        if (COEFF_IS_MPZ(*d))
            _fmpz_demote_val(d);
        if (COEFF_IS_MPZ(*a))
            _fmpz_demote_val(a);

        *d = n_gcdinv((mp_limb_t *) a, *f, *g);
    }
    else  /* g is large */
    {
        _fmpz_promote_val(d);
        _fmpz_promote_val(a);

        if (!COEFF_IS_MPZ(*f))  /* f is small */
        {
            mpz_t fptr;

            fptr->_mp_alloc = 1;
            fptr->_mp_size  = 1;
            fptr->_mp_d     = (mp_limb_t *) f;

            mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), NULL, 
                       fptr, COEFF_TO_PTR(*g));
        }
        else  /* f is large */
        {
            mpz_gcdext(COEFF_TO_PTR(*d), COEFF_TO_PTR(*a), NULL, 
                       COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
        }

        _fmpz_demote_val(d);
        _fmpz_demote_val(a);
    }
}

