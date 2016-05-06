/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)
{
    FLINT_ASSERT(fmpz_cmp(f, g) < 0);
    
    if (fmpz_is_zero(f))
    {
        fmpz_set(d, g);
        return;
    }

    if (!COEFF_IS_MPZ(*g))  /* g is small, hence f is small */
    {
        fmpz ff, gg;
        ff = *f;
        gg = *g;

        _fmpz_demote(d);
        _fmpz_demote(a);

        *d = n_gcdinv((mp_limb_t *) a, ff, gg);
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

