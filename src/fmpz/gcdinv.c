/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "ulong_extras.h"
#include "fmpz.h"

void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)
{
    FLINT_ASSERT(fmpz_cmp(f, g) < 0);

    if (fmpz_is_zero(f))
    {
        fmpz_set(d, g);
        fmpz_set_ui(a, 0);
        return;
    }

    if (!COEFF_IS_MPZ(*g))  /* g is small, hence f is small */
    {
        fmpz ff, gg;
        ff = *f;
        gg = *g;

        _fmpz_demote(d);
        _fmpz_demote(a);

        *d = n_gcdinv((ulong *) a, ff, gg);
    }
    else  /* g is large */
    {
        mpz_ptr atemp = _fmpz_new_mpz(), dtemp = _fmpz_new_mpz();

        _fmpz_promote_val(d);
        _fmpz_promote_val(a);

        if (!COEFF_IS_MPZ(*f))  /* f is small */
        {
            mpz_t fptr;

            fptr->_mp_alloc = 1;
            fptr->_mp_size  = 1;
            fptr->_mp_d     = (ulong *) f;

            mpz_gcdext(dtemp, atemp, NULL,
                    fptr, COEFF_TO_PTR(*g));
        }
        else  /* f is large */
        {
            mpz_gcdext(dtemp, atemp, NULL,
                    COEFF_TO_PTR(*f), COEFF_TO_PTR(*g));
        }

        if (mpz_cmp_ui(atemp, 0) < 0)
            mpz_add(atemp, atemp, COEFF_TO_PTR(*g));

        _fmpz_clear_mpz(*d);
        _fmpz_clear_mpz(*a);

        *d = PTR_TO_COEFF(dtemp);
        *a = PTR_TO_COEFF(atemp);

        _fmpz_demote_val(d);
        _fmpz_demote_val(a);
    }
}
