/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

__mpz_struct * _fmpz_new_mpz(void)
{
    __mpz_struct * mf = (__mpz_struct *) flint_malloc(sizeof(__mpz_struct));
    mpz_init2(mf, 2*FLINT_BITS);
    return mf;
}

void _fmpz_clear_mpz(fmpz f)
{
    mpz_clear(COEFF_TO_PTR(f));
    flint_free(COEFF_TO_PTR(f));  
}

void _fmpz_cleanup_mpz_content(void)
{
}

void _fmpz_cleanup(void)
{
}

__mpz_struct * _fmpz_promote(fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))  /* f is small so promote it first */
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        return mf;
    }
    else  /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

__mpz_struct * _fmpz_promote_val(fmpz_t f)
{
    fmpz c = *f;
    if (!COEFF_IS_MPZ(c))  /* f is small so promote it */
    {
        __mpz_struct * mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
        flint_mpz_set_si(mf, c);
        return mf;
    }
    else  /* f is large already, just return the pointer */
        return COEFF_TO_PTR(*f);
}

void _fmpz_demote_val(fmpz_t f)
{
    __mpz_struct * mf = COEFF_TO_PTR(*f);
    int size = mf->_mp_size;

    if (!(((unsigned int) size + 1U) & ~2U))  /* size +-1 */
    {
        ulong uval = mf->_mp_d[0];

        if (uval <= (ulong) COEFF_MAX)
        {
            _fmpz_clear_mpz(*f);
            *f = size * (fmpz) uval;
        }
    }
    else if (size == 0)  /* value is 0 */
    {
        _fmpz_clear_mpz(*f);
        *f = 0;
    }

    /* don't do anything if value has to be multi precision */
}

void _fmpz_init_readonly_mpz(fmpz_t f, const mpz_t z)
{
   __mpz_struct * mf = (__mpz_struct *) flint_malloc(sizeof(__mpz_struct));
    *f = PTR_TO_COEFF(mf);
    *mf = *z;
}

void _fmpz_clear_readonly_mpz(mpz_t z)
{
    if (((z->_mp_size == 1 || z->_mp_size == -1) && (z->_mp_d[0] <= COEFF_MAX))
        || (z->_mp_size == 0))
    {
        mpz_clear(z);
    }
}
