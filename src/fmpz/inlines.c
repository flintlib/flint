/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_INLINES_C

#include <gmp.h>
#include "fmpz.h"

/* NOTE: Herein we rely on that the fmpz collector never assigns an mpz with
   less than two limbs. */
void _fmpz_promote_set_ui(fmpz_t f, ulong v)
{
    mpz_ptr zf = _fmpz_promote(f);
    zf->_mp_d[0] = v;
    zf->_mp_size = 1;
}

void _fmpz_promote_neg_ui(fmpz_t f, ulong v)
{
    mpz_ptr zf = _fmpz_promote(f);
    zf->_mp_d[0] = v;
    zf->_mp_size = -1;
}

void _fmpz_promote_set_si(fmpz_t f, slong v)
{
    mpz_ptr zf = _fmpz_promote(f);
    zf->_mp_d[0] = FLINT_ABS(v);
    zf->_mp_size = (v < 0) ? -1 : 1;
}

void _fmpz_init_promote_set_ui(fmpz_t f, ulong v)
{
    mpz_ptr zf = _fmpz_new_mpz();
    *f = PTR_TO_COEFF(zf);
    zf->_mp_d[0] = v;
    zf->_mp_size = 1;
}

void _fmpz_init_promote_set_si(fmpz_t f, slong v)
{
    mpz_ptr zf = _fmpz_new_mpz();
    *f = PTR_TO_COEFF(zf);
    zf->_mp_d[0] = FLINT_ABS(v);
    zf->_mp_size = (v < 0) ? -1 : 1;
}

/* NOTE: We assume that hi != 0 for the following two functions. */
void _fmpz_promote_set_uiui(fmpz_t f, ulong hi, ulong lo)
{
    mpz_ptr zf = _fmpz_promote(f);
    zf->_mp_d[0] = lo;
    zf->_mp_d[1] = hi;
    zf->_mp_size = 2;
}

void _fmpz_promote_neg_uiui(fmpz_t f, ulong hi, ulong lo)
{
    mpz_ptr zf = _fmpz_promote(f);
    zf->_mp_d[0] = lo;
    zf->_mp_d[1] = hi;
    zf->_mp_size = -2;
}
