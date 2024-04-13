/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define FMPZ_INLINES_C

#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmissing-prototypes"
#endif

#include "fmpz.h"

#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif

#include "gmpcompat.h"

void _fmpz_promote_set_ui(fmpz_t f, ulong v)
{
    __mpz_struct * z = _fmpz_promote(f);
    flint_mpz_set_ui(z, v);
}

void _fmpz_promote_neg_ui(fmpz_t f, ulong v)
{
    __mpz_struct * z = _fmpz_promote(f);
    flint_mpz_set_ui(z, v);
    mpz_neg(z, z);
}

void _fmpz_promote_set_si(fmpz_t f, slong v)
{
    __mpz_struct * z = _fmpz_promote(f);
    flint_mpz_set_si(z, v);
}

void _fmpz_init_promote_set_ui(fmpz_t f, ulong v)
{
    __mpz_struct * z = _fmpz_new_mpz();
    *f = PTR_TO_COEFF(z);
    flint_mpz_set_ui(z, v);
}

void _fmpz_init_promote_set_si(fmpz_t f, slong v)
{
    __mpz_struct * z = _fmpz_new_mpz();
    *f = PTR_TO_COEFF(z);
    flint_mpz_set_si(z, v);
}
