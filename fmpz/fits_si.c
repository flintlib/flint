/* int mpz_fits_X_p (mpz_t z) -- test whether z fits signed type X.

Copyright 1997, 2000, 2001, 2002 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"

#if defined(_WIN64) || defined(__mips64)

#define FLINT_UI_MAX          ((mp_limb_t)(~(mp_limb_t)0))
#define FLINT_UI_HIBIT        (FLINT_UI_MAX ^ (FLINT_UI_MAX >> 1))
#define FLINT_SI_MAX          ((mp_limb_signed_t)(FLINT_UI_MAX ^ FLINT_UI_HIBIT))
#define FLINT_SI_MIN          ((mp_limb_signed_t)FLINT_UI_HIBIT)

int
flint_mpz_fits_si_p(mpz_srcptr z)
{
  mp_size_t n = z->_mp_size;
  mp_ptr p = z->_mp_d;
  mp_limb_t limb = p[0];

  if (n == 0)
    return 1;
  if (n == 1)
    return limb <= FLINT_SI_MAX;
  if (n == -1)
    return limb <= (mp_limb_t) FLINT_SI_MIN;
  return 0;
}

#endif

/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

int fmpz_fits_si(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
    {
        return 1;
    }
    else
    {
#if defined(_WIN64) || defined(__mips64)
       return flint_mpz_fits_si_p(COEFF_TO_PTR(*f));
#else
       return mpz_fits_slong_p(COEFF_TO_PTR(*f));
#endif
    }
}

