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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_set_mpf(fmpz_t f, const mpf_t x)
{
#if defined(__MPIR_VERSION)
    if (mpf_fits_si_p(x))
#else
    if (flint_mpf_fits_slong_p(x))
#endif
    {
        slong cx = flint_mpf_get_si(x);
        fmpz_set_si(f, cx);
    }
    else
    {
        __mpz_struct *z = _fmpz_promote(f);
        mpz_set_f(z, x);
    }
}
