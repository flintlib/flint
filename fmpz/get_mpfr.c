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

    Copyright (C) 2009 William Hart
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#if defined(_WIN64) || defined(__mips64)
#include <stdint.h> /* to enable mpfr_set_sj in mpfr.h */
#endif
#include <gmp.h>
#if defined( _WIN64) && defined( _MSC_MPIR_VERSION ) && __MPIR_RELEASE >= 20700
#  if defined( _MSC_VER ) && _MSC_VER >= 1600
#    include <stdint.h>
#    include <mpfr.h>
#    define mpfr_set_si mpfr_set_sj
#  endif
#endif
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*f))
#if defined(_WIN64) || defined(__mips64)
        mpfr_set_sj(x, *f, rnd);
#else
        mpfr_set_si(x, *f, rnd);    /* set x to small value */
#endif
    else
        mpfr_set_z(x, COEFF_TO_PTR(*f), rnd);   /* set x to large value */
}
