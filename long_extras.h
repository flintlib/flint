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

#ifndef LONG_EXTRAS_H
#define LONG_EXTRAS_H

#include <gmp.h>
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* Properties ****************************************************************/

size_t z_sizeinbase(len_t n, int b);

/* Randomisation  ************************************************************/

mp_limb_signed_t z_randtest(flint_rand_t state);

mp_limb_signed_t z_randtest_not_zero(flint_rand_t state);

mp_limb_signed_t z_randint(flint_rand_t state, mp_limb_t limit);

#ifdef __cplusplus
}
#endif

#endif

