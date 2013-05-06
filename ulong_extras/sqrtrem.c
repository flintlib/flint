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

******************************************************************************/

#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_sqrtrem(mp_limb_t * r, mp_limb_t a)
{
    mp_limb_t is;

    is = (mp_limb_t) sqrt((double) a);

    is -= (is*is > a);
#if FLINT64
    if (is == 4294967296UL) is--;
#endif
    (*r) = a - is*is;

    return is;
}
