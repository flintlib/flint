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
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_mulmod_precomp(mp_limb_t a, mp_limb_t b, mp_limb_t n, double npre)
{
    mp_limb_t quot;
    mp_limb_signed_t rem;

    quot = (mp_limb_t) ((double) a * (double) b * npre);
    rem  = a * b - quot * n;
    if (rem < 0) 
    {
        rem += n;
        if (rem < 0) return rem + n;
    }
    else if (rem >= n) return rem - n;
    return rem;
}
