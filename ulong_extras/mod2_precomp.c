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

mp_limb_t
n_mod2_precomp(mp_limb_t a, mp_limb_t n, double npre)
{
    mp_limb_t quot;
    long rem;

    if (a < n)
        return a;
    if ((mp_limb_signed_t) n < 0L)
        return a - n;

    if (n == 1)
    {
        quot = a;
        rem = 0;
    } else
    {
        quot = (mp_limb_t) ((double) a * npre);
        rem  = a - quot * n;
    }
    
    if (rem < (mp_limb_signed_t) (-n))
        quot -= (mp_limb_t) ((double) (-rem) * npre);
    else if (rem >= (long) n)
        quot += (mp_limb_t) ((double) rem * npre);
    else if (rem < 0L)
        return rem + n;
    else
        return rem;
    
    rem = a - quot * n;
    if (rem >= (long) n)
        return rem - n;
    else if (rem < 0L)
        return rem + n;
    else
        return rem;
}
