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

    Copyright (C) 2009, 2013 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t
n_powmod2_ui_preinv(mp_limb_t a, mp_limb_t exp, mp_limb_t n, mp_limb_t ninv)
{
    mp_limb_t x;
   
    if (n == UWORD(1) || (a == 0 && exp != 0)) return UWORD(0);

    x = UWORD(1);
    
    if (exp)
    {
       while ((exp & 1) == 0)
       {
          a = n_mulmod2_preinv(a, a, n, ninv);
          exp >>= 1;
       }

       if (a >= n)
          x = n_mod2_preinv(a, n, ninv);
       else
          x = a;
       
       while (exp >>= 1)
       {
          a = n_mulmod2_preinv(a, a, n, ninv);
          if (exp & 1) x = n_mulmod2_preinv(x, a, n, ninv);
       }
    }

    return x;
}

mp_limb_t
n_powmod2_preinv(mp_limb_t a, mp_limb_signed_t exp, mp_limb_t n, mp_limb_t ninv)
{
    ulong norm;
    
    if (exp < WORD(0))
    {
        a = n_invmod(a, n);
        exp = -exp;
    }

    count_leading_zeros(norm, n);

    return n_powmod_ui_preinv(a << norm, exp, n << norm, ninv, norm) >> norm;
}

