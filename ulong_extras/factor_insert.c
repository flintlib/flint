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

void n_factor_insert(n_factor_t * factors, mp_limb_t p, ulong exp)
{
   ulong i;
   
   for (i = 0; i < factors->num; i++)
   {
      if (factors->p[i] == p) break;
   }
   
   if (i != factors->num) 
   {
      factors->exp[i] += exp;
   } else
   {
      factors->p[i] = p;
      factors->exp[i] = exp;
      factors->num++;
   }
}
