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

#undef ulong /* prevent clash with standard library */
#include <math.h>
#define ulong unsigned long
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int mod64[64] = {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,
                 0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,
                 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0}; 

int mod65[65] = {1,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,
                 0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,
                 0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,1};

int mod63[63] = {1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,
                 0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,
                 0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0};
   
int n_is_square(mp_limb_t x)
{
    mp_limb_t sq;

    if (!mod64[x % 64UL]) return 0;
    if (!mod63[x % 63UL]) return 0;
    if (!mod65[x % 65UL]) return 0;

    sq = (mp_limb_t) (sqrt((double) x) + 0.5);
    
    return (x == sq*sq);
}
