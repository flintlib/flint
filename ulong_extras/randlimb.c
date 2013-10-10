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

#if FLINT64

mp_limb_t n_randlimb(flint_rand_t state) 
{   
    state->__randval = (state->__randval*UWORD(13282407956253574709) + UWORD(286824421));
    state->__randval2 = (state->__randval2*UWORD(7557322358563246341) + UWORD(286824421));

    return (state->__randval>>32) + ((state->__randval2>>32) << 32);
}

#else

mp_limb_t n_randlimb(flint_rand_t state) 
{   
    state->__randval = (state->__randval*UWORD(1543932465) +  UWORD(1626832771));
    state->__randval2 = (state->__randval2*UWORD(2495927737) +  UWORD(1626832771));

    return (state->__randval>>16) + ((state->__randval2>>16) << 16);
}

#endif
