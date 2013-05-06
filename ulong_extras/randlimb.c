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
    state->__randval = (state->__randval*13282407956253574709UL + 286824421UL);
    state->__randval2 = (state->__randval2*7557322358563246341UL + 286824421UL);

    return (state->__randval>>32) + ((state->__randval2>>32) << 32);
}

#else

mp_limb_t n_randlimb(flint_rand_t state) 
{   
    state->__randval = (state->__randval*1543932465UL +  1626832771UL);
    state->__randval2 = (state->__randval2*2495927737UL +  1626832771UL);

    return (state->__randval>>16) + ((state->__randval2>>16) << 16);
}

#endif
