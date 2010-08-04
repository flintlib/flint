/*============================================================================

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

===============================================================================*/
/****************************************************************************

   Copyright (C) 2010 William Hart
   Copyright (C) 2010 Daniel Woodhouse
   
*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mpoly.h"

void nmod_mpoly_init_preinv(nmod_mpoly_t poly, mp_limb_t n, mp_limb_t ninv, long vars, ulong ebits)
{
   poly->coeffs = NULL;
   poly->exps = NULL;
   
   poly->alloc = 0;
   poly->length = 0;

   poly->vars = vars;
   poly->ebits = ebits;

   poly->mod.n = n;
   poly->mod.ninv = ninv;
   count_leading_zeros(poly->mod.norm, n);
}

void nmod_mpoly_init(nmod_mpoly_t poly, mp_limb_t n, long vars, ulong ebits)
{
   nmod_mpoly_init_preinv(poly, n, n_preinvert_limb(n), vars, ebits);
}
