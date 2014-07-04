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

    Copyright (C) 2014 William Hart

******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, 
                                         const fmpz_t m, const fmpz_t n)
{
    fmpz_t t;
    slong i, B = fmpz_sizeinbase(m, 2);

    fmpz_init(t);
    fmpz_set_ui(Vm, 2);
    fmpz_set(Vm1, A);

    for (i = B - 1; i >= 0; i--)
    {
       if (fmpz_tstbit(m, i)) /* 1 in binary repn */
       {
          fmpz_mul(t, Vm, Vm1);
          fmpz_sub(t, t, A);
          fmpz_mod(Vm, t, n);

          fmpz_mul(t, Vm1, Vm1);
          fmpz_sub_ui(t, t, 2);
          fmpz_mod(Vm1, t, n);
       } else /* 0 in binary repn */
       {
          fmpz_mul(t, Vm, Vm1);
          fmpz_sub(t, t, A);
          fmpz_mod(Vm1, t, n);

          fmpz_mul(t, Vm, Vm);
          fmpz_sub_ui(t, t, 2);
          fmpz_mod(Vm, t, n);
       }
    }

    fmpz_clear(t);
}
