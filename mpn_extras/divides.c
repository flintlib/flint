/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

int flint_mpn_divides(mp_ptr q, mp_srcptr array1, 
      mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp)
{
   mpn_tdiv_qr(q, temp, 0, array1, limbs1, arrayg, limbsg);
   while ((limbsg) && temp[limbsg - 1] == 0) limbsg--;

   return (limbsg == 0);
}
