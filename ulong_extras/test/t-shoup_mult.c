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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   mp_limb_t w, t, p, w_pr, correct_w_prep;
   FLINT_TEST_INIT(state);
   
   flint_printf("shoup_mult....");
   fflush(stdout);  

   w = 1000;
   t = 14;
   p = 1301;
   correct_w_prep = 14178896290322483947;
   
   w_pr = w_prep(w, p);
   if (w_pr != correct_w_prep)
   {
      flint_printf("FAIL:\n");
      flint_printf("w = %wu, t = %wu, p = %wu\n", w, t, p);
      flint_printf("w' = %wu, should be %wu\n", w_pr, correct_w_prep);
      abort();
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}

