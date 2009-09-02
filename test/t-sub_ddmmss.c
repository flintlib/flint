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

   Copyright (C) 2009 William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("sub_ddmmss....");
   fflush(stdout);

   for (ulong i = 0; i < 1000000UL; i++)
   {
      mp_limb_t dh1, dl1, dh2, dl2, mh, ml, sh, sl;

      mh = n_randtest();
      ml = n_randtest();
      sh = n_randtest();
      sl = n_randtest();
      
      sub_ddmmss(dh1, dl1, mh, ml, sh, sl);
      
      dl2 = ml - sl;
      dh2 = -(sl > ml);
      dh2 += mh;
      dh2 -= sh;

      result = ((dh2 == dh1) && (dl2 == dl1));

      if (!result)
      {
         printf("FAIL\n");
         printf("mh = %lu, ml = %lu, sh = %lu, sl = %lu\n", mh, ml, sh, sl); 
         printf("dh2 = %lu, dh1 = %lu, dl2 = %lu, dl1 = %lu\n", dh2, dh1, dl2, dl1);
         abort();
      }
   }

   printf("PASS\n");
   return 0;
}
