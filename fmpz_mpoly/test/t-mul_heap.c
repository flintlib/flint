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

*****************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int main(void)
{
   int result;
   printf("mul_heap....");
   fflush(stdout);
   
   fmpz_mpoly_t a, a2, a3, a6, a12, b, b2, b3, b6, b12, c;
   ulong i;

   fmpz_mpoly_init2(a, 6, 5, 12);
   fmpz_mpoly_init2(b, 6, 5, 12);
   fmpz_mpoly_init(a2, 5, 12);
   fmpz_mpoly_init(a3, 5, 12);
   fmpz_mpoly_init(a6, 5, 12);
   fmpz_mpoly_init(a12, 5, 12);
   fmpz_mpoly_init(b2, 5, 12);
   fmpz_mpoly_init(b3, 5, 12);
   fmpz_mpoly_init(b6, 5, 12);
   fmpz_mpoly_init(b12, 5, 12);
   fmpz_mpoly_init(c, 5, 12);
   
   for (i = 0; i < 3; i++)
	  fmpz_set_ui(a->coeffs + i, 1);
   fmpz_set_ui(a->coeffs + 3, 2);
   fmpz_set_ui(a->coeffs + 4, 4);
   fmpz_set_ui(a->coeffs + 5, 5);
   
   fmpz_set_ui(b->coeffs, 1);
   fmpz_set_ui(b->coeffs + 1, 5);
   fmpz_set_ui(b->coeffs + 2, 3);
   fmpz_set_ui(b->coeffs + 3, 2);
   fmpz_set_ui(b->coeffs + 4, 1);
   fmpz_set_ui(b->coeffs + 5, 1);

   fmpz_set_ui(a->exps, 0);
   fmpz_set_ui(a->exps + 1, 1);
   fmpz_set_ui(a->exps + 2, 1UL<<12);
   fmpz_set_ui(a->exps + 3, 2UL<<24);
   fmpz_set_ui(a->exps + 4, 3UL<<36);
   fmpz_set_ui(a->exps + 5, 5UL<<48);

   fmpz_set_ui(b->exps, 0);
   fmpz_set_ui(b->exps + 1, 5);
   fmpz_set_ui(b->exps + 2, 3UL<<12);
   fmpz_set_ui(b->exps + 3, 2UL<<24);
   fmpz_set_ui(b->exps + 4, 1UL<<36);
   fmpz_set_ui(b->exps + 5, 1UL<<48);

   a->length = 6;
   b->length = 6;
   
   fmpz_mpoly_mul_heap(a2, a, a);
   fmpz_mpoly_mul_heap(a3, a2, a);
   fmpz_mpoly_mul_heap(a6, a3, a3);
   fmpz_mpoly_mul_heap(a12, a6, a6);
   
   fmpz_mpoly_mul_heap(b2, b, b);
   fmpz_mpoly_mul_heap(b3, b2, b);
   fmpz_mpoly_mul_heap(b6, b3, b3);
   fmpz_mpoly_mul_heap(b12, b6, b6);
   
   printf("a12->length = %ld\n", a12->length);
   printf("b12->length = %ld\n", b12->length);

   fmpz_mpoly_mul_heap(c, a12, b12);

   printf("c->length = %ld\n", c->length);

   fmpz_mpoly_clear(a);      
   fmpz_mpoly_clear(a2);      
   fmpz_mpoly_clear(a3);      
   fmpz_mpoly_clear(a6);      
   fmpz_mpoly_clear(a12);      
      
   fmpz_mpoly_clear(b);      
   fmpz_mpoly_clear(b2);      
   fmpz_mpoly_clear(b3);      
   fmpz_mpoly_clear(b6);      
   fmpz_mpoly_clear(b12);      

   fmpz_mpoly_clear(c);

   /*fmpz_mpoly_t c, c2, c3, c5, c10, c20, c30, c60;
   ulong i;

   fmpz_mpoly_init2(c, 5, 4, 15);
   fmpz_mpoly_init(c2, 4, 15);
   fmpz_mpoly_init(c3, 4, 15);
   fmpz_mpoly_init(c5, 4, 15);
   fmpz_mpoly_init(c10, 4, 15);
   fmpz_mpoly_init(c20, 4, 15);
   fmpz_mpoly_init(c30, 4, 15);
   fmpz_mpoly_init(c60, 4, 15);
   
   for (i = 0; i < 5; i++)
      fmpz_set_ui(c->coeffs + i, 1);

   fmpz_set_ui(c->exps, 0);
   for (i = 1; i < 5; i++)
	  fmpz_set_ui(c->exps + i, 1UL<<((i - 1)*15));
  
   c->length = 5;
   
   fmpz_mpoly_mul_heap(c2, c, c);
   fmpz_mpoly_mul_heap(c3, c2, c);
   fmpz_mpoly_mul_heap(c5, c3, c2);
   fmpz_mpoly_mul_heap(c10, c5, c5);
   fmpz_mpoly_mul_heap(c20, c10, c10);
   fmpz_mpoly_mul_heap(c30, c20, c10);
   printf("c30->length = %ld\n", c30->length);

   fmpz_mpoly_mul_heap(c60, c30, c30);

   printf("c60->length = %ld\n", c60->length);

   fmpz_mpoly_clear(c);      
   fmpz_mpoly_clear(c2);      
   fmpz_mpoly_clear(c3);      
   fmpz_mpoly_clear(c5);      
   fmpz_mpoly_clear(c10);      
   fmpz_mpoly_clear(c20);      
   fmpz_mpoly_clear(c30);      
   fmpz_mpoly_clear(c60); */  

   _fmpz_cleanup();
   printf("PASS\n");
   return 0;
}
