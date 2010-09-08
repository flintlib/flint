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
#include "longlong.h"
#include "nmod_vec.h"


int main(void)
{
   int result;
   printf("mul_heap....");
   fflush(stdout);
   
   // square (1 + xy + x^2)
   nmod_mpoly_t p;
   nmod_mpoly_t q;
   nmod_mpoly_t r;
   nmod_mpoly_t copyTest;
   mp_limb_t n = n_randtest_not_zero();

   nmod_mpoly_init2(p, n, 3, (long) 2, (ulong) 5);
   nmod_mpoly_init2(copyTest, n, 3, (long) 2, (ulong) 5);
   nmod_mpoly_init2(q, n, 3, (long) 2, (ulong) 5);
   nmod_mpoly_init2(r, n, 3, (long) 2, (ulong) 5);
   //nmod_mpoly_realloc(a, 3);
   
   nmod_mpoly_set(copyTest, p);
   if(!(nmod_mpoly_equal(copyTest,p)))
	printf("FAIL:\n");

   p->coeffs[0] = (mp_limb_t) 1;
   p->exps[0] = (mp_limb_t) 0;
   p->coeffs[1] = (mp_limb_t) 1;
   p->exps[1] = (mp_limb_t) 1 + ((mp_limb_t)1 <<5);  
   p->coeffs[2] = (mp_limb_t) 18446744073709551611;
   p->exps[2] = (mp_limb_t) 2 << 5;

   q->coeffs[0] = (mp_limb_t) 1;
   q->exps[0] = (mp_limb_t) 0;
   q->coeffs[1] = (mp_limb_t) 1;
   q->exps[1] = (mp_limb_t) 1 + ((mp_limb_t)1 <<5);  
   q->coeffs[2] = (mp_limb_t) 18446744073709551611;
   q->exps[2] = (mp_limb_t) 2 << 5;

   p->length = 3;
   q->length = 3;

   nmod_mpoly_mul_heap(r, p, q);
  

   if(r->length != (mp_limb_t) 6){
      printf("\nsquare (1 + xy + x^2)\n");
      printf("\n\nr->length = %u \n\n", r->length);
      int i;
      for(i=0; i < r->length; i++){
	   printf("%u %u + ", r->coeffs[i], r->exps[i]);
      } 
      printf("\nFAIL:\n");
   }

   if(nmod_mpoly_equal(copyTest,p))
	printf("FAIL:\n");

   if(nmod_mpoly_get_coeff(r, ((ulong)4 << 5)) != nmod_mpoly_get_coeff_of_product(p,q,((ulong)4 << 5))){
      printf("\n%lu", nmod_mpoly_get_coeff(r, (4 << 5)));
      printf("\n%lu", nmod_mpoly_get_coeff_of_product(p,q,((ulong)4 << 5)));
      printf("\nFAAIL:\n"); 
   }
/*
   if(nmod_mpoly_get_coeff_of_product(p,q,((ulong)4 << 5)) != (ulong)1){
      printf("\n%lu", nmod_mpoly_get_coeff_of_product(p,q,((ulong)4 << 5)));
      printf("\nFAAAIL:\n"); 
   }
*/
   nmod_mpoly_clear(p);     
   nmod_mpoly_clear(q);
   nmod_mpoly_clear(r);
   nmod_mpoly_clear(copyTest);

   //multiply (w + 2*x + 3*y + 4*z)^7
  
   nmod_mpoly_init2(p, 7, 4, (long) 4, (ulong) 6);
   nmod_mpoly_init2(q, 7, 4, (long) 4, (ulong) 6);
   nmod_mpoly_init2(r, 7, 4, (long) 4, (ulong) 6);

   p->coeffs[0] = (mp_limb_t) 4;
   p->exps[0] = (mp_limb_t) 1;
   p->coeffs[1] = (mp_limb_t) 3;
   p->exps[1] = (mp_limb_t) 1 << 6;
   p->coeffs[2] = (mp_limb_t) 2;
   p->exps[2] = (mp_limb_t) 1 << 12;
   p->coeffs[3] = (mp_limb_t) 1;
   p->exps[3] = (mp_limb_t) 1 << 18;

   p->length = 4;

   nmod_mpoly_mul_heap(q, p, p);
   nmod_mpoly_mul_heap(r, p, q);
/*   nmod_mpoly_mul_heap(r, p, r);
   nmod_mpoly_mul_heap(r, p, r);
   nmod_mpoly_mul_heap(r, p, r);
   nmod_mpoly_mul_heap(r, p, r);*/

   if(r->length != (mp_limb_t) 20){
      printf("\n (w + 2*x + 3*y + 4*z)^7 ");
      printf("\n\nr->length = %u \n\n", r->length);
      int i;
      for(i=0; i < r->length; i++){
	   printf("%u %u + ", r->coeffs[i], r->exps[i]);
      } 
      printf("\nFAIL:\n");
   }

   nmod_mpoly_clear(p);     
   
   nmod_mpoly_clear(r);

   nmod_mpoly_t a, a2, a3, a6, a12;
   nmod_mpoly_t c, c2, c3, c5, c10, c20, c30;
   ulong i;   

   nmod_mpoly_init2(a, (ulong) 100003, 6,  (long) 5, (ulong) 12);
   
   nmod_mpoly_init(a2, (ulong) 100003, 5, 12);
   nmod_mpoly_init(a3, (ulong) 100003, 5, 12);
   nmod_mpoly_init(a6, (ulong) 100003, 5, 12);
   nmod_mpoly_init(a12, (ulong) 100003, 5, 12);
   
   a->coeffs[0] = 1;
   a->coeffs[1] = 5;
   a->coeffs[2] = 3;
   a->coeffs[3] = 2;
   a->coeffs[4] = 1;
   a->coeffs[5] = 1;
   
   a->exps[0] = 0;
   a->exps[1] = 5;
   a->exps[2] = 3UL<<12;
   a->exps[3] = 2UL<<24;
   a->exps[4] = 1UL<<36;
   a->exps[5] = 1UL<<48;

   a->length = 6;

   nmod_mpoly_mul_heap(a2, a, a);
   nmod_mpoly_mul_heap(a3, a2, a);
   nmod_mpoly_mul_heap(a6, a3, a3);
   nmod_mpoly_mul_heap(a12, a6, a6);

   if (a12->length != 6188)
   {
	   printf("FAIL:\n");
	   printf("length = %ld\n", a12->length);
   }

   if(nmod_mpoly_get_coeff(a12, (ulong)5 <<36 + (ulong)3 <<48) != nmod_mpoly_get_coeff_of_product(a6,a6,(ulong)5 <<36 + (ulong)3 <<48)){
	   printf("FAIL:\n");
	   printf("val1 = %ld\n", nmod_mpoly_get_coeff(a12, (ulong)5 <<36 + (ulong)3 <<48));
           printf("val2 = %ld\n", nmod_mpoly_get_coeff_of_product(a6,a6,(ulong)5 <<36 + (ulong)3 <<48));
   }

   for (i = 0; i < 3; i++)
	  a->coeffs[i] = 1;
   a->coeffs[3] = 2;
   a->coeffs[4] = 3;
   a->coeffs[5] = 5;

   a->exps[0] = 0;
   a->exps[1] = 1;
   a->exps[2] = 1UL<<12;
   a->exps[3] = 2UL<<24;
   a->exps[4] = 3UL<<36;
   a->exps[5] = 5UL<<48;

   nmod_mpoly_mul_heap(a2, a, a);
   nmod_mpoly_mul_heap(a3, a2, a);
   nmod_mpoly_mul_heap(a6, a3, a3);
   nmod_mpoly_mul_heap(a12, a6, a6);

   a->length = 6;
 
   if (a12->length != 6188)
   {
	   printf("FAIL:\n");
	   printf("length = %ld\n", a12->length);
   }   

   nmod_mpoly_init2(c, (ulong) 100003, 5, 4, 15);
   nmod_mpoly_init(c2, (ulong) 100003, 4, 15);
   nmod_mpoly_init(c3, (ulong) 100003, 4, 15);
   nmod_mpoly_init(c5, (ulong) 100003, 4, 15);
   nmod_mpoly_init(c10, (ulong) 100003, 4, 15);
   nmod_mpoly_init(c20, (ulong) 100003, 4, 15);
   nmod_mpoly_init(c30, (ulong) 100003, 4, 15);
 
   for (i = 0; i < 5; i++)
	   c->coeffs[i] = 1;

   c->exps[0] = 0;
   for (i = 1; i < 5; i++)
	   c->exps[i] = 1UL<<((i - 1)*15);
  
   c->length = 5;

   nmod_mpoly_mul_heap(c2, c, c);
   nmod_mpoly_mul_heap(c3, c2, c);
   nmod_mpoly_mul_heap(c5, c3, c2);
   nmod_mpoly_mul_heap(c10, c5, c5);
   nmod_mpoly_mul_heap(c20, c10, c10);
   nmod_mpoly_mul_heap(c30, c20, c10);

   if (c30->length != 46376)
   {
	   printf("FAIL:\n");
	   printf("length = %ld\n", c30->length);
   }   

   nmod_mpoly_clear(c);      
   nmod_mpoly_clear(c2);      
   nmod_mpoly_clear(c3);      
   nmod_mpoly_clear(c5);      
   nmod_mpoly_clear(c10);      
   nmod_mpoly_clear(c20);      
   nmod_mpoly_clear(c30);    

   nmod_mpoly_clear(a);      
   nmod_mpoly_clear(a2);      
   nmod_mpoly_clear(a3);      
   nmod_mpoly_clear(a6);      
   nmod_mpoly_clear(a12);     

   printf("PASS\n");
   return 0;
}
