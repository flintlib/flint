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
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "nmod_mpoly.h"

void fmpz_mpoly_mod_mul(fmpz_mpoly_t res, fmpz_mpoly_t poly1, fmpz_mpoly_t poly2)
{

   //find out how big the coefficients can be so we can generate an appropriate 
   //number of primes.
   mp_bitcnt_t bits1, bits2, maxBitLength;
   bits1 = _fmpz_vec_max_bits(poly1->coeffs, poly1->length);
   bits2 = _fmpz_vec_max_bits(poly2->coeffs, poly2->length);
   
   maxBitLength = bits1 + bits2 + FLINT_BIT_COUNT(poly1->length) + 1;
   
   //allocate an array to store the primes in
   //We don't know how many we are going to use so we
   //have to allocate as we go.
   ulong * primes = (ulong *)malloc(sizeof(ulong)*2);
   long length_of_prime_array = 2;
   long num_of_primes = 0; //The number of primes
   mp_limb_t i = (ulong)1 <<52;
   
   mp_bitcnt_t prime_bits = 0;

   //keep on adding primes to our array until the bit length
   //exceeds the maximum bit legnth.
   while(prime_bits <= maxBitLength){
      if(length_of_prime_array == num_of_primes){
         length_of_prime_array *= 2;
         primes = (ulong *)realloc(primes, length_of_prime_array * sizeof(ulong));
      }
      if(n_is_prime(i)){
         primes[num_of_primes] = i;
         prime_bits += FLINT_BIT_COUNT(i);
         num_of_primes++;
      }
      i = i - 1;
   }
   
   //these are the iterators we are going to use
   // k - to keep track of primes
   // l - to keep track of polynomial length
   // j - to go through the coefficients
   long k, l, j;;
   
   

   //This is recycled from get_period_sequence
   ulong * cInverses = (ulong *)malloc(sizeof(ulong)*(num_of_primes -1)); 
   fmpz_t temp;
   fmpz_t temp1;
   fmpz_t temp2;
   fmpz_init(temp);
   fmpz_init(temp1);
   fmpz_init(temp2);

   //variables to perform CRT with
   fmpz_t m1;
   fmpz_t residue;
   double pre;
   fmpz_t inv;
   fmpz_t mm2;
   fmpz_t out;

   //nmod_polys we will multiply with     
   nmod_mpoly_t ppoly1;
   nmod_mpoly_t ppoly2;
  
   nmod_t mod;
   long coeffTemp;
   ulong coeffTemp2;
    
   //Calculate the inverses of our co-prime pairs
   fmpz_set_ui(temp1, primes[0]);  
  
   for(k=0; k<num_of_primes-1; k++){
      fmpz_set_ui(temp2, primes[k+1]);
      fmpz_invmod(temp, temp1, temp2);
      cInverses[k] = fmpz_get_ui(temp);
      fmpz_mul(temp, temp1, temp2);
      fmpz_set(temp1, temp);      
   };
   

   //Need some nmod_mpolys to store my answers in.
   nmod_mpoly_t * result_nmod_mpolys = (nmod_mpoly_t *)malloc(sizeof(nmod_mpoly_t) * num_of_primes);

   
   
   ulong ebits;
   ebits = FLINT_BITS / poly1->vars; 
    
   //now actually create and multiply our nmod_mpolys
   for(k=0; k<num_of_primes; k++){
      nmod_mpoly_init2(result_nmod_mpolys[k], primes[k], 0, poly1->vars, ebits);
      
      nmod_mpoly_init2(ppoly1, primes[k], poly1->length, (long) poly1->vars, (ulong) ebits);
      nmod_mpoly_init2(ppoly2, primes[k], poly2->length, (long) poly2->vars, (ulong) ebits);
      
      
      nmod_init(&mod, primes[k]);
       
      l = 0;
      //convert our poly1 to a nmod_mpoly
      for(j = 0;  j<poly1->length; j++){
         coeffTemp = fmpz_get_si(poly1->coeffs +j);
         coeffTemp2 = fmpz_get_ui(poly1->coeffs +j);
         if(coeffTemp < 0){
            while(coeffTemp <0)
               coeffTemp = coeffTemp + (long)primes[k];
            if(coeffTemp != 0){ 
               ppoly1->coeffs[j] = (ulong) coeffTemp;
               l++;
               ppoly1->exps[j] = fmpz_get_si(poly1->exps + j);
            }
         }
         else{
            NMOD_RED(coeffTemp2, coeffTemp2 , mod);
            if(coeffTemp2 != 0){
               ppoly1->coeffs[j] = (ulong) coeffTemp2;
               l++;
               ppoly1->exps[j] = fmpz_get_si(poly1->exps +j);
            }
         }       
      }
      ppoly1->length = l;
       
      l = 0;
      //convert our poly2 to an nmod_mpoly
      for(j = 0;   j<poly2->length; j++){
          coeffTemp = fmpz_get_si(poly2->coeffs + j);
          if(coeffTemp < 0){
             while(coeffTemp <0)
                coeffTemp = coeffTemp + (long)primes[k];
             if(coeffTemp != 0){ 
                ppoly2->coeffs[j] = (ulong) coeffTemp;
                l++;
                ppoly1->exps[j] = fmpz_get_si(poly2->exps + j);
             }
         }
          else{
             NMOD_RED(coeffTemp, coeffTemp, mod);
             if(coeffTemp != 0){
                ppoly2->coeffs[j] = (ulong) coeffTemp;
                l++;
                ppoly2->exps[j] = fmpz_get_si(poly2->exps + j);
            }
         }        
      }
      ppoly2->length = l;
      
      //perform the multiplication and store the resulting polynomial
      nmod_mpoly_mul_heap(result_nmod_mpolys[k], ppoly1, ppoly2);
      
   }
   
   //clear away our ppolys now that we are finished with them.
   nmod_mpoly_clear(ppoly1);
   nmod_mpoly_clear(ppoly2);
   /*
   for(k = 0; k<num_of_primes; k++){
      printf("\nnmod_mpoly %ld", k);
      for(j=0; j<result_nmod_mpolys[k]->length; j++){
          printf("\n%lu  %lu", result_nmod_mpolys[k]->coeffs[j], result_nmod_mpolys[k]->exps[j]);
      }

   }
   */
   //now allocate an array on longs to find the least exponent in each
   //nmod_mpoly
   long * mark = (long *)malloc(sizeof(long) * num_of_primes);
   ulong nextExponent = result_nmod_mpolys[0]->exps[0];
   ulong prevExponent = 0;

   
   for(k = 0; k< num_of_primes; k++){
      mark[k] = 0;
      if(result_nmod_mpolys[k]->exps[mark[k]]< nextExponent){
         nextExponent = result_nmod_mpolys[k]->exps[mark[k]];
      }    
   }
   //printf("\nnextExponent = %lu", nextExponent);
   int fin = 1;
   
   //create an array of ulongs to perform CRT on
   ulong * residues = (ulong *)malloc(sizeof(ulong)*num_of_primes);
   l = 0;
   fmpz_init(residue);
   
   while(fin){
      fin = 0;
      
      //take the next exponent into the residues and
      //move the corresponding pointers on one.
      for(k=0;k<num_of_primes;k++){
         
         if(mark[k] < result_nmod_mpolys[k]->length){
            fin = 1;
            residues[k] = nmod_mpoly_get_coeff(result_nmod_mpolys[k], nextExponent);
            //move the pointing thing on if there was an entry there.
            if(residues[k] !=0){
            mark[k]++;
            }
         }            
      } 

      //Now use CRT on the residues and add to the result polynomial




     fmpz_set_ui(residue, residues[0]);
     fmpz_set_ui(m1, primes[0]);

     for(k = 1; k<num_of_primes; k++){
                 
         pre = n_precompute_inverse(primes[k]);

         fmpz_init(out);
         fmpz_init(inv);
         fmpz_init(mm2);
         fmpz_set_ui(mm2, primes[k]);
         fmpz_invmod(inv, m1, mm2);
              
         fmpz_CRT_ui_precomp(out, residue, m1, residues[k], primes[k], cInverses[k-1], pre );  
         fmpz_set(temp, residue);
         fmpz_set(residue, out);
            
         fmpz_clear(out);
         fmpz_clear(inv);
         fmpz_clear(mm2);              
         
         fmpz_mul_ui(m1, m1, primes[k]);
                
      }
      fmpz_set_ui(temp1, (ulong) 2);
      //fmpz_divexact_ui(m1, m1, primes[numOfPrimes -1]);
      fmpz_fdiv_q(temp, m1, temp1);
      
      if(fmpz_cmpabs(residue, temp)>0){
         fmpz_sub(residue, residue, m1);
      }

     //need to allocate space for coefficient
     
     if(prevExponent != nextExponent || prevExponent == 0){
        //printf("\n%lu", l);
        l++;
        fmpz_mpoly_fit_length(res, l);
        fmpz_set(res->coeffs + (l-1), residue);
        //printf("\n");
        //fmpz_print(residue);
        fmpz_set_ui(res->exps + (l-1), nextExponent);   
        //printf(" %lu", nextExponent);  
     }
     prevExponent = nextExponent;
      
      //Now find the next Exponent
      for(k = 0; k< num_of_primes; k++){
         if(mark[k] < result_nmod_mpolys[k]->length){
            if((result_nmod_mpolys[k]->exps[mark[k]]< nextExponent) || (prevExponent == nextExponent)){
               nextExponent = result_nmod_mpolys[k]->exps[mark[k]];
            }
         }
      }
      
   } 
   
   res->length = l;

   fmpz_clear(m1);
   fmpz_clear(temp);
   fmpz_clear(temp1);
   fmpz_clear(temp2);
   fmpz_clear(residue);  

   free(mark);
   free(primes);
   free(residues);
}
