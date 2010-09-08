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

   Copyright (C) 2008, 2009 William Hart
   
*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void fmpz_mpoly_set_coeff_fmpz(fmpz_mpoly_t poly, ulong exponent, const fmpz_t x){

   //ensure that there is enough space
   fmpz_mpoly_fit_length(poly, poly->length + 1);

   fmpz_t fmpz_exp;
   fmpz_init(fmpz_exp);
   fmpz_set_ui(fmpz_exp, exponent);

   //temporary variables to store
   fmpz_t tempCoefficient;
   fmpz_t tempExponent;

   fmpz_init(tempCoefficient);
   fmpz_init(tempExponent);

   int i;
   int have_inserted = 0;
   int have_replaced = 0;

   //if x does not equal zero
   if( fmpz_is_zero(x) !=1 ){
   
      for(i=0; i<poly->length; i++){

         //if we are replacing coefficients
         if(fmpz_equal(poly->exps + i, fmpz_exp)){
            fmpz_set(poly->coeffs + i, x);
            have_replaced = 1;
            have_inserted = 1;
            break;
            }  
         if(have_inserted == 1){
            fmpz_swap(tempCoefficient, poly->coeffs + i);
            fmpz_swap(tempExponent, poly->exps + i);
         }
         //if we now have to insert
         if( fmpz_cmp(poly->exps + i, fmpz_exp) > 0 && have_inserted == 0){
            //ensure that there is enough space
            fmpz_mpoly_fit_length(poly, poly->length + 1);
            fmpz_set(tempCoefficient, poly->coeffs + i);
            fmpz_set(tempExponent, poly->exps + i);
            fmpz_set( poly->coeffs + i, x);
            fmpz_set( poly->exps + i, fmpz_exp);
            have_inserted = 1;            
         }          
      }

   //tag on the last entry
   if(have_inserted == 1 && have_replaced == 0){    
      fmpz_set(poly->coeffs + poly->length, tempCoefficient);
      fmpz_set(poly->exps + poly->length, tempExponent);
      poly->length += 1;
      }
   else if(have_inserted == 0){
      fmpz_mpoly_fit_length(poly, poly->length + 1);
      fmpz_set(poly->exps + poly->length, fmpz_exp);
      fmpz_set(poly->coeffs + poly->length, x);
      poly->length += 1;
      }
 
   }


   //if x does equal zero
   if( fmpz_is_zero(x) == 1){
      have_replaced = 0;
      
      for(i=0;i<poly->length;i++){
         
         if(have_replaced == 1){
            
            fmpz_set(poly->exps+(i-1), poly->exps+i);
            fmpz_set(poly->coeffs+(i-1), poly->coeffs+i);       

         }

         if(have_replaced == 0 && fmpz_equal(poly->exps + i, fmpz_exp)){ 
            have_replaced = 1;
         }
      }
 
     if(have_replaced == 1){
        poly->length -= 1;
        _fmpz_mpoly_truncate(poly, poly->length);   
     } 
   }

   fmpz_clear(fmpz_exp);
   fmpz_clear(tempCoefficient);
   fmpz_clear(tempExponent);

}
