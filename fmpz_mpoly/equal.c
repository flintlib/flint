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
   Copyright (C) 2010, Daniel Woodhouse
   
*****************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

int fmpz_mpoly_equal(fmpz_mpoly_t poly1, fmpz_mpoly_t poly2){

     if(poly1->length != poly2->length){
        return 0;
     }
     else if(poly1->vars != poly2->vars || poly1->ebits != poly2->ebits){
        return 0;  
     }

    int i;
  
    for(i=0; i < poly1->length; i++){
       
   
       if(fmpz_equal(poly1->exps + i, poly2->exps +i) == 0 || fmpz_equal(poly1->coeffs + i, poly2->coeffs + i) == 0)
          return 0;

    }

    return 1;
}


