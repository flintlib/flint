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

******************************************************************************/

#include "flint.h"
#include "fmpq.h"
#include "fmpq_mat.h"

void fmpq_mat_similarity(fmpq_mat_t A, slong r, fmpq_t d)
{
   slong n = A->r, i, j;
   
   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
         fmpq_addmul(fmpq_mat_entry(A, i, j), fmpq_mat_entry(A, i, r), d);
      
      for (j = r + 1; j < n; j++)
         fmpq_addmul(fmpq_mat_entry(A, i, j), fmpq_mat_entry(A, i, r), d); 
   }

   for (i = 0; i < n; i++)
   {
      for (j = 0; j < r - 1; j++)
         fmpq_submul(fmpq_mat_entry(A, r, i), fmpq_mat_entry(A, j, i), d);

      for (j = r + 1; j < n; j++)
         fmpq_submul(fmpq_mat_entry(A, r, i), fmpq_mat_entry(A, j, i), d);      
   }
}
