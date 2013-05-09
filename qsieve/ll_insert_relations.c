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
/******************************************************************************

 linear_algebra.c
 
 Routines for dealing with building the final F_2 matrix

 (C) 2006, 2011 William Hart

******************************************************************************/

#undef ulong /* avoid clash with stdlib */
#include <stdlib.h>
#include <stdio.h>
#define ulong unsigned long 

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"

/*=========================================================================
   
   Compare relations:
 
   Function: Compare two relations; used by qsort
 
==========================================================================*/

int qsieve_ll_relations_cmp(const void * a, const void * b)
{
  la_col_t * ra = *((la_col_t **) a);
  la_col_t * rb = *((la_col_t **) b);
  len_t point;

  if (ra->weight > rb->weight) return 1;
  else if (ra->weight < rb->weight) return -1;
  
  for (point = ra->weight - 1; point >= 0 
                            && (ra->data[point] == rb->data[point]); point--)
  {
      ;
  }
  
  if (point == -1) return 0;

  if (ra->data[point] > rb->data[point]) return 1;
  else if (ra->data[point] < rb->data[point]) return -1;

  return 0;
}

int qsieve_ll_relations_cmp2(const void * a, const void * b)
{
  la_col_t * ra = (la_col_t *) a;
  la_col_t * rb = (la_col_t *) b;
  len_t point;

  if (ra->weight > rb->weight) return 1;
  else if (ra->weight < rb->weight) return -1;
  
  for (point = ra->weight - 1; point >= 0 
                            && (ra->data[point] == rb->data[point]); point--)
  {
      ;
  }
  
  if (point == -1) return 0;

  if (ra->data[point] > rb->data[point]) return 1;
  else if (ra->data[point] < rb->data[point]) return -1;

  return 0;
}
  
/*==========================================================================
   Merge sort:

   Function: Merge a list of sorted new relations into a list of existing
             sorted relations. Sort is done using a merge sort algorithm
             with a short stack.
   
===========================================================================*/

len_t qsieve_ll_merge_sort(qs_t qs_inf)
{
   la_col_t * matrix = qs_inf->matrix;
   len_t columns = qs_inf->columns;
   la_col_t ** qsort_arr = qs_inf->qsort_arr;
   len_t num_unmerged = qs_inf->num_unmerged;
   len_t dups = 0;
   int comp;
   len_t i;

   for (i = columns + num_unmerged - 1L; i >= dups; i--) 
   {
      if (!columns) comp = -1;
      else if (!num_unmerged) comp = 1;
      else 
         comp = qsieve_ll_relations_cmp2(matrix + columns - 1L, qsort_arr[num_unmerged - 1L]);
      
      switch (comp)
      {
         case -1: 
         {
            copy_col(matrix + i, qsort_arr[num_unmerged - 1L]);
            clear_col(qsort_arr[num_unmerged - 1L]);
            num_unmerged--;
            break;
         }
         case 1: 
         {
            copy_col(matrix + i, matrix + columns - 1L);
            columns--;
            break;
         }
         case 0: 
         {
            free_col(qsort_arr[num_unmerged - 1L]);
            clear_col(qsort_arr[num_unmerged - 1L]);
            num_unmerged--;
            copy_col(matrix + i, matrix + columns - 1L);
            columns--;
            dups++;
            break;
         }
      } 
   }
   
   columns = qs_inf->columns + qs_inf->num_unmerged - dups;
   
   if (dups)
   {
      len_t i;
      for (i = 0; i < columns; i++)
          copy_col(matrix + i, matrix + i + dups);
      
      for ( ; i < columns + dups; i++)
          clear_col(matrix + i);
   }
            
   qs_inf->columns = columns;
   columns = qs_inf->num_unmerged - dups;
   qs_inf->num_unmerged = 0;

#if (QS_DEBUG & 64)
   printf("%ld new, %ld dups\n", columns, dups);
#endif   

   return columns;
}

/*==========================================================================
   Merge relations:

   Function: Merge unmerged relations into the matrix
   
===========================================================================*/

len_t qsieve_ll_merge_relations(qs_t qs_inf)
{
   const len_t num_unmerged = qs_inf->num_unmerged;
   la_col_t * unmerged = qs_inf->unmerged;
   la_col_t ** qsort_arr = qs_inf->qsort_arr;
   
   if (num_unmerged)
   {
      len_t i;

      for (i = 0; i < num_unmerged; i++)
         qsort_arr[i] = unmerged + i;
      
      qsort(qsort_arr, num_unmerged, sizeof(la_col_t *), qsieve_ll_relations_cmp);

      return qsieve_ll_merge_sort(qs_inf);
   }
   
   return 0;
}

/*==========================================================================
   Insert relation:

   Function: Insert the relation into the matrix and store the Y value
   
===========================================================================*/

len_t qsieve_ll_insert_relation(qs_t qs_inf, fmpz_t Y)
{
   la_col_t * unmerged = qs_inf->unmerged;
   len_t num_unmerged = qs_inf->num_unmerged;
   len_t * small = qs_inf->small;
   len_t num_factors = qs_inf->num_factors; 
   fac_t * factor = qs_inf->factor; 
   len_t * curr_rel = qs_inf->curr_rel;
   len_t fac_num = 0; 
   len_t i;

   clear_col(unmerged + num_unmerged);
   
   for (i = 0; i < qs_inf->small_primes; i++)
   {
       if (small[i] & 1) insert_col_entry(unmerged + num_unmerged, i);
       if (small[i]) 
       {
          curr_rel[2*fac_num + 1] = i;
          curr_rel[2*fac_num + 2] = small[i];
          fac_num++;
       }
   }
   
   for (i = 0; i < num_factors; i++)
   {
       if (factor[i].exp & 1) insert_col_entry(unmerged + num_unmerged, factor[i].ind);
       curr_rel[2*fac_num + 1] = factor[i].ind;
       curr_rel[2*fac_num + 2] = factor[i].exp;
       fac_num++;
   }

   curr_rel[0] = fac_num;
   
   unmerged[num_unmerged].orig = qs_inf->num_relations;
   
   fmpz_set(qs_inf->Y_arr + qs_inf->num_relations, Y); 
   
   qs_inf->curr_rel += qs_inf->max_factors*2;
   qs_inf->num_unmerged++;
   qs_inf->num_relations++;
   
   if (qs_inf->num_unmerged == qs_inf->qsort_rels)
      return qsieve_ll_merge_relations(qs_inf);
   
   return 0;
}

