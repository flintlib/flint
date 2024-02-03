/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

int _mpoly_heap_insert(mpoly_heap_s * heap, ulong * exp, void * x,
       slong * next_loc, slong * heap_len, slong N, const ulong * cmpmask)
{
   slong i = *heap_len, j, n = *heap_len;

   if (i != 1 && mpoly_monomial_equal(exp, heap[1].exp, N))
   {
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[1].next;
      heap[1].next = x;

      return 0;
   }

   if (*next_loc < *heap_len)
   {
      if (mpoly_monomial_equal(exp, heap[*next_loc].exp, N))
      {
         ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[*next_loc].next;
         heap[*next_loc].next = x;
         return 0;
      }
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (!mpoly_monomial_gt(exp, heap[j].exp, N, cmpmask))
         break;

      i = j;
   }

   if (j >= 1 && mpoly_monomial_equal(exp, heap[j].exp, N))
   {
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[j].next;
      heap[j].next = x;
      *next_loc = j;

      return 0;
   }

   (*heap_len)++;

   while (n > i)
   {
      heap[n] = heap[HEAP_PARENT(n)];
      n = HEAP_PARENT(n);
   }

   HEAP_ASSIGN(heap[i], exp, x);

   return 1;
}
