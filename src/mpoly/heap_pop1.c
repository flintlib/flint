/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void * _mpoly_heap_pop1(mpoly_heap1_s * heap, slong * heap_len, ulong maskhi)
{
   ulong exp;
   slong i, j, s = --(*heap_len);
   void * x = heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if ((heap[j].exp^maskhi) <= (heap[j + 1].exp^maskhi))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   exp = heap[s].exp;
   j = HEAP_PARENT(i);

   while (i > 1 && (exp^maskhi) > (heap[j].exp^maskhi))
   {
     heap[i] = heap[j];
     i = j;
     j = HEAP_PARENT(j);
   }

   heap[i] = heap[s];

   return x;
}
