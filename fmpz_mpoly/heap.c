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

    Copyright (C) 2010 William Hart

    This implementation was inspired by the reference implementation of
    Roman Pearce.

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

void
fmpz_mpoly_reheapify(fmpz_mpoly_heap_t * heap, ulong * n)
{
    ulong i = 1, j = 2;
    const ulong _n = (*n);
    const fmpz top_exp = heap[_n].exp;;

    while (j < _n)
    {
        if (heap[j].exp > heap[j + 1].exp)
            j++;                /* move right if smaller */

        if (heap[j].exp < top_exp)
        {
            heap[i] = heap[j];  /* move heap entry up into empty spot */
            i = j;
            j = 2 * i;          /* move down */
        }
        else
            break;              /* quit if heap bottom can go here */
    }

    heap[i] = heap[_n];         /* move heap bottom into empty spot */
    (*n)--;                     /* heap is now one smaller */
}

void
fmpz_mpoly_heap_insert(fmpz_mpoly_heap_t * heap, ulong * n,
                       fmpz_mpoly_entry_t * entry, fmpz exp)
{
    if ((*n) != 0)              /* make sure the heap has something in it */
    {
        ulong i, j, s;

        /* 
           first see if the new entry can be chained at the top 
           this happens often enough to optimise this case
         */
        if (exp == heap[1].exp)
        {
            entry->next = heap[1].entry;
            heap[1].entry = entry;
            return;
        }

        i = (*n) + 1;
        j = i / 2;

        /* find where to put the new entry */
        do
        {
            const mp_limb_signed_t diff = exp - heap[j].exp;

            if (diff < (mp_limb_signed_t) 0)
                j /= 2;
            else if (diff == (mp_limb_signed_t) 0)  /* exps are equal, chain here */
            {
                entry->next = heap[j].entry;
                heap[j].entry = entry;

                return;
            }
            else
                break;
        } while (j > 0);

        /* 
           s is the slot we want to put our entry _below_
           move entries down to create a space
         */
        s = j;

        while (i / 2 > s)
        {
            heap[i] = heap[i / 2];
            i /= 2;
        }

        /* place entry in spot */
        entry->next = NULL;
        heap[i].exp = exp;
        heap[i].entry = entry;
        (*n)++;

    }
    else                        /* heap is empty, insert entry straight in */
    {
        entry->next = NULL;
        heap[1].exp = exp;
        heap[1].entry = entry;
        (*n) = 1;
    }
}
