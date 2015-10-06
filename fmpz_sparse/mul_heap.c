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

    Authored 2015 by A. Whitman Groves; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"
#include "fmpz_vec.h"

typedef struct heap
{
  fmpz * expons;
  fmpz * coeffs;
  slong * p1;
  slong * p2;
  slong size;
  slong alloc;
} heap;

typedef struct heap heap_t[1];

void
fmpz_heap_init(heap * h, const slong size)
{
  h->expons = flint_calloc(size, sizeof(fmpz));
  h->coeffs = flint_calloc(size, sizeof(fmpz));
  h->p1 = flint_calloc(size, sizeof(slong));
  h->p2 = flint_calloc(size, sizeof(slong));
  h->size = 0;
  h->alloc = size;
}

void
fmpz_heap_free(heap_t heap)
{
  slong i;
  for(i = 0; i < heap->alloc; i++)
  {
    fmpz_clear(heap->coeffs + i);
    fmpz_clear(heap->expons + i);
  }
  flint_free(heap->expons);
  flint_free(heap->coeffs);
  flint_free(heap->p1);
  flint_free(heap->p2);
}

void
fmpz_heap_insert(heap_t heap, const fmpz_t expon, const fmpz_t coeff,
    const slong term1, const slong term2)
{

  if(heap->size == 0)
  {
    fmpz_set(heap->expons + 1, expon);
    fmpz_set(heap->coeffs + 1, coeff);
    heap->p1[1] = term1;
    heap->p2[1] = term2;
  }
  else
  {
    slong place = heap->size + 1;

    while(place != 1)
    {
      if(fmpz_cmp(heap->expons + place/2, expon) < 0)
      {
        fmpz_set(heap->expons + place, heap->expons + place/2);
        fmpz_set(heap->coeffs + place, heap->coeffs + place/2);
        heap->p1[place] = heap->p1[place/2];
        heap->p2[place] = heap->p2[place/2];
        place = place/2;
      }
      else
        break;
    }

    fmpz_set(heap->expons + place, expon);
    fmpz_set(heap->coeffs + place, coeff);
    heap->p1[place] = term1;
    heap->p2[place] = term2;
  }
  heap->size += 1;
}

void
fmpz_heap_pop(heap_t heap)
{
  slong place = 1;
  slong temp1 = heap->p1[heap->size];
  slong temp2 = heap->p2[heap->size];
  slong temp3;
  fmpz_t temp_e, temp_c;
  
  fmpz_init(temp_e);
  fmpz_init(temp_c);
  
  if(heap->size > 1)
  {
    fmpz_set(temp_e, heap->expons + heap->size);
    fmpz_set(temp_c, heap->coeffs + heap->size);
  }
  
  fmpz_clear(heap->expons + heap->size);
  fmpz_clear(heap->coeffs + heap->size);
  
  heap->size -= 1;
  while(place*2 <= heap->size)
  {
    if(place*2 + 1 <= heap->size)
      temp3 = place*2 + 1;
    else
      temp3 = place*2;
    
    if(fmpz_cmp(heap->expons + place*2, heap->expons + temp3) >= 0)
      temp3 = place*2;
    else
      temp3 = place*2 + 1;
    
    if(fmpz_cmp(temp_e, heap->expons + temp3) < 0)
    {
      fmpz_set(heap->expons + place, heap->expons + temp3);
      fmpz_set(heap->coeffs + place, heap->coeffs + temp3);
      heap->p1[place] = heap->p1[temp3];
      heap->p2[place] = heap->p2[temp3];
      place = temp3;
    }
    else
      break;
  }
  
  fmpz_set(heap->expons + place, temp_e);
  fmpz_set(heap->coeffs + place, temp_c);
  heap->p1[place] = temp1;
  heap->p2[place] = temp2;
  
  fmpz_clear(temp_c);
  fmpz_clear(temp_e);
}

void 
fmpz_sparse_mul_heaps(fmpz_sparse_t res, const fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2)
{
  slong i, j;
  if (res == poly1 || res == poly2)
  {
    fmpz_sparse_t temp;
    fmpz_sparse_init(temp);
    fmpz_sparse_mul_heaps(temp, poly1, poly2);
    fmpz_sparse_set(res, temp);
    fmpz_sparse_clear(temp);
    return;
  }
  
  if((poly1->length == 0) || (poly2->length == 0))
  {
    fmpz_sparse_zero(res);
    return;
  }

  if((poly1->length == 1) && (poly2->length == 1))
  {
    res->length = 1;
    fmpz_init(res->coeffs);
    fmpz_init(res->expons);
    fmpz_mul(res->coeffs, poly1->coeffs, poly2->coeffs);
    fmpz_add(res->expons, poly1->expons, poly2->expons);
    return;
  }

  if(poly2->length < poly1->length)
  {
    fmpz_sparse_mul_heaps(res, poly2, poly1);
  }
  else
  {
    heap_t heap;

    fmpz_t temp_e, temp_c;
   
    slong k;

    _fmpz_sparse_reserve(res, poly1->length*poly2->length);

    fmpz_heap_init(heap, poly1->length + 1);
   

    for(i = 0; i < poly1->length; ++i)
    {
      fmpz_init(temp_e);
      fmpz_init(temp_c);
      fmpz_mul(temp_c, poly1->coeffs + i, poly2->coeffs);
      fmpz_add(temp_e, poly1->expons + i, poly2->expons);

      fmpz_heap_insert(heap, temp_e, temp_c, i, 0);
      fmpz_clear(temp_e);
      fmpz_clear(temp_c);
    }

    k = 0;

    while(heap->size > 0)
    {
      fmpz_init(temp_e);
      fmpz_init(temp_c);
      i = heap->p1[1];
      j = heap->p2[1];

      if(k == 0)
      {
        fmpz_init(res->expons + k);
        fmpz_init(res->coeffs + k);
        fmpz_set(res->expons + k, heap->expons + 1);
        fmpz_set(res->coeffs + k, heap->coeffs + 1);
        k++;
      }
      else if(fmpz_is_zero(res->coeffs + k - 1))
      {
        fmpz_set(res->expons + k - 1, heap->expons + 1);
        fmpz_set(res->coeffs + k - 1, heap->coeffs + 1);
      }
      else if(fmpz_equal(res->expons + k - 1, heap->expons + 1) == 1)
      {
        fmpz_add(res->coeffs + k - 1, heap->coeffs + 1, res->coeffs + k - 1);
      }
      else
      {
        fmpz_init(res->expons + k);
        fmpz_init(res->coeffs + k);
        fmpz_set(res->expons + k, heap->expons + 1);
        fmpz_set(res->coeffs + k, heap->coeffs + 1);
        k++;
      }

      fmpz_heap_pop(heap);

      if(j < poly2->length - 1)
      {
        fmpz_mul(temp_c, poly1->coeffs + i, poly2->coeffs + j + 1);
        fmpz_add(temp_e, poly1->expons + i, poly2->expons + j + 1);
        fmpz_heap_insert(heap, temp_e, temp_c, i, j + 1);
      }

      fmpz_clear(temp_e);
      fmpz_clear(temp_c);
    }
   
    res->length = k;

    fmpz_heap_free(heap);
  }
}
