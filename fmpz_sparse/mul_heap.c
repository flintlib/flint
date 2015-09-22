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

    Authored 2015 by Daniel S. Roche; US Government work in the public domain. 

******************************************************************************/

#include "fmpz_sparse.h"
#include "fmpz_vec.h"

typedef struct heap
{
  fmpz_t * expons;
  fmpz_t * coeffs;
  slong size;
} heap;

typedef struct heap heap_t[1];

void
fmpz_heap_init(heap * h, const slong size)
{
  h->expons = flint_calloc(size, sizeof(fmpz));
  h->coeffs = flint_calloc(size, sizeof(fmpz));
  h->size = 0;
}

void
fmpz_heap_insert(heap_t heap, const fmpz_t expon, const fmpz_t coeff)
{
  if(heap->size == 0)
  {
    fmpz_set(heap->expons[1], expon);
    fmpz_set(heap->coeffs[1], coeff);
  }
  else
  {
    slong place = heap->size+1;
 
    fmpz_set(heap->expons[place], expon);
    fmpz_set(heap->coeffs[place], coeff);
    
    while(fmpz_cmp(heap->expons[place/2], expon) < 0)
    {
      fmpz_swap(heap->expons[place/2], heap->expons[place]);
      fmpz_swap(heap->coeffs[place/2], heap->coeffs[place]);
      place = place/2;
    }
  }
  heap->size += 1;
}

void
fmpz_heap_pop(fmpz_sparse_t res, heap_t heap)
{
  slong place = 1, k = res->length;

  if(res->coeffs + k - 1 == 0)
  {
    k--;
    fmpz_set(res->expons + k, heap->expons[place]);
    fmpz_set(res->coeffs + k, heap->coeffs[place]);
    k++;
  }
  else if(fmpz_equal(res->expons + k - 1, heap->expons[1]) == 0)
  {
    k--;
    fmpz_add(res->coeffs + k, heap->coeffs[place], res->coeffs + k);
    k++;
  }
  else
  {
    fmpz_init(res->expons + k);
    fmpz_init(res->coeffs + k);
    fmpz_set(res->expons + k, heap->expons[place]);
    fmpz_set(res->coeffs + k, heap->coeffs[place]);
    k++;
  }

 
  res->length = k;

  fmpz_swap(heap->expons[place], heap->expons[heap->size]);
  fmpz_swap(heap->coeffs[place], heap->coeffs[heap->size]);
  
  while((fmpz_cmp(heap->expons[place], heap->expons[place*2]) < 0) || 
      (fmpz_cmp(heap->expons[place], heap->expons[place*2+1]) < 0))
  {
    if(fmpz_cmp(heap->expons[place*2], heap->expons[place*2+1]) > 0)
    {
      fmpz_swap(heap->expons[place], heap->expons[place*2]);
      fmpz_swap(heap->expons[place], heap->expons[place*2]);
      place = place*2;
    }
    else
    {
      fmpz_swap(heap->expons[place], heap->expons[place*2+1]);
      fmpz_swap(heap->expons[place], heap->expons[place*2+1]);
      place = place*2 + 1;
    }
  }
  heap->size -= 1;
}

void 
fmpz_sparse_mul_heaps(fmpz_sparse_t res, const fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2)
{
  slong i, j, max = poly1->length * poly2->length;
  if((poly1->length == 0) || (poly2->length == 0))
  {
    fmpz_sparse_zero(res);
    return;
  }

  if((poly1->length == 1) && (poly2->length == 1))
  {
    res->length = 1;
    fmpz_mul(res->coeffs, poly1->coeffs, poly2->coeffs);
    fmpz_add(res->expons, poly1->expons, poly2->expons);
    return;
  }

  if (res == poly1 || res == poly2) 
  {
    fmpz_sparse_t temp;
    fmpz_sparse_init2(temp, max);
    fmpz_sparse_mul_heaps(temp, poly1, poly2);
    fmpz_sparse_init2(res, temp->length);
    fmpz_sparse_set(res, temp);
    fmpz_sparse_clear(temp);
  }
  else if(poly2->length < poly1->length)
  {
    fmpz_sparse_mul_heaps(res, poly2, poly1);
  }
  else
  {
    heap_t heap;

    fmpz_sparse_print(poly1), flint_printf("\n\n");
    fmpz_sparse_print(poly2), flint_printf("\n\n");
    fmpz_heap_init(heap, poly1->length + 1);

    fmpz_t temp_e, temp_c;
    
    fmpz_init(temp_e);
    fmpz_init(temp_c);
    
    for(j = 0; j < poly2->length; ++j)
    {
      for(i = 0; i < poly1->length; ++i)
      {
        fmpz_mul(temp_c, poly1->coeffs + i, poly2->coeffs + j);
        fmpz_add(temp_e, poly1->expons + i, poly2->expons + j);
      
        if(heap->size == poly1->length)
          fmpz_heap_pop(res, heap);
        
        fmpz_heap_insert(heap, temp_e, temp_c);
      }
    }
   
    fmpz_sparse_print(res), flint_printf("\n\n");
    
    fmpz_clear(temp_e);
    fmpz_clear(temp_c);
  }
}
