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

typedef struct node
{
  fmpz_t expons;
  fmpz_t coeffs;
  slong p1;
  slong p2;
} node;

typedef struct node node_t[1];

typedef struct new_heap
{
  node * nodes;
  slong size;
  slong alloc;
} heap;

typedef struct new_heap new_heap_t[1];

void
fmpz_node_init(node * n)
{
  fmpz_init(n->expons);
  fmpz_init(n->coeffs);
}

void
fmpz_new_heap_init(new_heap_t  h, const slong size)
{
  int i;
  h->nodes = flint_calloc(size, sizeof(node_t));
  for(i = 0; i < size; i++)
  {
    fmpz_node_init(h->nodes + i);
  }
  h->size = 0;
  h->alloc = size;
}

void
fmpz_node_free(node * n)
{
  fmpz_clear(n->expons);
  fmpz_clear(n->coeffs);
}

void
fmpz_new_heap_free(new_heap_t heap)
{
  slong i;
  for(i = 0; i < heap->alloc; i++)
  {
    fmpz_node_free(heap->nodes + i);
  }
}

void
fmpz_node_swap(node_t n1, node_t n2)
{
  node_t temp;
  *temp = *n1;
  *n1 = *n2;
  *n2 = *temp;
  fmpz_node_free(temp);
}

void
fmpz_heap_print(new_heap_t h)
{
  slong i;
  flint_printf("\n");
  for(i = 0; i < h->alloc - 1; ++i)
  {
    flint_printf("(");
    fmpz_print((h->nodes + i)->expons);
    flint_printf(")");
  }
  flint_printf("\n");
  for(i = 0; i < h->alloc - 1; ++i)
  {
    flint_printf("(%lu,%lu)", (h->nodes + i)->p1, (h->nodes + i)->p2);
  }
  flint_printf("\n");
}

void
fmpz_heap_replace(new_heap_t heap)
{
  slong place = 0;
  slong temp3;

  if(heap->size == 1)
    return;
 
  while(place*2 + 1 <= heap->size - 1)
  {
    if(place*2 + 2 <= heap->size - 1)
      temp3 = place*2 + 2;
    else
      temp3 = place*2 + 1;
    
    if(fmpz_cmp((heap->nodes + place*2 + 1)->expons, (heap->nodes + temp3)->expons) >= 0)
      temp3 = place*2 + 1;
    else
      temp3 = place*2 + 2;
    
    if(fmpz_cmp((heap->nodes + place)->expons, (heap->nodes + temp3)->expons) < 0)
    {
      fmpz_node_swap(heap->nodes + place, heap->nodes + temp3);
      place = temp3;
    }
    else
      break;
  }
}

void
fmpz_new_heap_pop(new_heap_t heap)
{
  slong place = 0;
  slong temp3;
 
  fmpz_node_swap(heap->nodes + place, heap->nodes + heap->size - 1);

  heap->size -= 1;
  while(place*2 + 1 <= heap->size - 1)
  {
    if(place*2 + 2 <= heap->size - 1)
      temp3 = place*2 + 2;
    else
      temp3 = place*2 + 1;
    
    if(fmpz_cmp((heap->nodes + place*2 + 1)->expons, (heap->nodes + temp3)->expons) >= 0)
      temp3 = place*2 + 1;
    else
      temp3 = place*2 + 2;
    
    if(fmpz_cmp((heap->nodes + place)->expons, (heap->nodes + temp3)->expons) < 0)
    {
      fmpz_node_swap(heap->nodes + place, heap->nodes + temp3);
      place = temp3;
    }
    else
      break;
  }
 
}

void
fmpz_sparse_new_heaps(fmpz_sparse_t res, const fmpz_sparse_t poly1, const fmpz_sparse_t poly2)
{
  if(res == poly1 || res == poly2)
  {
    fmpz_sparse_t temp;
    fmpz_sparse_init(temp);
    fmpz_sparse_new_heaps(temp, poly1, poly2);
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
    fmpz_sparse_new_heaps(res, poly2, poly1);
  }
  else
  {
    new_heap_t heap;
    
    slong k;
    
    slong i;
    
    _fmpz_sparse_reserve(res, poly1->length*poly2->length);
    
    fmpz_new_heap_init(heap, poly1->length + 1);
    
    
    for(i = 0; i < poly1->length; ++i)
    {
      fmpz_mul((heap->nodes + i)->coeffs, poly1->coeffs + i, poly2->coeffs);
      fmpz_add((heap->nodes + i)->expons, poly1->expons + i, poly2->expons);
      (heap->nodes + i)->p1 = i;
      (heap->nodes + i)->p2 = 0;
      heap->size += 1;
    }
    fmpz_heap_print(heap);
    flint_printf("\nORIGINAL HEAP SIZE: %lu", heap->size); 
    k = 0;
    
    while(heap->size > 0)
    {
      if(k == 0)
      {
        flint_printf("\nbegin number: (");
        fmpz_print((heap->nodes)->coeffs);
        flint_printf(",");
        fmpz_print((heap->nodes)->expons);
        flint_printf(") (%lu, %lu) size: %lu\n", (heap->nodes)->p1, (heap->nodes)->p2, heap->size);
        fmpz_init(res->expons);
        fmpz_init(res->coeffs);
        fmpz_set(res->expons, (heap->nodes)->expons);
        fmpz_set(res->coeffs, (heap->nodes)->coeffs);
        k++;
      }
      else if(fmpz_is_zero(res->coeffs + k - 1))
      {
        flint_printf("\ncancel zero: (");
        fmpz_print((heap->nodes)->coeffs);
        flint_printf(",");
        fmpz_print((heap->nodes)->expons);
        flint_printf(") (%lu, %lu) size: %lu\n", (heap->nodes)->p1, (heap->nodes)->p2, heap->size);
        fmpz_set(res->expons + k - 1, (heap->nodes)->expons);
        fmpz_set(res->coeffs + k - 1, (heap->nodes)->coeffs);
      }
      else if(fmpz_equal(res->expons + k - 1, (heap->nodes)->expons))
      {
        flint_printf("\nappend coeffs: (");
        fmpz_print((heap->nodes)->coeffs);
        flint_printf(",");
        fmpz_print((heap->nodes)->expons);
        flint_printf(") (%lu, %lu) size: %lu\n", (heap->nodes)->p1, (heap->nodes)->p2, heap->size);
        fmpz_add(res->coeffs + k - 1, (heap->nodes)->coeffs, res->coeffs + k - 1);
      }
      else
      {
        flint_printf("\nadd length: (");
        fmpz_print((heap->nodes)->coeffs);
        flint_printf(",");
        fmpz_print((heap->nodes)->expons);
        flint_printf(") (%lu, %lu) size: %lu\n", (heap->nodes)->p1, (heap->nodes)->p2, heap->size);
        fmpz_init(res->expons + k);
        fmpz_init(res->coeffs + k);
        fmpz_set(res->expons + k, (heap->nodes)->expons);
        fmpz_set(res->coeffs + k, (heap->nodes)->coeffs);
        k++;
      }
 
      /*if((heap->nodes)->p2 == poly2->length - 1)
      {*/
        flint_printf("\nPOP (%lu, %lu)\n", (heap->nodes)->p1, (heap->nodes)->p2);
        fmpz_new_heap_pop(heap);
        fmpz_heap_print(heap);
      /*}
      else
      {
        (heap->nodes)->p2 += 1;
        fmpz_mul((heap->nodes)->coeffs, poly1->coeffs + (heap->nodes)->p1, poly2->coeffs + (heap->nodes)->p2);
        fmpz_add((heap->nodes)->expons, poly1->expons + (heap->nodes)->p1, poly2->expons + (heap->nodes)->p2);
        fmpz_heap_replace(heap);
      }*/
    }
    
    res->length = k;
    
    fmpz_new_heap_free(heap);
  }
}
