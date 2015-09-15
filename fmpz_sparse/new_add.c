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

    Copyright (C) 2008, 2009 William Hart

******************************************************************************/

#include "fmpz_sparse.h"

void 
_fmpz_sparse_new_add(fmpz * res_c, fmpz * res_e, slong * res_len, const fmpz * poly1_c, 
    const fmpz * poly1_e, slong len1, const fmpz * poly2_c, const fmpz * poly2_e,
    slong len2)
{
    slong i = 0, j = 0, k = 0;

    while(i < len1 && j < len2)
    {
      fmpz_init(res_c + k);
      fmpz_init(res_e + k);

      if((poly1_e + i) > (poly2_e + j))
      {
        fmpz_set(res_c + k, poly1_c + i);
        fmpz_set(res_e + k, poly1_e + i);
        i++;
      }
      else if((poly1_e + i) < (poly2_e + j))
      {
        fmpz_set(res_c + k, poly2_c + j);
        fmpz_set(res_e + k, poly2_e + j);
        j++;
      }
      else
      {
        fmpz_add(res_c + k, poly1_c + i, poly2_c + j); 
        fmpz_set(res_e + k, poly1_e + i);
        i++;
        j++;

        if(fmpz_is_zero(res_c + k))
        {
          k--;
        }
      }

      k++;
    }

    if(i < len1)
      for(; i < len1; i++)
      {
        fmpz_set(res_c + k, poly1_c + i);
        fmpz_set(res_e + k, poly1_e + i);
        k++;
      }
    
    if(j < len2)
      for(; j < len2; j++)
      {
        fmpz_set(res_c + k, poly2_c + j);
        fmpz_set(res_e + k, poly2_e + j);
        k++;
      }

    *res_len = k;
}

void
fmpz_sparse_new_add(fmpz_sparse_t res, const fmpz_sparse_t poly1,
    const fmpz_sparse_t poly2)
{
  slong max_length = poly1->length + poly2->length;
  if (poly1 == res) 
  {
    fmpz_sparse_t temp;
    fmpz_sparse_init2(temp, max_length);

    _fmpz_sparse_new_add(temp->coeffs, temp->expons, &temp->length, poly1->coeffs,
        poly1->expons, poly1->length, poly2->coeffs, poly2->expons, poly2->length);

    fmpz_sparse_swap(res, temp);
    fmpz_sparse_clear(temp);
  }
  else if (poly2 == res)
  {
    fmpz_sparse_t temp;
    fmpz_sparse_init2(temp, max_length);

    _fmpz_sparse_new_add(temp->coeffs, temp->expons, &temp->length, poly1->coeffs,
         poly1->expons, poly1->length, poly2->coeffs, poly2->expons, poly2->length);
    
    fmpz_sparse_swap(res, temp);
    fmpz_sparse_clear(temp);
  }
  /*else if (poly1 == poly2) fmpz_sparse_scalar_mul(res, poly1, (fmpz) 2);*/
  else
  {
    fmpz_sparse_init2(res, max_length);
    
    _fmpz_sparse_new_add(res->coeffs, res->expons, &res->length, poly1->coeffs, 
        poly1->expons, poly1->length, poly2->coeffs, poly2->expons, poly2->length);
  }

  _fmpz_sparse_normalise(res);
}
