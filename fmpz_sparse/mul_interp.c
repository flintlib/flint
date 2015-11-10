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
void 
fmpz_sparse_mul_interp(fmpz_sparse_t res, flint_rand_t state, const fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2)
{
  /*
   * (substitutions for kronecker substitution are ignored for this implementation)
   * make copies for f, g, f_k, g_k, f_s, and g_s
   * -----------f_s and g_s have coefficients of 1
   * estimate structural sparsity
   * -----------choose primes p and p' to modulate the polynomials by ((f_s*g_S)^mod p)mod p' until half dense
   * compute structural support
   * -----------choose same prime p and make h_1 = ((f_s*g_s)^mod p)mod p
   * -----------f_2 = the sum of (e*L + 1)*z^(e mod p)
   * -----------g_2 = the sum of (e*L + 1)*z^(e mod p)
   * -----------h_2 = ((f_2*g_2)^mod p)mod L^2
   * -----------compute exponents of S with ratio ((c_2/c_1)-1)/L
   * compute arithmetic support
   * -----------f_k and g_k mod p and q ((f_k*g_k)^mod p) mod q where p|(q-1)
   * compute the coefficients
   * -----------multiple snapshots of ((f_k*g_k)^mod p) mod q with ascending q and the same p
   * -----------group the like terms of the snapshots and then compute CRT
   */
  fmpz_sparse_t f, g, f_s, g_s, h_1, h_2, f_2, g_2, temp_1, temp_2;
  fmpz_t p, q;
  slong i;

  fmpz_sparse_init(f);
  fmpz_sparse_init(g);
  fmpz_sparse_init(f_s);
  fmpz_sparse_init(g_s);
  fmpz_sparse_init(f_2);
  fmpz_sparse_init(g_2);
  fmpz_sparse_init(h_1);
  fmpz_sparse_init(h_2);
  fmpz_sparse_init(temp_1);
  fmpz_sparse_init(temp_2);
  
  fmpz_init(p);
  fmpz_one(q);

  fmpz_sparse_set(f, poly1);
  fmpz_sparse_set(g, poly2);

  flint_printf("BEFORE:\n");
  fmpz_sparse_print(f);
  flint_printf("\n");
  fmpz_sparse_print(g);
  flint_printf("\n");

  /* make f_s and g_s*/
  fmpz_sparse_set(f_s, poly1);
  fmpz_sparse_set(g_s, poly2);

  for(i = 0; i < poly1->length; i++)
  {
    fmpz_one(f_s->coeffs + i);
  }

  for(i = 0; i < poly2->length; i++)
  {
    fmpz_one(g_s->coeffs + i);
  }
  
  /*for(i = 0; i < f_s->length; i++)
  {
    fmpz_set(f_s->expons + i, f->expons + i);
  }

  for(i = 0; i < g_s->length; i++)
  {
    fmpz_set(g_s->expons + i, g->expons + i);
  }*/

  flint_printf("poly1: %d, f: %d, f_s: %d\n", poly1->length, f->length, f_s->length);
  /* estimate structural sparsity */
  do
  {
    fmpz_randm(p, state, f_s->expons);
    while(fmpz_is_zero(p))
    {
      fmpz_randm(p, state, f_s->expons);
    }

    fmpz_randm(q, state, f_s->expons);
    while(fmpz_is_zero(q))
    {
      fmpz_randm(q, state, f_s->expons);
    }

    fmpz_sparse_mul_heaps(temp_1, f_s, g_s);
    _fmpz_vec_scalar_mod_fmpz(temp_1->expons, temp_1->expons, temp_1->length, p);
    _fmpz_vec_scalar_mod_fmpz(temp_1->expons, temp_1->expons, temp_1->length, q);
    _fmpz_sparse_normalise(temp_1);
  } while(temp_1->length >= (f_s->length+g_s->length)/2);
  
  flint_printf("AFTER:\n");
  fmpz_sparse_print(f_s);
  flint_printf("\n");
  fmpz_sparse_print(g_s);
  flint_printf("\n");


  /* compute structural support */
  /*fmpz_sparse_set(temp_1, f_s);
  fmpz_sparse_set(temp_2, g_s);*/

  /* h_1 */
  /*_fmpz_vec_scalar_mod_fmpz(temp_1->expons, f_s->expons, f_s->length, p);
  _fmpz_vec_scalar_mod_fmpz(temp_2->expons, g_s->expons, g_s->length, p);
  _fmpz_sparse_normalise(temp_1);
  _fmpz_sparse_normalise(temp_2);
 fmpz_sparse_mul_heaps(h_1, temp_1, temp_2);
  _fmpz_vec_scalar_mod_fmpz(h_1->expons, h_1->expons, h_1->length, p);
  _fmpz_sparse_normalise(h_1);*/

  /* f_2 */
  /*fmpz_add(q, f_s->expons, g_s->expons);
  fmpz_mul_ui(q, q, 100);
 
  for(i = 0; i < f_s->length; i++)
  {
    fmpz_mul(f_s->coeffs + i, f_s->expons + i, q);
    fmpz_add_ui(f_s->coeffs + i, f_s->coeffs + i, 1);
  }

  _fmpz_vec_scalar_mod_fmpz(f_s->expons, f_s->expons, f_s->length, p);
  _fmpz_sparse_normalise(f_s);*/

  /* g_2 */
  /*for(i = 0; i < g_s->length; i++)
  {
    fmpz_mul(g_s->coeffs + i, g_s->expons + i, q);
    fmpz_add_ui(g_s->coeffs + i, g_s->coeffs + i, 1);
  }

  _fmpz_vec_scalar_mod_fmpz(g_s->expons, g_s->expons, g_s->length, p);
  _fmpz_sparse_normalise(g_s);*/

  /* h_2 */
  /*fmpz_sparse_mul_heaps(h_2, f_2, g_2);
  _fmpz_vec_scalar_mod_fmpz(h_2->expons, h_2->expons, h_2->length, p);

  fmpz_mul(p, q, q);
  _fmpz_vec_scalar_mod_fmpz(h_2->coeffs, h_2->coeffs, h_2->length, p);

  _fmpz_sparse_normalise(h_2);*/

  /* take coefficient ratios */
  /*for(i = 0; i < h_2->length; i++)
  {
    fmpz_cdiv_q(h_2->coeffs + i, h_2->coeffs + i, h_1->coeffs + i);
    fmpz_sub_ui(h_2->coeffs + i, h_2->coeffs + i, 1);
    fmpz_cdiv_q(h_2->coeffs + i, h_2->coeffs + i, q);
  }*/

  i = 2;
  fmpz_randtest(q, state, 100);
  fmpz_sparse_randtest(res, state, i, q, 100);
}
