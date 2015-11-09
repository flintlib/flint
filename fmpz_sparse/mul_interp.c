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
fmpz_sparse_mul_interp(fmpz_sparse_t res, const fmpz_sparse_t poly1, 
    const fmpz_sparse_t poly2)
{
  /*
   * (substitutions for kronecker substitution)
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

}
