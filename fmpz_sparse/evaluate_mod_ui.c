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

#include "ulong_extras.h"
#include "fmpz_sparse.h"

ulong fmpz_sparse_evaluate_mod_ui(const fmpz_sparse_t poly, ulong a, ulong m)
{
  ulong temp;
  ulong redco;
  slong i;
  slong res = UWORD(0);
  ulong minv;
  ulong redpow;
  ulong mphi;

  if (poly->length == 0) return res;

  minv = n_preinvert_limb(m);
  mphi = n_euler_phi(m);

  for (i=0; i<poly->length; ++i)
  {
    redpow = fmpz_fdiv_ui(poly->expons+i, mphi);
    temp = n_powmod2_ui_preinv(a, redpow, m, minv);
    redco = fmpz_fdiv_ui(poly->coeffs+i, m);
    temp = n_mulmod2_preinv(temp, redco, m, minv);
    res = n_addmod(res, temp, m);
  }

  return res;
}
