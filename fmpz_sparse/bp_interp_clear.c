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

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_sparse.h"

void fmpz_sparse_bp_interp_clear(fmpz_sparse_bp_interp_t res)
{
  fmpz_clear(res->q);
  fmpz_clear(res->order);
  _fmpz_vec_clear(res->sample_points, res->length);
  _fmpz_vec_clear(res->evaluations, res->length);
}
