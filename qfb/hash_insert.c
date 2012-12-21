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

    Copyright (C) 2012 William Hart

******************************************************************************/

#undef ulong /* prevent clash with stdlib */
#include <stdlib.h>
#define ulong unsigned long
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

void qfb_hash_insert(qfb_hash_t * qhash, qfb_t q, qfb_t q2, long iter, long depth)
{
   long size = (1L<<depth), i;
   fmpz_t r;

   fmpz_init(r);

   fmpz_fdiv_r_2exp(r, q->a, depth);
   i = fmpz_get_ui(r);

   while (!fmpz_is_zero(qhash[i].q->a))
   {
      i++;
      if (i == size)
         i = 0;
   }

   qfb_set(qhash[i].q, q);
   qhash[i].iter = iter;
   if (q2 != NULL)
      qfb_set(qhash[i].q2, q2);

   fmpz_clear(r);
}
