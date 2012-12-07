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

void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p)
{
   fmpz_t s, t;

   fmpz_init(s);
   fmpz_init(t);
   
   fmpz_mod(t, D, p);
   fmpz_sqrtmod(t, t, p);
   fmpz_sub(s, D, t);

   if (fmpz_is_odd(s))
      fmpz_sub(t, p, t);

   fmpz_set(r->a, p);
   fmpz_set(r->b, t);
   fmpz_mul(r->c, r->b, r->b);
   fmpz_sub(r->c, r->c, D);
   fmpz_divexact(r->c, r->c, r->a);
   fmpz_fdiv_q_2exp(r->c, r->c, 2);
   
   fmpz_clear(s);
   fmpz_clear(t);
}
