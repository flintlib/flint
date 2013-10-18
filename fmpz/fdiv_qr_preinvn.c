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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"

void
fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g, 
                         const fmpz_t h, const fmpz_preinvn_t inv)
{
    fmpz c1 = *g;
    fmpz c2 = *h;

    if (fmpz_is_zero(h))
    {
        flint_printf("Exception (fmpz_fdiv_q). Division by zero.\n");
        abort();
    }

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* h is also small */
        {
            mp_limb_t hi, lo, p[2], cy, dinv, q;

            /* set up dividend and divisor, normalising */
            mp_limb_t a = FLINT_ABS(c1);
            mp_limb_t d = (FLINT_ABS(c2) << inv->norm);
            dinv = inv->dinv[0];

            hi = r_shift(a, FLINT_BITS - inv->norm);
            lo = (a << inv->norm);

            /* precomputed divrem */
            if (hi)
            {
               umul_ppmm(p[1], p[0], dinv, hi);
               add_ssaaaa(cy, q, 0, p[1], 0, hi);

               umul_ppmm(p[1], p[0], d, q);
               
               cy = hi - p[1];
               sub_ddmmss(cy, lo, cy, lo, 0, p[0]);

               while (cy)
               {
                  sub_ddmmss(cy, lo, cy, lo, 0, d);
                  q++;
               }
            } else
               q = 0;

            if (lo >= d)
            {
               lo -= d;
               q++;
            }
           
            /* deal with signed inputs */
            if (((c2 ^ c1) < 0))
            {
               if (lo != 0)
               {
                  q++;
                  lo = d - lo;
               }

               q = -q;
            }

            fmpz_set_si(f, q);
            if (c2 < 0)
               fmpz_set_si(s, -(lo >> inv->norm)); /* shift remainder */
            else
               fmpz_set_si(s, lo >> inv->norm);
        }
        else                    /* h is large and g is small */
        {
            if (c1 == WORD(0))
            {
                fmpz_set_ui(f, WORD(0)); /* g is zero */
                fmpz_set_si(s, c1);
            }
            else if ((c1 < WORD(0) && fmpz_sgn(h) < 0) || (c1 > WORD(0) && fmpz_sgn(h) > 0))  /* signs are the same */
            {
                fmpz_zero(f);   /* quotient is positive, round down to zero */
                fmpz_set_si(s, c1);
            }
            else
            {
                fmpz_add(s, g, h);
                fmpz_set_si(f, WORD(-1));    /* quotient is negative, round down to minus one */
            }
        }
    }
    else /* g is large */
    {
        __mpz_struct *mpz_ptr, *mpz_ptr2;
        slong size1, size2 = FLINT_ABS(fmpz_size(h));
        mpz_t qt, rt;

        /* set up quotient and remainder for division */
        mpz_init(rt);
        mpz_mul_2exp(rt, COEFF_TO_PTR(c1), inv->norm);
        size1 = FLINT_ABS(rt->_mp_size);

        if (size1 < size2)
        {
           mpz_init(qt);
           mpz_abs(rt, rt);
        } else
        {
           mpz_init2(qt, (size1 - size2 + 1)*FLINT_BITS);

		     if (!COEFF_IS_MPZ(c2))  /* h is small */
           {
              ulong c2u = (FLINT_ABS(c2) << inv->norm);
              qt->_mp_d[size1 - size2] 
                 = flint_mpn_divrem_preinvn(qt->_mp_d, rt->_mp_d, size1, &c2u, 1, inv->dinv); 
              qt->_mp_size = size1 - size2 + (qt->_mp_d[size1 - size2] != 0);
              rt->_mp_size = 1 - (rt->_mp_d[0] == 0);
           }
           else /* both are large */
           {
              mp_ptr d;
            
              /* compute shifted d */
              if (inv->norm)
              {
                 d = flint_malloc(FLINT_ABS(COEFF_TO_PTR(c2)->_mp_size)*sizeof(mp_limb_t));
                 mpn_lshift(d, COEFF_TO_PTR(c2)->_mp_d, size2, inv->norm);
              } else
                 d = COEFF_TO_PTR(c2)->_mp_d;
            
              qt->_mp_d[size1 - size2] 
                 = flint_mpn_divrem_preinvn(qt->_mp_d, rt->_mp_d, size1, d, size2, inv->dinv); 
            
              /* normalise */
              qt->_mp_size = size1 - size2 + (qt->_mp_d[size1 - size2] != 0);
              while (size2 && rt->_mp_d[size2 - 1] == 0) size2--;
              rt->_mp_size = size2;

              if (inv->norm)
                 flint_free(d);
           }
        }
           
        /* shift remainer */
        mpz_fdiv_q_2exp(rt, rt, inv->norm);

        /* deal with signed inputs */
        if ((fmpz_sgn(g) ^ fmpz_sgn(h)) < 0)
        {
           if (mpz_sgn(rt) != 0)
           {
              mpz_add_ui(qt, qt, 1);
              if (!COEFF_IS_MPZ(c2))
                 mpz_ui_sub(rt, FLINT_ABS(c2), rt);
              else if (fmpz_sgn(h) < 0)
              {
                 mpz_add(rt, rt, COEFF_TO_PTR(c2));
                 mpz_neg(rt, rt);
              } else
                 mpz_sub(rt, COEFF_TO_PTR(c2), rt);
           }

           mpz_neg(qt, qt);
        }

        if (fmpz_sgn(h) < 0)
           mpz_neg(rt, rt);

        /* set results */
        _fmpz_promote(f); /* must not hang on to ptr whilst promoting s */
        mpz_ptr2 = _fmpz_promote(s);
		  mpz_ptr  = COEFF_TO_PTR(*f);

        mpz_swap(mpz_ptr, qt);
        mpz_swap(mpz_ptr2, rt);

        /* clean up */
        mpz_clear(qt);
        mpz_clear(rt);
            
        _fmpz_demote_val(f);    /* division by h may result in small value */
        _fmpz_demote_val(s);    /* division by h may result in small value */
    }
}
