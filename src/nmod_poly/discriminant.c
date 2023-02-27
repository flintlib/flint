/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

mp_limb_t 
_nmod_poly_discriminant(mp_srcptr poly, slong len, nmod_t mod)
{
   mp_ptr der = _nmod_vec_init(len - 1);
   slong dlen = len - 1;
   mp_limb_t res, pow;

   _nmod_poly_derivative(der, poly, len, mod);
   NMOD_VEC_NORM(der, dlen);

   if (dlen == 0)
   {
       _nmod_vec_clear(der);
       return 0;
   }

   res = _nmod_poly_resultant(poly, len, der, dlen, mod);
   pow = n_powmod2_preinv(poly[len - 1], len - dlen - 2, mod.n, mod.ninv);
   res = n_mulmod2_preinv(res, pow, mod.n, mod.ninv);

   if ((len & 3) == 0 || (len & 3) == 3) /* degree is not 0, 1 mod 4 */
      res = nmod_neg(res, mod);

   _nmod_vec_clear(der);

   return res;
}

mp_limb_t 
nmod_poly_discriminant(const nmod_poly_t f)
{
    const slong len = f->length;
    
    if (len <= 1)
        return 0;
    else
       return _nmod_poly_discriminant(f->coeffs, len, f->mod);
}

