/*
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nf_elem.h"


void _nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
{
   nf_elem_t cinv;

   nf_elem_init(cinv, nf);

   _nf_elem_inv(cinv, c, nf);
   _nf_elem_mul(a, b, cinv, nf);

   nf_elem_clear(cinv, nf);
}


void nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
{
   nf_elem_t t;

   if (a == b)
   {
      nf_elem_init(t, nf);

      _nf_elem_div(t, b, c, nf);
      nf_elem_swap(t, a, nf);

      nf_elem_clear(t, nf);
   }
   else
      _nf_elem_div(a, b, c, nf);

   nf_elem_canonicalise(a, nf);
}
