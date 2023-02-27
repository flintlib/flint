/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "qfb.h"

qfb_hash_t * qfb_hash_init(slong depth)
{
   slong i, size = (1L<<depth);
   qfb_hash_t * qhash = flint_malloc(size*sizeof(qfb_hash_t));

   for (i = 0; i < size; i++)
   {
      qfb_init(qhash[i].q);
      qfb_init(qhash[i].q2);
   }

   return qhash;
}
