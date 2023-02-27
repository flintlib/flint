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

void qfb_hash_clear(qfb_hash_t * qhash, slong depth)
{
   slong i, size = (1L<<depth);
   
   for (i = 0; i < size; i++)
   {
      qfb_clear(qhash[i].q);
      qfb_clear(qhash[i].q2);
   }

   flint_free(qhash);
}
