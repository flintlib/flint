/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2013 William Hart

******************************************************************************/

#include "nf.h"

void nf_clear(nf_t nf)
{
    fmpq_poly_clear(nf->pol);

    if (!(nf->flag & NF_MONIC))
       fmpz_preinvn_clear(nf->pinv.qq);

    if (nf->pol->length <= NF_POWERS_CUTOFF && nf->pol->length > 3)
    {
       if (nf->flag & NF_MONIC)
          _fmpz_poly_powers_clear(nf->powers.zz->powers, nf->powers.zz->len);
       else
          _fmpq_poly_powers_clear(nf->powers.qq->powers, nf->powers.qq->len);
    }

    fmpq_poly_clear(nf->traces);
}

