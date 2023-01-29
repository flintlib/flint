/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson
    Copyright (C) 2015 William Hart

******************************************************************************/

#include "nf_elem.h"

void nf_elem_print_pretty(const nf_elem_t a, const nf_t nf, const char * var)
{
    if (nf->flag & NF_LINEAR)
    {
        const fmpz * const den = LNF_ELEM_DENREF(a);
		fmpz_print(LNF_ELEM_NUMREF(a));
        if (!fmpz_is_one(den))
		{
		   flint_printf("/");
		   fmpz_print(LNF_ELEM_DENREF(a));
		}       
    } else if (nf->flag & NF_QUADRATIC)
    {
        const fmpz * const anum = QNF_ELEM_NUMREF(a);
        const fmpz * const aden = QNF_ELEM_DENREF(a);
        int den1 = fmpz_is_one(aden);
		int lead0 = fmpz_is_zero(anum + 1);
		
        if (!den1 && !lead0)
		   flint_printf("(");
        if (!lead0)
		{
		   fmpz_print(anum + 1);
		   flint_printf("*%s", var);
		   if (fmpz_sgn(anum) >= 0)
		   printf("+");
		}
        fmpz_print(anum);
        if (!den1 && !lead0)
		   flint_printf(")");
        if (!den1)
		{
		   flint_printf("/");
		   fmpz_print(aden);
		}
    } else
    {
        fmpq_poly_print_pretty(NF_ELEM(a), var);
    }
}

