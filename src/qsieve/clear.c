/*
    Copyright (C) 2006, 2011, 2016 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qsieve.h"

void qsieve_clear(qs_t qs_inf)
{
    fmpz_clear(qs_inf->n);
    fmpz_clear(qs_inf->kn);

    flint_free(qs_inf->factor_base);
    flint_free(qs_inf->sqrts);

    qs_inf->factor_base = NULL;
    qs_inf->sqrts       = NULL;

    flint_free(qs_inf->fname);
}
