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

void qsieve_init(qs_t qs_inf, const fmpz_t n)
{
    slong i;

    qs_inf->fname = (char *) flint_malloc(20); /* space for filename */

    /* store n in struct */
    fmpz_init_set(qs_inf->n, n);

    /* determine the number of bits of n */
    qs_inf->bits =  fmpz_bits(n);

    /* determine which index in the tuning table n corresponds to */
    for (i = 1; i < QS_TUNE_SIZE; i++)
    {
        if (qsieve_tune[i][0] > qs_inf->bits)
            break;
    }
    i--;

    qs_inf->ks_primes = qsieve_tune[i][1]; /* number of Knuth-Schroeppel primes */
    qs_inf->num_primes  = 0;
    qs_inf->num_relations = 0;
    qs_inf->full_relation = 0;
    qs_inf->num_cycles = 0;
    qs_inf->columns = 0;
    qs_inf->vertices = 0;
    qs_inf->components = 0;
    qs_inf->edges = 0;
#if QS_DEBUG
    qs_inf->poly_count = 0;
#endif

    fmpz_init(qs_inf->kn); /* initialise kn */

    qs_inf->factor_base = NULL;
    qs_inf->sqrts       = NULL;

    qs_inf->s = 0;
}
