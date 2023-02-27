/*
    Copyright (C) 2011, 2016 William Hart
    Copyright (C) 2011 Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  /* qsort */
#include "fmpz.h"
#include "fmpz_mat.h"

typedef struct
{
   ulong col;
   ulong hash;
} col_hash_struct;

typedef col_hash_struct * col_hash_ptr;
typedef const col_hash_struct * col_hash_srcptr;

/* cheap hash of all columns of M */
void
fmpz_mat_col_hash(col_hash_ptr col_h, fmpz_mat_t M)
{
    ulong i, j, hash;

    for (i = 0; i < M->c; i++)
    {
        col_h[i].col = i;
        hash = 0;

        for (j = 0; j < M->r; j++)
        {
            hash ^= fmpz_get_ui(M->rows[j] + i);
            hash = (hash << 1) + (hash >> (FLINT_BITS - 1));
        }

        col_h[i].hash = hash;
    }
}

int
fmpz_mat_col_hash_compare(const void * a, const void * b)
{
   col_hash_srcptr col_a = a;
   col_hash_srcptr col_b = b;

   if (col_a->hash == col_b->hash)
      return 0;

   return 2*(col_a->hash > col_b->hash) - 1;
}

int
fmpz_mat_col_partition(slong * part, fmpz_mat_t M, int short_circuit)
{
    slong start = 0, new_start = 0, upto = 1, p = 0, i, count;
    ulong hash;
    col_hash_ptr col_h;
    TMP_INIT;

    TMP_START;
    col_h = TMP_ALLOC(sizeof(col_hash_struct) * M->c);

    fmpz_mat_col_hash(col_h, M);

    qsort(col_h, M->c, sizeof(col_hash_struct), fmpz_mat_col_hash_compare);

    if (short_circuit)
    {
        hash = col_h[0].hash;
        count = 1;

        for (i = 1; i < M->c; i++)
        {
            if (col_h[i].hash != hash)
            {
                count++;
                hash = col_h[i].hash;
            }
        }

        if (count > M->r)
            goto cleanup;
    }

    for (i = 0; i < M->c; i++)
        part[i] = -WORD(1);

    while (start < M->c)
    {
        new_start = start;
        p++;

        if (short_circuit && p > M->r)
        {
            p = 0; /* too many partitions */
            goto cleanup;
        }

        part[col_h[start].col] = p;

        for (upto = start + 1; upto < M->c && col_h[upto].hash == col_h[start].hash; upto++)
        {
            if (part[col_h[upto].col] == -WORD(1))
            {
                if (!fmpz_mat_col_equal(M, col_h[start].col, col_h[upto].col))
                {
                    if (new_start == start)
                        new_start = upto;
                } else
                    part[col_h[upto].col] = p;
            }
        }

        start = start == new_start ? upto : new_start;
    }

cleanup:
    TMP_END;

    return p;
}
