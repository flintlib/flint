/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "bool_mat.h"

/*
 * Postorder traversal of a DAG follows
 * https://en.wikipedia.org/wiki/Topological_sorting#Depth-first_search
 */

typedef struct
{
    int *u; /* 'temporary mark' */
    int *v; /* 'permanent mark' */
    slong *post; /* postorder nodes */
    slong npost; /* number of postorder nodes so far */
    slong size;
} _toposort_s;

static void
_toposort_init(_toposort_s *s, slong size)
{
    s->size = size;
    s->u = flint_calloc(size, sizeof(int));
    s->v = flint_calloc(size, sizeof(int));
    s->post = flint_malloc(size * sizeof(slong));
    s->npost = 0;
}

static void
_toposort_clear(_toposort_s *s)
{
    flint_free(s->u);
    flint_free(s->v);
    flint_free(s->post);
}

static int
_toposort_visit(_toposort_s *s, const bool_mat_t A, slong n)
{
    if (s->u[n])
        return 1;
    if (!s->v[n])
    {
        slong m;
        s->u[n] = 1;
        for (m = 0; m < s->size; m++)
            if (bool_mat_get_entry(A, n, m))
                if (_toposort_visit(s, A, m))
                    return 1;
        s->v[n] = 1;
        s->u[n] = 0;
        s->post[s->npost++] = n;
    }
    return 0;
}

slong
bool_mat_nilpotency_degree(const bool_mat_t A)
{
    slong n;

    if (!bool_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "bool_mat_nilpotency_degree: a square matrix is required!\n");
    }

    if (bool_mat_is_empty(A))
        return 0;

    n = bool_mat_nrows(A);

    if (n == 1)
    {
        return bool_mat_get_entry(A, 0, 0) ? -1 : 1;
    }
    else
    {
        _toposort_s s;
        slong i;
        int has_cycle;
        int result;

        _toposort_init(&s, n);

        for (has_cycle = 0, i = 0; !has_cycle && i < n; i++)
            if (!s.v[i])
                has_cycle = _toposort_visit(&s, A, i);

        if (has_cycle)
        {
            result = -1;
        }
        else
        {
            /* Find the length of the longest path within the DAG */
            /* https://stackoverflow.com/a/10737524/4072759 */

            slong x, y, z;
            slong max_overall;
            fmpz_mat_t E;

            fmpz_mat_init(E, n, n);
            fmpz_mat_zero(E);
            max_overall = 0;
            for (i = n - 1; i >= 0; i--)
            {
                slong max_in = 0;
                y = s.post[i];
                for (x = 0; x < n; x++)
                {
                    max_in = FLINT_MAX(max_in,
                                       fmpz_get_si(fmpz_mat_entry(E, x, y)));
                }
                for (z = 0; z < n; z++)
                {
                    if (bool_mat_get_entry(A, y, z))
                    {
                        fmpz_set_si(fmpz_mat_entry(E, y, z), max_in + 1);
                        max_overall = FLINT_MAX(max_overall, max_in + 1);
                    }
                }
            }
            fmpz_mat_clear(E);
            result = max_overall + 1;
        }
        _toposort_clear(&s);
        return result;
    }
}
