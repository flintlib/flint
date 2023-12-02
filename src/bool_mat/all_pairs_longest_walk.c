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
 * Condensation of a matrix.
 * This is the directed acyclic graph of strongly connected components.
 */
typedef struct
{
    slong n; /* number of vertices in the original graph */
    slong k; /* number of strongly connected components (sccs) */
    bool_mat_t C; /* adjacency matrix of the sccs in the condensation */
    slong *partition; /* maps the vertex index to the scc index */
} _condensation_struct;

typedef _condensation_struct _condensation_t[1];

static void
_condensation_init(_condensation_t c, const bool_mat_t A)
{
    slong i, j, u, v;

    if (!bool_mat_is_square(A))
        flint_throw(FLINT_ERROR, "_condensation_init: A must be square\n");

    c->n = bool_mat_nrows(A);
    c->partition = flint_malloc(c->n * sizeof(slong));

    c->k = bool_mat_get_strongly_connected_components(c->partition, A);

    /*
     * Compute the adjacency matrix of the condensation.
     * This should be strict lower triangular, so that visiting the
     * vertices in increasing order corresponds to a postorder traversal.
     */
    bool_mat_init(c->C, c->k, c->k);
    bool_mat_zero(c->C);
    for (i = 0; i < c->n; i++)
    {
        for (j = 0; j < c->n; j++)
        {
            if (bool_mat_get_entry(A, i, j))
            {
                u = c->partition[i];
                v = c->partition[j];
                if (u != v)
                {
                    bool_mat_set_entry(c->C, u, v, 1);
                }
            }
        }
    }

    if (!bool_mat_is_lower_triangular(c->C) || bool_mat_trace(c->C))
        flint_throw(FLINT_ERROR, "_condensation_init: unexpected matrix structure\n");
}

static void
_condensation_clear(_condensation_t c)
{
    bool_mat_clear(c->C);
    flint_free(c->partition);
}



typedef struct
{
    _condensation_t con;
    bool_mat_t T; /* transitive closure of condensation */
    bool_mat_t P; /* is there a cycle in any component on a path from u to v */
    fmpz_mat_t Q; /* longest path, if any, from u to v */
    int *scc_has_cycle;
} _connectivity_struct;
typedef _connectivity_struct _connectivity_t[1];

static void
_connectivity_clear(_connectivity_t c)
{
    bool_mat_clear(c->T);
    bool_mat_clear(c->P);
    fmpz_mat_clear(c->Q);
    flint_free(c->scc_has_cycle);
    _condensation_clear(c->con);
}

static void
_connectivity_init_scc_has_cycle(_connectivity_t c, const bool_mat_t A)
{
    slong n, i, u;
    slong *scc_size;

    n = bool_mat_nrows(A);
    c->scc_has_cycle = flint_calloc(n, sizeof(int));

    /*
     * If a vertex of the original graph has a loop,
     * then the strongly connected component to which it belongs has a cycle.
     */
    for (i = 0; i < n; i++)
    {
        if (bool_mat_get_entry(A, i, i))
        {
            u = c->con->partition[i];
            c->scc_has_cycle[u] = 1;
        }
    }

    /*
     * If a strongly connected component contains more than one vertex,
     * then that component has a cycle.
     */
    scc_size = flint_calloc(c->con->k, sizeof(slong));
    for (i = 0; i < n; i++)
    {
        u = c->con->partition[i];
        scc_size[u]++;
    }
    for (u = 0; u < c->con->k; u++)
    {
        if (scc_size[u] > 1)
        {
            c->scc_has_cycle[u] = 1;
        }
    }
    flint_free(scc_size);
}

static void
_connectivity_init(_connectivity_t c, const bool_mat_t A)
{
    slong u, v, w;
    slong k;
    slong curr, rest;

    /* compute the condensation */
    _condensation_init(c->con, A);
    k = c->con->k;

    /* check whether each scc contains a cycle */
    _connectivity_init_scc_has_cycle(c, A);

    /* compute the transitive closure of the condensation */
    bool_mat_init(c->T, k, k);
    bool_mat_transitive_closure(c->T, c->con->C);

    /*
     * Is there a walk from u to v that passes through a cycle-containing scc?
     * Cycles in the components u and v themselves are not considered.
     * Remember that the condensation is a directed acyclic graph.
     */
    bool_mat_init(c->P, k, k);
    bool_mat_zero(c->P);
    for (w = 0; w < k; w++)
    {
        if (c->scc_has_cycle[w])
        {
            for (u = 0; u < k; u++)
            {
                for (v = 0; v < k; v++)
                {
                    if (bool_mat_get_entry(c->T, u, w) &&
                        bool_mat_get_entry(c->T, w, v))
                    {
                        bool_mat_set_entry(c->P, u, v, 1);
                    }
                }
            }
        }
    }

    /*
     * What is the max length path from u to v in the condensation graph?
     * If u==v or if v is unreachable from u then let this be zero.
     * Remember that the condensation is a directed acyclic graph,
     * and that the components are indexed in a post-order traversal.
     */
    fmpz_mat_init(c->Q, k, k);
    fmpz_mat_zero(c->Q);
    for (u = 0; u < k; u++)
    {
        for (w = 0; w < k; w++)
        {
            if (bool_mat_get_entry(c->con->C, u, w))
            {
                curr = fmpz_get_si(fmpz_mat_entry(c->Q, u, w));
                fmpz_set_si(
                        fmpz_mat_entry(c->Q, u, w),
                        FLINT_MAX(curr, 1));
                for (v = 0; v < k; v++)
                {
                    if (bool_mat_get_entry(c->T, w, v))
                    {
                        rest = fmpz_get_si(fmpz_mat_entry(c->Q, w, v));
                        curr = fmpz_get_si(fmpz_mat_entry(c->Q, u, v));
                        fmpz_set_si(
                                fmpz_mat_entry(c->Q, u, v),
                                FLINT_MAX(curr, rest + 1));
                    }
                }
            }
        }
    }
}


static void
_connectivity_entrywise_nilpotence_degree(
        fmpz_t N, _connectivity_t c, slong i, slong j)
{
    slong u, v;
    u = c->con->partition[i];
    v = c->con->partition[j];
    if (u == v)
    {
        if (c->scc_has_cycle[u])
        {
            fmpz_set_si(N, -1);
        }
        else
        {
            fmpz_one(N);
        }
    }
    else if (!bool_mat_get_entry(c->T, u, v))
    {
        fmpz_zero(N);
    }
    else if (
            c->scc_has_cycle[u] ||
            c->scc_has_cycle[v] ||
            bool_mat_get_entry(c->P, u, v))
    {
        fmpz_set_si(N, -1);
    }
    else
    {
        fmpz_add_ui(N, fmpz_mat_entry(c->Q, u, v), 1);
    }
}

slong
bool_mat_all_pairs_longest_walk(fmpz_mat_t B, const bool_mat_t A)
{
    slong n;

    if (!bool_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "bool_mat_all_pairs_longest_walk: "
                     "a square matrix is required!\n");
    }

    if (bool_mat_is_empty(A))
        return -1;

    n = bool_mat_nrows(A);

    if (n == 1)
    {
        if (bool_mat_get_entry(A, 0, 0))
        {
            fmpz_set_si(fmpz_mat_entry(B, 0, 0), -2);
            return -2;
        }
        else
        {
            fmpz_set_si(fmpz_mat_entry(B, 0, 0), 0);
            return 0;
        }
    }
    else
    {
        slong i, j, result;
        _connectivity_t c;

        _connectivity_init(c, A);

        result = -1;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                slong x;
                fmpz *p = fmpz_mat_entry(B, i, j);
                _connectivity_entrywise_nilpotence_degree(p, c, i, j);
                fmpz_sub_ui(p, p, 1);
                if (result != -2)
                {
                    x = fmpz_get_si(p);
                    if (x == -2)
                    {
                        result = -2;
                    }
                    else
                    {
                        result = FLINT_MAX(result, x);
                    }
                }
            }
        }

        _connectivity_clear(c);

        return result;
    }
}
