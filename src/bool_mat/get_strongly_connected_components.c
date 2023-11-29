/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "bool_mat.h"

/* fixed-capacity stack */
typedef struct
{
    slong *data;
    slong capacity;
    slong size;
} _si_stack_struct;

typedef _si_stack_struct _si_stack_t[1];

static void
_si_stack_init(_si_stack_t s, slong capacity)
{
    s->data = flint_malloc(capacity * sizeof(slong));
    s->capacity = capacity;
    s->size = 0;
}

static void
_si_stack_clear(_si_stack_t s)
{
    flint_free(s->data);
}

static void
_si_stack_push(_si_stack_t s, slong x)
{
    if (s->size >= s->capacity)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    s->data[s->size++] = x;
}

static slong
_si_stack_pop(_si_stack_t s)
{
    slong x;
    if (s->size <= 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    x = s->data[s->size - 1];
    s->size--;
    return x;
}


/* struct for Tarjan's strongly connected components algorithm */
typedef struct
{
    slong *index;
    slong *lowlink;
    int *onstack;
    _si_stack_t S;
    slong nsccs;
    slong dim;
    slong idx;
} _tarjan_struct;

typedef _tarjan_struct _tarjan_t[1];

static const slong _tarjan_UNDEFINED = -1;

static slong *
_tarjan_index(_tarjan_t t, slong i)
{
    return t->index + i;
}

static slong *
_tarjan_lowlink(_tarjan_t t, slong i)
{
    return t->lowlink + i;
}

static int
_tarjan_onstack(_tarjan_t t, slong i)
{
    return t->onstack[i];
}

static void
_tarjan_push(_tarjan_t t, slong v)
{
    _si_stack_push(t->S, v);
    t->onstack[v] = 1;
}

static slong
_tarjan_pop(_tarjan_t t)
{
    slong v;
    v = _si_stack_pop(t->S);
    t->onstack[v] = 0;
    return v;
}

static slong
_tarjan_next_scc(_tarjan_t t)
{
    return t->nsccs++;
}

static slong
_tarjan_next_idx(_tarjan_t t)
{
    return t->idx++;
}

static void
_tarjan_init(_tarjan_t t, slong dim)
{
    slong i;
    t->index = flint_calloc(dim, sizeof(slong));
    t->lowlink = flint_calloc(dim, sizeof(slong));
    t->onstack = flint_calloc(dim, sizeof(int));
    _si_stack_init(t->S, dim);
    t->dim = dim;
    t->nsccs = 0;
    t->idx = 0;
    for (i = 0; i < dim; i++)
    {
        t->index[i] = _tarjan_UNDEFINED;
    }
}

static void
_tarjan_clear(_tarjan_t t)
{
    flint_free(t->index);
    flint_free(t->lowlink);
    flint_free(t->onstack);
    _si_stack_clear(t->S);
}

static void
_tarjan_strongconnect(slong *sccs, _tarjan_t t, const bool_mat_t A, slong v)
{
    slong idx, w, scc;

    idx = _tarjan_next_idx(t);
    *_tarjan_index(t, v) = idx;
    *_tarjan_lowlink(t, v) = idx;
    _tarjan_push(t, v);
    for (w = 0; w < t->dim; w++)
    {
        if (bool_mat_get_entry(A, v, w))
        {
            if (*_tarjan_index(t, w) == _tarjan_UNDEFINED)
            {
                _tarjan_strongconnect(sccs, t, A, w);
                *_tarjan_lowlink(t, v) = FLINT_MIN(
                        *_tarjan_lowlink(t, v), *_tarjan_lowlink(t, w));
            }
            else if (_tarjan_onstack(t, w))
            {
                *_tarjan_lowlink(t, v) = FLINT_MIN(
                        *_tarjan_lowlink(t, v), *_tarjan_index(t, w));
            }
        }
    }
    if (*_tarjan_lowlink(t, v) == *_tarjan_index(t, v))
    {
        scc = _tarjan_next_scc(t);
        while (w != v)
        {
            w = _tarjan_pop(t);
            if (sccs[w] != _tarjan_UNDEFINED)
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);
            sccs[w] = scc;
        }
    }
}


/* following Tarjan */
slong
bool_mat_get_strongly_connected_components(slong *partition, const bool_mat_t A)
{
    slong i, n, result;
    _tarjan_t t;

    if (!bool_mat_is_square(A))
    {
        flint_throw(FLINT_ERROR, "bool_mat_get_strongly_connected_components: "
                     "a square matrix is required!\n");
    }

    if (bool_mat_is_empty(A))
        return 0;

    n = bool_mat_nrows(A);

    if (n == 1)
    {
        partition[0] = 0;
        return 1;
    }

    _tarjan_init(t, n);

    for (i = 0; i < n; i++)
    {
        partition[i] = _tarjan_UNDEFINED;
    }
    for (i = 0; i < n; i++)
    {
        if (*_tarjan_index(t, i) == _tarjan_UNDEFINED)
        {
            _tarjan_strongconnect(partition, t, A, i);
        }
    }

    result = t->nsccs;

    _tarjan_clear(t);

    return result;
}
