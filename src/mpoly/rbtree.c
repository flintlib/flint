/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpoly.h"

/* red-black tree with ui keys ***********************************************/

void mpoly_rbtree_ui_init(mpoly_rbtree_ui_t T, slong data_size)
{
    mpoly_rbnode_ui_struct * nodes;

    T->length = 0;
    T->node_alloc = 2;
    T->nodes = FLINT_ARRAY_ALLOC(T->node_alloc, mpoly_rbnode_ui_struct);
    T->data = NULL;
    T->data_alloc = 0;
    T->data_size = data_size;

    nodes = T->nodes + 2;

    /* pointer to head */
    nodes[-1].up = -2;
    nodes[-1].left = -2;    /* head */
    nodes[-1].right = -2;
    nodes[-1].key = 0;
    nodes[-1].color = 0;

    /* null node */
    nodes[-2].up = -2;
    nodes[-2].left = -2;
    nodes[-2].right = -2;
    nodes[-2].key = 0;
    nodes[-2].color = 0;
}

/* the memory management of the data chunks must be done prior */
void mpoly_rbtree_ui_clear(mpoly_rbtree_ui_t T)
{
    FLINT_ASSERT(T->node_alloc > 0);
    flint_free(T->nodes);
    flint_free(T->data);
}

static void mpoly_rbtree_ui_fit_length(mpoly_rbtree_ui_t T, slong len)
{
    slong dsize = T->data_size;

    if (len + 2 > T->node_alloc)
    {
        slong new_alloc = FLINT_MAX(len + 2, 2*T->node_alloc);
        T->nodes = FLINT_ARRAY_REALLOC(T->nodes, new_alloc, mpoly_rbnode_ui_struct);
        T->node_alloc = new_alloc;
    }

    if (dsize*len > T->data_alloc)
    {
        slong new_alloc = FLINT_MAX(dsize*len, 2*T->data_alloc);
        T->data = FLINT_ARRAY_REALLOC(T->data, new_alloc, char);
        T->data_alloc = new_alloc;
    }
}

/*
    If the key rcx exists, a pointer to its data will be returned and new
    will be set to 0. Otherwise, a pointer to a new data chunk is returned
    and new is set to 1.
*/
void * mpoly_rbtree_ui_lookup(mpoly_rbtree_ui_t T, int * new, ulong rcx)
{
    slong dsize = T->data_size;
    mpoly_rbnode_ui_struct * nodes = T->nodes + 2;
    slong rax, rdx, r8, r9, r10, r11;
    slong n = T->length;

    r10 = nodes[-1].left;

    if (n < 1)
    {
        mpoly_rbtree_ui_fit_length(T, 1);
        nodes = T->nodes + 2;
        rax = 0;
        nodes[rax].up = -1;
        nodes[rax].left = -2;
        nodes[rax].right = -2;
        nodes[rax].color = 0;
        nodes[rax].key = rcx;
        T->length = 1;
        *new = 1;
        nodes[-1].left = rax;
        return T->data + dsize*rax;
    }

Compare:
    FLINT_ASSERT(r10 >= 0);
    r8 = nodes[r10].left;
    r9 = nodes[r10].right;
    if (rcx < nodes[r10].key)
        goto GoLeft;
    if (rcx > nodes[r10].key)
        goto GoRight;

    rax = r10;
    * new = 0;
    return T->data + dsize*rax;

GoLeft:
    if (r8 < 0)
    {
        FLINT_ASSERT(r8 == -2);
        goto MakeNewLeft;
    }
    r10 = r8;
    goto Compare;

GoRight:
    if (r9 < 0)
    {
        FLINT_ASSERT(r9 == -2);
        goto MakeNewRight;
    }
    r10 = r9;
    goto Compare;

MakeNewLeft:
    mpoly_rbtree_ui_fit_length(T, n + 1);
    nodes = T->nodes + 2;
    rdx = n;
    nodes[r10].left = rdx;
    goto FixTree;

MakeNewRight:
    mpoly_rbtree_ui_fit_length(T, n + 1);
    nodes = T->nodes + 2;
    rdx = n;
    nodes[r10].right = rdx;

FixTree:
    nodes[rdx].up = r10;
    nodes[rdx].left = -2;
    nodes[rdx].right = -2;
    nodes[rdx].color = 1;
    nodes[rdx].key = rcx;
    T->length = n + 1;
    *new = 1;
    rax = rdx;

FixNode:

/*Case1:*/
    r8 = nodes[rdx].up;
    if (r8 < 0)
    {
        FLINT_ASSERT(r8 == -1);
        FLINT_ASSERT(nodes[-1].left == rdx);
        nodes[rdx].color = 0;
        return T->data + dsize*rax;
    }

/*Case2:*/
    if (nodes[r8].color == 0)
        return T->data + dsize*rax;

/*Case3:*/
    r9 = nodes[r8].up;
    r10 = nodes[r9].left;
    r11 = nodes[r9].right;
    if (r8 == r10)
        r10 = r11;
    if (r10 < 0 || nodes[r10].color == 0)
        goto Case4;
    rdx = r9;
    nodes[r8].color = 0;
    nodes[r9].color = 1;
    nodes[r10].color = 0;
    goto FixNode;

Case4:
    r10 = nodes[r9].up;
/*Case4A:*/
    if (rdx != nodes[r8].right || r8 != nodes[r9].left)
        goto Case4B;
    r11 = nodes[rdx].left;
    nodes[r9].left = rdx;
    nodes[rdx].left = r8;
    nodes[r8].right = r11;
    nodes[r8].up = rdx;
    nodes[rdx].up = r9;
    nodes[r11].up = r8;
    goto Case4Done;

Case4B:
    if (rdx != nodes[r8].left || r8 != nodes[r9].right)
        goto Case5;
    r11 = nodes[rdx].right;
    nodes[r9].right = rdx;
    nodes[rdx].right = r8;
    nodes[r8].left = r11;
    nodes[r8].up = rdx;
    nodes[rdx].up = r9;
    nodes[r11].up = r8;

Case4Done:
    FLINT_SWAP(slong, rdx, r8);

Case5:

    if (nodes[r10].right == r9)
        nodes[r10].right = r8;
    if (nodes[r10].left == r9)
        nodes[r10].left = r8;

    nodes[r8].up = r10;
    nodes[r8].color = 0;
    nodes[r9].up = r8;
    nodes[r9].color = 1;

    r11 = nodes[r8].right;
    r10 = nodes[r8].left;
    if (rdx == r10)
    {
        nodes[r8].right = r9;
        nodes[r9].left = r11;
        nodes[r11].up = r9;
    }
    else
    {
        nodes[r8].left = r9;
        nodes[r9].right = r10;
        nodes[r10].up = r9;
    }

    return T->data + dsize*rax;
}


/* red-black tree with fmpz keys *********************************************/

void mpoly_rbtree_fmpz_init(mpoly_rbtree_fmpz_t T, slong data_size)
{
    mpoly_rbnode_fmpz_struct * nodes;

    T->length = 0;
    T->node_alloc = 2;
    T->nodes = FLINT_ARRAY_ALLOC(T->node_alloc, mpoly_rbnode_fmpz_struct);
    T->data = NULL;
    T->data_alloc = 0;
    T->data_size = data_size;

    nodes = T->nodes + 2;

    /* pointer to head */
    nodes[-1].up = -2;
    nodes[-1].left = -2;    /* head */
    nodes[-1].right = -2;
    fmpz_init(nodes[-1].key);
    nodes[-1].color = 0;

    /* null node */
    nodes[-2].up = -2;
    nodes[-2].left = -2;
    nodes[-2].right = -2;
    fmpz_init(nodes[-2].key);
    nodes[-2].color = 0;
}

/* the memory management of the data chunks must be done prior */
void mpoly_rbtree_fmpz_clear(mpoly_rbtree_fmpz_t T)
{
    slong i;

    FLINT_ASSERT(T->node_alloc > 0);

    for (i = 0; i < T->node_alloc; i++)
        fmpz_clear(T->nodes[i].key);

    flint_free(T->nodes);
    flint_free(T->data);
}

static void mpoly_rbtree_fmpz_fit_length(mpoly_rbtree_fmpz_t T, slong len)
{
    slong i, dsize = T->data_size;

    if (len + 2 > T->node_alloc)
    {
        slong new_alloc = FLINT_MAX(len + 2, 2*T->node_alloc);
        T->nodes = FLINT_ARRAY_REALLOC(T->nodes, new_alloc, mpoly_rbnode_fmpz_struct);

        for (i = T->node_alloc; i < new_alloc; i++)
            fmpz_init(T->nodes[i].key);

        T->node_alloc = new_alloc;
    }

    if (dsize*len > T->data_alloc)
    {
        slong new_alloc = FLINT_MAX(dsize*len, 2*T->data_alloc);
        T->data = FLINT_ARRAY_REALLOC(T->data, new_alloc, char);
        T->data_alloc = new_alloc;
    }
}

/*
    If the key rcx exists, a pointer to its data will be returned and new
    will be set to 0. Otherwise, a pointer to a new data chunk is returned
    and new is set to 1.
*/
void * mpoly_rbtree_fmpz_lookup(mpoly_rbtree_fmpz_t T, int * new, const fmpz_t rcx)
{
    slong dsize = T->data_size;
    mpoly_rbnode_fmpz_struct * nodes = T->nodes + 2;
    slong rax, rdx, r8, r9, r10, r11;
    slong n = T->length;
    int cmp;

    r10 = nodes[-1].left;

    if (n < 1)
    {
        mpoly_rbtree_fmpz_fit_length(T, 1);
        nodes = T->nodes + 2;
        rax = 0;
        nodes[rax].up = -1;
        nodes[rax].left = -2;
        nodes[rax].right = -2;
        nodes[rax].color = 0;
        fmpz_set(nodes[rax].key, rcx);
        T->length = 1;
        *new = 1;
        nodes[-1].left = rax;
        return T->data + dsize*rax;
    }

Compare:
    FLINT_ASSERT(r10 >= 0);
    r8 = nodes[r10].left;
    r9 = nodes[r10].right;
    cmp = fmpz_cmp(rcx, nodes[r10].key);
    if (cmp < 0)
        goto GoLeft;
    if (cmp > 0)
        goto GoRight;

    rax = r10;
    * new = 0;
    return T->data + dsize*rax;

GoLeft:
    if (r8 < 0)
    {
        FLINT_ASSERT(r8 == -2);
        goto MakeNewLeft;
    }
    r10 = r8;
    goto Compare;

GoRight:
    if (r9 < 0)
    {
        FLINT_ASSERT(r9 == -2);
        goto MakeNewRight;
    }
    r10 = r9;
    goto Compare;

MakeNewLeft:
    mpoly_rbtree_fmpz_fit_length(T, n + 1);
    nodes = T->nodes + 2;
    rdx = n;
    nodes[r10].left = rdx;
    goto FixTree;

MakeNewRight:
    mpoly_rbtree_fmpz_fit_length(T, n + 1);
    nodes = T->nodes + 2;
    rdx = n;
    nodes[r10].right = rdx;

FixTree:
    nodes[rdx].up = r10;
    nodes[rdx].left = -2;
    nodes[rdx].right = -2;
    nodes[rdx].color = 1;
    fmpz_set(nodes[rdx].key, rcx);
    T->length = n + 1;
    *new = 1;
    rax = rdx;

FixNode:

/*Case1:*/
    r8 = nodes[rdx].up;
    if (r8 < 0)
    {
        FLINT_ASSERT(r8 == -1);
        FLINT_ASSERT(nodes[-1].left == rdx);
        nodes[rdx].color = 0;
        return T->data + dsize*rax;
    }

/*Case2:*/
    if (nodes[r8].color == 0)
        return T->data + dsize*rax;

/*Case3:*/
    r9 = nodes[r8].up;
    r10 = nodes[r9].left;
    r11 = nodes[r9].right;
    if (r8 == r10)
        r10 = r11;
    if (r10 < 0 || nodes[r10].color == 0)
        goto Case4;
    rdx = r9;
    nodes[r8].color = 0;
    nodes[r9].color = 1;
    nodes[r10].color = 0;
    goto FixNode;

Case4:
    r10 = nodes[r9].up;
/*Case4A:*/
    if (rdx != nodes[r8].right || r8 != nodes[r9].left)
        goto Case4B;
    r11 = nodes[rdx].left;
    nodes[r9].left = rdx;
    nodes[rdx].left = r8;
    nodes[r8].right = r11;
    nodes[r8].up = rdx;
    nodes[rdx].up = r9;
    nodes[r11].up = r8;
    goto Case4Done;

Case4B:
    if (rdx != nodes[r8].left || r8 != nodes[r9].right)
        goto Case5;
    r11 = nodes[rdx].right;
    nodes[r9].right = rdx;
    nodes[rdx].right = r8;
    nodes[r8].left = r11;
    nodes[r8].up = rdx;
    nodes[rdx].up = r9;
    nodes[r11].up = r8;

Case4Done:
    FLINT_SWAP(slong, rdx, r8);

Case5:

    if (nodes[r10].right == r9)
        nodes[r10].right = r8;
    if (nodes[r10].left == r9)
        nodes[r10].left = r8;

    nodes[r8].up = r10;
    nodes[r8].color = 0;
    nodes[r9].up = r8;
    nodes[r9].color = 1;

    r11 = nodes[r8].right;
    r10 = nodes[r8].left;
    if (rdx == r10)
    {
        nodes[r8].right = r9;
        nodes[r9].left = r11;
        nodes[r11].up = r9;
    }
    else
    {
        nodes[r8].left = r9;
        nodes[r9].right = r10;
        nodes[r10].up = r9;
    }

    return T->data + dsize*rax;
}

