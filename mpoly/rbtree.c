/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "mpoly.h"


/* old bad red-black tree: the keys are slong and fmpz's are hacked in */

void mpoly_rbtree_init(mpoly_rbtree_t tree)
{
    tree->size = 0;
    tree->head->up = tree->null;
    tree->head->left = tree->null;
    tree->head->right = tree->null;
    tree->head->key = 0;
    tree->head->col = 0;
    tree->null->up = NULL;
    tree->null->left = NULL;
    tree->null->right = NULL;
    tree->null->key = 0;
    tree->null->col = 0;
}


void _mpoly_rbnode_clear(mpoly_rbtree_t tree, mpoly_rbnode_t node,
                                 void ** dataout, slong * keysout, slong * idx)
{
    if (node->right != tree->null)
        _mpoly_rbnode_clear(tree, node->right, dataout, keysout, idx);
    dataout[*idx] = node->data;
    keysout[*idx] = node->key;
    (*idx)++;
    if (node->left != tree->null)
        _mpoly_rbnode_clear(tree, node->left, dataout, keysout, idx);
    flint_free(node);
}


void mpoly_rbtree_clear(mpoly_rbtree_t tree, void ** dataout, slong * keysout)
{
    slong idx = 0;
    if (tree->head->left != tree->null)
        _mpoly_rbnode_clear(tree, tree->head->left, dataout, keysout, &idx);
}

/*
    get the node with key "rcx"
    if such a node doesn't exist, one is created and "new" is set to 1
*/
mpoly_rbnode_struct * mpoly_rbtree_get(int * new, struct mpoly_rbtree * tree,
                                                                     slong rcx)
{
    mpoly_rbnode_struct * t, * head, * null;
    mpoly_rbnode_struct * rax, * rdx, * r8, * r9, * r10, * r11;

    * new = 0;
    head = tree->head;
    null = tree->null;

    r10 = head->left;

    if (tree->size == 0)
    {
        rax = flint_malloc(sizeof(struct mpoly_rbnode));
        rax->up = head;
        rax->left = null;
        rax->right = null;
        rax->data = NULL;
        rax->col = 0;
        rax->key = rcx;
        tree->size++;
        * new = 1;
        head->left = rax;
        return rax;
    }

Compare:
    r8 = r10->left;
    r9 = r10->right;
    if (rcx < r10->key)
        goto GoLeft;
    if (rcx > r10->key)
        goto GoRight;

    rax = r10;
    return rax;

GoLeft:
    if (r8 == null)
        goto MakeNewLeft;
    r10 = r8;
    goto Compare;

GoRight:
    if (r9 == null)
        goto MakeNewRight;
    r10 = r9;
    goto Compare;

MakeNewLeft:
    rdx = flint_malloc(sizeof(mpoly_rbnode_struct));
    r10->left = rdx;
    goto FixTree;

MakeNewRight:
    rdx = flint_malloc(sizeof(mpoly_rbnode_struct));
    r10->right = rdx;

FixTree:
    rdx->up = r10;
    rdx->left = null;
    rdx->right = null;
    rdx->data = NULL;
    rdx->col = 1;
    rdx->key = rcx;
    tree->size++;
    * new = 1;
    rax = rdx;

FixNode:

/*Case1:*/
    r8 = rdx->up;
    if (r8 == head)
    {
        rdx->col = 0;
        return rax;
    }


/*Case2:*/
    if (r8->col == 0)
        return rax;

/*Case3:*/
    r9 = r8->up;
    r10 = r9->left;
    r11 = r9-> right;
    if (r8 == r10)
        r10 = r11;
    if (r10 == null || r10->col == 0)
        goto Case4;
    rdx = r9;
    r8->col = 0;
    r9->col = 1;
    r10->col = 0;
    goto FixNode;

Case4:
    r10 = r9->up;
/*Case4A:*/
    if (rdx != r8->right || r8 != r9->left)
        goto Case4B;
    r11 = rdx->left;
    r9->left = rdx;
    rdx->left = r8;
    r8->right = r11;
    r8->up = rdx;
    rdx->up = r9;
    r11->up = r8;
    goto Case4Done;

Case4B:
    if (rdx != r8->left || r8 != r9->right)
        goto Case5;
    r11 = rdx->right;
    r9->right = rdx;
    rdx->right = r8;
    r8->left = r11;
    r8->up = rdx;
    rdx->up = r9;
    r11->up = r8;

Case4Done:
    t = rdx;
    rdx = r8;
    r8 = t;

Case5:
    if (r10->right == r9)
        r10->right = r8;
    if (r10->left == r9)
        r10->left = r8;

    r8->up = r10;
    r8->col = 0;
    r9->up = r8;
    r9->col = 1;

    r11 = r8->right;
    r10 = r8->left;
    if (rdx == r10)
    {
        r8->right = r9;
        r9->left = r11;
        r11->up = r9;
    } else
    {
        r8->left = r9;
        r9->right = r10;
        r10->up = r9;
    }

    return rax;
}

/*
    get the node with key "rcx"
    if such a node doesn't exist, one is created and "new" is set to 1
*/
mpoly_rbnode_struct * mpoly_rbtree_get_fmpz(int * new,
                                         struct mpoly_rbtree * tree, fmpz_t rcx)
{
    int cmp;
    mpoly_rbnode_struct * t, * head, * null;
    mpoly_rbnode_struct * rax, * rdx, * r8, * r9, * r10, * r11;
    * new = 0;
    head = tree->head;
    null = tree->null;
    r10 = head->left;
    if (tree->size == 0)
    {
        rax = flint_malloc(sizeof(struct mpoly_rbnode));
        rax->up = head;
        rax->left = null;
        rax->right = null;
        rax->data = NULL;
        rax->col = 0;
        fmpz_init_set(&rax->key, rcx);
        tree->size++;
        * new = 1;
        head->left = rax;
        return rax;
    }
Compare:
    r8 = r10->left;
    r9 = r10->right;
    cmp = fmpz_cmp(rcx, &r10->key);
    if (cmp < 0)
        goto GoLeft;
    if (cmp > 0)
        goto GoRight;
    rax = r10;
    return rax;
GoLeft:
    if (r8 == null)
        goto MakeNewLeft;
    r10 = r8;
    goto Compare;
GoRight:
    if (r9 == null)
        goto MakeNewRight;
    r10 = r9;
    goto Compare;
MakeNewLeft:
    rdx = flint_malloc(sizeof(mpoly_rbnode_struct));
    r10->left = rdx;
    goto FixTree;
MakeNewRight:
    rdx = flint_malloc(sizeof(mpoly_rbnode_struct));
    r10->right = rdx;
FixTree:
    rdx->up = r10;
    rdx->left = null;
    rdx->right = null;
    rdx->data = NULL;
    rdx->col = 1;
    fmpz_init_set(&rdx->key, rcx);
    tree->size++;
    * new = 1;
    rax = rdx;
FixNode:
/*Case1:*/
    r8 = rdx->up;
    if (r8 == head)
    {
        rdx->col = 0;
        return rax;
    }
/*Case2:*/
    if (r8->col == 0)
        return rax;
/*Case3:*/
    r9 = r8->up;
    r10 = r9->left;
    r11 = r9-> right;
    if (r8 == r10)
        r10 = r11;
    if (r10 == null || r10->col == 0)
        goto Case4;
    rdx = r9;
    r8->col = 0;
    r9->col = 1;
    r10->col = 0;
    goto FixNode;
Case4:
    r10 = r9->up;
/*Case4A:*/
    if (rdx != r8->right || r8 != r9->left)
        goto Case4B;
    r11 = rdx->left;
    r9->left = rdx;
    rdx->left = r8;
    r8->right = r11;
    r8->up = rdx;
    rdx->up = r9;
    r11->up = r8;
    goto Case4Done;
Case4B:
    if (rdx != r8->left || r8 != r9->right)
        goto Case5;
    r11 = rdx->right;
    r9->right = rdx;
    rdx->right = r8;
    r8->left = r11;
    r8->up = rdx;
    rdx->up = r9;
    r11->up = r8;
Case4Done:
    t = rdx;
    rdx = r8;
    r8 = t;
Case5:
    if (r10->right == r9)
        r10->right = r8;
    if (r10->left == r9)
        r10->left = r8;
    r8->up = r10;
    r8->col = 0;
    r9->up = r8;
    r9->col = 1;
    r11 = r8->right;
    r10 = r8->left;
    if (rdx == r10)
    {
        r8->right = r9;
        r9->left = r11;
        r11->up = r9;
    } else
    {
        r8->left = r9;
        r9->right = r10;
        r10->up = r9;
    }
    return rax;
}


/* red-black tree with dedicated types for the keys */

void mpoly_rbtree_ui_init(mpoly_rbtree_ui_t T)
{
    mpoly_rbnode_ui_struct * nodes;

    T->length = 0;
    T->nodes = (mpoly_rbnode_ui_struct *) flint_malloc(2*sizeof(mpoly_rbnode_ui_struct));
    T->node_alloc = 2;
    T->data = NULL;
    T->data_alloc = 0;

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

    if (T->data_alloc > 0)
        flint_free(T->data);
}

static void mpoly_rbtree_ui_fit_length(
    mpoly_rbtree_ui_t T,
    slong len,
    slong dsize)
{
    if (len + 2 > T->node_alloc)
    {
        slong new_alloc = FLINT_MAX(len + 2, 2*T->node_alloc);
        if (T->node_alloc > 0)
            T->nodes = (mpoly_rbnode_ui_struct *) flint_realloc(T->nodes,
                                   new_alloc * sizeof(mpoly_rbnode_ui_struct));
        else
            T->nodes = (mpoly_rbnode_ui_struct *) flint_malloc(
                                   new_alloc * sizeof(mpoly_rbnode_ui_struct));
        T->node_alloc = new_alloc;
    }

    if (dsize*len > T->data_alloc)
    {
        slong new_alloc = FLINT_MAX(dsize*len, 2*T->data_alloc);
        if (T->data_alloc > 0)
            T->data = (char *) flint_realloc(T->data, new_alloc);
        else
            T->data = (char *) flint_malloc(new_alloc);
        T->data_alloc = new_alloc;
    }
}

/*
    Tf the key rcx exists, a pointer to its data will be returned and new
    will be set to 0. Otherwise, a pointer to a new data chunk is returned
    and new is set to 1.

    Pass in the size of the data chunks in dsize. TODO move this to init.
*/
void * mpoly_rbtree_ui_lookup(
    mpoly_rbtree_ui_t T,
    int * new,
    ulong rcx,
    slong dsize)
{
    mpoly_rbnode_ui_struct * nodes = T->nodes + 2;
    slong rax, rdx, r8, r9, r10, r11;
    slong n = T->length;
    r10 = nodes[-1].left;

    if (n < 1)
    {
        mpoly_rbtree_ui_fit_length(T, 1, dsize);
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
    mpoly_rbtree_ui_fit_length(T, n + 1, dsize);
    nodes = T->nodes + 2;
    rdx = n;
    nodes[r10].left = rdx;
    goto FixTree;

MakeNewRight:
    mpoly_rbtree_ui_fit_length(T, n + 1, dsize);
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
    SLONG_SWAP(rdx, r8);

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

