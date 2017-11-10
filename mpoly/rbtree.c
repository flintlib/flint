/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "mpoly.h"


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

