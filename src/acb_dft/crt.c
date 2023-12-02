/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "acb_dft.h"

void
crt_init(crt_t c, ulong n)
{
    int k;
    n_factor_t fac;

    n_factor_init(&fac);
    if (n)
        n_factor(&fac, n, 1);
    else
        fac.num = 0;

    nmod_init(&c->n, FLINT_MAX(n, 1));

    c->num = fac.num;
    for (k = 0; k < fac.num; k++)
    {
        c->m[k] = n_pow(fac.p[k], fac.exp[k]);
        c->M[k] = n / c->m[k];
        c->vM[k] = c->M[k] * n_invmod(c->M[k] % c->m[k], c->m[k]);
        /*
        flint_printf("m[%ld]=%ld, M[%ld]=%wu, vM[%ld]=%wu\n",
                k, c->m[k], k, c->M[k], k, c->vM[k]);
                */
    }
}

void
crt_print(const crt_t c)
{
    slong k;
    if (c->num == 0)
    {
        flint_throw(FLINT_ERROR, "trivial group\n");
    }
    for (k = 0; k < c->num; k++)
        flint_printf("Z/%wuZ ", c->m[k]);
    flint_printf("\n");
}

#if 0
/* lexicographic index of crt elt j */
static ulong
index_crt(crt_t c, ulong j)
{
    int k;
    ulong res = 0;
    for (k = 0; k < c->num; k ++)
        res = res * c->m[k] + (j % c->m[k]);
    return res;
}

/* crt elt of lexicographic index i */
static ulong
crt_index(crt_t c, ulong i)
{
    int k;
    ulong j, res = 0;
    for (k = 0; k < c->num; k ++)
    {
        j = i % c->m[k];
        i = i / c->m[k];
        res = nmod_add(res,  j * c->vM[k], c->n);
    }
    return res;
}
/* for all elements can do fast conrey-like loop just adding vM[k] */
static acb_ptr
reorder_to_crt(acb_srcptr v, crt_t c, ulong len)
{
    ulong k;
    acb_ptr res;
    res = flint_malloc(len * sizeof(acb_struct));
    for (k = 0; k < len; k++)
        res[index_crt(c, k)] = v[k];
    return res;
}
static acb_ptr
reorder_from_crt(acb_srcptr v, crt_t c, ulong len)
{
    ulong k;
    acb_ptr res;
    res = flint_malloc(len * sizeof(acb_struct));
    for (k = 0; k < len; k++)
        res[k] = v[index_crt(c, k)];
    return res;
}
#endif

void
crt_decomp(acb_ptr y, acb_srcptr x, slong dv, const crt_t c, ulong len)
{
    int j, e[CRT_MAX];
    ulong k, l;

    for (j = 0; j < c->num; j++)
        e[j] = 0;

    l = 0;
    for(k = 0; k < len; k++)
    {
        acb_set(y + k, x + l * dv);
        for (j = c->num - 1; j >= 0; e[j] = 0, j--)
        {
            e[j]++; l = nmod_add(l, c->vM[j], c->n);
            if (e[j] < c->m[j])
                break;
        }
    }
}

void
crt_recomp(acb_ptr y, acb_srcptr x, const crt_t c, ulong len)
{
    int j, e[CRT_MAX];
    ulong k, l;

    for (j = 0; j < c->num; j++)
        e[j] = 0;

    l = 0;
    for(k = 0; k < len; k++)
    {
        acb_set(y + l, x + k);
        for (j = c->num - 1; j >= 0; e[j] = 0, j--)
        {
            e[j]++; l = nmod_add(l, c->M[j], c->n);
            if (e[j] < c->m[j])
                break;
        }
    }
}

void
_acb_dft_crt_init(acb_dft_crt_t crt, slong dv, slong len, slong prec)
{
    crt->n = len;
    crt_init(crt->c, len);
    crt->dv = dv;
    crt->cyc = _acb_dft_steps_prod(crt->c->m, crt->c->num, prec);
}

void
acb_dft_crt_init(acb_dft_crt_t crt, slong len, slong prec)
{
    crt->n = len;
    crt_init(crt->c, len);
    crt->dv = 1;
    crt->cyc = _acb_dft_steps_prod(crt->c->m, crt->c->num, prec);
}

void
acb_dft_crt_clear(acb_dft_crt_t crt)
{
    slong i;
    for (i = 0; i < crt->c->num; i++)
        acb_dft_precomp_clear(crt->cyc[i].pre);
    flint_free(crt->cyc);
}

void
acb_dft_crt_precomp(acb_ptr w, acb_srcptr v, const acb_dft_crt_t crt, slong prec)
{
    if (crt->n <= 1)
    {
        if (crt->n == 1)
            acb_set(w, v);
    }
    else
    {
        acb_ptr t;
        t = _acb_vec_init(crt->n);
        if (w == v)
        {
            _acb_vec_set(t, v, crt->n);
            v = t;
        }
        crt_decomp(w, v, crt->dv, crt->c, crt->n);
        acb_dft_step(t, w, crt->cyc, crt->c->num, prec);
        crt_recomp(w, t, crt->c, crt->n);
        _acb_vec_clear(t, crt->n);
    }
}

void
acb_dft_crt(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    if (len <= 1)
    {
        if (len == 1)
            acb_set(w, v);
    }
    else
    {
        crt_t c;
        acb_ptr t;
        t = _acb_vec_init(len);
        if (w == v)
        {
            _acb_vec_set(t, v, len);
            v = t;
        }
        crt_init(c, len);
        crt_decomp(w, v, 1, c, len);
        acb_dft_prod(t, w, c->m, c->num, prec);
        crt_recomp(w, t, c, len);
        _acb_vec_clear(t, len);
    }
}
