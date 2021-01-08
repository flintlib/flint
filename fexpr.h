/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FEXPR_H
#define FEXPR_H

#ifdef FEXPR_INLINES_C
#define FEXPR_INLINE
#else
#define FEXPR_INLINE static __inline__
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "calcium.h"

#define FEXPR_TYPE_SMALL_INT     UWORD(0)
#define FEXPR_TYPE_SMALL_SYMBOL  UWORD(1)
#define FEXPR_TYPE_SMALL_STRING  UWORD(2)
#define FEXPR_TYPE_BIG_INT_POS   UWORD(3)
#define FEXPR_TYPE_BIG_INT_NEG   UWORD(4)
#define FEXPR_TYPE_BIG_SYMBOL    UWORD(5)
#define FEXPR_TYPE_BIG_STRING    UWORD(6)
#define FEXPR_TYPE_CALL0         UWORD(7)
#define FEXPR_TYPE_CALL1         UWORD(8)
#define FEXPR_TYPE_CALL2         UWORD(9)
#define FEXPR_TYPE_CALL3         UWORD(10)
#define FEXPR_TYPE_CALL4         UWORD(11)
#define FEXPR_TYPE_CALLN         UWORD(12)

#define FEXPR_TYPE_BITS 4
#define FEXPR_TYPE_MASK ((UWORD(1) << FEXPR_TYPE_BITS) - 1)

#define FEXPR_COEFF_MAX ((WORD(1) << (FLINT_BITS - FEXPR_TYPE_BITS - 1)) - 1)
#define FEXPR_COEFF_MIN (-FEXPR_COEFF_MAX)

#define FEXPR_TYPE(head)  ((head) & FEXPR_TYPE_MASK)
#define FEXPR_SIZE(head)  ((slong) ((FEXPR_TYPE(head) <= FEXPR_TYPE_SMALL_STRING) ? 1 : (head) >> FEXPR_TYPE_BITS))

#define FEXPR_HEADER_SIZE        WORD(1)
#define FEXPR_SMALL_SYMBOL_LEN   ((FLINT_BITS / 8) - 1)

#define _FEXPR_SYMBOL_3(x,y,z) (FEXPR_TYPE_SMALL_SYMBOL | ((ulong)(x) << 8) | ((ulong)(y) << 16) | ((ulong)(z) << 24))

#define FEXPR_SYMBOL_Neg  _FEXPR_SYMBOL_3('N', 'e', 'g')
#define FEXPR_SYMBOL_Add  _FEXPR_SYMBOL_3('A', 'd', 'd')
#define FEXPR_SYMBOL_Sub  _FEXPR_SYMBOL_3('S', 'u', 'b')
#define FEXPR_SYMBOL_Mul  _FEXPR_SYMBOL_3('M', 'u', 'l')
#define FEXPR_SYMBOL_Div  _FEXPR_SYMBOL_3('D', 'i', 'v')
#define FEXPR_SYMBOL_Pow  _FEXPR_SYMBOL_3('P', 'o', 'w')


typedef struct
{
    ulong * data;
    slong alloc;
}
fexpr_struct;

typedef fexpr_struct fexpr_t[1];
typedef fexpr_struct * fexpr_ptr;
typedef const fexpr_struct * fexpr_srcptr;

FEXPR_INLINE void
fexpr_init(fexpr_t expr)
{
    expr->data = flint_malloc(sizeof(ulong));
    expr->data[0] = 0;
    expr->alloc = 1;
}

FEXPR_INLINE void
fexpr_clear(fexpr_t expr)
{
    flint_free(expr->data);
}

FEXPR_INLINE fexpr_ptr
_fexpr_vec_init(slong len)
{
    slong i;
    fexpr_ptr vec = flint_malloc(sizeof(fexpr_struct) * len);
    for (i = 0; i < len; i++)
        fexpr_init(vec + i);
    return vec;
}

FEXPR_INLINE void
_fexpr_vec_clear(fexpr_ptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        fexpr_clear(vec + i);
    flint_free(vec);
}

FEXPR_INLINE void
fexpr_fit_size(fexpr_t expr, slong size)
{
    if (expr->alloc < size)
    {
        size = FLINT_MAX(size, 2 * expr->alloc);
        expr->data = flint_realloc(expr->data, size * sizeof(ulong));
        expr->alloc = size;
    }
}

FEXPR_INLINE slong
_fexpr_size(const ulong * expr)
{
    ulong head = expr[0];
    return FEXPR_SIZE(head);
}

FEXPR_INLINE slong
fexpr_size(const fexpr_t expr)
{
    return _fexpr_size(expr->data);
}

FEXPR_INLINE void
fexpr_set(fexpr_t res, const fexpr_t expr)
{
    if (res != expr)
    {
        slong size = fexpr_size(expr);
        fexpr_fit_size(res, size);
        flint_mpn_copyi(res->data, expr->data, size);
    }
}

FEXPR_INLINE void
fexpr_swap(fexpr_t a, fexpr_t b)
{
    fexpr_struct tmp = *a;
    *a = *b;
    *b = tmp;
}


FEXPR_INLINE int
_mpn_equal(mp_srcptr a, mp_srcptr b, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (a[i] != b[i])
            return 0;

    return 1;
}

FEXPR_INLINE int
fexpr_equal(const fexpr_t a, const fexpr_t b)
{
    ulong ha, hb;
    slong sa, sb;

    ha = a->data[0];
    hb = b->data[0];

    if (ha != hb)
        return 0;

    sa = FEXPR_SIZE(ha);
    sb = FEXPR_SIZE(hb);

    if (sa != sb)
        return 0;

    return _mpn_equal(a->data + 1, b->data + 1, sa - 1);
}

FEXPR_INLINE int
_fexpr_is_integer(const ulong * expr)
{
    ulong type = FEXPR_TYPE(expr[0]);
    return (type == FEXPR_TYPE_SMALL_INT) || (type == FEXPR_TYPE_BIG_INT_POS) || (type == FEXPR_TYPE_BIG_INT_NEG);
}

FEXPR_INLINE int
fexpr_is_integer(const fexpr_t expr)
{
    return _fexpr_is_integer(expr->data);
}

FEXPR_INLINE int
_fexpr_is_symbol(const ulong * expr)
{
    ulong type = FEXPR_TYPE(expr[0]);
    return (type == FEXPR_TYPE_SMALL_SYMBOL) || (type == FEXPR_TYPE_BIG_SYMBOL);
}

FEXPR_INLINE int
fexpr_is_symbol(const fexpr_t expr)
{
    return _fexpr_is_symbol(expr->data);
}

FEXPR_INLINE int
_fexpr_is_string(const ulong * expr)
{
    ulong type = FEXPR_TYPE(expr[0]);
    return (type == FEXPR_TYPE_SMALL_STRING) || (type == FEXPR_TYPE_BIG_STRING);
}

FEXPR_INLINE int
fexpr_is_string(const fexpr_t expr)
{
    return _fexpr_is_string(expr->data);
}

FEXPR_INLINE int
_fexpr_is_atom(const ulong * expr)
{
    return FEXPR_TYPE(expr[0]) <= FEXPR_TYPE_BIG_STRING;
}

FEXPR_INLINE int
fexpr_is_atom(const fexpr_t expr)
{
    return _fexpr_is_atom(expr->data);
}

void fexpr_set_si(fexpr_t res, slong c);
void fexpr_set_fmpz(fexpr_t res, const fmpz_t c);
void fexpr_get_fmpz(fmpz_t c, const fexpr_t x);
void fexpr_set_fmpq(fexpr_t res, const fmpq_t x);

void fexpr_set_symbol_str(fexpr_t res, const char * s);
char * fexpr_get_symbol_str(const fexpr_t expr);

FEXPR_INLINE slong
fexpr_nargs(const fexpr_t expr)
{
    ulong type = FEXPR_TYPE(expr->data[0]);

    if (FEXPR_TYPE_CALL0 <= type && type <= FEXPR_TYPE_CALL4)
    {
        return type - FEXPR_TYPE_CALL0;
    }
    else if (type == FEXPR_TYPE_CALLN)
    {
        return expr->data[1];
    }
    else
    {
        return -1;
    }
}

void fexpr_func(fexpr_t res, const fexpr_t expr);
void fexpr_view_func(fexpr_t res, const fexpr_t expr);

void fexpr_arg(fexpr_t res, const fexpr_t expr, slong i);
void fexpr_view_arg(fexpr_t res, const fexpr_t expr, slong i);

FEXPR_INLINE void
fexpr_view_next(fexpr_t view)
{
    view->data += fexpr_size(view);
}

/* todo: handle aliasing! */
void fexpr_call0(fexpr_t res, const fexpr_t f);
void fexpr_call1(fexpr_t res, const fexpr_t f, const fexpr_t x1);
void fexpr_call2(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2);
void fexpr_call3(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3);
void fexpr_call4(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3, const fexpr_t x4);

void fexpr_call_vec(fexpr_t res, const fexpr_t f, fexpr_srcptr args, slong len);

void fexpr_print(const fexpr_t expr);

FEXPR_INLINE void
fexpr_neg(fexpr_t res, const fexpr_t a)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Neg;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call1(res, tmp, a);
}

FEXPR_INLINE void
fexpr_add(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Add;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

FEXPR_INLINE void
fexpr_sub(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Sub;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

FEXPR_INLINE void
fexpr_mul(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Mul;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

FEXPR_INLINE void
fexpr_div(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Div;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

FEXPR_INLINE void
fexpr_pow(fexpr_t res, const fexpr_t a, const fexpr_t b)
{
    fexpr_t tmp;
    ulong tmp_head = FEXPR_SYMBOL_Pow;
    tmp->data = &tmp_head;
    tmp->alloc = 1;
    fexpr_call2(res, tmp, a, b);
}

#ifdef __cplusplus
}
#endif

#endif

