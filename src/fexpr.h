/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FEXPR_H
#define FEXPR_H

#ifdef FEXPR_INLINES_C
#define FEXPR_INLINE
#else
#define FEXPR_INLINE static inline
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "mpoly_types.h"
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

#define FEXPR_BUILTIN_ID(head) ((head) >> 16)

typedef struct
{
    ulong * data;
    slong alloc;
}
fexpr_struct;

typedef fexpr_struct fexpr_t[1];
typedef fexpr_struct * fexpr_ptr;
typedef const fexpr_struct * fexpr_srcptr;

typedef struct
{
    fexpr_struct * entries;
    slong alloc;
    slong length;
}
fexpr_vec_struct;

typedef fexpr_vec_struct fexpr_vec_t[1];

#define fexpr_vec_entry(vec, i) ((vec)->entries + (i))

FEXPR_INLINE void
fexpr_init(fexpr_t expr)
{
    expr->data = (ulong *) flint_malloc(sizeof(ulong));
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
    fexpr_ptr vec = (fexpr_ptr) flint_malloc(sizeof(fexpr_struct) * len);
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
        expr->data = (ulong *) flint_realloc(expr->data, size * sizeof(ulong));
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
    FLINT_SWAP(fexpr_struct, *a, *b);
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

/* todo: document, test */
int fexpr_equal_si(const fexpr_t expr, slong c);
int fexpr_equal_ui(const fexpr_t expr, ulong c);
ulong fexpr_hash(const fexpr_t expr);

int fexpr_cmp_fast(const fexpr_t a, const fexpr_t b);

void
_fexpr_vec_sort_fast(fexpr_ptr vec, slong len);

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

/* todo: document, test */
int fexpr_is_neg_integer(const fexpr_t expr);

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

FEXPR_INLINE void
fexpr_zero(fexpr_t res)
{
    res->data[0] = 0;
}

FEXPR_INLINE int
fexpr_is_zero(const fexpr_t expr)
{
    return expr->data[0] == 0;
}

void fexpr_set_si(fexpr_t res, slong c);
void fexpr_set_ui(fexpr_t res, ulong c);
void fexpr_set_fmpz(fexpr_t res, const fmpz_t c);
int fexpr_get_fmpz(fmpz_t c, const fexpr_t x);
void fexpr_set_fmpq(fexpr_t res, const fmpq_t x);

FEXPR_INLINE void
fexpr_set_symbol_builtin(fexpr_t res, slong id)
{
    res->data[0] = FEXPR_TYPE_SMALL_SYMBOL | (id << 16);
}

void fexpr_set_symbol_str(fexpr_t res, const char * s);
char * fexpr_get_symbol_str(const fexpr_t expr);

void fexpr_set_string(fexpr_t res, const char * s);
char * fexpr_get_string(const fexpr_t expr);

slong fexpr_depth(const fexpr_t expr);
slong fexpr_num_leaves(const fexpr_t expr);

FEXPR_INLINE slong
fexpr_size_bytes(const fexpr_t expr)
{
    return fexpr_size(expr) * sizeof(ulong);
}

FEXPR_INLINE slong
fexpr_allocated_bytes(const fexpr_t expr)
{
    return expr->alloc * sizeof(ulong);
}

/* todo: document */
FEXPR_INLINE int
fexpr_is_any_builtin_symbol(const fexpr_t expr)
{
    ulong head;
    head = expr->data[0];
    return (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_SYMBOL) && (((head >> 8) & 0xff) == 0);
}

/* todo: document */
FEXPR_INLINE int
fexpr_is_builtin_symbol(const fexpr_t expr, slong i)
{
    ulong head;
    head = expr->data[0];
    return (FEXPR_TYPE(head) == FEXPR_TYPE_SMALL_SYMBOL) && (((head >> 8) & 0xff) == 0) && (FEXPR_BUILTIN_ID(head) == i);
}

/* todo: document, test */
int fexpr_is_builtin_call(const fexpr_t expr, slong i);
int fexpr_is_any_builtin_call(const fexpr_t expr);

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

void fexpr_call_builtin1(fexpr_t res, slong f, const fexpr_t x);
void fexpr_call_builtin2(fexpr_t res, slong f, const fexpr_t x, const fexpr_t y);

/* Subexpressions and replacement */

int fexpr_contains(const fexpr_t expr, const fexpr_t x);
int fexpr_replace(fexpr_t res, const fexpr_t expr, const fexpr_t x, const fexpr_t y);
int fexpr_replace2(fexpr_t res, const fexpr_t expr, const fexpr_t x1, const fexpr_t y1, const fexpr_t x2, const fexpr_t y2);
int fexpr_replace_vec(fexpr_t res, const fexpr_t expr, const fexpr_vec_t xs, const fexpr_vec_t ys);


/* Input/output */

void fexpr_write(calcium_stream_t stream, const fexpr_t expr);
void fexpr_print(const fexpr_t expr);
char * fexpr_get_str(const fexpr_t expr);

/* LaTeX output */

#define FEXPR_LATEX_SMALL 1
#define FEXPR_LATEX_LOGIC 2

void fexpr_write_latex(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_print_latex(const fexpr_t expr, ulong flags);
char * fexpr_get_str_latex(const fexpr_t expr, ulong flags);

void fexpr_write_latex_call(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_subscript(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_subscript_call(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_infix(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_mul(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_div(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_neg_pos(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_add(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_sub(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_pow(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_exp(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_factorial(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_integral(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_sum_product(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_divsum(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_limit(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_logic(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_collection(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_matrix(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_simple(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_simple2(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_simple2_small(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_alg_structure(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_setop(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_cases(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_where(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_show_form(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_range(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_decimal(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_call1_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_call2_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_sub1_call1_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_sub1_call2_optional_derivative(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_misc_special(calcium_stream_t out, const fexpr_t expr, ulong flags);
void fexpr_write_latex_residue(calcium_stream_t out, const fexpr_t expr, ulong flags);

void fexpr_set_arf(fexpr_t res, const arf_t x);
void fexpr_set_d(fexpr_t res, double x);
void fexpr_set_re_im_d(fexpr_t res, double x, double y);

void fexpr_neg(fexpr_t res, const fexpr_t a);
void fexpr_add(fexpr_t res, const fexpr_t a, const fexpr_t b);
void fexpr_sub(fexpr_t res, const fexpr_t a, const fexpr_t b);
void fexpr_mul(fexpr_t res, const fexpr_t a, const fexpr_t b);
void fexpr_div(fexpr_t res, const fexpr_t a, const fexpr_t b);
void fexpr_pow(fexpr_t res, const fexpr_t a, const fexpr_t b);

int fexpr_is_arithmetic_operation(const fexpr_t expr);

void fexpr_arithmetic_nodes(fexpr_vec_t nodes, const fexpr_t expr);

int fexpr_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fexpr_t expr, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx);

void fexpr_set_fmpz_mpoly(fexpr_t res, const fmpz_mpoly_t poly, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx);
void fexpr_set_fmpz_mpoly_q(fexpr_t res, const fmpz_mpoly_q_t frac, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx);

int fexpr_expanded_normal_form(fexpr_t res, const fexpr_t expr, ulong flags);

/* Vectors */

FEXPR_INLINE void
fexpr_vec_init(fexpr_vec_t vec, slong len)
{
    if (len == 0)
    {
        vec->entries = NULL;
        vec->length = 0;
        vec->alloc = 0;
    }
    else
    {
        slong i;
        vec->entries = (fexpr_ptr) flint_malloc(sizeof(fexpr_struct) * len);
        for (i = 0; i < len; i++)
            fexpr_init(vec->entries + i);
        vec->length = vec->alloc = len;
    }
}

FEXPR_INLINE void
fexpr_vec_print(const fexpr_vec_t F)
{
    slong i;

    flint_printf("[");
    for (i = 0; i < F->length; i++)
    {
        fexpr_print(F->entries + i);
        if (i < F->length - 1)
            flint_printf(", ");
    }
    flint_printf("]");
}

FEXPR_INLINE void
fexpr_vec_swap(fexpr_vec_t x, fexpr_vec_t y)
{
    fexpr_vec_t tmp;
    *tmp = *x;
    *x = *y;
    *y = *tmp;
}

FEXPR_INLINE void
fexpr_vec_fit_length(fexpr_vec_t vec, slong len)
{
    if (len > vec->alloc)
    {
        slong i;

        if (len < 2 * vec->alloc)
            len = 2 * vec->alloc;

        vec->entries = (fexpr_struct *) flint_realloc(vec->entries, len * sizeof(fexpr_struct));

        for (i = vec->alloc; i < len; i++)
            fexpr_init(vec->entries + i);

        vec->alloc = len;
    }
}

FEXPR_INLINE void
fexpr_vec_clear(fexpr_vec_t vec)
{
    slong i;

    for (i = 0; i < vec->alloc; i++)
        fexpr_clear(vec->entries + i);

    flint_free(vec->entries);
}

FEXPR_INLINE void
fexpr_vec_set(fexpr_vec_t dest, const fexpr_vec_t src)
{
    if (dest != src)
    {
        slong i;

        fexpr_vec_fit_length(dest, src->length);

        for (i = 0; i < src->length; i++)
            fexpr_set(dest->entries + i, src->entries + i);

        dest->length = src->length;
    }
}

FEXPR_INLINE void
fexpr_vec_append(fexpr_vec_t vec, const fexpr_t f)
{
    fexpr_vec_fit_length(vec, vec->length + 1);
    fexpr_set(vec->entries + vec->length, f);
    vec->length++;
}

FEXPR_INLINE slong
fexpr_vec_insert_unique(fexpr_vec_t vec, const fexpr_t f)
{
    slong i;

    for (i = 0; i < vec->length; i++)
    {
        if (fexpr_equal(vec->entries + i, f))
            return i;
    }

    fexpr_vec_append(vec, f);
    return vec->length - 1;
}

FEXPR_INLINE void
fexpr_vec_set_length(fexpr_vec_t vec, slong len)
{
    slong i;

    if (len > vec->length)
    {
        fexpr_vec_fit_length(vec, len);
        for (i = vec->length; i < len; i++)
            fexpr_zero(vec->entries + i);
    }
    else if (len < vec->length)
    {
        for (i = len; i < vec->length; i++)
           fexpr_zero(vec->entries + i);
    }

    vec->length = len;
}


#ifdef __cplusplus
}
#endif

#endif

