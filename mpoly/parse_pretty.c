/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "mpoly.h"

#define PREC_LOWEST    0
#define PREC_PLUS      1
#define PREC_MINUS     1
#define PREC_TIMES     2
#define PREC_DIVIDES   2
#define PREC_UPLUS     3
#define PREC_UMINUS    3
#define PREC_POWER     4
#define PREC_HIGHEST 255

#define OP_TIMES    0
#define OP_PLUS     1
#define OP_MINUS    2
#define OP_DIVIDES  3
#define OP_LROUND   4
#define FIX_INFIX       0
#define FIX_PREFIX      1
#define FIX_POSTFIX     2
#define FIX_MATCHFIX    3

static int _is_op(slong a)
{
    return a >= 0;
}

static slong _op_make(slong name, slong fix, slong prec)
{
    return (prec << 10) + (fix << 8) + (name << 0);
}

static slong _op_prec(slong a)
{
    return (ulong)(a) >> 10;
}

static slong _op_fix(slong a)
{
    return ((ulong)(a) >> 8) & 3;
}

static slong _op_name(slong a)
{
    return a&255;
}

void fparse_init(
    fparse_t E,
    void (*init_fxn)(void *, const void *),
    slong sz,   /* size of the element struct */
    const void * ctx)
{
    slong i;

    E->ctx = ctx;
    E->sz = sz;
    E->init_fxn = init_fxn;

    E->stack_len = 0;
    E->stack_alloc = 20;
    E->stack = FLINT_ARRAY_ALLOC(E->stack_alloc, slong);

    E->estore_len = 0;
    E->estore_alloc = 10;
    E->estore = flint_malloc(E->estore_alloc*E->sz);

    for (i = 0; i < E->estore_alloc; i++)
        E->init_fxn(E->estore + i*E->sz, E->ctx);

    E->terminals_len = 0;
    E->terminals_alloc = 5;
    E->terminal_strings = FLINT_ARRAY_ALLOC(E->terminals_alloc, string_with_length_struct);
    E->terminal_values = FLINT_ARRAY_ALLOC(E->terminals_alloc*E->sz, char);
    for (i = 0; i < E->terminals_alloc; i++)
    {
        E->terminal_strings[i].str = NULL;
        E->terminal_strings[i].str_len = 0;
        E->init_fxn(E->terminal_values + E->sz*i, ctx);
    }
}

void fparse_clear(fparse_t E)
{
    slong i;

    flint_free(E->stack);

    for (i = 0; i < E->estore_alloc; i++)
        E->clear_fxn(E->estore + E->sz*i, E->ctx);
    flint_free(E->estore);

    for (i = 0; i < E->terminals_alloc; i++)
    {
        flint_free(E->terminal_strings[i].str);
        E->clear_fxn(E->terminal_values + E->sz*i, E->ctx);
    }
    flint_free(E->terminal_strings);
    flint_free(E->terminal_values);
}

#define PTR_SWAP(T, a, b)   \
    do {                    \
        T * _tmp_ = a;      \
        a = b;              \
        b = _tmp_;          \
    } while (0);

void fparse_add_terminal(fparse_t E, const char * s, const void * val)
{
    slong l, n = E->terminals_len;

    if (n + 1 > E->terminals_alloc)
    {
        slong i = E->terminals_alloc;
        slong new_alloc = FLINT_MAX(n + 1, i + i/2);

        E->terminal_strings = (string_with_length_struct *) flint_realloc(
                                            E->terminal_strings, new_alloc*
                                            sizeof(string_with_length_struct));

        E->terminal_values = (char *) flint_realloc(E->terminal_values,
                                                              E->sz*new_alloc);
        for ( ; i < new_alloc; i++)
        {
            E->terminal_strings[i].str = NULL;
            E->terminal_strings[i].str_len = 0;
            E->init_fxn(E->terminal_values + E->sz*i, E->ctx);
        }

        E->terminals_alloc = new_alloc;
    }

    l = strlen(s);
    E->terminal_strings[n].str_len = l;

    E->terminal_strings[n].str = (char *) flint_realloc(E->terminal_strings[n].str, l + 1);
    memcpy(E->terminal_strings[n].str, s, l + 1);

    E->set_fxn(E->terminal_values + E->sz*n, val, E->ctx);

    E->terminals_len = n + 1;

    while (n > 0 && E->terminal_strings[n-1].str_len < E->terminal_strings[n].str_len)
    {
        PTR_SWAP(char, E->terminal_strings[n-1].str, E->terminal_strings[n].str);
        SLONG_SWAP(E->terminal_strings[n-1].str_len, E->terminal_strings[n].str_len);
        E->swap_fxn(E->terminal_values + E->sz*(n-1), E->terminal_values + E->sz*n, E->ctx);
        n--;
    }
}

static int fparse_top_is_expr(const fparse_t E)
{
    return E->stack_len > 0 && !_is_op(E->stack[E->stack_len - 1]);
}

static void * fparse_top_expr(fparse_t E)
{
    FLINT_ASSERT(E->stack_len > 0);
    FLINT_ASSERT(E->stack[E->stack_len - 1] < 0);
    return E->estore + E->sz*(-1 - E->stack[E->stack_len - 1]);
}

static void fparse_push_op(fparse_t E, slong op)
{
    FLINT_ASSERT(_is_op(op));
    _slong_array_fit_length(&E->stack, &E->stack_alloc, E->stack_len + 1);
    E->stack[E->stack_len] = op;
    E->stack_len++;
}

/* if the top is not an expr, push the tmp, otherwise fail */
static int fparse_push_expr(fparse_t E)
{
    if (fparse_top_is_expr(E))
        return -1;

    if (E->estore_len + 1 > E->estore_alloc)
    {
        slong i = E->estore_alloc;
        slong new_alloc = FLINT_MAX(E->estore_len + 1, i + i/2);
        E->estore = flint_realloc(E->estore, new_alloc*E->sz);
        for ( ; i < new_alloc; i++)
            E->init_fxn(E->estore + E->sz*i, E->ctx);
        E->estore_alloc = new_alloc;
    }

    _slong_array_fit_length(&E->stack, &E->stack_alloc, E->stack_len + 1);
    E->stack[E->stack_len] = -1 - E->estore_len;
    E->stack_len++;
    E->swap_fxn(E->estore + E->sz*E->estore_len, E->tmp, E->ctx);
    E->estore_len++;
    return 0;
}

/* if the top is an expr, pop it, otherwise fail */
static int fparse_pop_expr(fparse_t E)
{
    if (!fparse_top_is_expr(E))
        return -1;

    E->swap_fxn(E->tmp, E->estore + E->sz*(-1 - E->stack[E->stack_len - 1]), E->ctx);
    E->estore_len--;
    E->stack_len--;
    return 0;
}

/* if the top is an operation op, pop it, otherwise fail */
static int fparse_pop_op(fparse_t E, slong op)
{
    slong n = E->stack_len - 1;

    if (n < 0 || !_is_op(E->stack[n]) || _op_name(E->stack[n]) != op)
        return -1;

    E->stack_len = n;
    return 0;
}

/* pop ops with precedence > prec */
static int fparse_pop_prec(fparse_t E, slong prec)
{
    slong n, n1, n2, n3, p, l1, l3;

    if (E->stack_len < 1)
        return -1;

again:

    n = E->stack_len;
    if (n < 2)
        return 0;

    n1 = E->stack[n-1];
    n2 = E->stack[n-2];

    if (_is_op(n1) || !_is_op(n2))
        return 0;

    n1 = -1-n1;

    p = _op_prec(n2);

    if (p < prec)
        return 0;

    if (_op_fix(n2) == FIX_INFIX)
    {
        n3 = E->stack[n-3];
        FLINT_ASSERT(!_is_op(n3));
        n3 = -1 - n3;

        FLINT_ASSERT(n1 == n3 + 1);

        if (_op_name(n2) == OP_TIMES)
        {
            E->mul_fxn(E->tmp, E->estore + E->sz*n3, E->estore + E->sz*n1, E->ctx);
            E->swap_fxn(E->estore + E->sz*n3, E->tmp, E->ctx);
            E->estore_len -= 1;
            E->stack_len -= 2;
        }
        else if (_op_name(n2) == OP_PLUS)
        {
            l1 = E->length_fxn(E->estore + E->sz*n1, E->ctx);
            l3 = E->length_fxn(E->estore + E->sz*n3, E->ctx);

        do_plus:

            if (l1 > l3)
            {
                SLONG_SWAP(l3, l1);
                E->swap_fxn(E->estore + E->sz*n3, E->estore + E->sz*n1, E->ctx);
            }

            if (p > prec || 2*l1 >= l3)
            {
                E->add_fxn(E->estore + E->sz*n3, E->estore + E->sz*n3,
                                                 E->estore + E->sz*n1, E->ctx);
                E->estore_len -= 1;
                E->stack_len -= 2;
            }
            else
            {
                return 0;
            }
        }
        else if (_op_name(n2) == OP_MINUS)
        {
            l1 = E->length_fxn(E->estore + E->sz*n1, E->ctx);
            l3 = E->length_fxn(E->estore + E->sz*n3, E->ctx);

            if (4*l1 >= l3 || 4*l3 >= l1)
            {
                E->sub_fxn(E->estore + E->sz*n3, E->estore + E->sz*n3,
                                                 E->estore + E->sz*n1, E->ctx);
                E->estore_len -= 1;
                E->stack_len -= 2;
            }
            else
            {
                E->neg_fxn(E->estore + E->sz*n1, E->estore + E->sz*n1, E->ctx);
                E->stack[n-2] = _op_make(OP_PLUS, FIX_INFIX, PREC_PLUS);
                goto do_plus;
            }
        }
        else if (_op_name(n2) == OP_DIVIDES)
        {
            /*
            NOTE: if divides and times have the same precedence and the
                  multiplications were to be delayed as the addition are,
                  then there would have to be more shenenigans here.
            */
            if (!E->div_fxn(E->tmp, E->estore + E->sz*n3,
                                    E->estore + E->sz*n1, E->ctx))
                return -1;

            E->swap_fxn(E->estore + E->sz*n3, E->tmp, E->ctx);
            E->estore_len -= 1;
            E->stack_len -= 2;
        }
        else
        {
            flint_throw(FLINT_ERROR, "_pop_stack: internal error");
        }

        goto again;
    }
    else if (_op_fix(n2) == FIX_PREFIX)
    {
        if (_op_name(n2) == OP_MINUS)
            E->neg_fxn(E->estore + E->sz*n1, E->estore + E->sz*n1, E->ctx);

        E->stack[n-2] = -1-n1;
        E->stack_len -= 1;
        goto again;
    }
    else
    {
        return 0;
    }
}

static const char * _parse_int(fmpz_t c, const char * s, const char * end)
{
    char * buffer, * v;
    const char * send = s + 1;
    TMP_INIT;
    while (send < end && '0' <= *send && *send <= '9')
        send++;
    TMP_START;
    v = buffer = (char *) TMP_ALLOC((send - s + 1)*sizeof(char));
    while (s < send)
        *v++ = *s++;
    *v++ = '\0';
    fmpz_set_str(c, buffer, 10);
    TMP_END;
    return s;
}

int fparse_parse(fparse_t E, void * poly, const char * s, slong slen)
{
    const char * send = s + slen;
    fmpz_t c;
    int ret;

    fmpz_init(c);
    E->tmp = poly;

    while (s < send)
    {
        if ('0' <= *s && *s <= '9')
        {
            s = _parse_int(c, s, send);

            E->set_fmpz_fxn(E->tmp, c, E->ctx);
            if (fparse_push_expr(E))
                goto failed;
        }
        else if (*s == '^')
        {
            if (++s >= send || !('0' <= *s && *s <= '9'))
                goto failed;

            s = _parse_int(c, s, send);
            
            if (fparse_pop_prec(E, PREC_POWER))
                goto failed;

            if (!fparse_top_is_expr(E))
                goto failed;

            if (!E->pow_fmpz_fxn(fparse_top_expr(E), fparse_top_expr(E), c, E->ctx))
                goto failed;
        }
        else if (*s == '*')
        {
            if (!fparse_top_is_expr(E))
                goto failed;

            if (fparse_pop_prec(E, PREC_TIMES))
                goto failed;

            fparse_push_op(E, _op_make(OP_TIMES, FIX_INFIX, PREC_TIMES));
            s++;
        }
        else if (*s == '+')
        {
            if (!fparse_top_is_expr(E))
            {
                fparse_push_op(E, _op_make(OP_PLUS, FIX_PREFIX, PREC_UPLUS));
            }
            else
            {
                if (fparse_pop_prec(E, PREC_PLUS))
                    goto failed;

                fparse_push_op(E, _op_make(OP_PLUS, FIX_INFIX, PREC_PLUS));
            }
            s++;
        }
        else if (*s == '-')
        {
            if (!fparse_top_is_expr(E))
            {
                fparse_push_op(E, _op_make(OP_MINUS, FIX_PREFIX, PREC_UMINUS));
            }
            else
            {
                if (fparse_pop_prec(E, PREC_MINUS))
                    goto failed;

                fparse_push_op(E, _op_make(OP_MINUS, FIX_INFIX, PREC_MINUS));
            }
            s++;
        }
        else if (*s == '/')
        {
            if (!fparse_top_is_expr(E))
                goto failed;

            if (fparse_pop_prec(E, PREC_DIVIDES))
                goto failed;

            fparse_push_op(E, _op_make(OP_DIVIDES, FIX_INFIX, PREC_DIVIDES));
            s++;
        }
        else if (*s == ' ')
        {
            s++;
        }
        else if (*s == '(')
        {
            if (fparse_top_is_expr(E))
                goto failed;

            fparse_push_op(E, _op_make(OP_LROUND, FIX_MATCHFIX, PREC_LOWEST));
            s++;
        }
        else if (*s == ')')
        {
            if (fparse_pop_prec(E, PREC_LOWEST))
                goto failed;

            if (fparse_pop_expr(E))
                goto failed;

            if (fparse_pop_op(E, OP_LROUND))
                goto failed;

            if (fparse_push_expr(E))
                goto failed;

            s++;
        }
        else
        {
            slong k;
            for (k = 0; k < E->terminals_len; k++)
            {
                slong l = E->terminal_strings[k].str_len;
                if (0 == strncmp(s, E->terminal_strings[k].str, l))
                {
                    E->set_fxn(E->tmp, E->terminal_values + E->sz*k, E->ctx);
                    if (fparse_push_expr(E))
                        goto failed;

                    s += l;
                    goto continue_outer;
                }
            }

            goto failed;
        }
    continue_outer:;
    }

    if (fparse_pop_prec(E, PREC_LOWEST))
        goto failed;

    if (fparse_pop_expr(E))
        goto failed;

    if (E->stack_len != 0)
        goto failed;

    ret = 0;

done:

    fmpz_clear(c);
    return ret;

failed:

    ret = -1;
    goto done;
}
