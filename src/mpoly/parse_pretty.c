/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

/* initialize the R member first */
void mpoly_parse_init(mpoly_parse_t E)
{
    slong i;

    E->stack_len = 0;
    E->stack_alloc = 20;
    E->stack = FLINT_ARRAY_ALLOC(E->stack_alloc, slong);

    E->estore_len = 0;
    E->estore_alloc = 10;
    E->estore = flint_malloc(E->estore_alloc*E->R->elem_size);

    for (i = 0; i < E->estore_alloc; i++)
        E->R->init(E->estore + i*E->R->elem_size, E->R->ctx);

    E->terminals_len = 0;
    E->terminals_alloc = 5;
    E->terminal_strings = FLINT_ARRAY_ALLOC(E->terminals_alloc, string_with_length_struct);
    E->terminal_values = FLINT_ARRAY_ALLOC(E->terminals_alloc*E->R->elem_size, char);
    for (i = 0; i < E->terminals_alloc; i++)
    {
        E->terminal_strings[i].str = NULL;
        E->terminal_strings[i].str_len = 0;
        E->R->init(E->terminal_values + E->R->elem_size*i, E->R->ctx);
    }
}

void mpoly_parse_clear(mpoly_parse_t E)
{
    slong i;

    flint_free(E->stack);

    for (i = 0; i < E->estore_alloc; i++)
        E->R->clear(E->estore + E->R->elem_size*i, E->R->ctx);
    flint_free(E->estore);

    for (i = 0; i < E->terminals_alloc; i++)
    {
        flint_free(E->terminal_strings[i].str);
        E->R->clear(E->terminal_values + E->R->elem_size*i, E->R->ctx);
    }
    flint_free(E->terminal_strings);
    flint_free(E->terminal_values);
}

void mpoly_parse_add_terminal(mpoly_parse_t E, const char * s, const void * val)
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
                                                              E->R->elem_size*new_alloc);
        for ( ; i < new_alloc; i++)
        {
            E->terminal_strings[i].str = NULL;
            E->terminal_strings[i].str_len = 0;
            E->R->init(E->terminal_values + E->R->elem_size*i, E->R->ctx);
        }

        E->terminals_alloc = new_alloc;
    }

    l = strlen(s);
    E->terminal_strings[n].str_len = l;

    E->terminal_strings[n].str = (char *) flint_realloc(E->terminal_strings[n].str, l + 1);
    memcpy(E->terminal_strings[n].str, s, l + 1);

    E->R->set(E->terminal_values + E->R->elem_size*n, val, E->R->ctx);

    E->terminals_len = n + 1;

    while (n > 0 && E->terminal_strings[n-1].str_len < E->terminal_strings[n].str_len)
    {
        FLINT_SWAP(char *, E->terminal_strings[n-1].str, E->terminal_strings[n].str);
        FLINT_SWAP(slong, E->terminal_strings[n-1].str_len, E->terminal_strings[n].str_len);
        E->R->swap(E->terminal_values + E->R->elem_size*(n-1), E->terminal_values + E->R->elem_size*n, E->R->ctx);
        n--;
    }
}

static int mpoly_parse_top_is_expr(const mpoly_parse_t E)
{
    return E->stack_len > 0 && !_is_op(E->stack[E->stack_len - 1]);
}

static void * mpoly_parse_top_expr(mpoly_parse_t E)
{
    FLINT_ASSERT(E->stack_len > 0);
    FLINT_ASSERT(E->stack[E->stack_len - 1] < 0);
    return E->estore + E->R->elem_size*(-1 - E->stack[E->stack_len - 1]);
}

static void mpoly_parse_push_op(mpoly_parse_t E, slong op)
{
    FLINT_ASSERT(_is_op(op));
    _slong_array_fit_length(&E->stack, &E->stack_alloc, E->stack_len + 1);
    E->stack[E->stack_len] = op;
    E->stack_len++;
}

/* if the top is not an expr, push the tmp, otherwise fail */
static int mpoly_parse_push_expr(mpoly_parse_t E)
{
    if (mpoly_parse_top_is_expr(E))
        return -1;

    if (E->estore_len + 1 > E->estore_alloc)
    {
        slong i = E->estore_alloc;
        slong new_alloc = FLINT_MAX(E->estore_len + 1, i + i/2);
        E->estore = flint_realloc(E->estore, new_alloc*E->R->elem_size);
        for ( ; i < new_alloc; i++)
            E->R->init(E->estore + E->R->elem_size*i, E->R->ctx);
        E->estore_alloc = new_alloc;
    }

    _slong_array_fit_length(&E->stack, &E->stack_alloc, E->stack_len + 1);
    E->stack[E->stack_len] = -1 - E->estore_len;
    E->stack_len++;
    E->R->swap(E->estore + E->R->elem_size*E->estore_len, E->tmp, E->R->ctx);
    E->estore_len++;
    return 0;
}

/* if the top is an expr, pop it, otherwise fail */
static int mpoly_parse_pop_expr(mpoly_parse_t E)
{
    if (!mpoly_parse_top_is_expr(E))
        return -1;

    E->R->swap(E->tmp, E->estore + E->R->elem_size*(-1 - E->stack[E->stack_len - 1]), E->R->ctx);
    E->estore_len--;
    E->stack_len--;
    return 0;
}

/* if the top is an operation op, pop it, otherwise fail */
static int mpoly_parse_pop_op(mpoly_parse_t E, slong op)
{
    slong n = E->stack_len - 1;

    if (n < 0 || !_is_op(E->stack[n]) || _op_name(E->stack[n]) != op)
        return -1;

    E->stack_len = n;
    return 0;
}

/* pop ops with precedence > prec */
static int mpoly_parse_pop_prec(mpoly_parse_t E, slong prec)
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
            E->R->mul(E->tmp, E->estore + E->R->elem_size*n3, E->estore + E->R->elem_size*n1, E->R->ctx);
            E->R->swap(E->estore + E->R->elem_size*n3, E->tmp, E->R->ctx);
            E->estore_len -= 1;
            E->stack_len -= 2;
        }
        else if (_op_name(n2) == OP_PLUS)
        {
            l1 = E->R->length(E->estore + E->R->elem_size*n1, E->R->ctx);
            l3 = E->R->length(E->estore + E->R->elem_size*n3, E->R->ctx);

        do_plus:

            if (l1 > l3)
            {
                FLINT_SWAP(slong, l3, l1);
                E->R->swap(E->estore + E->R->elem_size*n3, E->estore + E->R->elem_size*n1, E->R->ctx);
            }

            if (p > prec || 2*l1 >= l3)
            {
                E->R->add(E->estore + E->R->elem_size*n3, E->estore + E->R->elem_size*n3,
                                                 E->estore + E->R->elem_size*n1, E->R->ctx);
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
            l1 = E->R->length(E->estore + E->R->elem_size*n1, E->R->ctx);
            l3 = E->R->length(E->estore + E->R->elem_size*n3, E->R->ctx);

            if (4*l1 >= l3 || 4*l3 >= l1)
            {
                E->R->sub(E->estore + E->R->elem_size*n3, E->estore + E->R->elem_size*n3,
                                                 E->estore + E->R->elem_size*n1, E->R->ctx);
                E->estore_len -= 1;
                E->stack_len -= 2;
            }
            else
            {
                E->R->neg(E->estore + E->R->elem_size*n1, E->estore + E->R->elem_size*n1, E->R->ctx);
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
            if (!E->R->divides(E->tmp, E->estore + E->R->elem_size*n3,
                                    E->estore + E->R->elem_size*n1, E->R->ctx))
            {
                return -1;
            }

            E->R->swap(E->estore + E->R->elem_size*n3, E->tmp, E->R->ctx);
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
            E->R->neg(E->estore + E->R->elem_size*n1, E->estore + E->R->elem_size*n1, E->R->ctx);

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

int mpoly_parse_parse(mpoly_parse_t E, void * poly, const char * s, slong slen)
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

            E->R->set_fmpz(E->tmp, c, E->R->ctx);
            if (mpoly_parse_push_expr(E))
                goto failed;
        }
        else if (*s == '^')
        {
            if (++s >= send || !('0' <= *s && *s <= '9'))
                goto failed;

            s = _parse_int(c, s, send);

            if (mpoly_parse_pop_prec(E, PREC_POWER))
                goto failed;

            if (!mpoly_parse_top_is_expr(E))
                goto failed;

            if (!E->R->pow_fmpz(mpoly_parse_top_expr(E), mpoly_parse_top_expr(E), c, E->R->ctx))
                goto failed;
        }
        else if (*s == '*')
        {
            if (!mpoly_parse_top_is_expr(E))
                goto failed;

            if (mpoly_parse_pop_prec(E, PREC_TIMES))
                goto failed;

            mpoly_parse_push_op(E, _op_make(OP_TIMES, FIX_INFIX, PREC_TIMES));
            s++;
        }
        else if (*s == '+')
        {
            if (!mpoly_parse_top_is_expr(E))
            {
                mpoly_parse_push_op(E, _op_make(OP_PLUS, FIX_PREFIX, PREC_UPLUS));
            }
            else
            {
                if (mpoly_parse_pop_prec(E, PREC_PLUS))
                    goto failed;

                mpoly_parse_push_op(E, _op_make(OP_PLUS, FIX_INFIX, PREC_PLUS));
            }
            s++;
        }
        else if (*s == '-')
        {
            if (!mpoly_parse_top_is_expr(E))
            {
                mpoly_parse_push_op(E, _op_make(OP_MINUS, FIX_PREFIX, PREC_UMINUS));
            }
            else
            {
                if (mpoly_parse_pop_prec(E, PREC_MINUS))
                    goto failed;

                mpoly_parse_push_op(E, _op_make(OP_MINUS, FIX_INFIX, PREC_MINUS));
            }
            s++;
        }
        else if (*s == '/')
        {
            if (!mpoly_parse_top_is_expr(E))
                goto failed;

            if (mpoly_parse_pop_prec(E, PREC_DIVIDES))
                goto failed;

            mpoly_parse_push_op(E, _op_make(OP_DIVIDES, FIX_INFIX, PREC_DIVIDES));
            s++;
        }
        else if (*s == ' ')
        {
            s++;
        }
        else if (*s == '(')
        {
            if (mpoly_parse_top_is_expr(E))
                goto failed;

            mpoly_parse_push_op(E, _op_make(OP_LROUND, FIX_MATCHFIX, PREC_LOWEST));
            s++;
        }
        else if (*s == ')')
        {
            if (mpoly_parse_pop_prec(E, PREC_LOWEST))
                goto failed;

            if (mpoly_parse_pop_expr(E))
                goto failed;

            if (mpoly_parse_pop_op(E, OP_LROUND))
                goto failed;

            if (mpoly_parse_push_expr(E))
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
                    E->R->set(E->tmp, E->terminal_values + E->R->elem_size*k, E->R->ctx);
                    if (mpoly_parse_push_expr(E))
                        goto failed;

                    s += l;
                    goto continue_outer;
                }
            }

            goto failed;
        }
    continue_outer:;
    }

    if (mpoly_parse_pop_prec(E, PREC_LOWEST))
        goto failed;

    if (mpoly_parse_pop_expr(E))
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
