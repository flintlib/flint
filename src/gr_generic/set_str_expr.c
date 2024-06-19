/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz.h"
#include "mpoly.h"
#include "gr.h"
#include "gr_special.h"
#include "gr_vec.h"
#include "gr_generic.h"

#define PREC_LOWEST      0
#define PREC_PLUS        1
#define PREC_MINUS       1
#define PREC_PLUSMINUS   1
#define PREC_TIMES       2
#define PREC_DIVIDES     2
#define PREC_UPLUS       3
#define PREC_UMINUS      3
#define PREC_UPLUSMINUS  3
#define PREC_POWER       4
#define PREC_HIGHEST   255

#define OP_TIMES       0
#define OP_PLUS        1
#define OP_MINUS       2
#define OP_DIVIDES     3
#define OP_LROUND      4
#define OP_POWER       5
#define OP_PLUSMINUS   6

#define FIX_INFIX      0
#define FIX_PREFIX     1
#define FIX_POSTFIX    2
#define FIX_MATCHFIX   3

typedef struct {
    gr_ctx_struct * R;
    slong * stack;
    slong stack_len;
    slong stack_alloc;
    char * estore;
    slong estore_len;
    slong estore_alloc;
    void * tmp;
    string_with_length_struct * terminal_strings;
    char * terminal_values;
    slong terminals_alloc;
    slong terminals_len;
    int flags;
    _gr_method_get_si_op size_func;
} gr_parse_struct;

typedef gr_parse_struct gr_parse_t[1];

void _gr_parse_init(gr_parse_t E);
void _gr_parse_clear(gr_parse_t E);
void _gr_parse_add_terminal(gr_parse_t E, const char * s, const void * v);
int _gr_parse_parse(gr_parse_t E, void * res, const char * s, slong len);

FLINT_FORCE_INLINE int _is_op(slong a)
{
    return a >= 0;
}

FLINT_FORCE_INLINE slong _op_make(slong name, slong fix, slong prec)
{
    return (prec << 10) + (fix << 8) + (name << 0);
}

FLINT_FORCE_INLINE slong _op_prec(slong a)
{
    return (ulong)(a) >> 10;
}

FLINT_FORCE_INLINE slong _op_fix(slong a)
{
    return ((ulong)(a) >> 8) & 3;
}

FLINT_FORCE_INLINE slong _op_name(slong a)
{
    return a&255;
}

/* initialize the R member first */
void _gr_parse_init(gr_parse_t E)
{
    slong i;

    E->flags = 0;
    E->size_func = (_gr_method_get_si_op) _gr_length;

    E->stack_len = 0;
    E->stack_alloc = 20;
    E->stack = FLINT_ARRAY_ALLOC(E->stack_alloc, slong);

    E->estore_len = 0;
    E->estore_alloc = 10;
    E->estore = gr_heap_init_vec(E->estore_alloc, E->R);

    E->terminals_len = 0;
    E->terminals_alloc = 5;
    E->terminal_values = gr_heap_init_vec(E->terminals_alloc, E->R);
    E->terminal_strings = FLINT_ARRAY_ALLOC(E->terminals_alloc, string_with_length_struct);

    for (i = 0; i < E->terminals_alloc; i++)
    {
        E->terminal_strings[i].str = NULL;
        E->terminal_strings[i].str_len = 0;
    }

}

void _gr_parse_clear(gr_parse_t E)
{
    slong i;

    flint_free(E->stack);

    gr_heap_clear_vec(E->estore, E->estore_alloc, E->R);
    gr_heap_clear_vec(E->terminal_values, E->terminals_alloc, E->R);

    for (i = 0; i < E->terminals_alloc; i++)
        flint_free(E->terminal_strings[i].str);

    flint_free(E->terminal_strings);
}

void _gr_parse_add_terminal(gr_parse_t E, const char * s, const void * val)
{
    slong l, n = E->terminals_len;
    slong sz = E->R->sizeof_elem;

    if (n + 1 > E->terminals_alloc)
    {
        slong i = E->terminals_alloc;
        slong new_alloc = FLINT_MAX(n + 1, i + i/2);

        E->terminal_strings = (string_with_length_struct *) flint_realloc(
                                            E->terminal_strings, new_alloc*
                                            sizeof(string_with_length_struct));

        E->terminal_values = (char *) flint_realloc(E->terminal_values, sz * new_alloc);
        for ( ; i < new_alloc; i++)
        {
            E->terminal_strings[i].str = NULL;
            E->terminal_strings[i].str_len = 0;
            gr_init(GR_ENTRY(E->terminal_values, i, sz), E->R);
        }

        E->terminals_alloc = new_alloc;
    }

    l = strlen(s);
    E->terminal_strings[n].str_len = l;

    E->terminal_strings[n].str = (char *) flint_realloc(E->terminal_strings[n].str, l + 1);
    memcpy(E->terminal_strings[n].str, s, l + 1);

    GR_MUST_SUCCEED(gr_set(GR_ENTRY(E->terminal_values, n, sz), val, E->R));

    E->terminals_len = n + 1;

    while (n > 0 && E->terminal_strings[n-1].str_len < E->terminal_strings[n].str_len)
    {
        FLINT_SWAP(char *, E->terminal_strings[n-1].str, E->terminal_strings[n].str);
        FLINT_SWAP(slong, E->terminal_strings[n-1].str_len, E->terminal_strings[n].str_len);
        gr_swap(GR_ENTRY(E->terminal_values, n - 1, sz), GR_ENTRY(E->terminal_values, n, sz), E->R);
        n--;
    }
}

static int gr_parse_top_is_expr(const gr_parse_t E)
{
    return E->stack_len > 0 && !_is_op(E->stack[E->stack_len - 1]);
}

static void * gr_parse_top_expr(gr_parse_t E)
{
    slong sz = E->R->sizeof_elem;

    FLINT_ASSERT(E->stack_len > 0);
    FLINT_ASSERT(E->stack[E->stack_len - 1] < 0);

    return GR_ENTRY(E->estore, -1 - E->stack[E->stack_len - 1], sz);
}

static void _gr_parse_push_op(gr_parse_t E, slong op)
{
    FLINT_ASSERT(_is_op(op));
    _slong_array_fit_length(&E->stack, &E->stack_alloc, E->stack_len + 1);
    E->stack[E->stack_len] = op;
    E->stack_len++;
}

/* if the top is not an expr, push the tmp, otherwise fail */
static int _gr_parse_push_expr(gr_parse_t E)
{
    slong sz = E->R->sizeof_elem;

    if (gr_parse_top_is_expr(E))
        return -1;

    if (E->estore_len + 1 > E->estore_alloc)
    {
        slong i = E->estore_alloc;
        slong new_alloc = FLINT_MAX(E->estore_len + 1, i + i/2);
        E->estore = flint_realloc(E->estore, new_alloc*sz);
        for ( ; i < new_alloc; i++)
            gr_init(GR_ENTRY(E->estore, i, sz), E->R);
        E->estore_alloc = new_alloc;
    }

    _slong_array_fit_length(&E->stack, &E->stack_alloc, E->stack_len + 1);
    E->stack[E->stack_len] = -1 - E->estore_len;
    E->stack_len++;
    gr_swap(GR_ENTRY(E->estore, E->estore_len, sz), E->tmp, E->R);
    E->estore_len++;
    return 0;
}

/* if the top is an expr, pop it, otherwise fail */
static int _gr_parse_pop_expr(gr_parse_t E)
{
    slong sz = E->R->sizeof_elem;

    if (!gr_parse_top_is_expr(E))
        return -1;

    gr_swap(E->tmp, GR_ENTRY(E->estore, -1 - E->stack[E->stack_len - 1], sz), E->R);
    E->estore_len--;
    E->stack_len--;
    return 0;
}

/* if the top is an operation op, pop it, otherwise fail */
static int _gr_parse_pop_op(gr_parse_t E, slong op)
{
    slong n = E->stack_len - 1;

    if (n < 0 || !_is_op(E->stack[n]) || _op_name(E->stack[n]) != op)
        return -1;

    E->stack_len = n;
    return 0;
}

/* pop ops with precedence > prec */
static int _gr_parse_pop_prec(gr_parse_t E, slong prec)
{
    slong n, n1, n2, n3, p, l1, l3;
    slong sz = E->R->sizeof_elem;

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
            if (GR_SUCCESS != gr_mul(E->tmp, GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n1, sz), E->R))
            {
                return -1;
            }

            gr_swap(GR_ENTRY(E->estore, n3, sz), E->tmp, E->R);
            E->estore_len -= 1;
            E->stack_len -= 2;
        }
        else if (_op_name(n2) == OP_PLUS)
        {
            if (E->flags & GR_PARSE_BALANCE_ADDITIONS)
            {
                l1 = E->size_func(GR_ENTRY(E->estore, n1, sz), E->R);
                l3 = E->size_func(GR_ENTRY(E->estore, n3, sz), E->R);
            }
            else
            {
                l1 = l3 = 0;
            }

        do_plus:

            if (l1 > l3)
            {
                FLINT_SWAP(slong, l3, l1);
                gr_swap(GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n1, sz), E->R);
            }

            if (p > prec || 2*l1 >= l3)
            {
                if (GR_SUCCESS != gr_add(GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n1, sz), E->R))
                {
                    return -1;
                }

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
            if (E->flags & GR_PARSE_BALANCE_ADDITIONS)
            {
                l1 = E->size_func(GR_ENTRY(E->estore, n1, sz), E->R);
                l3 = E->size_func(GR_ENTRY(E->estore, n3, sz), E->R);
            }
            else
            {
                l1 = l3 = 0;
            }

            if (4*l1 >= l3 || 4*l3 >= l1)
            {
                if (GR_SUCCESS != gr_sub(GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n3, sz),
                                                 GR_ENTRY(E->estore, n1, sz), E->R))
                {
                    return -1;
                }

                E->estore_len -= 1;
                E->stack_len -= 2;
            }
            else
            {
                if (GR_SUCCESS != gr_neg(GR_ENTRY(E->estore, n1, sz), GR_ENTRY(E->estore, n1, sz), E->R))
                {
                    return -1;
                }

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
            if (GR_SUCCESS != gr_div(E->tmp, GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n1, sz), E->R))
            {
                return -1;
            }

            gr_swap(GR_ENTRY(E->estore, n3, sz), E->tmp, E->R);
            E->estore_len -= 1;
            E->stack_len -= 2;
        }
        else if (_op_name(n2) == OP_POWER)
        {
            if (GR_SUCCESS != gr_pow(E->tmp, GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n1, sz), E->R))
            {
                return -1;
            }

            gr_swap(GR_ENTRY(E->estore, n3, sz), E->tmp, E->R);
            E->estore_len -= 1;
            E->stack_len -= 2;
        }
        else if (_op_name(n2) == OP_PLUSMINUS)
        {
            if (GR_SUCCESS != gr_set_interval_mid_rad(E->tmp, GR_ENTRY(E->estore, n3, sz), GR_ENTRY(E->estore, n1, sz), E->R))
            {
                return -1;
            }

            gr_swap(GR_ENTRY(E->estore, n3, sz), E->tmp, E->R);
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
        {
            if (GR_SUCCESS != gr_neg(GR_ENTRY(E->estore, n1, sz), GR_ENTRY(E->estore, n1, sz), E->R))
                return -1;
        }
        else if (_op_name(n2) == OP_PLUSMINUS)
        {
            gr_ptr zero;
            GR_TMP_INIT(zero, E->R);
            if (GR_SUCCESS != gr_set_interval_mid_rad(GR_ENTRY(E->estore, n1, sz), zero, GR_ENTRY(E->estore, n1, sz), E->R))
            {
                GR_TMP_CLEAR(zero, E->R);
                return -1;
            }
            GR_TMP_CLEAR(zero, E->R);
        }

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

    switch (send - s)
    {
        case 1:
            fmpz_set_ui(c, s[0] - '0'); 
            s += 1;
            break;
        case 2:
            fmpz_set_ui(c, (s[0] - '0') * 10 + (s[1] - '0'));
            s += 2;
            break;
        case 3:
            fmpz_set_ui(c, (s[0] - '0') * 100 + (s[1] - '0') * 10 + (s[2] - '0'));
            s += 3;
            break;
        default:
            TMP_START;
            v = buffer = (char *) TMP_ALLOC((send - s + 1)*sizeof(char));
            while (s < send)
                *v++ = *s++;
            *v++ = '\0';
            fmpz_set_str(c, buffer, 10);
            TMP_END;
    }

    return s;
}

FLINT_FORCE_INLINE int is_digit(int c)
{
    return '0' <= c && c <= '9';
}

static const char * _parse_decimal(fmpz_t c, fmpz_t d, const char * s, const char * end)
{
    char * buffer;
    slong int_digits = 1;
    slong frac_digits = 0;
    slong exp_digits = 0;
    slong i;
    int exp_minus = 0;
    const char * s_frac = s;
    const char * s_exp = s;

    TMP_INIT;

    TMP_START;

    while (s + int_digits < end && is_digit(s[int_digits]))
        int_digits++;

    s_frac = s + int_digits;

    if (s_frac < end && s_frac[0] == '.')
    {
        /* skip the . */
        s_frac++;
        frac_digits = 0;

        while (s_frac + frac_digits < end && is_digit(s_frac[frac_digits]))
            frac_digits++;

        s_exp = s_frac + frac_digits;
    }
    else
    {
        s_exp = s + int_digits;
    }

    if (s_exp + 1 < end && (s_exp[0] == 'e' || s_exp[0] == 'E') &&
                        (is_digit(s_exp[1]) || (s_exp + 2 < end && (s_exp[1] == '+' || s_exp[1] == '-') && is_digit(s_exp[2]))))
    {
        /* skip the e or E */
        s_exp++;

        if (s_exp[0] == '-')
        {
            exp_minus = 1;
            s_exp++;
        }
        else if (s_exp[0] == '+')
        {
            s_exp++;
        }

        exp_digits = 1;

        while (s_exp + exp_digits < end && is_digit(s_exp[exp_digits]))
            exp_digits++;
    }

    buffer = TMP_ALLOC((FLINT_MAX(int_digits + frac_digits, exp_digits) + 1) * sizeof(char));

    if (exp_digits)
    {
        for (i = 0; i < exp_digits; i++)
            buffer[i] = s_exp[i];
        buffer[exp_digits] = '\0';

        fmpz_set_str(d, buffer, 10);
        if (exp_minus)
            fmpz_neg(d, d);
    }
    else
    {
        fmpz_zero(d);
    }

    for (i = 0; i < int_digits; i++)
        buffer[i] = s[i];

    if (frac_digits)
    {
        for (i = 0; i < frac_digits; i++)
            buffer[int_digits + i] = s_frac[i];

        fmpz_sub_ui(d, d, frac_digits);
    }

    buffer[int_digits + frac_digits] = '\0';

    fmpz_set_str(c, buffer, 10);

    TMP_END;

    return s_exp + exp_digits;
}

int _gr_parse_parse(gr_parse_t E, void * poly, const char * s, slong slen)
{
    const char * send = s + slen;
    fmpz_t c, d;
    int ret;

    fmpz_init(c);
    fmpz_init(d);
    E->tmp = poly;

    while (s < send)
    {
        if ('0' <= *s && *s <= '9')
        {
#if 1
            s = _parse_decimal(c, d, s, send);

            if (fmpz_is_zero(d))
            {
                if (GR_SUCCESS != gr_set_fmpz(E->tmp, c, E->R))
                    goto failed;

                if (_gr_parse_push_expr(E))
                    goto failed;
            }
            else
            {
                if (GR_SUCCESS != gr_set_fmpz_10exp_fmpz(E->tmp, c, d, E->R))
                    goto failed;

                if (_gr_parse_push_expr(E))
                    goto failed;
            }

#else
            s = _parse_int(c, s, send);

            if (GR_SUCCESS != gr_set_fmpz(E->tmp, c, E->R))
                goto failed;

            if (_gr_parse_push_expr(E))
                goto failed;
#endif
        }
        else if (*s == '^')
        {
            if (E->flags & GR_PARSE_RING_EXPONENTS)
            {
                if (!gr_parse_top_is_expr(E))
                    goto failed;

                if (_gr_parse_pop_prec(E, PREC_POWER))
                    goto failed;

                _gr_parse_push_op(E, _op_make(OP_POWER, FIX_INFIX, PREC_POWER));
                s++;
            }
            else
            {
                if (++s >= send || !('0' <= *s && *s <= '9'))
                    goto failed;

                s = _parse_int(c, s, send);

                if (_gr_parse_pop_prec(E, PREC_POWER))
                    goto failed;

                if (!gr_parse_top_is_expr(E))
                    goto failed;

                if (GR_SUCCESS != gr_pow_fmpz(gr_parse_top_expr(E), gr_parse_top_expr(E), c, E->R))
                    goto failed;
            }
        }
        else if (*s == '*')
        {
            if (!gr_parse_top_is_expr(E))
                goto failed;

            if (_gr_parse_pop_prec(E, PREC_TIMES))
                goto failed;

            _gr_parse_push_op(E, _op_make(OP_TIMES, FIX_INFIX, PREC_TIMES));
            s++;
        }
        else if (*s == '+')
        {
            if (s + 2 < send && s[1] == '/' && s[2] == '-')
            {
                if (!gr_parse_top_is_expr(E))
                {
                    _gr_parse_push_op(E, _op_make(OP_PLUSMINUS, FIX_PREFIX, PREC_UPLUSMINUS));
                }
                else
                {
                    if (_gr_parse_pop_prec(E, PREC_PLUSMINUS))
                        goto failed;

                    _gr_parse_push_op(E, _op_make(OP_PLUSMINUS, FIX_INFIX, PREC_PLUSMINUS));
                }
                s += 3;
            }
            else
            {
                if (!gr_parse_top_is_expr(E))
                {
                    _gr_parse_push_op(E, _op_make(OP_PLUS, FIX_PREFIX, PREC_UPLUS));
                }
                else
                {
                    if (_gr_parse_pop_prec(E, PREC_PLUS))
                        goto failed;

                    _gr_parse_push_op(E, _op_make(OP_PLUS, FIX_INFIX, PREC_PLUS));
                }
                s++;
            }
        }
        else if (*s == '-')
        {
            if (!gr_parse_top_is_expr(E))
            {
                _gr_parse_push_op(E, _op_make(OP_MINUS, FIX_PREFIX, PREC_UMINUS));
            }
            else
            {
                if (_gr_parse_pop_prec(E, PREC_MINUS))
                    goto failed;

                _gr_parse_push_op(E, _op_make(OP_MINUS, FIX_INFIX, PREC_MINUS));
            }
            s++;
        }
        else if (*s == '/')
        {
            if (!gr_parse_top_is_expr(E))
                goto failed;

            if (_gr_parse_pop_prec(E, PREC_DIVIDES))
                goto failed;

            _gr_parse_push_op(E, _op_make(OP_DIVIDES, FIX_INFIX, PREC_DIVIDES));
            s++;
        }
        else if (*s == ' ')
        {
            s++;
        }
        else if (*s == '(' || *s == '[')
        {
            if (gr_parse_top_is_expr(E))
                goto failed;

            _gr_parse_push_op(E, _op_make(OP_LROUND, FIX_MATCHFIX, PREC_LOWEST));
            s++;
        }
        else if (*s == ')' || *s == ']')
        {
            if (_gr_parse_pop_prec(E, PREC_LOWEST))
                goto failed;

            if (_gr_parse_pop_expr(E))
                goto failed;

            if (_gr_parse_pop_op(E, OP_LROUND))
                goto failed;

            if (_gr_parse_push_expr(E))
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
                    GR_MUST_SUCCEED(gr_set(E->tmp, GR_ENTRY(E->terminal_values, k, E->R->sizeof_elem), E->R));
                    if (_gr_parse_push_expr(E))
                        goto failed;

                    s += l;
                    goto continue_outer;
                }
            }

            /* some builtin constants (for R and C) */
            /* todo: builtin functions */

            if (0 == strncmp(s, "pi", 2) || 0 == strncmp(s, "Pi", 2))
            {
                if (GR_SUCCESS != gr_pi(E->tmp, E->R))
                    goto failed;
                if (_gr_parse_push_expr(E))
                    goto failed;
                s += 2;
                goto continue_outer;
            }

            if (0 == strncmp(s, "i", 1) || 0 == strncmp(s, "I", 1))
            {
                if (GR_SUCCESS != gr_i(E->tmp, E->R))
                    goto failed;
                if (_gr_parse_push_expr(E))
                    goto failed;
                s += 1;
                goto continue_outer;
            }

            goto failed;
        }
    continue_outer:;
    }

    if (_gr_parse_pop_prec(E, PREC_LOWEST))
        goto failed;

    if (_gr_parse_pop_expr(E))
        goto failed;

    if (E->stack_len != 0)
        goto failed;

    ret = 0;

done:

    fmpz_clear(c);
    fmpz_clear(d);
    return ret;

failed:

    ret = -1;
    goto done;
}

int
gr_generic_set_str_expr(gr_ptr res, const char * s, int flags, gr_ctx_t ctx)
{
    gr_parse_t parse;
    gr_vec_t gens;
    slong i;
    char * g;
    int status;

    /* Quickly see if we simply have an integer literal, e.g. 0 */
    fmpz_t c;
    fmpz_init(c);
    if (!fmpz_set_str(c, s, 10))
    {
        status = gr_set_fmpz(res, c, ctx);
    }
    else
    {
        parse->R = ctx;
        _gr_parse_init(parse);
        parse->flags = flags;

        gr_vec_init(gens, 0, ctx);
        if (gr_gens_recursive(gens, ctx) == GR_SUCCESS)
        {
            for (i = 0; i < gens->length; i++)
            {
                GR_MUST_SUCCEED(gr_get_str(&g, gr_vec_entry_srcptr(gens, i, ctx), ctx));
                /* todo: version that consumes s and x */
                _gr_parse_add_terminal(parse, g, gr_vec_entry_srcptr(gens, i, ctx));
                flint_free(g);
            }
        }

        gr_vec_clear(gens, ctx);

        status = _gr_parse_parse(parse, res, s, strlen(s)) ? GR_UNABLE : GR_SUCCESS;

        _gr_parse_clear(parse);
    }

    fmpz_clear(c);

    return status;
}

int
gr_generic_set_str(gr_ptr res, const char * s, gr_ctx_t ctx)
{
    return gr_generic_set_str_expr(res, s, 0, ctx);
}

int
gr_generic_set_str_balance_additions(gr_ptr res, const char * s, gr_ctx_t ctx)
{
    return gr_generic_set_str_expr(res, s, GR_PARSE_BALANCE_ADDITIONS, ctx);
}

int
gr_generic_set_str_ring_exponents(gr_ptr res, const char * s, gr_ctx_t ctx)
{
    return gr_generic_set_str_expr(res, s, GR_PARSE_RING_EXPONENTS, ctx);
}
