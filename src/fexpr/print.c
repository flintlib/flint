/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fexpr.h"
#include "fexpr_builtin.h"

void
fexpr_write(calcium_stream_t stream, const fexpr_t expr)
{
    ulong type = FEXPR_TYPE(expr->data[0]);

    switch (type)
    {
        case FEXPR_TYPE_SMALL_INT:
            {
                slong c = ((slong) (expr->data[0])) >> FEXPR_TYPE_BITS;
                calcium_write_si(stream, c);
            }
            break;
        case FEXPR_TYPE_SMALL_SYMBOL:
            {
                slong i;

                if (((expr->data[0] >> 8) & 0xff) == 0)
                {
                    calcium_write(stream, fexpr_builtin_table[expr->data[0] >> 16].string);
                }
                else
                {
                    char tmp[FEXPR_SMALL_SYMBOL_LEN + 1];

                    tmp[FEXPR_SMALL_SYMBOL_LEN] = '\0';
                    for (i = 0; i < FEXPR_SMALL_SYMBOL_LEN; i++)
                    {
                        char c = expr->data[0] >> ((i + 1) * 8);
                        tmp[i] = c;
                        if (c == '\0')
                            break;
                    }

                    calcium_write(stream, tmp);
                }
            }
            break;
        case FEXPR_TYPE_BIG_SYMBOL:
            {
                calcium_write(stream, (const char *) (expr->data + 1));
            }
            break;
        case FEXPR_TYPE_SMALL_STRING:
            {
                slong i;
                char tmp[FEXPR_SMALL_SYMBOL_LEN + 3];

                /* todo: escape string */
                tmp[FEXPR_SMALL_SYMBOL_LEN] = '\0';
                for (i = 0; i < FEXPR_SMALL_SYMBOL_LEN; i++)
                {
                    char c = expr->data[0] >> ((i + 1) * 8);
                    tmp[i] = c;
                    if (c == '\0')
                        break;
                }

                calcium_write(stream, "\"");
                calcium_write(stream, tmp);
                calcium_write(stream, "\"");
            }
            break;
        case FEXPR_TYPE_BIG_STRING:
            {
                calcium_write(stream, "\"");
                /* todo: escape string */
                calcium_write(stream, (const char *) (expr->data + 1));
                calcium_write(stream, "\"");
            }
            break;
        case FEXPR_TYPE_BIG_INT_NEG:
        case FEXPR_TYPE_BIG_INT_POS: /* todo: print without copying */
            {
                fmpz_t c;
                fmpz_init(c);
                fexpr_get_fmpz(c, expr);
                calcium_write_fmpz(stream, c);
                fmpz_clear(c);
            }
            break;
        case FEXPR_TYPE_CALL0:
        case FEXPR_TYPE_CALL1:
        case FEXPR_TYPE_CALL2:
        case FEXPR_TYPE_CALL3:
        case FEXPR_TYPE_CALL4:
        case FEXPR_TYPE_CALLN:
            {
                fexpr_t t;
                ulong * ptr;
                slong i, num;

                if (type == FEXPR_TYPE_CALLN)
                {
                    num = expr->data[1];
                    ptr = expr->data + expr->data[2];
                }
                else
                {
                    ptr = expr->data + FEXPR_HEADER_SIZE;
                    num = type - FEXPR_TYPE_CALL0;
                }

                t->data = ptr;
                t->alloc = 1;
                fexpr_write(stream, t);
                t->data += fexpr_size(t);

                calcium_write(stream, "(");

                for (i = 0; i < num; i++)
                {
                    fexpr_write(stream, t);
                    t->data += fexpr_size(t);

                    if (i < num - 1)
                        calcium_write(stream, ", ");
                }

                calcium_write(stream, ")");
            }
            break;
        default:
            calcium_write(stream, "?UNKNOWN?");
    }
}

void
fexpr_print(const fexpr_t expr)
{
    calcium_stream_t t;
    calcium_stream_init_file(t, stdout);
    fexpr_write(t, expr);
}

char *
fexpr_get_str(const fexpr_t expr)
{
    calcium_stream_t t;
    calcium_stream_init_str(t);
    fexpr_write(t, expr);
    return t->s;
}
