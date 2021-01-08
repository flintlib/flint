/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fexpr.h"

void
fexpr_print(const fexpr_t expr)
{
    ulong type = FEXPR_TYPE(expr->data[0]);

    switch (type)
    {
        case FEXPR_TYPE_SMALL_INT:
            {
                slong c = ((slong) (expr->data[0])) >> FEXPR_TYPE_BITS;
                flint_printf("%ld", c);
            }
            break;
        case FEXPR_TYPE_SMALL_SYMBOL:
            {
                slong i;
                for (i = 0; i < FEXPR_SMALL_SYMBOL_LEN; i++)
                {
                    char c = expr->data[0] >> ((i + 1) * 8);

                    if (c == '\0')
                        break;
                    else
                        flint_printf("%c", c);
                }
            }
            break;
        case FEXPR_TYPE_BIG_SYMBOL:
            {
                flint_printf("%s", (const char *) (expr->data + 1));
            }
            break;
        case FEXPR_TYPE_BIG_INT_NEG:
        case FEXPR_TYPE_BIG_INT_POS: /* todo: print without copying */
            {
                fmpz_t c;
                fmpz_init(c);
                fexpr_get_fmpz(c, expr);
                fmpz_print(c);
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
                fexpr_print(t);
                t->data += fexpr_size(t);

                flint_printf("(");

                for (i = 0; i < num; i++)
                {
                    fexpr_print(t);
                    t->data += fexpr_size(t);

                    if (i < num - 1)
                        printf(", ");
                }

                flint_printf(")");
            }
            break;
        default:
            flint_printf("?UNKNOWN?");
    }
}
