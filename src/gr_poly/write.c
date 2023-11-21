/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <ctype.h>
#include "gr_poly.h"

#ifdef __GNUC__
# define strcmp __builtin_strcmp
#else
# include <string.h>
#endif

static int
want_parens(const char * s)
{
    if (s[0] == '(' || s[0] == '[' || s[0] == '{')
        return 0;

    if (s[0] == '-')
        s++;

    while (s[0] != '\0')
    {
        if (!isalnum(s[0]) && s[0] != '.')
            return 1;

        s++;
    }

    return 0;
}

int
gr_poly_write(gr_stream_t out, const gr_poly_t poly, const char * x, gr_ctx_t ctx)
{
    int status;
    slong i, n;
    slong sz;
    char * s;
    int printed_previously = 0;

    sz = ctx->sizeof_elem;
    n = gr_poly_length(poly, ctx);
    status = GR_SUCCESS;

    if (n == 0)
    {
        gr_stream_write(out, "0");
        return status;
    }

    for (i = 0; i < n; i++)
    {
        if (gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx) == T_TRUE)
            continue;

        gr_get_str(&s, GR_ENTRY(poly->coeffs, i, sz), ctx);

        if (i >= 1 && !strcmp(s, "1"))
        {
            if (printed_previously)
                gr_stream_write(out, " + ");

            gr_stream_write(out, x);

            if (i >= 2)
            {
                gr_stream_write(out, "^");
                gr_stream_write_si(out, i);
            }
        }
        else if (i >= 1 && !strcmp(s, "-1"))
        {
            if (printed_previously)
                gr_stream_write(out, " - ");
            else
                gr_stream_write(out, "-");

            gr_stream_write(out, x);

            if (i >= 2)
            {
                gr_stream_write(out, "^");
                gr_stream_write_si(out, i);
            }
        }

        else
        {
            if (want_parens(s))
            {
                if (printed_previously)
                    gr_stream_write(out, " + ");

                gr_stream_write(out, "(");
                gr_stream_write_free(out, s);
                gr_stream_write(out, ")");
            }
            else
            {
                if (printed_previously && s[0] == '-')
                {
                    gr_stream_write(out, " - ");
                    gr_stream_write(out, s + 1);
                    flint_free(s);
                }
                else
                {
                    if (printed_previously)
                        gr_stream_write(out, " + ");

                    gr_stream_write_free(out, s);
                }
            }

            if (i == 1)
            {
                gr_stream_write(out, "*");
                gr_stream_write(out, x);
            }
            else if (i >= 2)
            {
                gr_stream_write(out, "*");
                gr_stream_write(out, x);
                gr_stream_write(out, "^");
                gr_stream_write_si(out, i);
            }
        }

        printed_previously = 1;

/*
        gr_stream_write(out, "(");
        status |= gr_write(out, GR_ENTRY(poly->coeffs, i, sz), ctx);
        gr_stream_write(out, ")");

        if (i == 1)
        {
            gr_stream_write(out, "*x");
        }
        else if (i >= 2)
        {
            gr_stream_write(out, "*x^");
            gr_stream_write_si(out, i);
        }

        if (i < n - 1)
            gr_stream_write(out, " + ");
*/
    }
    return status;
}

int
gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    return gr_poly_write(out, poly, "x", ctx);
}
