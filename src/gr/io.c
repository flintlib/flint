/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz.h"
#include "gr.h"

#ifdef __GNUC__
# define memcpy __builtin_memcpy
# define strlen __builtin_strlen
#else
# include <math.h>
#endif

/* todo: error handling */

void gr_stream_init_file(gr_stream_t out, FILE * fp)
{
    out->fp = (FLINT_FILE *) fp;
    out->s = NULL;
}

void gr_stream_init_str(gr_stream_t out)
{
    out->fp = NULL;
    out->s = flint_malloc(16);
    out->s[0] = '\0';
    out->len = 0;
    out->alloc = 16;
}

void gr_stream_write(gr_stream_t out, const char * s)
{
    if (out->fp != NULL)
    {
        fprintf((FILE *) out->fp, "%s", s);
    }
    else
    {
        slong len, alloc;

        len = strlen(s);
        alloc = out->len + len + 1;

        if (alloc > out->alloc)
        {
            alloc = FLINT_MAX(alloc, out->alloc * 2);
            out->s = flint_realloc(out->s, alloc);
            out->alloc = alloc;
        }

        memcpy(out->s + out->len, s, len + 1);
        out->len += len;
    }
}

void
gr_stream_write_si(gr_stream_t out, slong x)
{
    if (out->fp != NULL)
    {
        flint_fprintf((FILE *) out->fp, "%wd", x);
    }
    else
    {
        char tmp[22];
        sprintf(tmp, WORD_FMT "d", x);
        gr_stream_write(out, tmp);
    }
}

void
gr_stream_write_ui(gr_stream_t out, ulong x)
{
    if (out->fp != NULL)
    {
        flint_fprintf((FILE *) out->fp, "%wu", x);
    }
    else
    {
        char tmp[22];
        sprintf(tmp, WORD_FMT "u", x);
        gr_stream_write(out, tmp);
    }
}

void
gr_stream_write_free(gr_stream_t out, char * s)
{
    gr_stream_write(out, s);
    flint_free(s);
}

void
gr_stream_write_fmpz(gr_stream_t out, const fmpz_t x)
{
    gr_stream_write_free(out, fmpz_get_str(NULL, 10, x));
}

int
gr_ctx_print(gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_ctx_write(out, ctx);
    return GR_SUCCESS;
}

int
gr_ctx_println(gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_ctx_write(out, ctx);
    gr_stream_write(out, "\n");
    return GR_SUCCESS;
}

int
gr_print(gr_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_write(out, x, ctx);
    return GR_SUCCESS;
}

int
gr_println(gr_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_write(out, x, ctx);
    gr_stream_write(out, "\n");
    return GR_SUCCESS;
}

int
gr_ctx_get_str(char ** s, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_str(out);
    gr_ctx_write(out, ctx);
    *s = out->s;
    return GR_SUCCESS;
}

int
gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_str(out);
    gr_write(out, x, ctx);
    *s = out->s;
    return GR_SUCCESS;
}

int
gr_get_str_n(char ** s, gr_srcptr x, slong n, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_str(out);
    gr_write_n(out, x, n, ctx);
    *s = out->s;
    return GR_SUCCESS;
}
