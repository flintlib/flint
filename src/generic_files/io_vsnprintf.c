/*
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <ctype.h> /* isdigit */
#include <stdint.h> /* intmax_t */
#include <stdio.h>
#include <string.h> /* memcpy, strncmp and strchr */
#include <stdarg.h>
#include <wchar.h> /* wchar_t and wint_t */
#include "nmod_types.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpq_types.h"
#include "fmpq.h"
#include "arf_types.h"
#include "arb.h"
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"
#include "gr_mat.h"
#include "gr_mat/impl.h"

typedef struct
{
    char * buf;
    size_t size;
    size_t used;
} flint_vsnprintf_out;

static void flint_vsnprintf_init(flint_vsnprintf_out * out, char * buf, size_t size)
{
    out->buf = buf;
    out->size = size;
    out->used = 0;

    if (out->buf != NULL && out->size > 0)
        out->buf[0] = '\0';
}

static size_t flint_vsnprintf_write(const void * buf, size_t len, flint_vsnprintf_out * out)
{
    if (out->buf != NULL && out->size > 0 && out->used < out->size - 1)
    {
        size_t avail = out->size - 1 - out->used;
        size_t copy = len < avail ? len : avail;

        memcpy(out->buf + out->used, buf, copy);
        out->used += copy;
        out->buf[out->used] = '\0';
    }

    return len;
}

static int flint_vsnprintf_putc(int ch, flint_vsnprintf_out * out)
{
    if (out->buf != NULL && out->size > 0 && out->used < out->size - 1)
    {
        out->buf[out->used] = (char) ch;
        out->used++;
        out->buf[out->used] = '\0';
    }
    return (unsigned char) ch;
}

static int flint_vsnprintf_vprintf(flint_vsnprintf_out * out, const char * fmt, va_list ap)
{
    va_list ap_copy;
    int res;
    char dummy[1];
    char * dst;
    size_t avail;

    if (out->buf != NULL && out->size > 0 && out->used < out->size)
    {
        dst = out->buf + out->used;
        avail = out->size - out->used;
    }
    else
    {
        dst = dummy;
        avail = 1;
    }

    va_copy(ap_copy, ap);
    res = vsnprintf(dst, avail, fmt, ap_copy);
    va_end(ap_copy);

    if (res < 0)
        return res;

    if (out->buf != NULL && out->size > 0 && out->used < out->size)
    {
        if ((size_t) res < avail)
            out->used += (size_t) res;
        else
            out->used = out->size - 1;
    }

    return res;
}

static int flint_vsnprintf_printf(flint_vsnprintf_out * out, const char * fmt, ...)
{
    va_list ap;
    int res;

    va_start(ap, fmt);
    res = flint_vsnprintf_vprintf(out, fmt, ap);
    va_end(ap);

    return res;
}

#define FLINT_VPRINTF_FUNCTION flint_vsnprintf
#define FLINT_VPRINTF_FUNCTION_ARGS char * s, size_t n
#define FLINT_VPRINTF_OUT_T flint_vsnprintf_out
#define FLINT_VPRINTF_INIT(out) flint_vsnprintf_init((out), s, n)
#define FLINT_VPRINTF_PRINTF(out, ...) flint_vsnprintf_printf((out), __VA_ARGS__)
#define FLINT_VPRINTF_VPRINTF(out, fmt, ap) flint_vsnprintf_vprintf((out), (fmt), (ap))
#define FLINT_VPRINTF_WRITE(buf, len, out) flint_vsnprintf_write((buf), (len), (out))
#define FLINT_VPRINTF_PUTC(ch, out) flint_vsnprintf_putc((ch), (out))
#define FLINT_VPRINTF_PUTC_ERRVAL (-1)
#define FLINT_VPRINTF_GR_STREAM_INIT(gr_out, out) gr_stream_init_str((gr_out))
#define FLINT_VPRINTF_GR_STREAM_FLUSH(gr_out, out) \
    do { \
        FLINT_VPRINTF_WRITE((gr_out)->s, (gr_out)->len, (out)); \
        flint_free((gr_out)->s); \
    } while (0)

#include "io_vprintf_impl.h"

#undef FLINT_VPRINTF_FUNCTION
#undef FLINT_VPRINTF_FUNCTION_ARGS
#undef FLINT_VPRINTF_OUT_T
#undef FLINT_VPRINTF_INIT
#undef FLINT_VPRINTF_PRINTF
#undef FLINT_VPRINTF_VPRINTF
#undef FLINT_VPRINTF_WRITE
#undef FLINT_VPRINTF_PUTC
#undef FLINT_VPRINTF_PUTC_ERRVAL
#undef FLINT_VPRINTF_GR_STREAM_INIT
#undef FLINT_VPRINTF_GR_STREAM_FLUSH

int flint_snprintf(char * s, size_t n, const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vsnprintf(s, n, str, vlist);
   va_end(vlist);

   return ret;
}

int flint_vsprintf(char * s, const char * str, va_list vlist)
{
    return flint_vsnprintf(s, INT_MAX, str, vlist);
}

int flint_sprintf(char * s, const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vsprintf(s, str, vlist);
   va_end(vlist);

   return ret;
}
