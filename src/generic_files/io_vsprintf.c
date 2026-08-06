/*
    Copyright (C) 2026 Lars Göttgens
    Copyright (C) 2026 Edgar Costa

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

/*
    Sink for flint_sprintf-style functions: writes into a caller-supplied
    buffer with no length bound. Uses system vsprintf rather than vsnprintf
    because some 32-bit glibc versions corrupt output when vsnprintf is
    called with a very large size_t (issue #2646).
*/
typedef struct
{
    char * buf;
    size_t used;
} flint_vsprintf_out;

static void flint_vsprintf_init(flint_vsprintf_out * out, char * buf)
{
    out->buf = buf;
    out->used = 0;
    out->buf[0] = '\0';
}

static int flint_vsprintf_vprintf(flint_vsprintf_out * out, const char * fmt, va_list ap)
{
    va_list ap_copy;
    int res;

    va_copy(ap_copy, ap);
    res = vsprintf(out->buf + out->used, fmt, ap_copy);
    va_end(ap_copy);

    if (res > 0)
        out->used += (size_t) res;

    return res;
}

static int flint_vsprintf_printf(flint_vsprintf_out * out, const char * fmt, ...)
{
    va_list ap;
    int res;

    va_start(ap, fmt);
    res = flint_vsprintf_vprintf(out, fmt, ap);
    va_end(ap);

    return res;
}

static size_t flint_vsprintf_write(const void * buf, size_t len, flint_vsprintf_out * out)
{
    memcpy(out->buf + out->used, buf, len);
    out->used += len;
    out->buf[out->used] = '\0';
    return len;
}

static int flint_vsprintf_putc(int ch, flint_vsprintf_out * out)
{
    out->buf[out->used] = (char) ch;
    out->used++;
    out->buf[out->used] = '\0';
    return (unsigned char) ch;
}

#define FLINT_VPRINTF_FUNCTION flint_vsprintf
#define FLINT_VPRINTF_FUNCTION_ARGS char * s
#define FLINT_VPRINTF_OUT_T flint_vsprintf_out
#define FLINT_VPRINTF_INIT(out) flint_vsprintf_init((out), s)
#define FLINT_VPRINTF_PRINTF(out, ...) flint_vsprintf_printf((out), __VA_ARGS__)
#define FLINT_VPRINTF_VPRINTF(out, fmt, ap) flint_vsprintf_vprintf((out), (fmt), (ap))
#define FLINT_VPRINTF_WRITE(buf, len, out) flint_vsprintf_write((buf), (len), (out))
#define FLINT_VPRINTF_PUTC(ch, out) flint_vsprintf_putc((ch), (out))
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

int flint_sprintf(char * s, const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vsprintf(s, str, vlist);
   va_end(vlist);

   return ret;
}
