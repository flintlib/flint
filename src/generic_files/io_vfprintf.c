/*
    Copyright (C) 2023 Albin Ahlbäck
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


typedef struct {
    FILE * fs;
} flint_vfprintf_out;

static void flint_vfprintf_init(flint_vfprintf_out * out, FILE * fs)
{
    out->fs = fs;
}

static int flint_vfprintf_vprintf(flint_vfprintf_out * out, const char * fmt, va_list ap)
{
    return vfprintf(out->fs, fmt, ap);
}

static int flint_vfprintf_printf(flint_vfprintf_out * out, const char * fmt, ...)
{
    va_list ap;
    int res;

    va_start(ap, fmt);
    res = flint_vfprintf_vprintf(out, fmt, ap);
    va_end(ap);

    return res;
}

static size_t flint_vfprintf_write(const void * buf, size_t len, flint_vfprintf_out * out)
{
    return fwrite(buf, sizeof(char), len, out->fs);
}

static int flint_vfprintf_putc(int ch, flint_vfprintf_out * out)
{
    return fputc(ch, out->fs);
}

#define FLINT_VPRINTF_FUNCTION flint_vfprintf
#define FLINT_VPRINTF_FUNCTION_ARGS FILE * fs
#define FLINT_VPRINTF_OUT_T flint_vfprintf_out
#define FLINT_VPRINTF_INIT(out) flint_vfprintf_init((out), fs)
#define FLINT_VPRINTF_PRINTF(out, ...) flint_vfprintf_printf((out), __VA_ARGS__)
#define FLINT_VPRINTF_VPRINTF(out, fmt, ap) flint_vfprintf_vprintf((out), (fmt), (ap))
#define FLINT_VPRINTF_WRITE(buf, len, out) flint_vfprintf_write((buf), (len), (out))
#define FLINT_VPRINTF_PUTC(ch, out) flint_vfprintf_putc((ch), (out))
#define FLINT_VPRINTF_PUTC_ERRVAL (EOF)
#define FLINT_VPRINTF_GR_STREAM_INIT(gr_out, out) gr_stream_init_file((gr_out), (out)->fs)
#define FLINT_VPRINTF_GR_STREAM_FLUSH(gr_out, out) do { (void) (gr_out); (void) (out); } while (0)

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


int flint_fprintf(FILE * fs, const char * str, ...)
{
   va_list vlist;
   int ret;

   va_start(vlist, str);
   ret = flint_vfprintf(fs, str, vlist);
   va_end(vlist);

   return ret;
}
