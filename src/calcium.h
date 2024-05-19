/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CALCIUM_H
#define CALCIUM_H

#ifdef CALCIUM_INLINES_C
#define CALCIUM_INLINE
#else
#define CALCIUM_INLINE static inline
#endif

#include "acb_types.h"
#include "ca_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Input and output */

#ifdef FLINT_HAVE_FILE
void calcium_stream_init_file(calcium_stream_t out, FILE * fp);
#endif

CALCIUM_INLINE
void calcium_stream_init_str(calcium_stream_t out)
{
    out->fp = NULL;
    out->s = (char *) flint_malloc(16);
    out->s[0] = '\0';
    out->len = 0;
    out->alloc = 16;
}

void calcium_write(calcium_stream_t out, const char * s);
void calcium_write_si(calcium_stream_t out, slong x);
void calcium_write_fmpz(calcium_stream_t out, const fmpz_t c);
void calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags);

CALCIUM_INLINE
void calcium_write_free(calcium_stream_t out, char * s)
{
    calcium_write(out, s);
    flint_free(s);
}

/* Triple-valued logic */

/* TODO: Either remove this one or thruth_println in gr.h */
CALCIUM_INLINE void truth_print(truth_t t)
{
    if (t == T_TRUE) flint_printf("T_TRUE");
    if (t == T_FALSE) flint_printf("T_FALSE");
    if (t == T_UNKNOWN) flint_printf("T_UNKNOWN");
}

const char * calcium_func_name(calcium_func_code func);

/* Flint extras */

ulong calcium_fmpz_hash(const fmpz_t x);

#ifdef __cplusplus
}
#endif

#endif
