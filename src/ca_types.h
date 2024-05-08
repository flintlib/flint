/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_TYPES_H
#define CA_TYPES_H

#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* streams *******************************************************************/

/* TODO: Deprecate these and simply replace with gr_stream_struct */
#define calcium_stream_struct gr_stream_struct
#define calcium_stream_t gr_stream_t

#ifdef __cplusplus
}
#endif

#endif /* CA_TYPES_H */
