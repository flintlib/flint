/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>
#include "test_helpers.h"
#include "flint.h"

TEST_FUNCTION_START(memory_manager, state)
{
    int result;
    size_t alignment;
    void * (* fmalloc)(size_t);
    void * (* fcalloc)(size_t, size_t);
    void * (* frealloc)(void *, size_t);
    void (* ffree)(void *);
    void * (* faligned_alloc)(size_t, size_t);
    void (* faligned_free)(void *);

    /* Standard FLINT align functions */
    for (alignment = sizeof(ulong); alignment < (WORD(1) << 14); alignment <<= 1)
    {
        slong ix;

        for (ix = 0; ix < 10 * flint_test_multiplier(); ix++)
        {
            size_t size;
            void * ptr;

            /* Allocate at most around 1 MB */
            size = alignment * (1 + n_randint(state, 100));

            ptr = flint_aligned_alloc(alignment, size);

            /* Set memory */
            memset(ptr, 0, size);

            flint_aligned_free(ptr);
        }
    }

    /* Use non-aligned backend */
    __flint_set_memory_functions(malloc, calloc, realloc, free);

    __flint_get_memory_functions(&fmalloc, &fcalloc, &frealloc, &ffree);
    result = (fmalloc == malloc
            && fcalloc == calloc
            && frealloc == realloc
            && ffree == free);
    if (!result)
        TEST_FUNCTION_FAIL("Non-standard non-aligned memory functions are different.");

    for (alignment = sizeof(ulong); alignment < (WORD(1) << 14); alignment <<= 1)
    {
        slong ix;

        for (ix = 0; ix < 10 * flint_test_multiplier(); ix++)
        {
            size_t size;
            void * ptr;

            /* Allocate at most around 1 MB */
            size = alignment * (1 + n_randint(state, 100));

            ptr = flint_aligned_alloc(alignment, size);

            /* Set memory */
            memset(ptr, 0, size);

            flint_aligned_free(ptr);
        }
    }

#if defined(_MSC_VER) || defined(__MINGW32__) || defined(__MINGW64__)
# define aligned_alloc _aligned_malloc
# define aligned_free _aligned_free
#else
# define aligned_free free
#endif
    __flint_set_all_memory_functions(malloc, calloc, realloc, free, aligned_alloc, aligned_free);
    __flint_get_all_memory_functions(&fmalloc, &fcalloc, &frealloc, &ffree, &faligned_alloc, &faligned_free);

    result = (fmalloc == malloc
            && fcalloc == calloc
            && frealloc == realloc
            && ffree == free
            && faligned_alloc == aligned_alloc
            && faligned_free == aligned_free);

    if (!result)
        TEST_FUNCTION_FAIL("Non-standard aligned memory functions are different.");

    TEST_FUNCTION_END(state);
}
