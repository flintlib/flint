.. _memory:

**Memory management**
===============================================================================

Memory management in FLINT
-------------------------------------------------------------------------------

The file ``flint.h`` defines functions ``flint_malloc``,
``flint_realloc``, ``flint_calloc`` and ``flint_free``.
They have the same interface as the standard library functions, but
may perform additional error checking.

FLINT may cache some data (such as allocated integers
and tables of prime numbers) to speed up various computations.
If FLINT is built in threadsafe mode, cached data is kept in thread-local
storage by default (unless configured otherwise). Cached data can be freed
by calling the ``flint_cleanup()`` function. It is recommended to call
``flint_cleanup()`` right before exiting a thread, and at the end of the
main program.

The user can register additional cleanup functions to be invoked
by ``flint_cleanup()`` by passing a pointer
to a function with signature ``void cleanup_function(void)``
to ``flint_register_cleanup_function()``.

Flint also makes use of ``malloc``, ``realloc``, ``calloc`` and
``free`` by default to allocate memory internally. The user is able to
override this behaviour by calling ``__flint_set_memory_functions``
passing the ``malloc``, ``realloc``, ``calloc`` and ``free`` function
pointers as parameters (see ``flint.h`` for the exact prototype).

The current memory functions can be returned in a similar manner by calling
``__flint_get_memory_functions`` passing the address of pointers in which
the function pointers can be stored.

Temporary allocation
-------------------------------------------------------------------------------

FLINT allows for temporary allocation of memory using ``alloca``
to allocate on the stack if the allocation is small enough.

The following program demonstrates how to use this facility to
allocate two different arrays.

.. code-block:: C

    #include <gmp.h>
    #include "flint.h"
    
    void myfun(void)
    {
       /* other variable declarations */
       mp_ptr a, b;
       TMP_INIT;
    
       /* arbitrary code*/
    
       TMP_START; /* we are about to do some allocation */
    
       /* arbitrary code */
    
       a = TMP_ALLOC(32*sizeof(mp_limb_t));
       b = TMP_ALLOC(64*sizeof(mp_limb_t));
    
       /* arbitrary code */
    
       TMP_END; /* cleans up a and b */
    
       /* arbitrary code */
    }

It is very important to note that temporary allocations should not be
made in recursive functions, as many small allocations on the stack
can exhaust the stack causing a stack overflow.

