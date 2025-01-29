.. _memory:

**Memory management**
===============================================================================

Memory allocation functions
-------------------------------------------------------------------------------

The file ``flint.h`` defines functions ``flint_malloc``,
``flint_realloc``, ``flint_calloc`` and ``flint_free``.
They have the same interface as the standard library functions, but
may perform additional error checking.

By default the memory allocation functions wrap the system's
``malloc``, ``realloc``, ``calloc`` and ``free``.
The user can override this behaviour by calling ``__flint_set_memory_functions``
passing the ``malloc``, ``realloc``, ``calloc`` and ``free`` function
pointers as parameters (see ``flint.h`` for the exact prototype).
The current memory functions can be returned in a similar manner by calling
``__flint_get_memory_functions`` passing the address of pointers in which
the function pointers can be stored.

Memory allocated with ``flint_malloc`` must be freed with
``flint_free`` and not with ``free``.

Global caches and cleanup
-------------------------------------------------------------------------------

FLINT may cache some data (such as allocated integers
and tables of prime numbers) to speed up various computations.
If FLINT is built in threadsafe mode, most caches are thread-local
(some are always global and shared among the threads).

Data cached by the current thread can be freed by calling the
``flint_cleanup()`` function. The user can register additional cleanup
functions to be invoked
by ``flint_cleanup()`` by passing a pointer
to a function with signature ``void cleanup_function(void)``
to ``flint_register_cleanup_function()``.

The user should call ``flint_cleanup_master()`` exactly once
right before exiting a program. This cleans up all caches in all threads
and should result in a clean output with tools like ``valgrind``
if there are no memory leaks.

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
       nn_ptr a, b;
       TMP_INIT;
    
       /* arbitrary code */
    
       TMP_START; /* we are about to do some allocation */
    
       /* arbitrary code */
    
       a = TMP_ALLOC(32*sizeof(ulong));
       b = TMP_ALLOC(64*sizeof(ulong));
    
       /* arbitrary code */
    
       TMP_END; /* cleans up a and b */
    
       /* arbitrary code */
    }

It is very important to note that temporary allocations should not be
made in recursive functions or in loop bodies, as many small allocations
on the stack can exhaust the stack causing a stack overflow.

