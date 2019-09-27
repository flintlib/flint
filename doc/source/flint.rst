.. _flint:

**flint.h** -- global definitions
===============================================================================

.. function:: flint_rand_s * flint_rand_alloc()

    Allocates a ``flint_rand_t`` object to be used like a heap-allocated
    ``flint_rand_t`` in external libraries.
    The random state is not initialised.

.. function:: void flint_rand_free(flint_rand_s * state)
   
    Frees a random state object as allocated using ``flint_rand_alloc``.

