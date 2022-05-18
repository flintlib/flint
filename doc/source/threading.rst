.. _threading:

**Threading**
===============================================================================

Multithreaded FLINT
-------------------------------------------------------------------------------

FLINT provides a number of multithreaded functions, which use multiple threads
by default if FLINT was built with at least pthreads. (This functionality works
best when thread local storage is also available on the system.)

By default, FLINT will just use one thread. To control the maximum number of
threads FLINT uses, one can call the function ``flint_set_num_threads(n)``,
where `n` is the maximum number of threads to use.

One can also query the current thread limit by calling
``flint_get_num_threads()``.

Each version of FLINT brings new functions that are threaded by default.

Many core algorithms such as the FFT (for large integer and polynomial
operations, including some factoring algorithms), integer factoring and
multivariate polynomial algorithms are threaded in FLINT.

Writing threaded functions in FLINT
-----------------------------------

Flint uses a custom thread pool for threading. This involves creating a worker
function, requesting threads from the thread pool, starting the threads,
waiting for them to finish, then giving the threads back to the pool. Simple
examples of this include ``fmpz_mod_mat_mul_classical_threaded`` and
``fmpz_poly_taylor_shift_multi_mod``.

The user should not have to run specialised versions of functions to get
threading. This means that user facing functions should generally not have
``_threaded`` appended to their name. Either there is a single function that does
the job, and it happens to be threaded, or there is a best-of-breed function
that calls the appropriate threaded function when this is the best strategy.

There are some instances where it may be desirable (e.g. for testing purposes,
or because naming proves difficult) where one wants a _threaded in the name.
But these cases should be rare.

In some cases, one does not want functions to request threads from the pool
themselves, but to accept threads from another function which has already
obtained them. Such functions will accept an array of thread pool handles
and a number of threads. The naming convention for such functions is to append
``_threaded_pool`` to the end of their name. However, the usual distinctions
between underscore and non-underscore functions should still apply.

Functions should request ``flint_get_num_threads()`` threads from the thread pool.
The function should not exceed this number of threads in total. In general a
thread that is woken should start zero additional workers. However, if this is
not the desired behaviour, an option exists to the function for waking worker
threads to alter how many threads it can start. In some cases it is also
necessary to temporarily restrict the number of worker threads a given function
can start. This is accomplished by calling flint_set_num_workers() and then
once the function is called, calling flint_reset_num_workers(). Any function
threaded function which calls flint_get_num_threads() to determine how many
threads to request from the thread pool will be appropriately restricted by
such calls.

Note that if ``flint_get_num_threads()`` returns ``n`` then the number of workers that
can be started is ``n - 1`` (in addition to the thread the function is already
running in). For this reason our documentation often distinguishes number of
workers and number of threads. Please refer to the thread pool interface and
Flint threading interface documentation to see the exact specification.

Functional parallel programming helpers
---------------------------------------

The following convenience function are defined in ``thread_support.h``.
They are currently experimental, and
the interfaces might change in a future version.

.. type:: void (* do_func_t)(slong i, void * args)

.. function:: void flint_parallel_do(do_func_t f, void * args, slong n, int thread_limit, int flags)

    Evaluate ``f(i, args)`` for `0 \le i < n - 1` in parallel using
    up to ``thread_limits`` threads (including the master thread).
    If *thread_limit* is nonpositive, the number of threads defaults to
    ``flint_get_num_threads()``.

    The following ``flags`` are supported:

    ``FLINT_PARALLEL_UNIFORM`` - assumes that the cost of function
    calls is roughly constant, so that scheduling uniformly into
    blocks is efficient.

    ``FLINT_PARALLEL_STRIDED`` - assumes that the cost increases
    or decreases monotonically with ``i``, so that strided
    scheduling is efficient.

    ``FLINT_PARALLEL_DYNAMIC`` (not implemented) - use dynamic
    scheduling.

    ``FLINT_PARALLEL_VERBOSE`` - print information.

.. type:: void (* bsplit_merge_func_t)(void *, void *, void *, void *)

.. type:: void (* bsplit_basecase_func_t)(void *, slong, slong, void *)

.. type:: void (* bsplit_init_func_t)(void *, void *)

.. type:: void (* bsplit_clear_func_t)(void *, void *)

.. function:: void flint_parallel_binary_splitting(void * res, bsplit_basecase_func_t basecase, bsplit_merge_func_t merge, size_t sizeof_res, bsplit_init_func_t init, bsplit_clear_func_t clear, void * args, slong a, slong b, slong basecase_cutoff, int thread_limit, int flags)

    Sets ``res`` to `f(a) \circ f(a+1) \circ \cdots \circ f(b - 1)`
    computed using parallel binary splitting, using
    up to ``thread_limits`` threads (including the master thread).
    If *thread_limit* is nonpositive, the number of threads defaults to
    ``flint_get_num_threads()``.

    The function ``basecase(res, a, b, args)`` gets called
    when `b - a` does not exceed ``basecase_cutoff``, which
    must be at least 1.

    The function ``merge(res, x, y, args)`` implements the
    associative operation (`x \circ y`), writing the result to ``res``.
    If called with ``FLINT_PARALLEL_BSPLIT_LEFT_INPLACE`` in ``flags``,
    the same space will be used for ``res`` and ``x``.

    A result is assumed to fit in a structure of size ``sizeof_res``.
    The functions ``init(res, args)`` and ``clear(res, args)``
    initialize and clear intermediate result objects.

