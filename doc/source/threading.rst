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
where `n` is the maximum number of threads to use. Note that the limit only
applies to threads started directly by FLINT; external libraries such as
OpenBLAS may start more than `n` threads even when called from FLINT.

One can also query the current thread limit by calling
``flint_get_num_threads()``.

Each version of FLINT brings new functions that are threaded by default.

Many core algorithms such as the FFT (for large integer and polynomial
operations, including some factoring algorithms), integer factoring and
multivariate polynomial algorithms are threaded in FLINT.

Writing threaded functions in FLINT
-----------------------------------

FLINT uses a custom thread pool for threading. This involves creating a worker
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
once the function is called, calling flint_reset_num_workers(). Any
threaded function which calls flint_get_num_threads() to determine how
many threads to request from the thread pool will be appropriately
restricted by such calls.

Note that if ``flint_get_num_threads()`` returns ``n`` then the number of workers that
can be started is ``n - 1`` (in addition to the thread the function is already
running in). For this reason our documentation often distinguishes number of
workers and number of threads. Please refer to the thread pool interface and
FLINT threading interface documentation to see the exact specification.

Functional parallel programming helpers
---------------------------------------

The following convenience function are defined in ``thread_support.h``.
They are currently experimental, and
the interfaces might change in a future version.

.. function:: slong flint_get_num_available_threads()

    Returns the number of threads that are not currently in use.

.. type:: void (* do_func_t)(slong i, void * args)

.. function:: void flint_parallel_do(do_func_t f, void * args, slong n, int thread_limit, int flags)

    Evaluate ``f(i, args)`` for `0 \le i < n` in parallel using
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

Shared caches
-------------------------------------------------------------------------------

A shared cache, declared in ``thread_support.h``, holds one process-wide,
read-only *version* of some lazily computed data -- a table of
constants, say -- that any thread reads without locking: an acquire load
of the current version, a plain load on x86-64.  The cache grows by
replacement: a thread that finds the current version insufficient takes
the cache's lock, checks again (another thread may have grown the cache
meanwhile), builds a larger version, usually from the current one, and
publishes it with a release store, so that a reader seeing the new
version also sees its contents.  Superseded versions are retired, not
freed, since other threads may still be reading them; everything is
freed by :func:`flint_cleanup_master` (not by the per-thread
:func:`flint_cleanup`, which a worker thread may call while others are
still computing).  Versions should therefore grow geometrically, so
that the retired ones hold a bounded multiple of the current one's
memory.  The ``mp_real`` module keeps its reduction tables in such
caches.

.. type:: flint_shared_cache_struct

    A cache: the lock, the current version, the retired versions, the
    destructor of a version and the link into the cleanup registry.
    Caches are static objects initialized with
    ``FLINT_SHARED_CACHE_INIT``; the version objects are whatever the
    module chooses, handled as ``void *``.

.. macro:: FLINT_SHARED_CACHE_INIT(free_version)

    The static initializer of a cache whose versions are freed by
    ``void free_version(void * version)``, for example
    ``static flint_shared_cache_struct my_cache = FLINT_SHARED_CACHE_INIT(my_free);``.

.. type:: flint_shared_cache_free_t

    ``void (* f)(void * version)``: frees one version.

.. function:: void * flint_shared_cache_get(flint_shared_cache_struct * cache)

    Returns the current version, or ``NULL`` if none has been published.
    The version stays valid until :func:`flint_cleanup_master` (or
    :func:`flint_shared_cache_clear`), whatever other threads publish
    meanwhile.

.. type:: flint_shared_cache_sufficient_t

    ``int (* f)(const void * version, void * args)``: whether a version
    satisfies the request described by *args*.

.. type:: flint_shared_cache_build_t

    ``void * (* f)(const void * version, void * args)``: a new version
    satisfying the request, built from the current one (``NULL`` if
    none) -- which must not be modified -- and typically covering it, so
    that the new version satisfies every request the old one did.

.. function:: void * flint_shared_cache_ensure(flint_shared_cache_struct * cache, flint_shared_cache_sufficient_t sufficient, flint_shared_cache_build_t build, void * args)

    Returns a version satisfying the request *args*: the current one if
    ``sufficient(version, args)`` holds, else the version
    ``build(current, args)``, built and published under the cache's lock
    after a second check of the current version.  This is the whole
    protocol; the following primitives let a module run it by hand.

.. function:: void flint_shared_cache_lock(flint_shared_cache_struct * cache)
              void flint_shared_cache_unlock(flint_shared_cache_struct * cache)

    Takes and releases the cache's lock, which serializes the growth of
    the cache.  Readers do not take it.

.. function:: void flint_shared_cache_publish(flint_shared_cache_struct * cache, void * version)

    Under the lock: makes *version* the current version, retiring the
    previous one, and registers the cache for :func:`flint_cleanup_master`.

.. function:: void flint_shared_cache_clear(flint_shared_cache_struct * cache)

    Frees every version of the cache, leaving it empty.  No other
    thread may be reading the cache; this is for tuning and testing
    programs measuring cold builds.  (:func:`flint_cleanup_master` does
    this for every cache that has published a version.)
