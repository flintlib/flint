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

Please refer to `code_conventions.txt` for some preliminary guidance on how to
write threaded functions in FLINT.

