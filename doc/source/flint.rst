.. _flint:

**flint.h** -- global definitions
===============================================================================

Macros
-----------------------------------------------

The file ``flint.h`` contains various useful macros.

.. macro:: __FLINT_VERSION
           __FLINT_VERSION_MINOR
           __FLINT_VERSION_PATCHLEVEL

    The major, minor and patch for current version of FLINT.

.. macro:: __FLINT_RELEASE

    Equivalent to ``10000 * __FLINT_VERSION + 100 * __FLINT_VERSION_MINOR +
    __FLINT_VERSION_PATCHLEVEL``.

.. macro:: FLINT_VERSION

    A static text string giving the version number, e.g. ``3.1.0`` or ``3.2.0-dev``.

.. macro:: FLINT_BITS

    The constant defining how many bits per limb on the machine. We require this
    to be either 32 or 64. This constant is set during the configuration.

.. macro:: FLINT_D_BITS

    A constant set at compile time to be the number of bits per double on the
    machine or one less than the number of bits per limb, whichever is smaller.
    This will have the value *31* on 32-bit systems and *53* on 64-bit systems.
    Numerous internal functions using precomputed inverses only support operands
    up to ``FLINT_D_BITS`` bits, hence the macro.

.. macro:: FLINT_ABS(x)

    Returns the absolute value of *x* for primitive signed numerical types.  It
    might fail for least negative values such as *INT_MIN* and *LONG_MIN*.

.. macro:: FLINT_UABS(x)

    Returns the absolute value of *x* for primitive signed numerical types,
    casting the result to an *ulong*. The result is well-defined
    for least negative values.

.. macro:: FLINT_MIN(x, y)
           FLINT_MAX(x, y)

    Returns the minimum or maximum of *x* and *y* for primitive types. This
    macro is only safe to use when *x* and *y* are of the same type, to avoid
    problems with integer promotion.

.. macro:: FLINT_SWAP(T, x, y)

    Swaps *x* and *y*, both of types *T*. For instance, with *x* and *y* of type
    ``fmpz_poly_t``, one can write ``FLINT_SWAP(fmpz_poly_struct, *x, *y)`` to
    swap the content of *x* with the content of *y*.

.. macro:: FLINT_SGN(x)

    Returns the sign of `x` where `x` is interpreted as a :type:`slong`, that
    is, returns `-1` if `x < 0`, `0` if `x = 0` and `1` if `x > 0`.

Integer types
-----------------------------------------------

The *char*, *short* and *int* types are assumed to be two's complement types
with exactly 8, 16 and 32 bits. Although this is not guaranteed prior to C23, it
is true on all mainstream platforms prior to this.

Since the C types *long* and *unsigned long* do not have a standardised size in
practice, FLINT defines *slong* and *ulong* types which are guaranteed to be 32
bits on a 32-bit system and 64 bits on a 64-bit system. They are also guaranteed
to have the same size as GMP's *mp_limb_t*. GMP builds with a different limb
size configuration are not supported at all.

.. type:: ulong

    The *ulong* type is used for integer-valued coefficients that are known to
    be unsigned, and for values that require the full 32-bit or 64-bit range.
    In method names, a *ulong* parameter is denoted by *ui*, for example
    :func:`arb_add_ui`.

    The constant *UWORD_MAX* gives the range of this type.
    This type can be printed with *flint_printf* using the format string ``%wu``.

    This is equivalent to GMP's *mp_limb_t*.

.. type:: slong

    The *slong* type is used for precisions, loop indices, array sizes, and the
    like, even when those values are known to be nonnegative. It is also used
    for small integer-valued coefficients. In method names, an *slong* parameter
    is denoted by *si*, for example :func:`arb_add_si`.

    This type can be printed with *flint_printf* using the format string ``%wd``
    or ``%{slong}``.

    This is equivalent to GMP's *mp_limb_signed_t*. Furthermore, for UNIX-type
    systems it is also equivalent to *mp_size_t*.

.. macro:: UWORD_MIN
           UWORD_MAX
           WORD_MIN
           WORD_MAX

    The minimum and maximum values that a *ulong* and *slong* can hold,
    respectively.

.. type:: flint_bitcnt_t

    A bit offset within an array of limbs (always nonnegative).

.. type:: nn_ptr

    Pointer to a writable array of limbs.

    This is equivalent to GMP's *mp_ptr*.

.. type:: nn_srcptr

    Pointer to a read-only array of limbs.

    This is equivalent to GMP's *mp_srcptr*.


Allocation Functions
-----------------------------------------------

.. function:: void * flint_malloc(size_t size)

   Allocate *size* bytes of memory.

.. function:: void * flint_realloc(void * ptr, size_t size)

   Reallocate an area of memory previously allocated by :func:`flint_malloc`,
   :func:`flint_realloc`, or :func:`flint_calloc`.

.. function:: void * flint_calloc(size_t num, size_t size)

   Allocate *num* objects of *size* bytes each, and zero the allocated memory.

.. function:: void flint_free(void * ptr)

   Free a section of memory allocated by  :func:`flint_malloc`,
   :func:`flint_realloc`, or :func:`flint_calloc`.


Random Numbers
------------------

.. type:: flint_rand_struct

    A structure holding the state of the FLINT pseudo random number generator.

.. type:: flint_rand_t

    An array of length 1 of :type:`flint_rand_struct`.

.. function:: void flint_rand_init(flint_rand_t state)
              void flint_rand_clear(flint_rand_t state)

    Initialises or clears a :type:`flint_rand_t`:.


Thread functions
-----------------------

.. function:: void flint_set_num_threads(int num_threads)

    Set up a thread pool of ``num_threads - 1`` worker threads (in addition
    to the master thread) and set the maximum number of worker threads the
    master thread can start to ``num_threads - 1``.

    This function may only be called globally from the master thread. It can
    also be called at a global level to change the size of the thread pool, but
    an exception is raised if the thread pool is in use (threads have been
    woken but not given back). The function cannot be called from inside
    worker threads.

.. function:: int flint_get_num_threads(void)

    When called at the global level, this function returns one more than the
    number of worker threads in the FLINT thread pool, i.e. it returns the
    number of workers in the thread pool plus one for the master thread.

    In general, this function returns one more than the number of additional
    worker threads that can be started by the current thread.

    Use :func:`thread_pool_wake` to set this number for a given worker thread.

    See also: :func:`flint_get_num_available_threads`.

.. function:: int flint_set_num_workers(int num_workers)

    Restricts the number of worker threads that can be started by the current
    thread to ``num_workers``. This function can be called from any thread.

    Assumes that the FLINT thread pool is already set up.

    The function returns the old number of worker threads that can be started.

    The function can only be used to reduce the number of workers that can be
    started from a thread. It cannot be used to increase the number. If a
    higher number is passed, the function has no effect.

    The number of workers must be restored to the original value by a call to
    :func:`flint_reset_num_workers` before the thread is returned to the thread
    pool.

    The main use of this function and :func:`flint_reset_num_workers` is to cheaply
    and temporarily restrict the number of workers that can be started, e.g. by
    a function that one wishes to call from a thread, and cheaply restore the
    number of workers to its original value before exiting the current thread.

.. function:: void flint_reset_num_workers(int num_workers)

    After a call to :func:`flint_set_num_workers` this function must be called to
    set the number of workers that may be started by the current thread back to
    its original value.

Input/Output
-----------------

.. function:: int flint_printf(const char * format, ...)
              int flint_fprintf(FILE * fs, const char * format, ...)
              int flint_vprintf(const char * format, va_list vlist)
              int flint_vfprintf(FILE * fs, const char * format, va_list vlist)

    These functions are extensions of the C standard library functions
    ``printf``, ``fprintf``, ``vprintf``, and ``vfprintf``.

    The first extension is the addition of the length modifier ``w``, used for
    printing the types :type:`ulong`, :type:`slong` and :type:`ulong`. As
    these types are either defined as signed and unsigned ``long int`` or
    ``long long int``, this comes in handy. Just like ``long int`` and ``long
    long int``, the conversion format specifier are allowed to be ``d``, ``i``,
    ``o``, ``x``, ``X`` and ``u``.

    The second and final extension is printing of FLINT types. Currently
    supported types are the base types :type:`ulong`, :type:`slong`,
    :type:`fmpz_t`, :type:`fmpq_t`, :type:`mag_t`, :type:`arf_t`, :type:`arb_t`
    and :type:`acb_t` as well as the context structures for modulo arithmetic
    :type:`nmod_t` and :type:`fmpz_mod_ctx_t`. We also support the GMP types
    ``mpz_t`` and ``mpq_t``.

    We currently support printing vectors of pointers to the following base
    types: :type:`slong`, :type:`ulong`, :type:`fmpz`, :type:`fmpq`,
    :type:`mag_struct`, :type:`arf_struct`, :type:`arb_struct` and
    :type:`acb_struct`.

    We also support printing matrices of the following types:
    :type:`nmod_mat_t`, :type:`fmpz_mat_t`, :type:`fmpq_mat_t`,
    :type:`arb_mat_t` and :type:`acb_mat_t`.

    Finally, we currently support printing polynomial of the following types:
    :type:`nmod_poly_t`, :type:`fmpz_poly_t`, :type:`fmpq_poly_t`,
    :type:`arb_poly_t` and :type:`acb_poly_t`.

.. code-block:: c

    ulong bulong;
    slong bslong;
    fmpz_t bfmpz;
    fmpq_t bfmpq;
    mag_t bmag;
    arf_t barf;
    arb_t barb;
    acb_t bacb;
    nmod_t bnmod;
    fmpz_mod_ctx_t bfmpz_mod_ctx;
    mpz_t bmpz;
    mpq_t bmpq;

    /* Initialize and set variables */

    flint_printf(
        "ulong: %{ulong}\n"
        "slong: %{slong}\n"
        "fmpz: %{fmpz}\n"
        "fmpq: %{fmpq}\n"
        "mag: %{mag}\n"
        "arf: %{arf}\n"
        "arb: %{arb}\n"
        "acb: %{acb}\n"
        "nmod: %{nmod}\n"
        "fmpz_mod_ctx: %{fmpz_mod_ctx}\n"
        "mpz: %{mpz}\n"
        "mpq: %{mpq}\n",
        bulong,
        bslong,
        bfmpz,
        bfmpq,
        bmag,
        barf,
        barb,
        bacb,
        bnmod,
        bfmpz_mod_ctx,
        bmpz,
        bmpq);

.. code-block:: c

    slong * vslong; slong vslong_len;
    nn_ptr vnmod; slong vnmod_len; /* The base type for nmod is ulong */
    fmpz * vfmpz; slong vfmpz_len;
    /* fmpz_mod vectors are given by the type `fmpz *' */
    fmpq * vfmpq; slong vfmpq_len;
    mag_ptr vmag; slong vmag_len;
    arf_ptr varf; slong varf_len;
    arb_ptr varb; slong varb_len;
    acb_ptr vacb; slong vacb_len;

    /* Initialize and set variables */

    flint_printf(
        "slong vector: %{slong*}\n"
        "nmod vector: %{ulong*}\n"
        "fmpz vector: %{fmpz*}\n"
        "fmpq vector: %{fmpq*}\n"
        "mag vector: %{mag*}\n"
        "arf vector: %{arf*}\n"
        "arb vector: %{arb*}\n"
        "acb vector: %{acb*}\n"
        vslong, vslong_len, /* They require a vector length specifier */
        vnmod, vnmod_len,
        vfmpz, vfmpz_len,
        vfmpq, vfmpq_len,
        vmag, vmag_len,
        varf, varf_len,
        varb, varb_len,
        vacb, vacb_len);

.. code-block:: c

    nmod_mat_t mnmod;
    fmpz_mat_t mfmpz;
    fmpz_mod_mat_t mfmpz_mod;
    fmpq_mat_t mfmpq;
    arb_mat_t marb;
    acb_mat_t macb;

    /* Initialize and set variables */

    flint_printf(
        "nmod matrix: %{nmod_mat}\n"
        "fmpz matrix: %{fmpz_mat}\n"
        "fmpz_mod matrix: %{fmpz_mod_mat}\n"
        "fmpq matrix: %{fmpq_mat}\n"
        "arb vector: %{arb_mat}\n"
        "acb vector: %{acb_mat}\n"
        mnmod,
        mfmpz,
        mfmpz_mod,
        mfmpq,
        marb,
        macb);

.. code-block:: c

    nmod_poly_t pnmod;
    fmpz_poly_t pfmpz;
    fmpz_mod_poly_t pfmpz_mod;
    fmpq_poly_t pfmpq;
    arb_poly_t parb;
    acb_poly_t pacb;

    /* Initialize and set variables */

    flint_printf(
        "nmod polynomial: %{nmod_poly}\n"
        "fmpz polynomial: %{fmpz_poly}\n"
        "fmpz_mod polynomial: %{fmpz_mod_poly}\n"
        "fmpq polynomial: %{fmpq_poly}\n"
        "arb polynomial: %{arb_poly}\n"
        "acb polynomial: %{acb_poly}\n"
        pnmod,
        pfmpz,
        pfmpz_mod,
        pfmpq,
        parb,
        pacb);

.. note::

    Printing of FLINT types does not currently support any flags.

.. note::

    Any use of ``%n`` flags will be invalid, but will not generate any error.

.. note::

    Invalid formats using variable minimum field width and/or precision such as
    ``"%* p"`` may be wrongly parsed, and may result in a different result
    compared to the C standard library functions.

.. function:: int flint_sprintf(char * s, const char * str, ...)

    This functions is an extensions of the C standard library functions
    ``sprintf``. It is currently advised to not use this function as it is
    currently not coherent with :func:`flint_printf`.

.. function:: int flint_scanf(const char * str, ...)
              int flint_fscanf(FILE * f, const char * str, ...)
              int flint_sscanf(const char * s, const char * str, ...)

     These are equivalent to the standard library functions ``scanf``,
     ``fscanf``, and ``sscanf`` with an additional length modifier "w" for
     reading an :type:`ulong` type.

Exceptions
-----------------

.. function:: void flint_abort(void)

    FLINT version of the C standard function ``abort``.

.. function:: void flint_set_abort(void (* func)(void))

    Sets the :func:`flint_abort` function to call ``func`` instead of
    ``abort``.

.. enum:: flint_err_t

    An error code with one of the following values

    .. macro:: FLINT_ERROR

        Describes a generic error.

    .. macro:: FLINT_OVERFLOW

        Describes an overflow.

    .. macro:: FLINT_IMPINV

        Describes an impossible inversion.

    .. macro:: FLINT_DOMERR

        Describes a domain error.

    .. macro:: FLINT_DIVZERO

        Describes a division by zero.

    .. macro:: FLINT_EXPOF

        Describes a exponent overflow.

    .. macro:: FLINT_INEXACT

        Describes an inexact operation.

    .. macro:: FLINT_TEST_FAIL

        Describes a test fail.

.. function:: void flint_throw(flint_err_t exc, const char * msg, ...)

    Throws an error of type ``exc`` with message ``msg`` and aborts via
    :func:`flint_abort`. The printing back-end function is
    :func:`flint_fprintf`, and so it allows for printing of FLINT types as
    well.
