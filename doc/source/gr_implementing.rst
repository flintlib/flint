.. _gr-implementing:

**gr.h (continued)** -- implementing rings
===============================================================================

.. highlight:: c

Defining a ring requires putting appropriate data into a
:type:`gr_ctx_t` parent object, most importantly the method table
and the size of elements.

Example
--------------------------------------------------------------------------------

This is an extract from the ``fmpz`` wrapper in ``gr/fmpz.c``::

    /* Some methods */
    ...
    int
    _gr_fmpz_add(fmpz_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
    {
        fmpz_add(res, x, y);
        return GR_SUCCESS;
    }
    ...

    /* The method table */

    int _fmpz_methods_initialized = 0;

    gr_static_method_table _fmpz_methods;

    gr_method_tab_input _fmpz_methods_input[] =
    {
        {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
        ...
        {GR_METHOD_INIT,            (gr_funcptr) _gr_fmpz_init},
        {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fmpz_clear},
        ...
        {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_fmpz_add},
        ...
        {0,                         (gr_funcptr) NULL},
    };

    /* Context object initializer */

    void
    gr_ctx_init_fmpz(gr_ctx_t ctx)
    {
        ctx->which_ring = GR_CTX_FMPZ;
        ctx->sizeof_elem = sizeof(fmpz);
        ctx->size_limit = WORD_MAX;

        ctx->methods = _fmpz_methods;

        if (!_fmpz_methods_initialized)
        {
            gr_method_tab_init(_fmpz_methods, _fmpz_methods_input);
            _fmpz_methods_initialized = 1;
        }
    }

Note that the method table only has to be constructed once,
allowing new context objects for the same domain
to be initialized cheaply.

Method table
--------------------------------------------------------------------------------

.. type:: gr_funcptr

    Typedef for a pointer to a function with signature ``int func(void)``,
    used to represent method table entries.

.. type:: gr_method

    Enumeration type for indexing method tables. Enum values named
    ``GR_METHOD_INIT``,  ``GR_METHOD_ADD_UI``, etc.
    correspond to methods ``gr_init``, ``gr_add_ui``, etc.
    The number of methods is given by ``GR_METHOD_TAB_SIZE``,
    which can be used to declare static method tables.

.. type:: gr_static_method_table

    Typedef for an array of length ``GR_METHOD_TAB_SIZE``
    with :type:`gr_funcptr` entries.

.. type:: gr_method_tab_input

    Typedef representing a (index, function pointer) pair.

.. function:: void gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab)

    Initializes the method table *methods*. This first inserts
    default and generic methods in all slots, and then overwrites
    with the specialized methods listed in *tab*.

Placeholder and trivial methods
--------------------------------------------------------------------------------

.. function:: int gr_not_implemented(void)

    This function does nothing and returns ``GR_UNABLE``. It is used
    as a generic fallback method when no implementation is available.

.. function:: int gr_not_in_domain(void)

    This function does nothing and returns ``GR_DOMAIN``. It can be used
    for an operation that never makes sense in the present domain,
    e.g., for the constant `\pi` in the rational numbers.

.. function:: truth_t gr_generic_ctx_predicate(gr_ctx_t ctx)

    Does nothing and returns ``T_UNKNOWN``, used as a generic
    fallback for predicate methods.

.. function:: truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx)

    A predicate that does nothing and returns ``T_TRUE``.

.. function:: truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx)

    A predicate that does nothing and returns ``T_FALSE``.


Required methods
--------------------------------------------------------------------------------

A context object must at minimum define the following methods for a ring:

* ``init``
* ``clear``
* ``swap``
* ``randtest``
* ``write``
* ``zero``
* ``one``
* ``equal``
* ``set``
* ``set_si``
* ``set_ui``
* ``set_fmpz``
* ``neg``
* ``add``
* ``sub``
* ``mul``

Other methods have generic defaults which may be
overridden for performance or completeness.

Implementing context predicates (``ctx_is_integral_domain``, ``ctx_is_field``, etc.)
is strongly recommended so that the most appropriate algorithms
can be used in generic implementations.

Rings with cheap operations on single elements should also provide
non-generic versions of performance-critical vector operations
to minimize overhead. The most important vector operations include:

* ``vec_init``
* ``vec_clear``
* ``vec_swap``
* ``vec_zero``
* ``vec_neg``
* ``vec_add``
* ``vec_sub``
* ``vec_mul_scalar_ui``/``si``
* ``vec_addmul_scalar_ui``/``si``
* ``vec_dot``
* ``vec_dot_rev``

Dot products, for example, are the main building block for
classical polynomial multiplication and matrix multiplication.
The methods

* ``poly_mullow``
* ``matrix_mul``

should be overridden for rings where faster-than-classical polynomial and
matrix multiplication is possible.
Other higher-complexity generic algorithms will try to
reduce to polynomial and matrix multiplication automatically, but may
in turn need to be overridden to select accurate
cutoffs between different algorithms.

Testing rings
--------------------------------------------------------------------------------

.. function:: void gr_test_ring(gr_ctx_t R, slong iters, int test_flags)

    Test correctness of the ring *R*. This calls test functions for
    various methods, each being repeated up to *iters* times.



.. raw:: latex

    \newpage
