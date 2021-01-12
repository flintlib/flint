.. _fexpr:

**fexpr.h** -- flat-packed symbolic expressions
===============================================================================

This module supports working with symbolic expressions. Formally,
a symbolic expression is either:

* An atom, being one of the following:

  * An integer, for example 0 or -34.

  * A symbol, for example ``x``, ``Mul``, ``SomeUserNamedSymbol``.
    Symbols should be valid C identifiers (containing only the
    characters ``A-Z``, ``a-z``, ``0-9``, ``_``,
    and not starting with a digit).

  * A string, for example ``"Hello, world!"``. For the moment, we
    only consider ASCII strings, but there is no obstacle in
    principle to supporting UTF-8.

* A non-atomic expression, representing a function call
  `e_0(e_1, \ldots, e_n)` where `e_0, \ldots, e_n` are symbolic
  expressions.

The meaning of an expression depends on the interpretation
of symbols in a given context.
For example, with a standard intepretation (used within Calcium) of the symbols
``Mul``, ``Add`` and
``Neg``, the expression ``Mul(3, Add(Neg(x), y))``
encodes the formula `3 \cdot ((-x)+y)`
where ``x`` and ``y`` are symbolic variables.

Flat-packed representation
-----------------------------------------------------------------------

Symbolic expressions are often implemented using trees of pointers
(often together with hash tables for uniqueness),
requiring some form of memory management.
The :type:`fexpr_t` type, by contrast, stores a symbolic expression
using a "flat-packed" representation without internal pointers.
The expression data is just an array of words (``ulong``).
The first word is a header encoding type information (whether
the expression is a function call or an atom, and the type
of the atom) and the total number of words
in the expression.
For atoms, the data is stored either in the header word itself (small
integers and short symbols/strings) or in the following words.
For function calls, the header is followed by
the expressions `e_0`, ..., `e_n` packed contiguously in memory.

Pros:

* Memory management is trivial.
* Copying an expression is just copying an array of words.
* Comparing expressions for equality is just comparing arrays of words.
* Merging expressions is basically just concatenating arrays of words.
* Expression data can be shared freely in binary form between
  threads and even between machines (as long as all machines
  have the same word size and endianness).

Cons:

* Repeated instances of the same subexpression cannot share memory
  (a workaround is to introduce local dummy symbols for repeated
  subexpressions).
* Extracting a subexpression for modification generally
  requires making a complete
  copy of that subxepression (however, for read-only access
  to subexpressions, one can use "view" expressions which have
  zero overhead).
* Manipulating a part of an expression generally requires rebuilding
  the whole expression.
* Building an expression incrementally is typically `O(n^2)`.
  As a workaround, it is a good idea to work with balanced (low-depth)
  expressions and try to construct an expression in one go
  (for example, to create a sum, create a single ``Add`` expression
  with many arguments instead of chaining binary ``Add`` operations).

Types and macros
-------------------------------------------------------------------------------

.. type:: fexpr_struct

.. type:: fexpr_t

    An *fexpr_struct* consists of a pointer to an array of words along
    with a record of the number of allocated words.

    An *fexpr_t* is defined as an array of length one of type
    *fexpr_struct*, permitting an *fexpr_t* to be passed by
    reference.

.. type:: fexpr_ptr

   Alias for ``fexpr_struct *``, used for arrays of expressions.

.. type:: fexpr_srcptr

   Alias for ``const fexpr_struct *``, used for arrays of expressions
   when passed as constant input to functions.

.. type:: fexpr_vec_struct

.. type:: fexpr_vec_t

    A type representing a vector of expressions with managed length.
    The structure contains an :type:`fexpr_ptr` *entries* for
    the entries, an integer *length* (the size of the vector), and
    an integer *alloc* (the number of allocated entries).

.. macro:: fexpr_vec_entry(vec, i)

    Returns a pointer to entry *i* in the vector *vec*.

Memory management
-------------------------------------------------------------------------------

.. function:: void fexpr_init(fexpr_t expr)

    Initializes *expr* for use. Its value is set to the atomic
    integer 0.

.. function:: void fexpr_clear(fexpr_t expr)

    Clears *expr*, freeing its allocated memory.

.. function:: fexpr_ptr _fexpr_vec_init(slong len)

    Returns a heap-allocated vector of *len* initialized expressions.

.. function:: void _fexpr_vec_clear(fexpr_ptr vec, slong len)

    Clears the *len* expressions in *vec* and frees *vec* itself.

.. function:: void fexpr_fit_size(fexpr_t expr, slong size)

    Ensures that *expr* has room for *size* words.

.. function:: slong fexpr_size(const fexpr_t expr)

    Returns the number of words in the internal representation
    of *expr*.

.. function:: void fexpr_set(fexpr_t res, const fexpr_t expr)

    Sets *res* to the a copy of *expr*.

.. function:: void fexpr_swap(fexpr_t a, fexpr_t b)

    Swaps *a* and *b* efficiently.

Comparisons
-------------------------------------------------------------------------------

.. function:: int fexpr_equal(const fexpr_t a, const fexpr_t b)

    Checks if *a* and *b* are exactly equal as expressions.

.. function:: ulong fexpr_hash(const fexpr_t expr)

    Returns a hash of the expression *expr*.

.. function:: int fexpr_cmp_fast(const fexpr_t a, const fexpr_t b)

    Compares *a* and *b* using an ordering based on the internal
    representation, returning -1, 0 or 1. This can be used, for
    instance, to maintain sorted arrays of expressions for binary
    search; the sort order has no mathematical significance.


Atoms
-------------------------------------------------------------------------------

.. function:: int fexpr_is_integer(const fexpr_t expr)

    Returns whether *expr* is an atomic integer

.. function:: int fexpr_is_symbol(const fexpr_t expr)

    Returns whether *expr* is an atomic symbol.

.. function:: int fexpr_is_string(const fexpr_t expr)

    Returns whether *expr* is an atomic string.

.. function:: int fexpr_is_atom(const fexpr_t expr)

    Returns whether *expr* is any atom.

.. function:: void fexpr_zero(fexpr_t res)

    Sets *res* to the atomic integer 0.

.. function:: int fexpr_is_zero(const fexpr_t expr)

    Returns whether *expr* is the atomic integer 0.

.. function:: void fexpr_set_si(fexpr_t res, slong c)
              void fexpr_set_ui(fexpr_t res, ulong c)
              void fexpr_set_fmpz(fexpr_t res, const fmpz_t c)

    Sets *res* to the atomic integer *c*.

.. function:: void fexpr_get_fmpz(fmpz_t res, const fexpr_t expr)

    Sets *res* to the atomic integer in *expr*. This aborts
    if *expr* is not an atomic integer.

.. function:: void fexpr_set_symbol_str(fexpr_t res, const char * s)

    Sets *res* to the symbol given by *s*.

.. function:: char * fexpr_get_symbol_str(const fexpr_t expr)

    Returns the symbol in *expr* as a string. The string must
    be freed with :func:`flint_free`.
    This aborts if *expr* is not an atomic symbol.

Input and output
------------------------------------------------------------------------

.. function:: void fexpr_write(calcium_stream_t stream, const fexpr_t expr)

    Writes *expr* to *stream*.

.. function:: void fexpr_print(const fexpr_t expr)

    Prints *expr* to standard output.

.. function:: char * fexpr_get_str(const fexpr_t expr)

    Returns a string representation of *expr*. The string must
    be freed with :func:`flint_free`.

Function call structure
------------------------------------------------------------------------

.. function:: slong fexpr_nargs(const fexpr_t expr)

    Returns the number of arguments *n* in the function call
    `f(e_1,\ldots,e_n)` represented
    by *expr*. If *expr* is an atom, returns -1.

.. function:: void fexpr_func(fexpr_t res, const fexpr_t expr)

    Assuming that *expr* represents a function call
    `f(e_1,\ldots,e_n)`, sets *res* to the function expression *f*.

.. function:: void fexpr_view_func(fexpr_t view, const fexpr_t expr)

    As :func:`fexpr_func`, but sets *view* to a shallow view
    instead of copying the expression.
    The variable *view* must not be initialized before use or
    cleared after use, and *expr* must not be modified or cleared
    as long as *view* is in use.

.. function:: void fexpr_arg(fexpr_t res, const fexpr_t expr, slong i)

    Assuming that *expr* represents a function call
    `f(e_1,\ldots,e_n)`, sets *res* to the argument `e_{i+1}`.
    Note that indexing starts from 0.
    The index must be in bounds, with `0 \le i < n`.

.. function:: void fexpr_view_arg(fexpr_t view, const fexpr_t expr, slong i)

    As :func:`fexpr_arg`, but sets *view* to a shallow view
    instead of copying the expression.
    The variable *view* must not be initialized before use or
    cleared after use, and *expr* must not be modified or cleared
    as long as *view* is in use.

.. function:: void fexpr_view_next(fexpr_t view)

    Assuming that *view* is a shallow view of a function argument `e_i`
    in a function call `f(e_1,\ldots,e_n)`, sets *view* to
    a view of the next argument `e_{i+1}`.
    This function can be called when *view* refers to the last argument
    `e_n`, provided that *view* is not used afterwards.
    This function can also be called when *view* refers to the function *f*,
    in which case it will make *view* point to `e_1`.


Composition
------------------------------------------------------------------------

.. function:: void fexpr_call0(fexpr_t res, const fexpr_t f)
              void fexpr_call1(fexpr_t res, const fexpr_t f, const fexpr_t x1)
              void fexpr_call2(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2)
              void fexpr_call3(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3)
              void fexpr_call4(fexpr_t res, const fexpr_t f, const fexpr_t x1, const fexpr_t x2, const fexpr_t x3, const fexpr_t x4)
              void fexpr_call_vec(fexpr_t res, const fexpr_t f, fexpr_srcptr args, slong len)

    Creates the function call `f(x_1,\ldots,x_n)`.
    The *vec* version takes the arguments as an array *args*
    and *n* is given by *len*.

Arithmetic expressions
------------------------------------------------------------------------

.. function:: void fexpr_set_fmpq(fexpr_t res, const fmpq_t x)

    Sets *res* to the rational number *x*. This creates an atomic
    integer if the denominator of *x* is one, and otherwise creates a
    division expression.

.. function:: void fexpr_neg(fexpr_t res, const fexpr_t a)
              void fexpr_add(fexpr_t res, const fexpr_t a, const fexpr_t b)
              void fexpr_sub(fexpr_t res, const fexpr_t a, const fexpr_t b)
              void fexpr_mul(fexpr_t res, const fexpr_t a, const fexpr_t b)
              void fexpr_div(fexpr_t res, const fexpr_t a, const fexpr_t b)
              void fexpr_pow(fexpr_t res, const fexpr_t a, const fexpr_t b)

    Constructs an arithmetic expression with given arguments.
    No simplifications whatsoever are performed.

.. function:: int fexpr_is_arithmetic_operation(const fexpr_t expr)

    Returns whether *expr* is of the form `f(e_1,\ldots,e_n)`
    where *f* is one of the arithmetic operators ``Pos``, ``Neg``,
    ``Add``, ``Sub``, ``Mul``, ``Div``.

.. function:: void fexpr_arithmetic_nodes(fexpr_vec_t nodes, const fexpr_t expr)

    Sets *nodes* to a vector of subexpressions of *expr* such that *expr*
    is an arithmetic expression with *nodes* as leaves.
    More precisely, *expr* will be constructed out of nested application
    the arithmetic operators
    ``Pos``, ``Neg``, ``Add``, ``Sub``, ``Mul``, ``Div`` with
    integers and expressions in *nodes* as leaves.
    Powers ``Pow`` with an atomic integer exponent are also allowed.
    The nodes are output without repetition but are not automatically sorted in
    a canonical order.

.. function:: int fexpr_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fexpr_t expr, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to the expression *expr* as a formal rational
    function of the subexpressions in *vars*.
    The vector *vars* must have the same length as the number of
    variables specified in *ctx*.
    To build *vars* automatically for a given expression,
    :func:`fexpr_arithmetic_nodes` may be used.

    Returns 1 on success and 0 on failure. Failure can occur for the
    following reasons:

    * A subexpression is encountered that cannot be interpreted
      as an arithmetic operation and does not appear (exactly) in *vars*.
    * Overflow (too many terms or too large exponent).
    * Division by zero (a zero denominator is encountered).

    It is important to note that this function views *expr* as
    a formal rational function with *vars* as formal indeterminates.
    It does thus not check for algebraic relations between *vars*
    and can implicitly divide by zero if *vars* are not algebraically
    independent.

.. function:: void fexpr_set_fmpz_mpoly(fexpr_t res, const fmpz_mpoly_t poly, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)
              void fexpr_set_fmpz_mpoly_q(fexpr_t res, const fmpz_mpoly_q_t frac, const fexpr_vec_t vars, const fmpz_mpoly_ctx_t ctx)

    Sets *res* to an expression for the multivariate polynomial *poly*
    (or rational function *frac*),
    using the expressions in *vars* as the variables. The length
    of *vars* must agree with the number of variables in *ctx*.
    If *NULL* is passed for *vars*, a default choice of symbols
    is used.

.. function:: int fexpr_expanded_normal_form(fexpr_t res, const fexpr_t expr, ulong flags)

    Sets *res* to *expr* converted to expanded normal form viewed
    as a formal rational function with its non-arithmetic subexpressions
    as terminal nodes.
    This function first computes nodes with :func:`fexpr_arithmetic_nodes`,
    sorts the nodes, evaluates to a rational function with
    :func:`fexpr_get_fmpz_mpoly_q`, and then converts back to an
    expression with :func:`fexpr_set_fmpz_mpoly_q`.
    Optional *flags* are reserved for future use.


Vectors
------------------------------------------------------------------------

.. function:: void fexpr_vec_init(fexpr_vec_t vec)

    Initializes *vec* to a vector of length zero.

.. function:: void fexpr_vec_clear(fexpr_vec_t vec)

    Clears the vector *vec*.

.. function:: void fexpr_vec_print(const fexpr_vec_t vec)

    Prints *vec* to standard output.

.. function:: void fexpr_vec_swap(fexpr_vec_t x, fexpr_vec_t y)

    Swaps *x* and *y* efficiently.

.. function:: void fexpr_vec_fit_length(fexpr_vec_t vec, slong len)

    Ensures that *vec* has space for *len* entries.

.. function:: void fexpr_vec_set(fexpr_vec_t dest, const fexpr_vec_t src)

    Sets *dest* to a copy of *src*.

.. function:: void fexpr_vec_append(fexpr_vec_t vec, const fexpr_t expr)

    Appends *expr* to the end of the vector *vec*.

.. function:: slong fexpr_vec_insert_unique(fexpr_vec_t vec, const fexpr_t expr)

    Inserts *expr* without duplication into vec, returning its
    position. If this expression already exists, *vec* is unchanged.
    If this expression does not exist in *vec*, it is appended.

.. function:: void fexpr_vec_set_length(fexpr_vec_t vec, slong len)

    Sets the length of *vec* to *len*, truncating or zero-extending as needed.

.. function:: void _fexpr_vec_sort_fast(fexpr_ptr vec, slong len)

    Sorts the *len* entries in *vec* using 
    the comparison function :func:`fexpr_cmp_fast`.

.. raw:: latex

    \newpage
