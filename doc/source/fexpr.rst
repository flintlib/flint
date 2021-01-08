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

   Alias for ``fexpr_struct *``, used for vectors of expressions.

.. type:: fexpr_srcptr

   Alias for ``const fexpr_struct *``, used for vectors of expressions
   when passed as constant input to functions.

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

.. function:: void fexpr_set_si(fexpr_t res, slong c)
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

.. function:: void fexpr_print(const fexpr_t expr)

    Prints *expr* to standard output.

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


.. raw:: latex

    \newpage
