.. _contributing:

Contributing to Arb
===============================================================================

The Arb project welcomes new feature contributions in
addition to patches for bugs and performance problems.

What are appropriate new features?
Most of the numerical functionality that can be found in a general-purpose
computer algebra system is certainly within scope
(of course, the main restriction is that the algorithm
must have a proof of correctness).
However, if the functionality is easily accomplished by combining existing
functions in Arb, consider whether it is worth the increase in code size,
maintenance effort and test time.
Prospective contributors are recommended to discuss their ideas
on the mailing list or the issue tracker.

The process for actually submitting code is simple: anyone can submit a pull
request on GitHub. If the patch looks good to the maintainer and the test
code passes, the patch will get merged into the git master.

Code conventions
-------------------------------------------------------------------------------

Four steps are needed to add a new function:

* Add the function ``module_foo()`` in a new file ``module/foo.c``.
* Add a corresponding test program in a new file ``module/test/t-foo.c``.
* Add the function prototype to ``module.h``.
* Document the function in ``doc/source/module.rst``.

The build system takes care of everything else automatically.

Test code (see below)
can be omitted if ``module_foo()`` is a trivial helper function, but it should
at least be tested indirectly via another function in that case.
Auxiliary functions needed to implement ``module_foo()`` but which have no
use elsewhere should be declared as ``static`` in ``module/foo.c``.
If ``module_foo()`` is very short, it can be declared inline directly
in ``module.h`` with the ``MODULE_INLINE`` macro.

Use the following checklist regarding code style:

* Try to keep names and function arguments consistent with existing code.
* Follow the conventions regarding types, aliasing rules, etc. described
  in :ref:`issues` and in ``code_conventions.txt`` in FLINT (https://github.com/wbhart/flint2/blob/trunk/code_conventions.txt).
* Use basic FLINT constants, types and functions: ``FLINT_BITS``, ``flint_malloc``/``flint_free``, ``flint_abort``, ``flint_printf``, etc.
* Complex macros should be avoided.
* Indentation is four spaces.
* Curly braces normally go on a new line.
* Binary operators are surrounded by spaces (but parentheses and brackets are not).
* Logically distinct chunks of code (variable declarations, initialization,
  precomputations, the main loop, cleanup, etc.) should be separated by
  a single blank line.
* Lines are up to 79 characters long, but this rule can be broken if it helps readability.
* Add correct copyright notices at the top of each file.

Test code
-------------------------------------------------------------------------------

See :ref:`setup` for instructions on running test code.

The easiest way to write a test program for a new function
is to adapt the test code for an existing, similar function.

Most of the test code in Arb uses the strategy of computing the same
mathematical quantity in two or more different ways (for example, using
functional equations, interchanging the order of parameter, or varying
the precision and other algorithm parameters) and verifying that
the results are consistent.
It is also a good idea to test that aliasing works.
Input data is usually generated randomly, but in some cases
including precomputed reference values also makes sense.

Faster test code is better. A single test program should not take more
than 10 seconds to run, and preferably no more than 1 second. Most functions
can be tested effectively in less than 0.1 seconds. Think of what the corner
cases are and try to generate random input biased toward such cases.
The ``randtest()`` functions attempt to generate corner cases automatically, but
some thought may be needed to use them optimally. Try to ensure that the test
code fails if you deliberately break the tested function in any way. It is also
a good idea to run the test code once with ``ARB_TEST_MULTIPLIER=10.0`` or higher.
If a function's input space is too large to probe effectively for corner cases
with random input, that can be a hint that the function should be split into
smaller logical parts that can be tested separately.

The test code must complete without errors when run with ``valgrind``.
The most common mistake leading to memory corruption or memory leaks
is to miss or duplicate an ``init()`` or ``clear()`` call.
Check that the ``init()`` and ``clear()`` calls exactly match the variable
declarations in each code block, including the test code itself.

Profiling code is not needed in most cases, but it is often a good idea to
run some benchmarks at least during the initial development of a new feature.
The ``TIMEIT_START``/``TIMEIT_STOP`` and ``SHOW_MEMORY_USAGE`` macros
in FLINT are useful for quick measurements.


