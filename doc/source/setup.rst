.. _setup:

Setup
===============================================================================

Installation
-------------------------------------------------------------------------------

Calcium has the following dependencies:

* FLINT (http://www.flintlib.org) and its dependencies (GMP/MPIR and MPFR).
  Calcium will require FLINT 2.7 (unreleased) or later; currently
  a git checkout of https://github.com/wbhart/flint2 is necessary.
* Arb (http://arblib.org) version 2.18.1 or later.
* Antic (https://github.com/wbhart/antic/) - use a git checkout.

To compile, test and install Calcium from source assuming that
all dependencies have been installed, run the following
commands in the source directory::

    ./configure <options>
    make
    make check       (optional)
    make install

If GMP/MPIR, MPFR, FLINT, Arb or Antic are installed in some other
location than the default path ``/usr/local``, pass
``--with-gmp=...``, ``--with-mpfr=...``, ``--with-flint=...``,
``--with-arb=...``, ``--with-antic=...`` with
the correct path to configure (type ``./configure --help`` to show
more options).

After the installation, you may have to run ``ldconfig``
to make sure that the system's dynamic linker finds the library.

On a multicore system, ``make`` can be run with the ``-j`` flag to build
in parallel. For example, use ``make -j4`` on a quad-core machine.

Running tests
-------------------------------------------------------------------------------

After running ``make``, it is recommended to also run ``make check``
to verify that all unit tests pass.

By default, the unit tests run a large number of iterations to improve
the chances of detecting subtle problems.
The test suite will take several minutes on a single core
(``make -jN check`` if you have more cores to spare).
You can adjust the number of test iterations via
the ``CALCIUM_TEST_MULTIPLIER`` environment variable. For example, the following
will only run 10% of the default iterations::

    export CALCIUM_TEST_MULTIPLIER=0.1
    make check

It is also possible to run the unit tests for a single module, for instance::

    make check MOD=ca

Running code
-------------------------------------------------------------------------------

Here is an example program to get started using Calcium:

.. code-block:: c

    #include "calcium/ca.h"

    int main()
    {
        ca_ctx_t ctx;
        ca_t x;
        ca_ctx_init(ctx);
        ca_init(x, ctx);

        ca_pi(x, ctx);               /* x = pi */
        ca_sub_ui(x, x, 3, ctx);     /* x = x - 3 */
        ca_pow_ui(x, x, 2, ctx);     /* x = x^2 */
        ca_print(x, ctx); printf("\n");
        printf("Computed with calcium-%s\n", calcium_version());

        ca_clear(x, ctx);
        ca_ctx_clear(ctx);
        flint_cleanup();
        return EXIT_SUCCESS;
    }

Compile it with::

    gcc test.c -lcalcium

Depending on the environment, you may also have to pass
the flags ``-larb``, ``-lantic``, ``-lflint``, ``-lmpfr``, ``-lgmp``
to the compiler.
On some Debian based systems, ``-larb`` needs to be replaced
with ``-lflint-arb``.

If the header and library files are not in a standard location
(``/usr/local`` on most systems), you may also have to provide flags such as::

    -I/path/to/calcium -I/path/to/arb -I/path/to/flint -L/path/to/calcium -L/path/to/flint -L/path/to/arb

Finally, to run the program, make sure that the linker
can find the libraries. If they are installed in a
nonstandard location, you can for example add this path to the
``LD_LIBRARY_PATH`` environment variable.

The output of the example program should be something like the following::

    0.0200485 {a^2-6*a+9 where a = 3.14159 [Pi]}
    Computed with calcium-0.0.0


.. raw:: latex

    \newpage
