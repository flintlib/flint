Building FLINT2 with Microsoft Visual Studio 2015
-------------------------------------------------

To build FLINT2 with Microsoft Visual Studio 2015 you will
need Visual Studio 2015 Community (or a higher version) and:

   a) an installed version of Python 3
   b) an installed version of Python Tools for 
      Visual Studio (http://pytools.codeplex.com/)

1. Clone FLINT2 using GIT from:

     git@github.com:BrianGladman/flint2.git

or unpack a released distribution.
   
2. Open the solution file flint.sln in Visual Studio 2015. This
   includes the projects:

     dll_flint     - a Visual Studio build project for a
                     FLINT2 Dynamic Link Library
     lib_flint     - a Visual Studio build project for a
                     fLINT2 Static Library
     flint_config  - a Python program for creating the Visual 
                     Studio build files
     build_tests   - a Python program for building the FLINT2
                     tests (after they have been created)
     run_tests     - a Python program for running the FLINT2
                     tests (after they have been built)

3. If you have Python and Python Tools for Visual Studio (PTVS) 
   installed, set the project flint_config as the start-up project
   and run it.  If you don't have PTVS installed but you do have
   Python, you can run flint_config.py directly without Visual
   Studio.

   By default flint_config is set to create the FLINT2 tests and
   profiling but it can also create the Visual Studio 2015 build
   files for FLINT2 DLL and Static Libraries by changing these
   defines at the start of flint_config.py: 

      build_lib = False
      build_dll = False
      build_tests = True
      build_profiles = True

   Rebuilding the library build files will be necessary if FLINT2
   has been updated since it was first downloaded. 

   Without Python you can use the pre-installed lib_flint and 
   dll_flint projects but these may not reflect the latest changes
   in FLINT2 (note that the FLINT tests and profiling applications
   are not pre-installed so these won't be available)
   
4. The lib_flint and dll_flint projects build the static library 
   and dynamic library versions of FLINT2.  Before building one
   or both of these, you need to select the architecture <Win32|x64>
   and the build type <Release|Debug>.

5. After the FLINT2 tests have been created using flint_config.py,
   they can then be built by running the build_tests.py program.

6. There are also a number of Visual Studio solution files that
   provide an *alternative* way of building the FLINT2 tests and
   profiling.  Their use is not recommended however because each
   of the multiple solution files flint-tests<NN>.sln (where NN
   is a number) has to be loaded and built by Visual Studio (this
   approach is used because it takes Visual Studio too long to
   load the tests from a single solution file).

7. Once the tests have been built, the Python project run_tests can be
   set as the start-up project and started to run all the tests (or 
   the Python file run_tests.py can be run outside Visual Studio).
   
7. After building FLINT2, the libraries and the header files you need
   to use FLINT2 are placed in the directories:
   
      lib\<Win32|x64>\<Debug|Release>
      dll\<Win32|x64>\<Debug|Release>

   depending on the version(s) that have been built (<a|b> here means
   either 'a' or 'b' depending on which versions have been built.
   
      Brian Gladman
      7th August 2015
