Building FLINT2 with Microsoft Visual Studio 2013
-------------------------------------------------

1. Clone FLINT2 using GIT from:

     git@github.com:BrianGladman/flint2.git
   
2. Open the solution file flint2.sln in Visual Studio 2013

3. If you have Python and you have installed Python Tools for 
   Visual Studio (PTVS) installed from:
   
     http://pytools.codeplex.com/
     
   set the project flint_config as the start-up project and 
   run it.
   
   If you don't have PTVS installed but you have Python, run
   flint_config.py outside Visual Studio and then open flint2.sln
   in Visual Studio. 
   
   Without either of these tools you can simply use the pre-installed
   lib_flint and dll_flint projects but these may not reflect the 
   latest changes in FLINT2 (note that the FLINT tests and profiling 
   applications are not pre-installed so these won't be available)
   
4. The projects lib_flint and dll_flint build the static library 
   and dynamic library versions of FLINT2.  Before building one
   or both of these, you need to select the architecture (win32|x64)
   and the build type (Release|Debug).
   
5. There are 14 solution files for the FLINT2 tests and one for the
   FLINT2 profiles (there are 14 test solutions because it takes 
   Visual Studio too long to load the 1348 tests from a single
   solution file - this step needs to be recast as a makefile once I 
   have worked out how to do this).  To run the tests each of these
   14 solution files need to be loaded and built in Visual Studio.
   
6. Once the tests have been built, the Python project run_tests can be
   set as the start-up project and started to run all the tests (or 
   the Python file run_tests.py can be run outside Visual Studio
   
7. After building FLINT2, the libraries and the header files you need
   to use FLINT2 are placed in the directories:
   
      lib\<arch>\<config>\
      dll\<arch>\<config>\
   
   depending on the version(s) that have been built (<arch> is either
   'Win32' or 'x64' and  config is either 'Release' or 'Debug').
   
   
      Brian Gladman
      22nd August 2014
