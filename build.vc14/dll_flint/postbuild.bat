@echo off

if /i %2 EQU LIB (set two=lib)
if /i %2 EQU DLL (set two=dll)

if not exist ..\..\%two% (md ..\..\%two%)
call :copy_files %1 ..\..\%two%\%1 %two%
exit /b 0

:copy_files
call :copy_rename ..\..\flint.h %2 > nul 2>&1
call :copy_rename ..\config.h %2 > nul 2>&1

rem copy C/C++ header files
for %%i in ("..\..\*.h") do (call :copy_rename %%i %2 %%~nxi > nul 2>&1)

rem copy C++ headers in flintxx
rem for %%i in ("..\..\flintxx\*.h") do (call :copy_rename %%i %2 %%~nxi > nul 2>&1)

rem copy the FLINT static library and related files
if /i %3 EQU LIB (
    if exist ..\%1\lib_flint.lib (
        echo ..\%1\lib_flint.lib
        call :copy_rename ..\%1\lib_flint.lib %2 > nul 2>&1  	
        if exist ..\%1\lib_flint.pdb (
            call :copy_rename ..\%1\lib_flint.pdb %2 > nul 2>&1
            )
        )
)

rem copy the FLINT dynamic library and related files
if /i %3 EQU DLL (
    if exist ..\%1\dll_flint.dll (
        call :copy_rename ..\%1\dll_flint.dll %2 > nul 2>&1
        if exist ..\%1\dll_flint.pdb (
            call :copy_rename ..\%1\dll_flint.pdb %2 > nul 2>&1
            )
        if exist ..\%1\dll_flint.lib (
            call :copy_rename ..\%1\dll_flint.lib %2 > nul 2>&1
            )
        if exist ..\%1\dll_flint.exp (
            call :copy_rename ..\%1\dll_flint.exp %2 > nul 2>&1
            )
        )
)
exit /b 0

rem copy rename 'in_file_name directory out_file_name'
:copy_rename
@echo off
if not exist %1 goto cr_nofile

if "%3" EQU "" (set outname=%~nx1) else (set outname=%3)
if exist %2 goto cr_next

echo creating directory %2
md %2 > nul

:cr_next
rem strip quotes if present
set str=%2
for /f "useback tokens=*" %%a in ('%str%') do set str=%%~a

rem add a backslash if the output directory lacks one
set str=%str:~-1%

if "%str%" == "\" (set outf=%2%outname%) else (set outf=%2\%outname%)
if exist "%outf%" goto cr_check

echo copying %1 to %outf%
copy %1 %outf% > nul
goto cr_exit

:cr_check
rem don't overwrite if output exists and is not changed
fc %1 %outf% > nul && if not %errorlevel 1 goto cr_exit
echo overwriting...... %outf% with %1
copy %1 %outf% > nul
goto cr_exit

:cr_nofile
echo %1 not found

:cr_exit
exit /b 0
