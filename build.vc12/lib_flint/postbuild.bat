@echo off
rem input is the output directory

set outdir=%1
if "%outdir%" EQU "" (set outdir=..\..\lib\x64\release)

call :copy_rename ..\..\flint.h %outdir% 

for %%i in ("..\..\*.h") do (call :copy_rename %%i %outdir% %%~nxi)

rem for %%i in ("..\..\flintxx\*.h") do (call :copy_rename %%i %outdir% %%~nxi)

call :copy_rename ..\..\qadic\CPimport.txt %outdir%

exit /b 0

rem copy rename 'in_file_name directory out_file_name''
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
