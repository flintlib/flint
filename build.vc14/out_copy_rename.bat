@echo off
rem %1 = source file path
rem %2 = destination directory
rem %3 = destination file name

if not exist %1 goto nofile
if exist %2 goto next

echo creating directory %2
md %2 > nul

:next
rem strip quotes if present
set str=%2
for /f "useback tokens=*" %%a in ('%str%') do set str=%%~a

rem add a backslash if the output directory lacks one
set str=%str:~-1%
if "%str%" == "\" (set outf=%2%3) else (set outf=%2\%3)

set op=copying
if not exist "%outf%" goto copy

rem don't overwrite if output exists and is not changed
fc %1 %outf% > nul && if not %errorlevel 1 goto exit
set op=overwriting

:copy
if "%4" NEQ "" (echo %op% %outf% from %1)
copy %1 %outf% > nul
goto exit

:nofile
echo %1 not found

:exit
