@echo off
setlocal
set PROJECT=IBLBaker

if exist Build64 (rd /s /q Build64)

mkdir Build64
cd Build64

cmake.exe -Wno-dev -G "Visual Studio 15 2017 Win64" ../

endlocal

pause
