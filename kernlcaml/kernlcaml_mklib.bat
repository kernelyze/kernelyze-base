REM Copyright (c) 2017 by Kernelyze LLC
REM Author: Thomas A. Knox
REM This program is free software: you can redistribute it and/or modify
REM it under the terms of the GNU Affero General Public License as
REM published by the Free Software Foundation, either version 3 of the
REM License, or (at your option) any later version.
REM
REM This program is distributed in the hope that it will be useful,
REM but WITHOUT ANY WARRANTY; without even the implied warranty of
REM MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
REM GNU Affero General Public License for more details.
REM 
REM You should have received a copy of the GNU Affero General Public License
REM along with this program.  If not, see <http://www.gnu.org/licenses/>.

@echo off

REM This file is in C:\Users\thoma\Documents\OCamlProjects\kernlcaml

set intel_libs=C:\intel64libs
set kernelyze_c=C:\Users\thoma\Documents\Kernelyze\KernelyzeLibs\KernelyzeBaseBindC\KernelyzeBaseBindC\x64\Release\KernelyzeBaseBindC.lib

echo Compiling C files
cl /c /MD /I "C:\ocamlms64\lib" kernlcaml_stubs.c
cl /c /MD /I "C:\ocamlms64\lib" mainc.c

echo Compiling OCaml files to bytecode
ocamlc -c kernlcaml.ml bigarray.cma

echo Compiling OCaml files to native code 
ocamlopt -c kernlcaml.ml bigarray.cmxa


echo Making libraries
REM I had to be quite careful in two regards here: first, it is important to pass bigarray.lib and not any cmx or cma of bigarray, and second, it is 
REM critical to avoid spaces in paths to libraries -- quoting conventions seems to be handled differently by different tools, so even  with \" escaped quotes
REM a path with spaces will create problems downstream when ocamlopt is called to build an executable with the cmxa native-code library.  Enclosing
REM in multiple nested \" does not help.
ocamlmklib -verbose -o kernlcaml -L"%intel_libs%" -I "C:\ocamlms64\lib" -L"C:\ocamlms64\lib" kernlcaml.ml kernlcaml_stubs.obj mainc.obj C:\ocamlms64\lib\bigarray.lib %kernelyze_c%

echo Compiling custom OCaml toplevel
ocamlmktop -o kernlcaml.exe -custom -I "%intel_libs%" bigarray.cma kernlcaml.ml kernlcaml_stubs.obj mainc.obj %kernelyze_c%

REM To use in custom toplevel, run "kernlcaml"

echo Building test programs
ocamlc -o kernlcaml_test_noscript.byte bigarray.cma kernlcaml.cma kernlcaml_test_noscript.ml
ocamlopt -o kernlcaml_test_noscript.exe bigarray.cmxa kernlcaml.cmxa kernlcaml_test_noscript.ml mainc.obj

echo Running test programs
echo Bytecode test:
ocamlrun kernlcaml_test_noscript.byte
echo Native-code test:
kernlcaml_test_noscript