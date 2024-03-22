# Challenge 1

Luca Brambilla 10510718-919812

Change `PACS_ROOT` to have `muparser` work!

`make` compiles and executes.

## Folder structure

- `./include` contains all external libraries header files.

- `./lib` contains the shared objects, but probably not all of them... Self-contained use of `muparser` doesn't work.

- `./report` contains a short `.pdf` report and source LaTeX file.

- `main.cpp` source file.
- `optimization.hpp` utility containing minimization techniques.
- `muparser_fun.hpp` utility to evaluate parsed functions.
- `data.pot` contains the input parameters.
- `Makefile`