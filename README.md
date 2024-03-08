# Challenge 1

Luca Brambilla 10510718-919812

Change `PACS_ROOT` to have `muparser` work!

`make` compiles and executes.

In `Makefile` the we have the following definitions for `CPPFLAGS`:

- `GNUPLOT` to be able to use `gnuplot` directly from the source code.
- `PLOTSAVE` to save the plots to `.pdf` format.
- `GNUTEMP` to plot from a temporary file.
- `GNUPERM` to plot from an existing file.

## Folder structure

- `./include` contains all external libraries header files.

- `./output` contains the files resulting from each solution in the format `resuls_N###.dat.` (where `###` is the number of nodes of the discretization) and a file resulting from convergence analysis in `convergence.dat`.

- `./figures` contains the `.pdf` files with the plots from each solution in the format `plot_N###.pdf` (where `###` is the number of nodes of the discretization) and the ones resulting from convergence analyses in the format `conv_theta###.pdf` (where `###` is the value of `theta` used as input of the solver).

- `./lib` contains the shared objects, but probably not all of them... Self-contained use of `muparser` doesn't work.

- `./report` contains a short `.pdf` report and source LaTeX file.

- `main.cpp` source file.
- `solver_theta.hpp` utility containing the Newton solver.
- `muparser_fun.hpp` utility to evaluate parsed functions.
- `data.pot` contains the input parameters.
- `Makefile`