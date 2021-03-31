# PGMG: An MPI-based parallel geometric multigrid library in C++

PGMG is a software library for solving a class of [elliptic partial differential
equations](https://en.wikipedia.org/wiki/Elliptic_partial_differential_equation)
on three-dimensional rectangular grids using the [multigrid
method](https://en.wikipedia.org/wiki/Multigrid_method).

## Compiling the code examples
The code is written in C++ and has been tested on Linux, MacOS, and Windows
(via [Cygwin](https://www.cygwin.com)).The following documentation assumes you
are familiar with the Linux/Mac/Cygwin
[command-line interface](https://en.wikipedia.org/wiki/Command-line_interface).

To compile the code it is necessary to create a common configuration file
called **config.mk** in the parent directory where the **pgmg** repository is
stored. Several templates are provided in the **config** directory. To use,
copy one of the templates into the parent directory. From the **pgmg**
directory, on a Linux computer, type
```Shell
cp config/config.mk.linux ../config.mk
```
On a Mac using GCC 10 installed via [MacPorts](http://www.macports.org), type
```Shell
cp config/config.mk.mac_mp ../config.mk
```
On a Mac using GCC installed via [Homebrew](http://brew.sh), type
```Shell
cp config/config.mk.mac_hb ../config.mk
```
On a Windows computer with Cygwin installed, type
```Shell
cp config/config.mk.win_cw ../config.mk
```
After this, the code examples can be compiled by typing
```Shell
make
```
## Contact
For questions about the code, contact [Chris Rycroft](http://seas.harvard.edu/~chr/).

## Acknowledgments
This work has been partially supported by the National Science Foundation under
Grant Nos. DMR-1409560 and DMS-1753203, and by the Applied Mathematics Program
of the U.S. DOE Office of Science Advanced Scientific Computing Research under
contract number DE-AC02-05CH11231.

## Bibliography
1. James W. Demmel, *Applied Numerical Linear Algebra*, SIAM (1997).
   [doi:10.1137/1.9781611971446](https://doi.org/10.1137/1.9781611971446)

2. William L. Briggs, Van Emden Henson, and Steve F. McCormick, *A Multigrid
   Tutorial, Second Edition*, SIAM (2000).
   [doi:10.1137/1.9780898719505](https://doi.org/10.1137/1.9780898719505)

3. Nicholas M. Boffi and Chris H. Rycroft, *Parallel three-dimensional
   simulations of quasi-static elastoplastic solids*, Comput. Phys. Commun.
   **257**, 107254 (2020). [doi:10.1016/j.cpc.2020.107254](https://doi.org/10.1016/j.cpc.2020.107254)

4. Nicholas M. Boffi and Chris H. Rycroft, *Coordinate transformation
   methodology for simulating quasistatic elastoplastic solids*, Phys. Rev. E
   **101**, 053304 (2020). [doi:10.1103/PhysRevE.101.053304](https://doi.org/10.1103/PhysRevE.101.053304)
