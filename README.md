Solitons
========

This directory contains codes to solve the Korteveg De Vries equation using
finite difference methods. There are two independent codes: a C code and an
ipython notebook. Both solve the problem in the same way (central differences
for the spatial derivatives and a 4th order Runge Kutta for the time 
derivative).

Initial conditions are embedded in the code but are highlighted by comments and
should be easy to change.

C Code
------
The C code is all inside a single file, it uses malloc and math.h but otherwise
should be pretty portable - just remember to add "-lm" to link in the math
libraries. To compile it try...

gcc -o kdv kdv.c -lm

then to run it try

./kdv > kdv.dat

The resulting .dat file can be used to create an animated gif of your solution
using the kdv.gnu gnuplot script. As long as gnuplot is installed, just running
this script (./kdv.gnu) should output a file called kdv.gif with your solution
animated - any web browser should be able to display it.


iPython notebook
----------------
The ipython notebook should run on most versions of ipython but needs the numpy
and matplotlib libraries installed. Like the C code, it will try to animate the
solution it produces. The animation does not work with the plot embedded into
the notebook at the moment but should open in a separate window. To run it
change into the directory containing Solitons.ipynb and start ipython notebook

ipython notebook

Select the Solitons notebook and run all of the cells.
