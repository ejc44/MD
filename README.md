# MD
Molcular dyanmics implementation for CS 6.338 Final Project

This repo includes a serial and several parallel versions of a molecular dynamics simulation in the Julia programming language contained in Julia notebooks. Notable features include

1) The option to read parameters (like number of simulation steps) and/or initial data (starting configuration: position, velocity, acceleration, forces) from external files. These files must be in the same folder as this notebook and have the correct names (as specified in the code).

2) The option to directly specify the parameters in the notebook. Note: These parameters are all constants, so one must restart the kernel to redefine them.

3) The option to save the parameters and output data to external files.

4) The option to model finite and infinite systems.

5) The ability to make finite systems periodic or non-periodic.

6) If the initial data / configuration is not specified in a file, it can be generated in the code. For example, the starting positions are random within the specified box size. The user can modify any of the initial conditions by altering the initialize() function.

7) The user can also modify the form and strength of the forces by altering the find_forces() and gen_interaction() functions.

8) For 2 or 3 dimensional simulations, the system can be visually displayed in a plot. For non-Windows systems, the user can also use the Interact package to manipulate the plot to see the movement of the particles in the system over time. One of the parameters sets the frequency with which the program saves the particles' positions.

Note: The parameters are all constants, so one must restart the kernel to redefine them.

Note: If the parameters and/or initial data is read from files, the program assumes the files are compatible (i.e. the dimension in the dim.txt file matches the dimensions of the particles' positions, etc.). This is not a problem if the user is restarting from previously saved files but could be a problem if the user has made changes to the files.

Specific details regarding each version are included in each notebook, and separate timing notebooks are provided that use the BenchmarkTools.jl package but do not contain the plotting functionality.  There are also three Testing folders containing the parameter and data files used for the timing tests.  The last folder, called Incomplete, contains a not fully implemented GPU version of the code.  My raw timing data, presentation, and final report are all also included.

References:

For Julia help: The documentation at https://julialang.org/

For calculating forces: Article titled "Verlet integration" at https://www.saylor.org/site/wp-content/uploads/2011/06/MA221-6.1.pdf

For general layout of basic MD code: John Burkardt's website at https://people.sc.fsu.edu/~jburkardt/py_src/md/md.html
