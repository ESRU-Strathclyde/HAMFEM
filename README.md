# HAMFEM
Heat, Air, and Moisture 1D to 3D simulation for porous materials using the Finite Element Method

## Documentation
https://appdocs.esru.strath.ac.uk/

## Compile
Install dependencies:
	sudo apt-get install libblas-dev liblapack-dev

Use make to compile the application.

After compilation, move all files and subfolders to folder /opt/hamfem/
This location is hardcoded on Hamfem.f90. The lenght of the path is also hardcoded (character(12)::installpath).

Create a symbolic link to the executable in a folder already on the PATH:
	cd /usr/local/bin; sudo rm hamfem; sudo ln -s /opt/hamfem/Release/hamfem hamfem

## Execute
To run sample cases, copy the folder tester/to_run to the home folder.
Execute the command in the folder where the simulation configuration is located (extension .ham)
Example:
	hamfem  00basevers.ham

## Test
The tester compares current hamfem results against a set of precalculated results for selected cases.

Compile the tester:
  f95 -O3 -v  Hamfem-tester.f90  -o hamfem-tester
Copy the folder to_run to the home folder.
Move the tester executable to the to_run folder (in home), and execute it there.
