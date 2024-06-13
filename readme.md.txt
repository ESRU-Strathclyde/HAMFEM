Install dependencies:
	sudo apt-get install libblas-dev liblapack-dev


use make to compile the application. Inside the folder with the file "Makefile"
	make
after compilation, move all files and subfolders to the instalation folder (this is hardcoded on Hamfem.f90). the lenght of the path is also hardcoded (character(12)::installpath)
	/opt/hamfem/
create a symbolic link to the executable in a folder already on the PATH
	cd /usr/local/bin; sudo rm hamfem; sudo ln -s /opt/hamfem/Release/hamfem hamfem

to run sample cases, copy the folder tester/to_run
run an example, such as:
	hamfem  00basevers.ham



compile the tester:
f95 -O3 -v  Hamfem-tester.f90  -o hamfem-tester
move the executable to the to_run folder, and execute it rhere 



MSYS requirements (still to be completed, not working)
pacman -S mingw64/mingw-w64-x86_64-openblas
pacman -S mingw64/mingw-w64-x86_64-lapack
pacman -S mingw64/mingw-w64-x86_64-arpack