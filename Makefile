
FORTRAN  = f95
GCC      = gcc
OPTS     = -O3 -v



all: socketlib.o HAM_thermodynamics_library.o HAM_material_library.o HAM_boundary_conditions.o FEM_shape_functions.o Hamfem.o
	${FORTRAN}  ${OPTS} -lnsl -o ./Release/hamfem socketlib.o HAM_thermodynamics_library.o HAM_material_library.o HAM_boundary_conditions.o FEM_shape_functions.o Hamfem.o -L/usr/local/lib -llapack -lblas

socketlib.o: socketlib.c
	${GCC} -c socketlib.c		
HAM_thermodynamics_library.o: HAM_thermodynamics_library.f90
	${FORTRAN} -c ${OPTS} HAM_thermodynamics_library.f90
HAM_material_library.o: HAM_material_library.f90 HAM_thermodynamics_library.f90
	${FORTRAN} -c ${OPTS} HAM_material_library.f90
HAM_boundary_conditions.o: HAM_boundary_conditions.f90 HAM_thermodynamics_library.f90
	${FORTRAN} -c ${OPTS} HAM_boundary_conditions.f90
FEM_shape_functions.o: FEM_shape_functions.f90
	${FORTRAN} -c ${OPTS} FEM_shape_functions.f90
Hamfem.o: HAM_thermodynamics_library.f90 Hamfem.f90 Mainprog.f90
	${FORTRAN} -c ${OPTS} Hamfem.f90

clean:
	rm -f *.o
