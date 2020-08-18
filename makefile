program_name = run/heat_eq 

IFORT_COMPILE = mpiifort -c -O3
IFORT_LINK = mpiifort -O3

OBJ = main.o partition_mod.o mpicom_mod.o halo_mod.o grid_function_mod.o data_mod.o

$(program_name): $(OBJ)
	$(IFORT_LINK) $(OBJ) -o $(program_name)

data_mod.o: data_mod.f90 grid_function_mod.o
	$(IFORT_COMPILE) data_mod.f90

grid_function_mod.o: grid_function_mod.f90
	$(IFORT_COMPILE) grid_function_mod.f90

halo_mod.o: halo_mod.f90 
	$(IFORT_COMPILE) halo_mod.f90

mpicom_mod.o: mpicom_mod.f90 halo_mod.o partition_mod.o data_mod.o
	$(IFORT_COMPILE) mpicom_mod.f90

partition_mod.o: partition_mod.f90
	$(IFORT_COMPILE) partition_mod.f90

main.o: main.f90 mpicom_mod.o partition_mod.o grid_function_mod.o data_mod.o
	$(IFORT_COMPILE) main.f90

clear: 
	rm *.o* *.e* *.mod
