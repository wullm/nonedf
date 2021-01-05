#Compiler options
GCC = mpicc

#Libraries
INI_PARSER = parser/minIni.o
STD_LIBRARIES = -lm
FFTW_LIBRARIES = -lfftw3 -lfftw3_omp -lfftw3_mpi
HDF5_LIBRARIES = -lhdf5
GSL_LIBRARIES = -lgsl -lgslcblas

GSL_INCLUDES =

HDF5_INCLUDES += -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
HDF5_LIBRARIES += -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -I/usr/include/hdf5/openmpi

#Putting it together
INCLUDES = $(HDF5_INCLUDES) $(GSL_INCLUDES) $(FIREBOLT_INCLUDES)
LIBRARIES = $(INI_PARSER) $(STD_LIBRARIES) $(FFTW_LIBRARIES) $(HDF5_LIBRARIES) $(GSL_LIBRARIES) $(FIREBOLT_LIBRARIES)
CFLAGS = -Wall -Wshadow=global -fopenmp -march=native -O4
LDFLAGS =

OBJECTS = lib/*.o

all:
	make minIni
	mkdir -p lib
	$(GCC) src/input.c -c -o lib/input.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/random.c -c -o lib/random.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/cosmology.c -c -o lib/cosmology.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/fft.c -c -o lib/fft.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/mesh_grav.c -c -o lib/mesh_grav.o $(INCLUDES) $(CFLAGS)
	$(GCC) src/phase_space.c -c -o lib/phase_space.o $(INCLUDES) $(CFLAGS)

	$(GCC) src/nonedf.c -o nonedf $(INCLUDES) $(OBJECTS) $(LIBRARIES) $(CFLAGS) $(LDFLAGS)

minIni:
	cd parser && make
