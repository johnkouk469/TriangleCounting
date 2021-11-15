CC=gcc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3

default: all

sequential_masked_triangle_counting:
	$(CC) $(FLAGS) sequential_masked_triangle_counting.c mmio.c -o sequential_masked_triangle_counting

triangles_opencilk:
	$(CILKCC) $(FLAGS) triangles_opencilk.c mmio.c -o triangles_opencilk -fcilkplus

triangles_openmp:
	$(CC) $(FLAGS) triangles_openmp.c mmio.c -o triangles_openmp -fopenmp

all: sequential_masked_triangle_counting triangles_opencilk triangles_openmp

clean:
	rm -f sequential_masked_triangle_counting triangles_opencilk triangles_openmp