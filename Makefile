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

triangles_pthreads:
	$(CC) $(FLAGS) -pthread triangles_pthreads.c mmio.c -o triangles_pthreads

all: sequential_masked_triangle_counting triangles_opencilk triangles_openmp triangles_pthreads

test:
	@printf "\nPthreads: \n"
	./triangles_pthreads mtx/belgium_osm.mtx
	./triangles_pthreads mtx/com-Youtube.mtx
	./triangles_pthreads mtx/mycielskian13.mtx
	./triangles_pthreads mtx/dblp-2010.mtx
	./triangles_pthreads mtx/NACA0015.mtx
	@printf "\nOpenCilk: \n"
	./triangles_opencilk mtx/belgium_osm.mtx
	./triangles_opencilk mtx/com-Youtube.mtx
	./triangles_opencilk mtx/mycielskian13.mtx
	./triangles_opencilk mtx/dblp-2010.mtx
	./triangles_opencilk mtx/NACA0015.mtx
	@printf "\nOpenMP: \n"
	./triangles_openmp mtx/belgium_osm.mtx
	./triangles_openmp mtx/com-Youtube.mtx
	./triangles_openmp mtx/mycielskian13.mtx
	./triangles_openmp mtx/dblp-2010.mtx
	./triangles_openmp mtx/NACA0015.mtx
	@printf "\nSequential: \n"
	./sequential_masked_triangle_counting mtx/belgium_osm.mtx
	./sequential_masked_triangle_counting mtx/com-Youtube.mtx
	./sequential_masked_triangle_counting mtx/mycielskian13.mtx
	./sequential_masked_triangle_counting mtx/dblp-2010.mtx
	./sequential_masked_triangle_counting mtx/NACA0015.mtx

clean:
	rm -f sequential_masked_triangle_counting triangles_opencilk triangles_openmp triangles_pthreads