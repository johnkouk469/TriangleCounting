CCILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/c;ang
CFLAGS=-O3

default: all

v3_opencilk:
	$(CILKCC) $(FLAGS) v3_opencilk_1111.c mmio.c -o v3_opencilk_1111 -fcilkplus
