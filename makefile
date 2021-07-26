#
# Makefile for snarkhunter
#

SHELL = /bin/sh

# Compiling executing program with DWORDSIZE=32 is slightly faster, 
# but limits the order of the graphs to 32.
# -mpopcnt is necessary for __builtin_popcount
CC32 = gcc -DWORDSIZE=32 -DMAXN=WORDSIZE -mpopcnt
CC64 = gcc -DWORDSIZE=64 -DMAXN=WORDSIZE -mpopcnt
CFLAGS = -O4

all :
	rm -rf main 
	${CC32} $(CFLAGS) main.c nautyW1.a -o CubicGraphGen

64bit :
	rm -rf main-64	
	${CC64} $(CFLAGS) main.c nautyL1.a -o CubicGraphGen-64

profile :
	rm -rf main-profile
	${CC32} -pg -g main.c nautyW1.a -o CubicGraphGen-profile

