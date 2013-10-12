all : build run

run : ca-stripped
	time ./ca-stripped 01110110 20000000 20000 0

build : ca ca.lst

ca-stripped : ca
	strip -s ca -o ca-stripped

ca : ca.c Makefile
	gcc -Ofast -mtune=native -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mpopcnt -Wall -Wextra ca.c -o ca -lpthread

ca.lst : ca
	objdump -d ca > ca.lst

clean :
	rm -rf ca