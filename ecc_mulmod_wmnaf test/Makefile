CFLAGS= -O3 -march=native -mtune=native -mfpmath=sse -fomit-frame-pointer -funroll-loops -pipe -lgmp -I.
#CFLAGS= -O0 -march=native -mtune=native -mfpmath=sse -fomit-frame-pointer -funroll-loops -pipe -lgmp -Wall -I.
#GNUTLS_LOCAL_LIB_PATH= ./gnutls/lib/.libs/libgnutls.so
GNUTLS_LOCAL_LIB_PATH= ../../../external/gnutls-test/lib/.libs/libgnutls.so
CC= gcc $(CFLAGS)

all: test 

test: my_ecc.o test.o 
	$(CC) -o test $(GNUTLS_LOCAL_LIB_PATH) my_ecc.o test.o

my_ecc.o: my_ecc.c my_ecc.h
	$(CC) -c my_ecc.c

test.o: test.c my_ecc.h
	$(CC) -c test.c

clean:
	rm -f *.o test *\~
