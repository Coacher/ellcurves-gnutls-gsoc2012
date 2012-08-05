CFLAGS= -O3 -march=native -mtune=native -mfpmath=sse -fomit-frame-pointer -funroll-loops -pipe -lgmp -I.
#CFLAGS= -O0 -march=native -mtune=native -mfpmath=sse -fomit-frame-pointer -funroll-loops -pipe -lgmp -Wall -I.

# GNUTLS_LOCAL_LIB_DIR must be set for this to work
GNUTLS_LOCAL_LIB_PATH= ${GNUTLS_LOCAL_LIB_DIR}/libgnutls.so

CC= gcc $(CFLAGS)

all: test

run: test
	LD_LIBRARY_PATH="${GNUTLS_LOCAL_LIB_DIR}" ./test 

test: my_ecc.o test.o
	$(CC) -o test $(GNUTLS_LOCAL_LIB_PATH) my_ecc.o test.o

my_ecc.o: my_ecc.c my_ecc.h
	$(CC) -c my_ecc.c

test.o: test.c ecc.c my_ecc.h
	$(CC) -c test.c

clean:
	rm -f *.o test *\~
