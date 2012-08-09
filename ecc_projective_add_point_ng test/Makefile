CFLAGS= -O3 -march=native -mtune=native -mfpmath=sse -fomit-frame-pointer -funroll-loops -pipe -lgmp -I.
#CFLAGS= -O0 -march=native -mtune=native -mfpmath=sse -fomit-frame-pointer -funroll-loops -pipe -lgmp -Wall -I.

# GNUTLS_LOCAL_LIB_DIR must be set for this to work
# See README
GNUTLS_LOCAL_LIB_PATH= ${GNUTLS_LOCAL_LIB_DIR}/libgnutls.so

CC= gcc $(CFLAGS)

all: run

run: test
	LD_LIBRARY_PATH="${GNUTLS_LOCAL_LIB_DIR}" ./test

test: ecc_point.o test.o
	$(CC) -o test $(GNUTLS_LOCAL_LIB_PATH) ecc_point.o test.o

ecc_point.o: ecc_point.c ecc_point.h
	$(CC) -c ecc_point.c

test.o: test.c ecc.h ecc_point.h
	$(CC) -c test.c

clean:
	rm -f *.o test *\~
