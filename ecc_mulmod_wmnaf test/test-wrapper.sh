#!/bin/bash

#GNUTLS_LOCAL_LIB_PATH="./gnutls/lib/.libs/libgnutls.so"
GNUTLS_LOCAL_LIB_PATH="../../../external/gnutls-test/lib/.libs/libgnutls.so"

LD_LIBRARY_PATH="$GNUTLS_LOCAL_LIB_PATH" test
