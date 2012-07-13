#!/bin/bash

GNUTLS_LOCAL_LIB_PATH="./gnutls/lib/.libs/"
#GNUTLS_LOCAL_LIB_PATH="../../../external/gnutls-test/lib/.libs/"

LD_LIBRARY_PATH="$GNUTLS_LOCAL_LIB_PATH" ./test
