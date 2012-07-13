#!/bin/bash

# we are in "ecc_mulmod_wmnaf test"
DIRNAME="ecc_mulmod_wmnaf test"

# step up & create diff
cd ../
diff -Nur ./gnutls/lib ./lib > ./wmnaf-tmp.patch

# step into gnutls tree
cd ./gnutls/

# patch & clean
patch -Np1 -i "../wmnaf-tmp.patch"
rm "../wmnaf-tmp.patch"
