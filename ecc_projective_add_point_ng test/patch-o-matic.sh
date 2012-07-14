#!/bin/bash

# we are in "ecc_projective_add_point_ng test"
DIRNAME="ecc_projective_add_point_ng test"

# step up & create diff
cd ../
diff -Nur ./gnutls/lib ./lib > ./wmnaf-tmp.patch

# step into gnutls tree
cd ./gnutls/

# patch & clean
patch -Np1 -i "../wmnaf-tmp.patch"
rm "../wmnaf-tmp.patch"
