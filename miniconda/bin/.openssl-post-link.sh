#!/bin/bash

# break hard link
cp $PREFIX/lib/libcrypto.so.1.0.0 $PREFIX/lib/libcrypto.so.1.0.0-tmp
mv $PREFIX/lib/libcrypto.so.1.0.0-tmp $PREFIX/lib/libcrypto.so.1.0.0

$PREFIX/bin/.openssl-libcrypto-fix $PREFIX
