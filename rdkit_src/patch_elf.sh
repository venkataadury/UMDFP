#!/bin/bash
#
exe_file=$1
if test -z "$exe_file"
then
	echo "Please pick an executable to patch"
	exit 1
fi

mkdir -p libs
ldd $exe_file |grep -i rdkit| awk '{print $3}' | xargs -I {} cp {} ./libs/
ldd $exe_file |grep -i boost| awk '{print $3}' | xargs -I {} cp {} ./libs/
patchelf --set-rpath '$ORIGIN/libs' $exe_file
