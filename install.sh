#!/bin/bash

if [ "x$1" != "x" ]; then
    prefix=$1
else
    prefix=$HOME
fi
bin_path=${prefix}/bin
isis_path=${prefix}/lib/isis
python_path=${prefix}/lib/python

echo "bin: $bin_path"
echo "isis: $isis_path"
echo "python: $python_path"

install isis/*.sl $isis_path
install isis/klines/* $isis_path/klines/
install isis/xrsdemo/* $isis_path/xrsdemo/
install python/*.py $python_path
install bin/cxrs ${prefix}/bin
chmod a+x ${prefix}/bin/cxrs

echo "done"
