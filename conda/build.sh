#!/bin/bash
mkdir -p $PREFIX/bin
GOOS=linux GOARCH=amd64 go build -o $PREFIX/bin/
chmod a+x $PREFIX/bin/gopeaks
