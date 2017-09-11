#!/bin/sh
aclocal -I m4
autoreconf --install
automake --add-missing --copy >/dev/null 2>&1
