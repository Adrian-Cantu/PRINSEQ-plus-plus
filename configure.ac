#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([PRINSEQ++], [1.0.1], [garbanyo@gmail.com] ,[prinseq-plus-plus], [https://github.com/Adrian-Cantu/PRINSEQ-plus-plus])
##AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cpp])
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11()

AC_CONFIG_HEADERS([config.h])



# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AX_BOOST_BASE([1.41], [], AC_MSG_ERROR([Could not find a useful version of boost]))
AX_BOOST_FILESYSTEM
AX_BOOST_IOSTREAMS
AX_CHECK_ZLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL


# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 README.md])
AC_OUTPUT
