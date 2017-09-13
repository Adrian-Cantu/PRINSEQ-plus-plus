#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pprinseqc], [0.1], [garbanyo@gmail.com])

AC_CONFIG_SRCDIR([src/pprinseqc.cpp])
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11()

AC_CONFIG_HEADERS([config.h])



# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL


# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
