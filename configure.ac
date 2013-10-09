#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

m4_define([VERSION_H],m4_esyscmd([cat version.h]))
changequote({{,}})dnl
m4_define({{XX_VERSION}},m4_bregexp(VERSION_H,{{"\(.+\)\\0";}},{{\1}}))
changequote([,])dnl

AC_INIT([libsoxrate], [XX_VERSION], [honeycomb77@gmail.com])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AM_PROG_CC_C_O
LT_LANG([Windows Resource])
AM_CONDITIONAL([HAVE_WINDRES], [test x$ac_cv_prog_ac_ct_RC != x])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
