#pragma once

#include <cstdlib>
#include <stdio.h>
#include <iostream>

//! Prints the error message, the stack trace, and exits
#define FatalError(s) {                                             \
  printf("Fatal error '%s' at %s:%d\n",s,__FILE__,__LINE__);        \
  exit(1); }
