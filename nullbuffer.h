#ifndef NULLBUFFER_H
#define NULLBUFFER_H
#include <streambuf>

class NullBuffer : public std::streambuf{
 public:
  int overflow(int c) { return c; }
};
#endif
