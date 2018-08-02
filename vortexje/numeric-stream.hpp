//
// Vortexje -- Numeric stream supporting textual and binary streaming of numeric values
//
// Copyright (C) 2018 hrobeers
//
// Authors: Hans Robeers (hrobeers)
//   This file is distributed under the MIT License.
//   See notice at the end of this file.
//

#ifndef __NUMERIC_STREAM_HPP__
#define __NUMERIC_STREAM_HPP__

#include <iostream>

namespace nstream
{
  enum serial_mode_e {
    ASCII,
    BINARY
  };
  enum float_type_e {
    FLOAT,
    DOUBLE
  };
  enum byte_order_e {
    BIGENDIAN,
    LITTLEENDIAN
  };

  byte_order_e get_byte_order()
  {
    union {
      uint32_t i;
      char c[4];
    } bint = {0x01020304};
    return bint.c[0] == 1? BIGENDIAN : LITTLEENDIAN;
  }

  template <typename T>
  T bswap(T val) {
    T retVal;
    char *pVal = (char*)&val;
    char *pRetVal = (char*)&retVal;
    int size = sizeof(T);
    for(int i=0; i<size; i++)
      pRetVal[size-1-i] = pVal[i];

    return retVal;
  }

  template <typename T>
  std::ostream& write_bin(std::ostream &stream, T val, byte_order_e byte_order)
  {
    T ordered = (byte_order==get_byte_order())? val : bswap<T>(val);
    stream.write((char*)&ordered, sizeof(T));
    return stream;
  }

  // output numeric stream
  struct onstream
  {
    std::ostream &s;
    serial_mode_e mode;
    float_type_e float_type;
    byte_order_e byte_order;
    onstream(std::ostream &s, serial_mode_e m, float_type_e ft = DOUBLE, byte_order_e bo = BIGENDIAN)
      : s(s), mode(m), float_type(ft), byte_order(bo) {}
  };


  template <typename T>
  onstream& write_float(onstream &ns, T val) {
    if (ns.mode != BINARY) {
      ns.s << val;
      return ns;
    }

    switch (ns.float_type) {
      case FLOAT:
        write_bin<float>(ns.s, val, ns.byte_order);
        break;
      case DOUBLE:
        write_bin<double>(ns.s, val, ns.byte_order);
        break;
    }
    return ns;
  }
  inline onstream& operator<<(onstream &ns, float f) {
    return write_float(ns, f);
  }
  inline onstream& operator<<(onstream &ns, double d) {
    return write_float(ns, d);
  }
  inline onstream& operator<<(onstream &ns, char c) {
    // ignore char in binary mode (typically whitespace)
    if (ns.mode != BINARY)
      ns.s << c;
    return ns;
  }
  inline onstream& operator<<(onstream &ns, const char *cstr) {
    // ignore cstr in binary mode (typically whitespace)
    if (ns.mode != BINARY)
      ns.s << cstr;
    return ns;
  }
  inline onstream& operator<<(onstream &ns, std::ostream&(*f)(std::ostream&)) {
    // ignore manipulators like std::endl in binary mode
    if(ns.mode != BINARY)
      ns.s << f;
    return ns;
  }
};

#endif // __NUMERIC_STREAM_HPP__

/* ----------------------------------------------------------------------
 * Copyright (C) 2018 Hans Robeers. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */
