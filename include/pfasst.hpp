#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include "pfasst/interfaces.hpp"

#include <string>
using namespace std;

struct OUT
{
  public:
    static const string black;
    static const string red;
    static const string green;
    static const string yellow;
    static const string blue;
    static const string magenta;
    static const string cyan;
    static const string white;

    static const string bold;
    static const string underline;

    static const string reset;
};

const string OUT::black = "\033[30m";
const string OUT::red = "\033[31m";
const string OUT::green = "\033[32m";
const string OUT::yellow = "\033[33m";
const string OUT::blue = "\033[34m";
const string OUT::magenta = "\033[35m";
const string OUT::cyan = "\033[36m";
const string OUT::white = "\033[37m";
const string OUT::bold = "\033[1m";
const string OUT::underline = "\033[4m";
const string OUT::reset = "\033[0m";

#endif
