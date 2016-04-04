//---------------------------------------------------------------------------
// �Nikolai V. Shokhirev, 2004-2008  <nikolai@shokhirev.com>
// http://www.shokhirev.com/nikolai.html
// Reduced demo version
//---------------------------------------------------------------------------

#include "utility.h"
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

//---------------------------------------------------------------------------

bool Equal(real x1, real x2, real eps)
{
    return (fabs(x1 - x2) < eps);
}
bool Equal(uint x1, uint x2, real eps)
{
    return x1 == x2;
}
// FORTRAN abs
real abs(real x)
{
    return ((x >= 0.0) ? x : -x);
}

// FORTRAN iabs
int iabs(int x)
{
    return ((x >= 0.0) ? x : -x);
}
// more functions ...

bool Zero(real x, real eps)
{
    return (fabs(x) < eps);
}

bool CleanFile(const string& FileName)
{
    ofstream ofs(FileName, std::ofstream::out | std::ofstream::trunc);
    if (ofs.is_open()) {
        ofs.close();
        return true;
    }
    else {
        ofs.close();
        return false;
    }
}

bool DoesFileExist(const string& FileName)
{
    ofstream ofs(FileName, std::ofstream::in);
    if (ofs.is_open()) {
        ofs.close();
        return true;
    }
    else {
        ofs.close();
        return false;
    }
}

std::string ToString(const double x, const int width, const int decDigits)
{
    stringstream ss;
    ss << fixed << right;
    ss.fill(' '); // fill space around displayed #
    ss.width(width); // set  width around displayed #
    ss.precision(decDigits); // set # places after decimal
    ss << x;
    return ss.str();
}
std::string Center(const string s, const int w)
{
    stringstream ss, spaces;
    auto padding = w - s.size(); // count excess room to pad
    for (int i = 0; i < padding / 2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str(); // format with padding
    if (padding > 0 && padding % 2 != 0) // if odd #, add 1 space
        ss << " ";
    return ss.str();
}
string Right(const string s, const int w)
{
    stringstream ss, spaces;
    auto padding = w - s.size(); // count excess room to pad
    for (int i = 0; i < padding; ++i)
        spaces << " ";
    ss << spaces.str() << s; // format with padding
    return ss.str();
}

string Left(const string s, const int w)
{
    stringstream ss, spaces;
    auto padding = w - s.size(); // count excess room to pad
    for (int i = 0; i < padding; ++i)
        spaces << " ";
    ss << s << spaces.str(); // format with padding
    return ss.str();
}
