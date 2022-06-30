// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9

#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include "xcode_redirect.hpp"
#include <algorithm>
#include <getopt.h>
#include "TSP.h"
#include <cmath>
#include <algorithm>
using namespace std;

int main(int argc, char * argv[])
{
    ios_base::sync_with_stdio(false);
    xcode_redirect(argc, argv);
    //cerr << fixed << showpoint << setprecision(2); // for debuggging purpose delete later
    cout << std::setprecision(2);
    cout << std::fixed;
    TSP problem;
    problem.get_options(argc, argv);
    problem.read_file();
    problem.all_modes_implementation();
    return 0;
}
