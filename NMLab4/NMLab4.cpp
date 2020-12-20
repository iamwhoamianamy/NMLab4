#include <iostream>
#include <fstream>
#include <vector>
#include "Vector.h"
#include "Newton.h"

using namespace std;
typedef double real;

int main()
{
   Test2_2 test;
   Newton newtonSolver = Newton(test);

   cout << "k\tbeta\tx\ty\t\norm" << endl;
   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve(1);

   cout << "k\tbeta\tx\ty\t\norm" << endl;
   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve(2);

   cout << "k\tbeta\tx\ty\t\norm" << endl;
   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve(6);
   
   cout << newtonSolver.xk[0] << " " << newtonSolver.xk[1] << endl;

   cout << "bl";
}