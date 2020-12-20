#include <iostream>
#include <fstream>
#include <vector>
#include "Vector.h"
#include "Newton.h"

using namespace std;
typedef double real;

int main()
{
   Test1_3 test;
   Newton newtonSolver = Newton(test);

   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve(1, "results/result_1.txt");

   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve(2, "results/result_2.txt");

   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve(6, "results/result_6.txt");
   
   cout << newtonSolver.xk[0] << " " << newtonSolver.xk[1] << endl;

   cout << "bl";
}