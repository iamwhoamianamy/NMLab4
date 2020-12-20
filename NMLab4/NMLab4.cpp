#include <iostream>
#include <fstream>
#include <vector>
#include "Vector.h"
#include "Newton.h"

using namespace std;
typedef double real;

int main()
{
   //vector<vector<real>> A = { {1,2,5},{0,3,6} };
   //vector<real> res(2);

   //Newton ns;

   //ns.forward_gauss(A);
   //ns.backward_gauss(A, res);

   Test2_1 test;
   Newton newtonSolver = Newton(test);
   newtonSolver.read_info("tests/info1.txt");
   newtonSolver.solve();
   
   cout << newtonSolver.xk[0] << " " << newtonSolver.xk[1] << endl;

   cout << "bl";
}