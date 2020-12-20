#pragma once

using namespace std;
typedef double real;

real sq(real val)
{
   return val * val;
}

// Две окружности не пересекаются
class Test1_1
{
public:
   virtual int n_func() { return 2; };
   virtual int n_var() { return 2; };

   Test1_1()
   {

   }

   real F1(vector<real>& xs) { return sq(xs[0] - 2) + sq(xs[1] - 2) - 4; }

   real F2(vector<real>& xs) { return sq(xs[0] + 2) + sq(xs[1] + 2) - 4; }

   real F1dx(vector<real>& xs) { return 2 * xs[0] - 4; }

   real F1dy(vector<real>& xs) { return 2 * xs[1] - 4; }

   real F2dx(vector<real>& xs) { return 2 * xs[0] + 4; }

   real F2dy(vector<real>& xs) { return 2 * xs[1] + 4; }

   vector<vector<real>> Jacobi(vector<real>& xs) { return { { F1dx(xs), F1dy(xs) }, { F2dx(xs), F2dy(xs) } }; }

   vector<real> Fk(vector<real>& xs) { return { F1(xs), F2(xs) }; }
};

// Две окружности пересекаются в 1 точке
class Test1_2
{
public:
   virtual int n_func() { return 2; };
   virtual int n_var() { return 2; };

   Test1_2()
   {

   }

   real F1(vector<real>& xs) { return sq(xs[0] - 2) + sq(xs[1] - 4) - 4; }

   real F2(vector<real>& xs) { return sq(xs[0] + 2) + sq(xs[1] - 4) - 4; }

   real F1dx(vector<real>& xs) { return 2 * xs[0] - 4; }

   real F1dy(vector<real>& xs) { return 2 * xs[1] - 8; }

   real F2dx(vector<real>& xs) { return 2 * xs[0] + 4; }

   real F2dy(vector<real>& xs) { return 2 * xs[1] - 8; }

   vector<vector<real>> Jacobi(vector<real>& xs) { return { { F1dx(xs), F1dy(xs) }, { F2dx(xs), F2dy(xs) } }; }

   vector<real> Fk(vector<real>& xs) { return { F1(xs), F2(xs) }; }
};

// Две окружности пересекаются в 2 точках
class Test1_3
{
public:
   virtual int n_func() { return 2; };
   virtual int n_var() { return 2; };

   Test1_3()
   {

   }

   real F1(vector<real>& xs) { return sq(xs[0] - 2) + sq(xs[1] - 4) - 9; }

   real F2(vector<real>& xs) { return sq(xs[0] + 2) + sq(xs[1] - 4) - 9; }

   real F1dx(vector<real>& xs) { return 2 * xs[0] - 4; }

   real F1dy(vector<real>& xs) { return 2 * xs[1] - 8; }

   real F2dx(vector<real>& xs) { return 2 * xs[0] + 4; }

   real F2dy(vector<real>& xs) { return 2 * xs[1] - 8; }

   vector<vector<real>> Jacobi(vector<real>& xs) { return { { F1dx(xs), F1dy(xs) }, { F2dx(xs), F2dy(xs) } }; }

   vector<real> Fk(vector<real>& xs) { return { F1(xs), F2(xs) }; }
};

// Две окружности и прямая не пересекаются
class Test2_1
{
public:
   virtual int n_func() { return 3; };
   virtual int n_var() { return 2; };

   Test2_1()
   {

   }

   real F1(vector<real>& xs) { return sq(xs[0] - 2) + sq(xs[1] - 2) - 4; }

   real F2(vector<real>& xs) { return sq(xs[0] + 2) + sq(xs[1] + 2) - 4; }

   real F3(vector<real>& xs) { return xs[0] + xs[1]; }

   real F1dx(vector<real>& xs) { return 2 * xs[0] - 4; }

   real F1dy(vector<real>& xs) { return 2 * xs[1] - 4; }

   real F2dx(vector<real>& xs) { return 2 * xs[0] + 4; }

   real F2dy(vector<real>& xs) { return 2 * xs[1] + 4; }

   real F3dx(vector<real>& xs) { return 1; }

   real F3dy(vector<real>& xs) { return 1; }

   vector<vector<real>> Jacobi(vector<real>& xs) { return { { F1dx(xs), F1dy(xs) }, { F2dx(xs), F2dy(xs) }, { F3dx(xs), F3dy(xs) } }; }

   vector<real> Fk(vector<real>& xs) { return { F1(xs), F2(xs), F3(xs) }; }
};

// Две окружности и прямая пересекаются
class Test2_2
{
public:
   virtual int n_func() { return 3; };
   virtual int n_var() { return 2; };

   Test2_2()
   {

   }

   real F1(vector<real>& xs) { return sq(xs[0] - 2) + sq(xs[1] - 2) - 4; }

   real F2(vector<real>& xs) { return sq(xs[0] + 2) + sq(xs[1] - 2) - 4; }

   real F3(vector<real>& xs) { return xs[0]/* + xs[1] - 2*/; }

   real F1dx(vector<real>& xs) { return 2 * xs[0] - 4; }

   real F1dy(vector<real>& xs) { return 2 * xs[1] - 4; }

   real F2dx(vector<real>& xs) { return 2 * xs[0] + 4; }

   real F2dy(vector<real>& xs) { return 2 * xs[1] - 4; }

   real F3dx(vector<real>& xs) { return 1; }

   real F3dy(vector<real>& xs) { return 0; }

   vector<vector<real>> Jacobi(vector<real>& xs) { return { { F1dx(xs), F1dy(xs) }, { F2dx(xs), F2dy(xs) }, { F3dx(xs), F3dy(xs) } }; }

   vector<real> Fk(vector<real>& xs) { return { F1(xs), F2(xs), F3(xs) }; }
};
