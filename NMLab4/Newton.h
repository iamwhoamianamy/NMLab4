#pragma once
#include "Tests.h"

typedef double real;
using namespace std;

class Newton
{
public:
   Test2_1 test;

   int max_iter_k;
   int max_iter_v;
   real eps1, eps2;
   real norm_F0;

   vector<vector<real>> A;
   vector<real> Fk;
   vector<real> Fvk;
   vector<real> xk;
   vector<real> xk1;
   vector<real> dxk;

   Newton(Test2_1 _test)
   {
      test = _test;
      A.resize(test.n_func());

      for (int i = 0; i < test.n_func(); i++)
         A[i].resize(test.n_var());

      Fk.resize(test.n_var());
      Fvk.resize(test.n_var());

      xk.resize(test.n_var());
      xk1.resize(test.n_var());
      dxk.resize(test.n_var());
   }

   Newton() {};

   void read_info(string file_name)
   {
      ifstream fin;
      fin.open(file_name);

      string fict;

      fin >> fict;
      fin >> eps1;

      fin >> fict;
      fin >> eps2;

      fin >> fict;
      fin >> max_iter_k;

      fin >> fict;
      fin >> max_iter_v;

      fin >> fict;
      for(int i = 0; i < test.n_var(); i++)
         fin >> xk[i];

      get_Fk(xk, Fk);
      norm_F0 = norm(Fk);

      fin.close();
   }

   void get_Fk(vector<real>& xs, vector<real>& res)
   {
      res = test.Fk(xs);
   }

   void find_mins(vector<real>& vec, vector<real>& indices, int n)
   {
      int added = 0;
      for(int i = 0; i < test.n_func(); i++)
      {

      }

   }

   void build_jacobi_2(vector<real>& xs)
   {
      A = test.Jacobi(xs);
      
      get_Fk(xs, Fk);

      for(int i = 0; i < test.n_func(); i++)
         A[i].push_back(- Fk[i]);

      int min_index = 0;

      for(int i = 0; i < test.n_func(); i++)
         if(abs(A[i][test.n_var() - 1]) < abs(A[min_index][test.n_var() - 1]))
            min_index = i;

      A.erase(A.begin() + min_index);
   }

   void forward_gauss(vector<vector<real>>& mat)
   {
      int n = mat.size() + 1;
      int rowWithMaxElem = 0;
      vector<real> rowAdress(n);
      real maxElem = 0;

      for (int j = 0; j < n - 2; j++)
      {
         int rowNumber = j;
         for (int currentCol = j; currentCol < n - 1; currentCol++)
         {
            if (fabs(mat[currentCol][j]) > maxElem)
            {
               maxElem = fabs(mat[currentCol][j]);
               rowWithMaxElem = currentCol;
            }
         }

         maxElem = 0;
         rowAdress = mat[rowWithMaxElem];
         mat[rowWithMaxElem] = mat[rowNumber];
         mat[rowNumber] = rowAdress;

         for (int k = 1 + j; k < n - 1; k++)
         {
            real factor = mat[k][rowNumber] / mat[rowNumber][rowNumber];
            if (factor != 0)
            {
               for (int i = rowNumber; i < n; i++)
               {
                  real tmp = mat[rowNumber][i] * factor;
                  mat[k][i] -= tmp;
               }
            }
         }
      }
   }

   void backward_gauss(vector<vector<real>>& mat, vector<real>& res)
   {
      int n = mat.size() + 1;
      res.resize(n - 1);

      for (int i = n - 2; i >= 0; i--)
      {
         real sum = 0;
         for (int j = i + 1; j < n - 1; j++)
         {
            sum += res[j] * mat[i][j];
         }

         res[i] = (mat[i][n - 1] - sum) / mat[i][i];
      }
   }

   void solve()
   {
      for(int k = 0; k < max_iter_k && norm(Fk) / norm_F0 > eps2; k++)
      {
         build_jacobi_2(xk);
         forward_gauss(A);
         backward_gauss(A, dxk);

         for(int i = 0; i < test.n_var(); i++)
            if(abs(dxk[i]) == INFINITY)
            {
               cout << "Cant solve!" << endl;
               return;
            }

         real beta = 1;

         for(int v = 0; v < max_iter_v; v++)
         {
            xk1 = xk + beta * dxk;
            get_Fk(xk1, Fvk);

            if(norm(Fvk) < norm(Fk) || beta < eps1)
               break;
            else
               beta /= 2;
         }
         xk = xk1;

         cout << "k:       " << setw(5) << k << endl;
         cout << "beta:    " << setw(5) << beta << endl;
         cout << "xk:      " << setw(5);
         for(int i = 0; i < test.n_var(); i++)
            cout << xk[i] << " ";
         cout << endl;
         cout << "norm Fk: " << setw(5) << norm(Fk) << endl << endl;
      }
   }
};