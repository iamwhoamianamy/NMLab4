#pragma once
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "Tests.h"
#include "Vector.h"

typedef double real;
using namespace std;

class Newton
{
public:
   Test2_2 test;                      // Информация о СНУ

   int max_iter_k;                    // Максимальное число итераций цикла k
   int max_iter_v;                    // Максимальное число итераций цикла v
   real eps0 = 1e-3;                  // Шаг для численного вычисления производной
   real eps1;                         // Условие малости для выхода из цикла v
   real eps2;                         // Условие малости для выхода из цикла k
   real norm_F0;                      // Норма вектор-функции F при x=x0

   vector<vector<real>> A;            // Матрица Якоби
   vector<real> Fk;                   // Значение вектор-функции F при x=xk
   vector<real> Fvk;                  // Значение вектор-функции F при x=xk
                                      // на v-той итерации
   vector<real> xk;                   // Вектор x на k-той итерации
   vector<real> xk1;                  // Вектор x на k+1-той итерации
   vector<real> dxk;                  // Вектор направления поиска решения
                                      // на k-той итерации

   // Решаем систему, описанную уравнениями, заданными соответствующими тестами
   Newton(Test2_2 _test)
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

   // Считываем необходимую для решения информацию из файла file_name
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

   // Получаем значение вектор-функции F при x=xk
   void get_Fk(vector<real>& xs, vector<real>& res)
   {
      res = test.Fk(xs);
   }

   // Получаем индексы наибольших элементов матрицы Якоби на удаление
   vector<int> best_indices(vector<real>& vec, const int& n)
   {
      vector<int> indices(vec.size());
      iota(indices.begin(), indices.end(), 0);

      partial_sort(indices.begin(), indices.begin() + n, indices.end(),
         [&vec](int i, int j) {return abs(vec[i]) < abs(vec[j]); });

      return vector<int>(indices.begin(), indices.begin() + n);
   }

   // Получаем индексы наименьших элементов матрицы Якоби на удаление
   vector<int> worst_indices(vector<real>& vec, const int& n)
   {
      vector<int> indices(vec.size());
      iota(indices.begin(), indices.end(), 0);

      partial_sort(indices.begin(), indices.begin() + n, indices.end(),
         [&vec](int i, int j) {return abs(vec[i]) > abs(vec[j]); });

      return vector<int>(indices.begin(), indices.begin() + n);
   }

   // Построение матрицы Якоби согласно 1 варианту
   void build_jacobi_1(vector<real>& xs)
   {
      A = test.Jacobi(xs);

      get_Fk(xs, Fk);

      for(int i = 0; i < test.n_func(); i++)
         A[i].push_back(-Fk[i]);

      vector<real> max_abs(test.n_var() + 1);

      for(int i = 0; i < test.n_func(); i++)
      {
         max_abs[i] = 0;
         for(int j = 0; j < test.n_var(); j++)
            if(abs(A[i][j]) > abs(max_abs[i]))
               max_abs[i] = abs(A[i][j]);
      }

      vector<int> indexes = worst_indices(max_abs, test.n_func() - test.n_var());

      sort(indexes.begin(), indexes.end());

      for(int i = indexes.size() - 1; i >= 0; i--)
         A.erase(A.begin() + indexes[i]);
   }

   // Построение матрицы Якоби согласно 2 варианту
   void build_jacobi_2(vector<real>& xs)
   {
      A = test.Jacobi(xs);
      
      get_Fk(xs, Fk);

      for(int i = 0; i < test.n_func(); i++)
         A[i].push_back(- Fk[i]);

      vector<int> indexes = best_indices(Fk, test.n_func() - test.n_var());

      sort(indexes.begin(), indexes.end());

      for(int i = indexes.size() - 1; i >= 0; i--)
         A.erase(A.begin() + indexes[i]);
   }

   // Построение матрицы Якоби согласно 6 варианту
   void build_jacobi_6(vector<real>& xs)
   {
      A.resize(test.n_func());
    
      for(int i = 0; i < test.n_func(); i++)
      {
         A[i].resize(test.n_var());
         for(int j = 0; j < test.n_var(); j++)
         {
            vector<real> eps(test.n_var());
            eps[j] = xs[j] + eps0;

            A[i][j] = test.Fk(eps)[i];

            eps[j] = xs[j] - eps0;

            A[i][j] -= test.Fk(eps)[i];
            A[i][j] /= 2 * eps0;
         }
      }

      get_Fk(xs, Fk);

      for(int i = 0; i < test.n_func(); i++)
         A[i].push_back(-Fk[i]);

      vector<int> indexes = best_indices(Fk, test.n_func() - test.n_var());

      sort(indexes.begin(), indexes.end());

      for(int i = indexes.size() - 1; i >= 0; i--)
         A.erase(A.begin() + indexes[i]);
   }

   // Прямой ход решателя методом Гаусса
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

   // Обратный ход решателя методом Гаусса
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

   // Функция решения СНУ
   void solve(int var)
   {
      for(int k = 0; k < max_iter_k && norm(Fk) / norm_F0 > eps2; k++)
      {
         switch(var)
         {
         case 1:
         {
            build_jacobi_1(xk);
            break;
         }
         case 2:
         {
            build_jacobi_2(xk);
            break;
         }
         case 6:
         {
            build_jacobi_6(xk);
            break;
         }
         }

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

         // Блок вывода информации о текущей итерации в консоль
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