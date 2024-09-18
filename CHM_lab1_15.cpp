
#include <iostream>
#include <vector>
#include <cmath> // Для функции sqrt
#include <iomanip> // Для форматирования вывода
#include <fstream> // Для работы с файлами;
using namespace std;

// Функция для вывода матрицы
void printMatrix(const vector<vector<double>>& matrix) 
{
   int n = matrix.size();
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++) 
         cout << setw(10) << matrix[i][j] << " ";
      cout << endl;
   }
}

// Функция для вывода вектора
void printVector(const vector<double>& vec) 
{
   for (double v : vec) 
      cout << setw(10) << v << " ";
   cout << endl;
}

// Функция для выполнения LU(sq)-разложения
bool LU_SQ_Decomposition(const vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) 
{
   int n = A.size();
   L = vector<vector<double>>(n, vector<double>(n, 0));
   U = vector<vector<double>>(n, vector<double>(n, 0));

   for (int i = 0; i < n; i++) 
   {
      for (int j = 0; j <= i; j++) 
      {
         double sum = 0;

         // Вычисление элементов нижнетреугольной матрицы L
         if (i == j) {  // Диагональный элемент
            for (int k = 0; k < j; k++) 
            {
               sum += L[j][k] * L[j][k];
            }
            double value = A[j][j] - sum;
            if (value <= 0) 
            {
               cout << "Матрица не является положительно определенной." << endl;
               return false;
            }
            L[j][j] = sqrt(value);
            U[j][j] = L[j][j]; // L и U имеют одинаковые диагональные элементы
         }
         else 
         {  // Вне диагонали
            for (int k = 0; k < j; k++) 
            {
               sum += L[i][k] * L[j][k];
            }
            L[i][j] = (A[i][j] - sum) / L[j][j];
            U[j][i] = L[i][j]; // U является транспонированной копией L
         }
      }
   }

   return true;
}

// Прямой ход: решение системы Ly = F
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& F) {
   int n = L.size();
   vector<double> y(n);

   for (int i = 0; i < n; i++) 
   {
      double sum = 0;
      for (int j = 0; j < i; j++) 
      {
         sum += L[i][j] * y[j];
      }
      y[i] = (F[i] - sum) / L[i][i];
   }

   return y;
}

// Обратный ход: решение системы Ux = y
vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y) 
{
   int n = U.size();
   vector<double> x(n);

   for (int i = n - 1; i >= 0; i--) 
   {
      double sum = 0;
      for (int j = i + 1; j < n; j++) 
      {
         sum += U[i][j] * x[j];
      }
      x[i] = (y[i] - sum) / U[i][i];
   }

   return x;
}

// Функция для чтения матрицы A и вектора F из файла
bool readMatrixAndVector(const string& matrixFile, const string& vectorFile, vector<vector<double>>& A, vector<double>& F) {
   ifstream matrixStream(matrixFile);
   ifstream vectorStream(vectorFile);

   if (!matrixStream.is_open() || !vectorStream.is_open()) 
   {
      cerr << "Не удалось открыть файл." << endl;
      return false;
   }

   int n;
   matrixStream >> n; // Считываем размер матрицы
   A = vector<vector<double>>(n, vector<double>(n, 0));
   F = vector<double>(n, 0);

   // Считывание матрицы A
   for (int i = 0; i < n; i++) 
   {
      for (int j = 0; j < n; j++) 
      {
         matrixStream >> A[i][j];
      }
   }

   // Считывание вектора F
   for (int i = 0; i < n; i++) 
   {
      vectorStream >> F[i];
   }

   matrixStream.close();
   vectorStream.close();
   return true;
}

int main() {
   setlocale(LC_ALL, "Russian");

   vector<vector<double>> A, L, U;
   vector<double> F, y, x;

   // Считывание матрицы и вектора из файлов
   if (!readMatrixAndVector("matrix.txt", "vector.txt", A, F)) 
   {
      return 1;
   }

   // 1. LU(sq)-разложение
   if (LU_SQ_Decomposition(A, L, U)) 
   {
      cout << "Матрица L:" << endl;
      printMatrix(L);
      cout << endl << "Матрица U:" << endl;
      printMatrix(U);

      // 2. Прямой ход: решение системы Ly = F
      y = forwardSubstitution(L, F);
      cout << endl << "Вектор y (решение системы Ly = F):" << endl;
      printVector(y);

      // 3. Обратный ход: решение системы Ux = y
      x = backwardSubstitution(U, y);
      cout << endl << "Вектор x (решение системы Ux = y):" << endl;
      printVector(x);
   }
   else {
      cout << "Разложение невозможно." << endl;
   }

   return 0;
}
