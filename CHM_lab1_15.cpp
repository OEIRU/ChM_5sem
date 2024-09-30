#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>

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

// LU-разложение с диагональным преобладанием
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
            if (i == j) 
            { 
                for (int k = 0; k < j; k++) 
                    sum += L[j][k] * L[j][k];
                double value = A[j][j] - sum;
                if (value <= 0) 
                {
                    cout << "Матрица не является положительно определенной." << endl;
                    return false;
                }
                L[j][j] = sqrt(value);
                U[j][j] = L[j][j];
            } 
            else 
            {
                for (int k = 0; k < j; k++) 
                    sum += L[i][k] * L[j][k];
                L[i][j] = (A[i][j] - sum) / L[j][j];
                U[j][i] = L[i][j];
            }
        }
    }
    return true;
}

// Прямой ход
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& F) 
{
    int n = L.size();
    vector<double> y(n);
    for (int i = 0; i < n; i++) 
    {
        double sum = 0;
        for (int j = 0; j < i; j++) 
            sum += L[i][j] * y[j];
        y[i] = (F[i] - sum) / L[i][i];
    }
    return y;
}

// Обратный ход
vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y) 
{
    int n = U.size();
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) 
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++) 
            sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / U[i][i];
    }
    return x;
}

// Построение матрицы с регулируемым числом обусловленности
vector<vector<double>> buildMatrixWithConditioning(int n, double alpha) 
{
    vector<vector<double>> A(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            if (i == j)
                A[i][j] = alpha * (i + 1);
            else
                A[i][j] = 1.0 / (i + j + 1);
        }
    }
    return A;
}

// Вычисление ошибки
double computeError(const vector<double>& x, const vector<double>& x_exact) 
{
    double error = 0;
    for (int i = 0; i < x.size(); i++) 
        error += pow(x[i] - x_exact[i], 2);
    return sqrt(error);
}

// Построение матрицы Гильберта
vector<vector<double>> buildHilbertMatrix(int n) 
{
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            A[i][j] = 1.0 / (i + j + 1);
        }
    }
    return A;
}

// Метод Гаусса с выбором ведущего элемента
bool gaussianElimination(vector<vector<double>>& A, vector<double>& b) 
{
    int n = A.size();
    for (int i = 0; i < n; i++) 
    {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) 
        {
            if (abs(A[k][i]) > abs(A[maxRow][i])) 
                maxRow = k;
        }

        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; k++) 
        {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) 
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) 
    {
        x[i] = b[i] / A[i][i];
        for (int j = i - 1; j >= 0; j--) 
            b[j] -= A[j][i] * x[i];
    }
    return true;
}

// Подсчет операций для LU-разложения
int countOperationsLU(int n) 
{
    return (n * n * n) / 3;
}

// Подсчет операций для метода Гаусса
int countOperationsGauss(int n) 
{
    return (2 * n * n * n) / 3;
}

int main() 
{
    setlocale(LC_ALL, "Russian");

    // Размер матрицы и параметр для регулировки числа обусловленности
    int n = 5;
    double alpha = 1e3;

    vector<vector<double>> A_cond = buildMatrixWithConditioning(n, alpha);
    vector<vector<double>> L, U;
    vector<double> F(n, 1);  // Вектор правой части
    vector<double> x_exact(n, 1); // Точное решение

    // Решение с помощью LU-разложения
    cout << "Решение для матрицы с регулируемым числом обусловленности:" << endl;
    if (LU_SQ_Decomposition(A_cond, L, U)) 
    {
        vector<double> y = forwardSubstitution(L, F);
        vector<double> x = backwardSubstitution(U, y);
        double error = computeError(x, x_exact);
        cout << "Погрешность решения: " << error << endl;
    }

    // Решение для матрицы Гильберта
    vector<vector<double>> A_hilbert = buildHilbertMatrix(n);
    cout << "\nРешение для матрицы Гильберта:" << endl;
    if (LU_SQ_Decomposition(A_hilbert, L, U)) 
    {
        vector<double> y = forwardSubstitution(L, F);
        vector<double> x = backwardSubstitution(U, y);
        double error = computeError(x, x_exact);
        cout << "Погрешность решения: " << error << endl;
    }

    // Метод Гаусса
    vector<vector<double>> A_gauss = A_hilbert;
    vector<double> F_gauss = F;
    cout << "\nРешение с использованием метода Гаусса:" << endl;
    if (gaussianElimination(A_gauss, F_gauss)) 
    {
        cout << "Метод Гаусса успешно выполнен." << endl;
    }

    // Подсчет операций
    cout << "\nЧисло операций для LU-разложения: " << countOperationsLU(n) << endl;
    cout << "Число операций для метода Гаусса: " << countOperationsGauss(n) << endl;

    return 0;
}
