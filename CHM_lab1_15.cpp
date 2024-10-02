#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>
#include <string>
#include <limits>

// Выбор точности: раскомментируйте нужное определение
typedef double precision_t;
// typedef float precision_t;

using namespace std;

// Функция для вычисления ошибки между найденным и точным решениями
template <typename T>
double computeError(const vector<T>& x, const vector<T>& x_exact) {
    if (x.size() != x_exact.size()) {
        cerr << "Error: Векторы x и x_exact должны иметь одинаковый размер." << endl;
        return -1.0;
    }
    double error = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        error += pow(static_cast<double>(x[i] - x_exact[i]), 2);
    }
    return sqrt(error);
}

// Функция для вывода матрицы
template <typename T>
void printMatrix(const vector<vector<T>>& matrix) {
    for (const auto& row : matrix) {
        for (T val : row)
            cout << setw(10) << fixed << setprecision(4) << val << " ";
        cout << endl;
    }
}

// Функция для вывода вектора
template <typename T>
void printVector(const vector<T>& vec) {
    for (T v : vec)
        cout << setw(10) << fixed << setprecision(4) << v << " ";
    cout << endl;
}

// Чтение матрицы из файла
bool readMatrixFromFile(const string& filename, vector<vector<double>>& matrix) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Не удалось открыть файл матрицы: " << filename << endl;
        return false;
    }

    int n;
    file >> n;
    if (n <= 0) {
        cerr << "Некорректный размер матрицы: " << n << endl;
        return false;
    }

    matrix.assign(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(file >> matrix[i][j])) {
                cerr << "Ошибка при чтении элементов матрицы в строке " << i + 1 << ", столбце " << j + 1 << endl;
                return false;
            }
        }
    }

    file.close();
    return true;
}

// Чтение вектора из файла
bool readVectorFromFile(const string& filename, vector<double>& vec) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Не удалось открыть файл вектора: " << filename << endl;
        return false;
    }

    int n;
    file >> n;
    if (n <= 0) {
        cerr << "Некорректный размер вектора: " << n << endl;
        return false;
    }

    vec.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        if (!(file >> vec[i])) {
            cerr << "Ошибка при чтении элементов вектора в позиции " << i + 1 << endl;
            return false;
        }
    }

    file.close();
    return true;
}

// LU-разложение для симметричных положительно определённых матриц (Cholesky-подобное)
bool LU_SQ_Decomposition(const vector<vector<double>>& A_input, vector<vector<double>>& L, vector<vector<double>>& U, unsigned long long& operation_count) 
{
    int n = A_input.size();
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j <= i; j++) 
        {
            double sum = 0.0;
            if (i == j) 
            { 
                for (int k = 0; k < j; k++) 
                {
                    sum += L[j][k] * L[j][k];
                    operation_count += 2; // Умножение и сложение
                }
                double value = A_input[j][j] - sum;
                operation_count += 1; // Вычитание
                if (value <= 0) 
                {
                    cout << "Матрица не является положительно определенной." << endl;
                    // return false;
                }
                L[j][j] = sqrt(value);
                operation_count += 1; // Взятие квадратного корня
                U[j][j] = L[j][j];
                operation_count += 1; // Присваивание
            } 
            else 
            {
                for (int k = 0; k < j; k++) 
                {
                    sum += L[i][k] * L[j][k];
                    operation_count += 2; // Умножение и сложение
                }
                L[i][j] = (A_input[i][j] - sum) / L[j][j];
                operation_count += 2; // Вычитание и деление
                U[j][i] = L[i][j];
                operation_count += 1; // Присваивание
            }
        }
    }
    return true;
}

// Прямой ход (решение L*y = F)
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& F, unsigned long long& operation_count) 
{
    int n = L.size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++) 
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++) 
        {
            sum += L[i][j] * y[j];
            operation_count += 2; // Умножение и сложение
        }
        y[i] = (F[i] - sum) / L[i][i];
        operation_count += 2; // Вычитание и деление
    }
    return y;
}

// Обратный ход (решение U*x = y)
vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y, unsigned long long& operation_count) 
{
    int n = U.size();
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) 
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) 
        {
            sum += U[i][j] * x[j];
            operation_count += 2; // Умножение и сложение
        }
        x[i] = (y[i] - sum) / U[i][i];
        operation_count += 2; // Вычитание и деление
    }
    return x;
}

// Метод Гаусса с выбором главного элемента и подсчетом операций
bool GaussianElimination_func(vector<vector<precision_t>> A, vector<precision_t> b, vector<precision_t>& x, unsigned long long& operation_count) {
    int n = A.size();
    vector<int> pivot(n);
    for (int i = 0; i < n; ++i)
        pivot[i] = i;

    // Прямой ход
    for (int k = 0; k < n; ++k) {
        // Выбор главного элемента
        precision_t max_val = abs(A[k][k]);
        int maxRow = k;
        for (int i = k + 1; i < n; ++i) {
            if (abs(A[i][k]) > max_val) {
                max_val = abs(A[i][k]);
                maxRow = i;
            }
            operation_count++; // Сравнение
        }

        if (max_val == 0) {
            cerr << "Матрица вырождена!" << endl;
            return false;
        }

        // Перестановка строк
        if (maxRow != k) {
            swap(A[k], A[maxRow]);
            swap(b[k], b[maxRow]);
            swap(pivot[k], pivot[maxRow]);
            operation_count += 3 * n; // Перестановка строк
        }

        // Нули под диагональю
        for (int i = k + 1; i < n; ++i) {
            precision_t factor = A[i][k] / A[k][k];
            operation_count += 1; // Деление
            A[i][k] = 0;
            for (int j = k + 1; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
                operation_count += 2; // Умножение и вычитание
            }
            b[i] -= factor * b[k];
            operation_count += 1; // Вычитание
        }
    }

    // Обратный ход
    x.assign(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        if (A[i][i] == 0) {
            cerr << "Нулевой элемент на диагонали в строке " << i << endl;
            return false;
        }
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
            operation_count += 2; // Умножение и вычитание
        }
        x[i] /= A[i][i];
        operation_count += 1; // Деление
    }
    return true;
}

int main() {
    // Файлы ввода
    string matrixFile = "matrix.txt";
    string vectorFile = "vector.txt";

    // Чтение матрицы и вектора из файлов
    vector<vector<double>> A;
    if (!readMatrixFromFile(matrixFile, A)) {
        return 1;
    }

    vector<double> b;
    if (!readVectorFromFile(vectorFile, b)) {
        return 1;
    }

    int n = A.size();
    if (b.size() != n) {
        cerr << "Размер вектора правой части не соответствует размерности матрицы." << endl;
        return 1;
    }

    // Задание точного решения (предполагается, что оно известно)
    // В этом примере x_exact = (1, 2, ..., n)
    vector<double> x_exact(n, 0.0);
    for (int i = 0; i < n; ++i)
        x_exact[i] = i + 1;

    // ---------- LU_SQ_Decomposition ---------- //
    cout << "Выполнение LU_SQ_Decomposition..." << endl;
    vector<vector<double>> L, U;
    unsigned long long operations_LU_SQ = 0;

    auto start_LU_SQ = chrono::high_resolution_clock::now();
    bool success_LU_SQ = LU_SQ_Decomposition(A, L, U, operations_LU_SQ);
    auto end_LU_SQ = chrono::high_resolution_clock::now();
    double time_LU_SQ = chrono::duration<double>(end_LU_SQ - start_LU_SQ).count();

    // Объявляем x_LU_SQ вне блока if-else
    vector<double> x_LU_SQ;

    if (!success_LU_SQ) {
        cout << "ну ты и кринжа навалил!" << endl; 
        // cerr << "LU_SQ_Decomposition не удалось." << endl;
    } else {
        cout << "Матрица L:" << endl;
        printMatrix(L);
        cout << endl;
        cout << "Матрица U:" << endl;
        printMatrix(U);
        cout << endl;

        // Прямой ход
        unsigned long long operations_forward_LU_SQ = 0;
        vector<double> y = forwardSubstitution(L, b, operations_forward_LU_SQ);

        // Обратный ход
        unsigned long long operations_backward_LU_SQ = 0;
        x_LU_SQ = backwardSubstitution(U, y, operations_backward_LU_SQ);

        // Вычисление полной погрешности
        double error_LU_SQ = computeError(x_LU_SQ, x_exact);

        // Вывод результатов LU_SQ_Decomposition
        cout << "Результаты LU_SQ_Decomposition:" << endl;
        cout << "Время выполнения: " << fixed << setprecision(6) << time_LU_SQ << " секунд" << endl;
        cout << "Количество операций:" << endl;
        cout << "  LU_SQ_Decomposition: " << operations_LU_SQ << endl;
        cout << "  Forward Substitution: " << operations_forward_LU_SQ << endl;
        cout << "  Backward Substitution: " << operations_backward_LU_SQ << endl;
        cout << "  Всего: " << operations_LU_SQ + operations_forward_LU_SQ + operations_backward_LU_SQ << endl;
        cout << "Погрешность: " << fixed << setprecision(6) << error_LU_SQ << endl;
        cout << endl;
    }

    // ---------- Метод Гаусса ---------- //
    cout << "Выполнение метода Гаусса..." << endl;
    // Преобразуем A и b в тип precision_t, если это необходимо
    vector<vector<precision_t>> A_Gauss(n, vector<precision_t>(n, 0.0));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            A_Gauss[i][j] = static_cast<precision_t>(A[i][j]);

    vector<precision_t> b_Gauss(n, 0.0);
    for(int i=0;i<n;i++)
        b_Gauss[i] = static_cast<precision_t>(b[i]);

    vector<precision_t> x_Gauss;
    unsigned long long operations_Gauss = 0;

    auto start_Gauss = chrono::high_resolution_clock::now();
    bool success_Gauss = GaussianElimination_func(A_Gauss, b_Gauss, x_Gauss, operations_Gauss);
    auto end_Gauss = chrono::high_resolution_clock::now();
    double time_Gauss = chrono::duration<double>(end_Gauss - start_Gauss).count();

    // Объявляем x_Gauss_double вне блока if-else
    vector<double> x_Gauss_double;

    if (!success_Gauss) {
        cerr << "Метод Гаусса не удалось." << endl;
    } else {
        // Преобразование x_Gauss к типу double для совместимости
        x_Gauss_double.assign(x_Gauss.begin(), x_Gauss.end());

        // Вычисление полной погрешности для метода Гаусса
        double error_Gauss = computeError(x_Gauss_double, x_exact);

        // Вывод результатов метода Гаусса
        cout << "Результаты метода Гаусса:" << endl;
        cout << "Время выполнения: " << fixed << setprecision(6) << time_Gauss << " секунд" << endl;
        cout << "Количество операций: " << operations_Gauss << endl;
        cout << "Погрешность: " << fixed << setprecision(6) << error_Gauss << endl;
        cout << endl;
    }

    // ---------- Сравнение Результатов ---------- //
    cout << "Сравнение LU_SQ_Decomposition и метода Гаусса:" << endl;

    // Проверка успешности обоих методов
    if (success_LU_SQ && success_Gauss) {
        double error_LU_SQ = computeError(x_LU_SQ, x_exact);
        double error_Gauss = computeError(x_Gauss_double, x_exact);

        cout << "Погрешность LU_SQ_Decomposition: " << fixed << setprecision(6) << error_LU_SQ << endl;
        cout << "Погрешность метода Гаусса: " << fixed << setprecision(6) << error_Gauss << endl;
        cout << "Количество операций LU_SQ_Decomposition: " << (operations_LU_SQ + 0 /* Forward and Backward are printed earlier */) << endl;
        cout << "Количество операций метода Гаусса: " << operations_Gauss << endl;
    } else {
        cout << "Некоторые методы не были выполнены успешно." << endl;
    }

    return 0;
}