#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

// Общие переменные для работы с матрицей
int n, nn;

// Функция для нахождения максимального значения
int max_(int a, int b) {
    return (a > b) ? a : b;
}

// Функция LU-разложения матрицы
void calc_lusq(const std::vector<int>& ia, std::vector<float>& di, std::vector<float>& al, std::vector<float>& au) {
    for (int i = 0; i < n; ++i) {
        int num_elem_i = ia[i + 1] - ia[i];  // Количество элементов в i-ой строке
        int ind = ia[i];  // Индекс в массиве al
        float sum_ = 0.0;  // Сумма для диагонали

        for (int j = i - num_elem_i; j < i; ++j) {
            float sum_l = 0.0;
            float sum_u = 0.0;
            int num_elem_j = ia[j + 1] - ia[j];  // Количество элементов в j-ой строке

            int k = max_(i - num_elem_i, j - num_elem_j);
            int ind1 = ia[i] - (i - num_elem_i) + k;
            int ind2 = ia[j] - (j - num_elem_j) + k;
            int end_k = j - 1 - k + ind1;

            for (; ind1 <= end_k; ++ind1) {
                sum_l += al[ind1] * au[ind2];
                sum_u += au[ind1] * al[ind2];
                ind2++;
            }

            al[ind] = (al[ind] - sum_l) / di[j];
            au[ind] = (au[ind] - sum_u) / di[j];
            sum_ += al[ind] * au[ind];
            ind++;
        }

        if (di[i] - sum_ <= 0.0) {
            std::cerr << "Matrix is NOT LU(sq) decomposable!" << std::endl;
            exit(1);
        }

        di[i] = std::sqrt(di[i] - sum_);
    }
}

// Функция для вычисления вектора y
void calc_y(const std::vector<int>& ia, const std::vector<float>& di, const std::vector<float>& al, std::vector<float>& vect) {
    for (int i = 0; i < n; ++i) {
        float sum_ = 0.0;
        int num_elem_i = ia[i + 1] - ia[i];
        int ind = ia[i];

        for (int j = i - num_elem_i; j < i; ++j) {
            sum_ += al[ind] * vect[j];
            ind++;
        }

        vect[i] = (vect[i] - sum_) / di[i];
    }
}

// Функция для вычисления вектора x
void calc_x(const std::vector<int>& ia, const std::vector<float>& di, const std::vector<float>& au, std::vector<float>& vect) {
    for (int i = n - 1; i >= 0; --i) {
        float xi = vect[i];
        xi /= di[i];
        vect[i] = xi;

        int ind = ia[i + 1] - 1;
        int num_elem_i = ia[i + 1] - ia[i];

        for (int j = i - 1; j >= i - num_elem_i; --j) {
            vect[j] -= au[ind] * xi;
            ind--;
        }
    }
}

// Функция для чтения данных из файла
void read_(std::vector<float>& mem, std::vector<int>& ia, int& di, int& al, int& au, int& vect) {
    std::ifstream file_matrix("profile.txt");
    std::ifstream file_vector("vector.txt");

    if (!file_matrix.is_open() || !file_vector.is_open()) {
        std::cerr << "ERROR: Can't read profile.txt or vector.txt" << std::endl;
        exit(1);
    }

    file_matrix >> n;
    ia.resize(n + 1);
    for (int i = 0; i < n + 1; ++i) {
        file_matrix >> ia[i];
    }

    nn = ia[n] - 1;
    di = 0;
    mem.resize(nn * 2 + n + 1);

    for (int i = 0; i < n; ++i) {
        file_matrix >> mem[di + i];
    }

    al = di + n + 1;
    for (int i = 0; i < nn; ++i) {
        file_matrix >> mem[al + i];
    }

    au = al + nn;
    for (int i = 0; i < nn; ++i) {
        file_matrix >> mem[au + i];
    }

    vect = au + nn;
    for (int i = 0; i < n; ++i) {
        file_vector >> mem[vect + i];
    }
}

// Функция для записи данных в файл
void write_(const std::vector<float>& vect) {
    std::ofstream file_out("tmp3.txt");
    for (int i = 0; i < n; ++i) {
        file_out << vect[i] << std::endl;
    }
}

int main() {
    std::vector<int> ia(4097);
    int di, al, au, vect;
    std::vector<float> mem;

    // Чтение данных
    read_(mem, ia, di, al, au, vect);

    // Создаем временные векторы для передачи в функции
    std::vector<float> di_vec(mem.begin() + di + 1, mem.begin() + di + 1 + n);
    std::vector<float> al_vec(mem.begin() + al + 1, mem.begin() + al + 1 + nn);
    std::vector<float> au_vec(mem.begin() + au + 1, mem.begin() + au + 1 + nn);
    std::vector<float> vect_vec(mem.begin() + vect + 1, mem.begin() + vect + 1 + n);

    // Вызываем функции с исправленными аргументами
    calc_lusq(ia, di_vec, al_vec, au_vec);
    calc_y(ia, di_vec, al_vec, vect_vec);
    calc_x(ia, di_vec, au_vec, vect_vec);

    // Записываем результат в файл
    write_(vect_vec);

    // Выводим результат на экран
    for (int i = 0; i < n; ++i) {
        std::cout << vect_vec[i] << std::endl;
    }

    return 0;
}