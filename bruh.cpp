#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <random>
#include <iomanip>
#include <cmath>

// Размерность матрицы
const int size = 1000;

// Функция для генерации части матрицы
void generateMatrixSection(std::vector<std::vector<int>>& A, int startRow, int endRow, std::mt19937& gen, std::uniform_int_distribution<> distrib, std::uniform_int_distribution<> distrib_diag) {
    for (int i = startRow; i < endRow; ++i) {
        int rowSum = 0;
        for (int j = 0; j < size; ++j) {
            if (j < i) {
                // Копируем симметричное значение
                A[i][j] = A[j][i];
            } else if (j == i) {
                // Диагональный элемент будет определен после генерации остальных элементов
                A[i][j] = 0;
            } else {
                // Генерируем случайное значение для верхней треугольной части
                A[i][j] = distrib(gen);
                // Копируем в нижнюю треугольную часть для симметрии
                A[j][i] = A[i][j];
                rowSum += std::abs(A[i][j]);
            }
        }
        // Устанавливаем диагональный элемент
        // Обеспечиваем диагональную доминантность
        A[i][i] = rowSum + distrib_diag(gen);
    }
}

int main() {
    // Инициализация матрицы нулями
    std::vector<std::vector<int>> A(size, std::vector<int>(size, 0));
    
    // Настройка генераторов случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(1, 100);      // Для ненулевых элементов
    std::uniform_int_distribution<> distrib_diag(1, 10);  // Для диагональных элементов
    
    // Определение количества доступных потоков
    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) numThreads = 4; // Если определить не удалось, используем 4 потока по умолчанию
    
    // Определение диапазонов строк для каждого потока
    std::vector<std::thread> threads;
    int rowsPerThread = size / numThreads;
    int remainingRows = size % numThreads;
    int currentRow = 0;
    
    for (unsigned int t = 0; t < numThreads; ++t) {
        int startRow = currentRow;
        int endRow = startRow + rowsPerThread;
        if (t == numThreads - 1) {
            endRow += remainingRows; // Последнему потоку даем оставшиеся строки
        }
        threads.emplace_back(generateMatrixSection, std::ref(A), startRow, endRow, std::ref(gen), std::ref(distrib), std::ref(distrib_diag));
        currentRow = endRow;
    }
    
    // Ожидание завершения всех потоков
    for (auto& th : threads) {
        th.join();
    }
    
    // Запись матрицы в файл
    std::ofstream file("matrix.txt");
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла matrix.txt для записи!" << std::endl;
        return 1;
    }
    
    // Запись размерности матрицы
    file << size << "\n";
    
    // Запись элементов матрицы построчно
    for (const auto& row : A) {
        for (int j = 0; j < size; ++j) {
            file << row[j];
            if (j < size - 1) file << " ";
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "Симметричная положительная матрица записана в файл 'matrix.txt'" << std::endl;
    
    return 0;
}