#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>
#include <random>
#include <iomanip>

const int size = 1000; // Размер матрицы
std::mutex mtx; // Мьютекс для защиты записи в файл

// Функция для генерации строки матрицы
std::string generateRow(int rowIndex) {
    std::random_device rd; // Источник случайных чисел
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(1, 100); // Распределение от 1 до 100

    std::string row;

    for (int j = 0; j < size; ++j) {
        if (j < rowIndex) {
            // Используем значение из предыдущей строки для симметрии
            row += std::to_string(distrib(gen)) + (j < size - 1 ? " " : "");
        } else {
            int value = distrib(gen); // Генерируем случайное значение
            row += std::to_string(value) + (j < size - 1 ? " " : "");
        }
    }
    return row;
}

// Функция для записи строк в файл
void writeRowsToFile(std::ofstream &file, int start, int end) {
    for (int i = start; i < end; ++i) {
        std::string row = generateRow(i);
        std::lock_guard<std::mutex> lock(mtx); // Защита записи в файл
        file << row << "\n";
    }
}

int main() {
    std::ofstream file("matrix.txt");

    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла!" << std::endl;
        return 1;
    }

    const int numThreads = std::thread::hardware_concurrency(); // Количество доступных потоков
    std::vector<std::thread> threads;

    int rowsPerThread = size / numThreads; // Количество строк для каждого потока

    // Запуск потоков
    for (int i = 0; i < numThreads; ++i) {
        int start = i * rowsPerThread;
        int end = (i == numThreads - 1) ? size : start + rowsPerThread; // Последнему потоку даем оставшиеся строки
        threads.emplace_back(writeRowsToFile, std::ref(file), start, end);
    }

    // Ожидание завершения всех потоков
    for (auto &t : threads) {
        t.join();
    }

    file.close();
    std::cout << "Симметричная матрица записана в файл 'matrix.txt'" << std::endl;

    return 0;
}
