﻿using System;

class GaussMethod
{
    static void Main()
    {
        double[,] matrix = {
            { 2, 1, -1, 8 },
            { -3, -1, 2, -11 },
            { -2, 1, 2, -3 }
        };

        int n = matrix.GetLength(0); // Розмірність матриці

        Console.WriteLine("Початкова матриця:");
        PrintMatrix(matrix);

        // Застосовуємо метод Гаусса з вибором головного елемента по стовпцях
        for (int k = 0; k < n - 1; k++)
        {
            int maxRow = FindMaxRow(matrix, k);
            SwapRows(matrix, k, maxRow);

            for (int i = k + 1; i < n; i++)
            {
                double factor = matrix[i, k] / matrix[k, k];
                for (int j = k; j < n + 1; j++)
                {
                    matrix[i, j] -= factor * matrix[k, j];
                }
            }
        }

        Console.WriteLine("\nМатриця після застосування методу Гаусса:");
        PrintMatrix(matrix);
  
        // Знаходимо розв'язок системи
        double[] solution = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < n; j++)
            {
                sum += matrix[i, j] * solution[j];
            }
            solution[i] = (matrix[i, n] - sum) / matrix[i, i];
        }
        Console.WriteLine("\nРозв'язок системи:");
        PrintVector(solution);

        // Знаходимо вектор нев'язки
        double[] residuals = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < n; j++)
            {
                sum += matrix[i, j] * residuals[j];
            }
            residuals[i] = (matrix[i, n] - sum) / matrix[i, i];
        }
        Console.WriteLine("\nВектор нев'язки:");
        PrintVector(residuals);

        // Знаходимо число обумовленості матриці
        double conditionNumber = FindConditionNumber(matrix);
        Console.WriteLine("\nЧисло обумовленості матриці: " + conditionNumber);

        // Знаходимо визначник матриці
        double determinant = 1;
        for (int i = 0; i < n; i++)
        {
            determinant *= matrix[i, i];
        }

        Console.WriteLine("\nВизначник матриці: " + determinant);

        double[,] inverseMatrix = InverseMatrix(matrix);
        Console.WriteLine("\nОберерна матриця:");
        PrintMatrix(inverseMatrix);
        double[,] MultipliedMatrix = MultiplyInverseMatrix(inverseMatrix, matrix);
        Console.WriteLine("\nМноження матриць:");
        PrintMatrix(MultipliedMatrix);
    }



    // Функція для знаходження рядка з максимальним елементом у стовпці k
    static int FindMaxRow(double[,] matrix, int k)
    {
        int maxRow = k;
        double maxValue = Math.Abs(matrix[k, k]);

        for (int i = k + 1; i < matrix.GetLength(0); i++)
        {
            double value = Math.Abs(matrix[i, k]);
            if (value > maxValue)
            {
                maxValue = value;
                maxRow = i;
            }
        }

        return maxRow;
    }

    // Функція для обміну рядків матриці
    static void SwapRows(double[,] matrix, int row1, int row2)
    {
        int n = matrix.GetLength(1);
        for (int j = 0; j < n; j++)
        {
            double temp = matrix[row1, j];
            matrix[row1, j] = matrix[row2, j];
            matrix[row2, j] = temp;
        }
    }

    // Функція для виведення матриці на консоль
    static void PrintMatrix(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                Console.Write(matrix[i, j] + "\t\t");
            }
            Console.WriteLine();
        }
    }

    // Функція для виведення вектора на консоль
    static void PrintVector(double[] vector)
    {
        int n = vector.Length;

        for (int i = 0; i < n; i++)
        {
            Console.WriteLine(vector[i]);
        }
    }

    // Функція для знаходження числа обумовленості матриці
    static double FindConditionNumber(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        double maxNorm = FindMatrixNorm(matrix);
        double[,] inverseMatrix = InverseMatrix(matrix);
        double inverseNorm = FindMatrixNorm(inverseMatrix);
        double conditionNumber = maxNorm * inverseNorm;
        return conditionNumber;
    }

    // Функція для знаходження норми матриці
    static double FindMatrixNorm(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        double maxNorm = 0;

        for (int i = 0; i < n; i++)
        {
            double rowSum = 0;
            for (int j = 0; j < n; j++)
            {
                rowSum += Math.Abs(matrix[i, j]);
            }
            if (rowSum > maxNorm)
            {
                maxNorm = rowSum;
            }
        }

        return maxNorm;
    }

    // Функція для знаходження оберненої матриці
    static double[,] InverseMatrix(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        double[,] augmentedMatrix = new double[n, 2 * n];

        // Створюємо розширену матрицю, де права частина - одинична матриця
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmentedMatrix[i, j] = matrix[i, j];
            }
            augmentedMatrix[i, i + n] = 1;
        }

        // Застосовуємо метод Гаусса з вибором головного елемента по стовпцях до розширеної матриці
        for (int k = 0; k < n; k++)
        {
            int maxRow = FindMaxRow(augmentedMatrix, k);
            SwapRows(augmentedMatrix, k, maxRow);

            for (int i = k + 1; i < n; i++)
            {
                double factor = augmentedMatrix[i, k] / augmentedMatrix[k, k];
                for (int j = k; j < 2 * n; j++)
                {
                    augmentedMatrix[i, j] -= factor * augmentedMatrix[k, j];
                }
            }
        }

        // Зводимо ліву частину до одиничної матриці
        for (int i = n - 1; i > 0; i--)
        {
            for (int j = i - 1; j >= 0; j--)
            {
                double factor = augmentedMatrix[j, i] / augmentedMatrix[i, i];
                for (int k = i; k < 2 * n; k++)
                {
                    augmentedMatrix[j, k] -= factor * augmentedMatrix[i, k];
                }
            }
        }

        // Нормалізуємо рядки, щоб на головній діагоналі отримати одиниці
        for (int i = 0; i < n; i++)
        {
            double factor = augmentedMatrix[i, i];
            for (int j = i; j < 2 * n; j++)
            {
                augmentedMatrix[i, j] /= factor;
            }
        }

        // Виділяємо обернену матрицю з розширеної матриці
        double[,] inverseMatrix = new double[n, n];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                inverseMatrix[i, j] = augmentedMatrix[i, j + n];
            }
        }

        return inverseMatrix;
    }
    // Метод для множення оберненої матриці на матрицю
    static double[,] MultiplyInverseMatrix(double[,] inverseMatrix, double[,] matrix)
    {
        int n = inverseMatrix.GetLength(0);
        int m = matrix.GetLength(1);
        double[,] result = new double[n, m];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                {
                    sum += inverseMatrix[i, k] * matrix[k, j];
                }
                result[i, j] = sum;
            }
        }

        return result;
    }

}