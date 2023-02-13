import csv

import numpy as np
from scipy import linalg


def testCase1():
    def f(i, j):
        return np.where(i > j, 0, np.where(i == j, 21 - i, -1))
    return np.fromfunction(f, (20, 21))

def testCase2():
    def f(i, j):
        return np.where(i > j, 0, np.where(i == j, 1, -1))
    return np.fromfunction(f, (30, 30))


def testCases3():
    return np.random.randint(1, 11,size=(4, 3))


def testCase4():
    a = [1,6,3,6,4,6,2,6,6,10,7,2]
    a = np.array(a)
    a.resize((4, 3))
    return a

def testCasebidiag():
    a = [2.3 ,9.7 ,0.0
        ,0.0 ,5.5 ,5.4
        ,0.0 ,0.0 ,7.12
        ,0.0 ,0.0 ,0.0]
    a = np.array(a)
    a.resize((4,3))
    return a



def write_matrix_to_file(matrix, filename):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        for row in matrix:
            writer.writerow(row)

def write_vector_to_file(vec, filename):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(vec)


def createTestCases():
    m1 = testCase1()
    m2 = testCase2()
    m3 = testCases3()
    m4 = testCase4()
    m5 = testCasebidiag()

    matrices = [m1, m2, m3, m4, m5]

    for i in range(5):
        filename = "test" + str(i) + "matrix.csv"
        write_matrix_to_file(matrices[i], filename)
        v = solveCase(matrices[i])
        filename = "test" + str(i) + "svd.csv"
        write_vector_to_file(v, filename)


def solveCase(matrix):
    return linalg.svdvals(matrix)


def generate_random_matrices(start_size, size_increment, num_matrices):
    for i in range(num_matrices):
        matrix_size = start_size + i * size_increment
        matrix = np.random.randint(0, 21, size=(matrix_size, matrix_size))
        vector = linalg.svdvals(matrix)
        vector = np.atleast_2d(vector)
        filename_matrix = "testmatrix_" + str(matrix_size) + "x" + str(matrix_size) + ".csv"
        with open(filename_matrix, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(matrix)
        filename_svd = "testmatrix_" + str(matrix_size) + "x" + str(matrix_size) + "svd.csv"
        with open(filename_svd, 'w', newline='') as g:
            writer = csv.writer(g)
            writer.writerows(vector)




if __name__ == '__main__':
    createTestCases()
    generate_random_matrices(5, 5, 10)
