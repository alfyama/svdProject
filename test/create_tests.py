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


def testCases():
    return np.random.rand()


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

    write_matrix_to_file(m1, "test1matrix.csv")
    write_matrix_to_file(m2, "test2matrix.csv")

    v1 = solveCase(m1)
    v2 = solveCase(m2)

    write_vector_to_file(v1, "test1svd.csv")
    write_vector_to_file(v2, "test2svd.csv")


def solveCase(matrix):
    return linalg.svdvals(matrix)



if __name__ == '__main__':
    createTestCases()
