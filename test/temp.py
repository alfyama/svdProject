import csv

import numpy as np
from scipy import linalg


def testCase1(type):
    types = {
        "float": np.float32,
        "double": np.float64,
        "longdouble": np.longdouble
    }
    def f(i, j):
        return np.where(i > j, 0, np.where(i == j, 21 - i, -1))
    return np.fromfunction(f, (20, 21), dtype=types[type])

def testCase2(type):
    types = {
        "float": np.float32,
        "double": np.float64,
        "longdouble": np.longdouble
    }
    def f(i, j):
        return np.where(i > j, 0, np.where(i == j, 1, -1))
    return np.fromfunction(f, (30, 30), dtype=types[type])


def testCase3(type):
    types = {
        "float": np.float32,
        "double": np.float64,
        "longdouble": np.longdouble
    }
    return np.random.randint(1, 11, size=(4, 3), dtype=types[type])


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
    types = ["float", "double", "longdouble"]

    matrices = []

    for t in types:
        matrices.append(testCase1(t))
        matrices.append(testCase2(t))
        matrices.append(testCase3(t))

    for i in range(len(matrices)):
        filename = "test" + str(i) + "matrix.csv"
        write_matrix_to_file(matrices[i], filename)
        v = solveCase(matrices[i])
        filename = "test" + str(i) + "svd.csv"
        write_vector_to_file(v, filename)


def solveCase(matrix):
    return linalg.svdvals(matrix)



if __name__ == '__main__':
    createTestCases()
