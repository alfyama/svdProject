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
    return np.random.randint(1, 11,size = (4,3))


def testCase4():
    a = [1,6,3,6,4,6,2,6,6,10,7,2]
    a = np.array(a)
    a.resize((4,3))
    print(a)
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

    matrices = [m1, m2, m3, m4]

    for i in range(4):
        filename = "test" + str(i) + "matrix.csv"
        write_matrix_to_file(matrices[i], filename)
        v = solveCase(matrices[i])
        filename = "test" + str(i) + "svd.csv"
        write_vector_to_file(v, filename)


def solveCase(matrix):
    return linalg.svdvals(matrix)



if __name__ == '__main__':
    createTestCases()
