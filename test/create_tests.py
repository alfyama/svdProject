import csv
import os
import time
import re
import fnmatch
import pandas as pd
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


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
        start_time = time.time()
        vector = linalg.svdvals(matrix)
        vector = np.atleast_2d(vector)
        end_time = time.time()
        end_time = end_time - start_time
        end_time = np.atleast_2d(end_time)
        filename_matrix = "test" + str(matrix_size) + "x" + str(matrix_size) + "matrix.csv"
        with open(filename_matrix, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(matrix)
        filename_svd = "test" + str(matrix_size) + "x" + str(matrix_size) + "Nsvd.csv"
        with open(filename_svd, 'w', newline='') as g:
            writer = csv.writer(g)
            writer.writerows(vector)
            writer.writerow(end_time)


def create_test_results():
    file1="test0svd.csv"
    file2="test1svd.csv"

    resultfile0g="../results/test0DOgr.csv"
    resultfile1g="../results/test1DOgr.csv"
    resultfile0d="../results/test0DOdqds.csv"
    resultfile1d="../results/test1DOdqds.csv"

    svdf1 = 0
    svdf2 = 0
    svdgr0 = 0
    svdgr1 = 0
    svddqds0 = 0
    svddqds1 = 0

    # with open(resultfile0g, 'r',newline='') as file:
    #     reader = csv.reader(file)
    #     h = next(reader)
    #     v = next(reader)
    #     singvals = next(reader)
    #     singvals = [float(x) for x in singvals]
    #     svdgr0 = np.array(singvals)
        # intersection = np.dot(np.array(singvals), sv[])

    with open(resultfile1g, 'r',newline='') as file:
        reader = csv.reader(file)
        h = next(reader)
        v = next(reader)
        singvals = next(reader)
        singvals = [float(x) for x in singvals]
        svdgr1 = np.array(singvals)

    # with open(resultfile0d, 'r',newline='') as file:
    #     reader = csv.reader(file)
    #     h = next(reader)
    #     v = next(reader)
    #     singvals = next(reader)
    #     singvals = [float(x) for x in singvals]
    #     svddgds0 = np.array(singvals)

    with open(resultfile1d, 'r',newline='') as file:
        reader = csv.reader(file)
        h = next(reader)
        v = next(reader)
        singvals = next(reader)
        singvals = [float(x) for x in singvals]
        svddgds1 = np.array(singvals)


    # with open(file1, 'r',newline='') as file:
    #     reader = csv.reader(file1)
    #     singvals = next(reader)
    #     singvals = [float(x) for x in singvals]
    #     singvals = np.array(singvals)
    #     intersection = np.dot(singvals, svdgr0)
    #     union = np.sum(singvals) + np.sum(svdgr0) - intersection
    #     print("Similarity for test0svd with test0DOgr {:.2f}".format(intersection/union))

    #     intersection = np.dot(singvals, svdgr1)
    #     union = np.sum(singvals) + np.sum(svdgr1) - intersection
    #     print("Similarity for test0svd with test1DOgr {:.2f}".format(intersection/union))

    #     intersection = np.dot(singvals, svddgds0)
    #     union = np.sum(singvals) + np.sum(svddgds0) - intersection
    #     print("Similarity for test0svd with test0DOdgds {:.2f}".format(intersection/union))
    #     intersection = np.dot(singvals, svddgds1)
    #     union = np.sum(singvals) + np.sum(svddgds1) - intersection
    #     print("Similarity for test0svd with test1DOgdgds {:.2f}".format(intersection/union))

    with open(file2, 'r',newline='') as file:
        reader = csv.reader(file)
        singvals = next(reader)
        singvals = [float(x) for x in singvals]
        singvals = np.array(singvals)
        # intersection = np.dot(singvals, svdgr0)
        # union = np.sum(singvals) + np.sum(svdgr0) - intersection
        # print("Similarity for test1svd with test0DOgr {:.2f}".format(intersection/union))

        sim = np.dot(singvals, svdgr1) / (np.linalg.norm(singvals) * np.linalg.norm(svdgr1))
        union = 100 * (sim + 1) / 2
        print("Similarity for test1svd with test1DOgr {:.2f}%".format(union))

        # intersection = np.dot(singvals, svddgds0)
        # union = np.sum(singvals) + np.sum(svddgds0) - intersection
        # print("Similarity for test1svd with test0DOdgds {:.2f}".format(intersection/union))
        sim = np.dot(singvals, svdgr1) / (np.linalg.norm(singvals) * np.linalg.norm(svdgr1))
        union = 100 * (sim + 1) / 2

        print("Similarity for test1svd with test1DOgdgds {:.2f}".format(union))


def create_plots():
    path = "../results"
    path2 = "."
    timefgr = []
    sizefgr = []
    timedgr = []
    sizedgr = []
    timeldgr = []
    sizeldgr = []

    timefds = []
    sizefds = []
    timedds = []
    sizedds = []
    timeldds = []
    sizeldds = []

    for f in os.listdir(path):
        file_path = os.path.join(path, f)
        if os.path.isfile(file_path):
            with open(file_path, 'r', newline='') as file:
                reader = csv.reader(file)
                headers = next(reader)
                info = next(reader)
                if (str(f).find("gr") != -1):
                    if(str(f).find("F") != -1):
                        timefgr.append(info[0])
                        sizefgr.append(info[1])
                    elif(str(f).find("DO") != -1):
                        timedgr.append(info[0])
                        sizedgr.append(info[1])
                    elif(str(f).find("LD") != -1):
                        timeldgr.append(info[0])
                        sizeldgr.append(info[1])

                if (str(f).find("dqds") != -1):
                    if(str(f).find("F") != -1):
                        timefds.append(info[0])
                        sizefds.append(info[1])
                    elif(str(f).find("DO") != -1):
                        timedds.append(info[0])
                        sizedds.append(info[1])
                    elif(str(f).find("LD") != -1):
                        timeldds.append(info[0])
                        sizeldds.append(info[1])

    timeC = []
    sizeC = []

    # Scipy results
    for f in os.listdir(path2):
        if f.endswith("Nsvd.csv"):
            file_path = os.path.join(path2, f)
            with open(file_path, "r") as file:
                reader = csv.reader(file)
                header = next(reader)
                info = next(reader)
                sizeC.append(re.search("\d+",f).group(0))
                timeC.append(info[0][1:-1])

    fig, ax = plt.subplots()

    dfc = pd.DataFrame({'Size': sizeC, 'Time': timeC})
    dfc['Size'] =  dfc['Size'].astype(float)
    dfc['Time'] =  dfc['Time'].astype(float)
    dfc.sort_values(by="Size",inplace=True)


    df1 = pd.DataFrame({'Size': sizefgr, 'Time':timefgr})
    df1['Size'] =  df1['Size'].astype(float)
    df1['Time'] =  df1['Time'].astype(float)
    df1.sort_values(by="Size",inplace=True)

    df2 = pd.DataFrame({'Size': sizedgr, 'Time':timedgr})
    df2['Size'] =  df2['Size'].astype(float)
    df2['Time'] =  df2['Time'].astype(float)
    df2.sort_values(by="Size",inplace=True)

    df3 = pd.DataFrame({'Size': sizeldgr, 'Time':timeldgr})
    df3['Size'] =  df3['Size'].astype(float)
    df3['Time'] =  df3['Time'].astype(float)
    df3.sort_values(by="Size",inplace=True)

    df4 = pd.DataFrame({'Size': sizefds, 'Time':timefds})
    df4['Size'] =  df4['Size'].astype(float)
    df4['Time'] =  df4['Time'].astype(float)
    df4.sort_values(by="Size",inplace=True)

    df5 = pd.DataFrame({'Size': sizedds, 'Time':timedds})
    df5['Size'] =  df5['Size'].astype(float)
    df5['Time'] =  df5['Time'].astype(float)
    df5.sort_values(by="Size",inplace=True)

    df6 = pd.DataFrame({'Size': sizeldds, 'Time':timeldds})
    df6['Size'] =  df6['Size'].astype(float)
    df6['Time'] =  df6['Time'].astype(float)
    df6.sort_values(by="Size",inplace=True)


    dfc.plot(x='Size', y='Time', label='Scipy values',ax=ax)
    df1.plot(x="Size", y="Time", label="float-gr", ax=ax)
    df2.plot(x="Size", y="Time", label="double-gr", ax=ax)
    df3.plot(x="Size", y="Time", label="long double-gr", ax=ax)
    df4.plot(x="Size", y="Time", label="float-dqds", ax=ax)
    df5.plot(x="Size", y="Time", label="double-dqds", ax=ax)
    df6.plot(x="Size", y="Time", label="long double-dqds", ax=ax)

    ax.set_xlabel('Size of matrix N')
    ax.set_ylabel('Time seconds')
    ax.set_title("Performance comparison")
    plt.legend()

    plt.show()

if __name__ == '__main__':
    createTestCases()
    generate_random_matrices(5, 5, 20)
    create_plots()
    create_test_results()
