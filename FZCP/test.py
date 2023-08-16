import sys
import time

sys.path.append("..")
import pandas as pd
import solv_fzcp as sfzcp
import uvarprob as uvpr


def log_point(x, points_list):
    points_list.append(x)


def read_problems(fname):
    data = pd.read_csv(fname, index_col='name', comment='#')
    return data


points_db = {}
psl_lipint_points_list = []
psl_lip_points_list = []
psqe_lipint_points_list = []
psqe_lip_points_list = []
bnb2_lipint_points_list = []
bnb2_lip_points_list = []


def getMax(matrix):
    max_values = []
    for j in range(len(matrix[0])):
        max_value = matrix[0][j]  # 初始化为第一个元素
        for i in range(len(matrix)):
            if matrix[i][j] > max_value:
                max_value = matrix[i][j]
        max_values.append(max_value)
    return max_values


def getMaxRatio(matrix):
    max_values = []
    end = len(matrix[0]) - 1
    for j in range(len(matrix[0])):
        max_value = matrix[0][j] / matrix[0][end]  # 初始化为第一个元素
        for i in range(len(matrix)):
            if matrix[i][j] / matrix[i][end] > max_value:
                max_value = matrix[i][j] / matrix[i][end]
        max_values.append(max_value)
    return max_values


def getMin(matrix):
    min_values = []
    for j in range(len(matrix[0])):
        min_value = matrix[0][j]  # 初始化为第一个元素
        for i in range(len(matrix)):
            if matrix[i][j] < min_value:
                min_value = matrix[i][j]
        min_values.append(min_value)
    return min_values


def getMinRatio(matrix):
    min_values = []
    end = len(matrix[0]) - 1
    for j in range(len(matrix[0])):
        min_value = matrix[0][j] / matrix[0][end]  # 初始化为第一个元素
        for i in range(len(matrix)):
            if matrix[i][j] / matrix[i][end] < min_value:
                min_value = matrix[i][j] / matrix[i][end]
        min_values.append(min_value)
    return min_values


def getAvg(matrix):
    avg_values = []
    for j in range(len(matrix[0])):
        sum_value = 0
        for i in range(len(matrix)):
            sum_value += matrix[i][j]
        avg_value = sum_value / len(matrix)
        avg_values.append(avg_value)
    return avg_values


def getAvgRatio(matrix):
    avg_values = []
    end = len(matrix[0]) - 1
    for j in range(len(matrix[0])):
        sum_value = 0
        for i in range(len(matrix)):
            sum_value += matrix[i][j] / matrix[i][end]
        avg_value = sum_value / len(matrix)
        avg_values.append(avg_value)
    return avg_values


def print_row1(row, index):
    print('%s & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f \\\\' % (index,
                                                                               row[0], row[1],
                                                                               row[2], row[3],
                                                                               row[4], row[5],
                                                                               row[6], row[7]))


def print_row2(row, index):
    print('%s & %d & %d & %d & %d & %d & %d & %d & %d \\\\' % (index,
                                                               row[0], row[1],
                                                               row[2], row[3],
                                                               row[4], row[5],
                                                               row[6], row[7]))


def print_row3(row, index):
    print('%s' % index, end=' ')
    for roxi in row:
        print('& %.1f' % roxi, end=' ')

    print('\\\\')
    # print('%s & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f \\\\' % (index,
    #                                                                            row[0], row[1],
    #                                                                            row[2], row[3],
    #                                                                            row[4], row[5],
    #                                                                            row[6], row[7]))


def test_alg2(df, eps, global_lip, repeat):
    # print('test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R')
    it_list = []
    time_list = []

    for test in df.itertuples():
        it_list_row = [None] * 8
        time_list_row = [0.] * 8
        points_db[test.Index] = {'bnb2_pslint_points_list': []}
        prob = uvpr.UniVarProblem(test.Index, test.objective, test.a, test.b, test.min_f, test.min_x,
                                  lambda x: log_point(x, points_db[test.Index]['bnb2_pslint_points_list']))
        T1 = time.perf_counter()
        for num in range(0, repeat):
            PC_N = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=False).nsteps
        T2 = time.perf_counter()
        time_PC_N = (T2 - T1) / repeat
        it_list_row[0] = PC_N
        time_list_row[0] = time_PC_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            PI_N = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=False).nsteps
        T2 = time.perf_counter()
        time_PI_N = (T2 - T1) / repeat
        it_list_row[1] = PI_N
        time_list_row[1] = time_PI_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QC_N = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=False).nsteps
        T2 = time.perf_counter()
        time_QC_N = (T2 - T1) / repeat
        it_list_row[2] = QC_N
        time_list_row[2] = time_QC_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QI_N = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=False).nsteps
        T2 = time.perf_counter()
        time_QI_N = (T2 - T1) / repeat
        it_list_row[3] = QI_N
        time_list_row[3] = time_QI_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            PC_R = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=True).nsteps
        T2 = time.perf_counter()
        time_PC_R = (T2 - T1) / repeat
        it_list_row[4] = PC_R
        time_list_row[4] = time_PC_R

        T1 = time.perf_counter()
        for num in range(0, repeat):
            PI_R = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=True).nsteps
        T2 = time.perf_counter()
        time_PI_R = (T2 - T1) / repeat
        it_list_row[5] = PI_R
        time_list_row[5] = time_PI_R

        T1 = time.perf_counter()
        for num in range(0, repeat * 10):
            QC_R = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=True).nsteps
        T2 = time.perf_counter()
        time_QC_R = (T2 - T1) / repeat / 10
        it_list_row[6] = QC_R
        time_list_row[6] = time_QC_R

        T1 = time.perf_counter()
        for num in range(0, repeat * 10):
            QI_R = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=True).nsteps
        T2 = time.perf_counter()
        time_QI_R = (T2 - T1) / repeat / 10
        it_list_row[7] = QI_R
        time_list_row[7] = time_QI_R

        it_list.append(it_list_row)
        time_list.append(time_list_row)

        # print('%s & %d & %d & %d & %d & %d & %d & %d & %d \\\\' %(test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R))
        # print('%s & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f \\\\' % (
        # test.Index, time_PC_N, time_PI_N, time_QC_N, time_QI_N, time_PC_R, time_PI_R, time_QC_R, time_QI_R))
    index = 0
    print('test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R')
    for time_row in time_list:
        index = index + 1
        print('%s & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f \\\\' % (index,
                                                                                   time_row[0], time_row[1],
                                                                                   time_row[2], time_row[3],
                                                                                   time_row[4], time_row[5],
                                                                                   time_row[6], time_row[7]))
    min_list = getMinRatio(time_list)
    max_list = getMaxRatio(time_list)
    avg_list = getAvgRatio(time_list)
    print_row3(min_list, 'Min')
    print_row3(max_list, 'Max')
    print_row3(avg_list, 'Average')

    index = 0
    print('test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R')
    for it_row in it_list:
        index = index + 1
        print('%s & %d & %d & %d & %d & %d & %d & %d & %d \\\\' % (index,
                                                                   it_row[0], it_row[1],
                                                                   it_row[2], it_row[3],
                                                                   it_row[4], it_row[5],
                                                                   it_row[6], it_row[7]))
    min_list = getMinRatio(it_list)
    max_list = getMaxRatio(it_list)
    avg_list = getAvgRatio(it_list)
    print_row3(min_list, 'Min')
    print_row3(max_list, 'Max')
    print_row3(avg_list, 'Average')


def test_casado(df, eps, repeat):
    # print('test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R')
    it_list = []
    time_list = []

    for test in df.itertuples():
        it_list_row = [None] * 3
        time_list_row = [0.] * 3
        points_db[test.Index] = {'bnb2_pslint_points_list': []}
        prob = uvpr.UniVarProblem(test.Index, test.objective, test.a, test.b, test.min_f, test.min_x,
                                  lambda x: log_point(x, points_db[test.Index]['bnb2_pslint_points_list']))
        T1 = time.perf_counter()
        for num in range(0, repeat):
            Cas = sfzcp.cas(prob=prob, sym=False, epsilon=eps).nsteps

        T2 = time.perf_counter()
        time_Cas = (T2 - T1) / repeat
        it_list_row[0] = Cas
        time_list_row[0] = time_Cas

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QI_R = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=False, estimator=2,
                                    reduction=2).nsteps
        T2 = time.perf_counter()
        time_QI_R = (T2 - T1) / repeat
        it_list_row[1] = QI_R
        time_list_row[1] = time_QI_R

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QI_R = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=True, estimator=2,
                                    reduction=2).nsteps
        T2 = time.perf_counter()
        time_QI_R = (T2 - T1) / repeat
        it_list_row[2] = QI_R
        time_list_row[2] = time_QI_R

        it_list.append(it_list_row)
        time_list.append(time_list_row)

        # print('%s & %d & %d & %d & %d & %d & %d & %d & %d \\\\' %(test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R))
        # print('%s & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f \\\\' % (
        # test.Index, time_PC_N, time_PI_N, time_QC_N, time_QI_N, time_PC_R, time_PI_R, time_QC_R, time_QI_R))
    index = 0
    print('Index & IBB & local Lip & global Lip \\')
    for time_row in time_list:
        index = index + 1
        print('%s & %.6f & %.6f & %.6f \\\\' % (index,
                                                time_row[0], time_row[1], time_row[2]))
    min_list = getMinRatio(time_list)
    max_list = getMaxRatio(time_list)
    avg_list = getAvgRatio(time_list)
    print_row3(min_list, 'Min')
    print_row3(max_list, 'Max')
    print_row3(avg_list, 'Average')

    index = 0
    print('Index & IBB & local Lip & global Lip \\')
    for it_row in it_list:
        index = index + 1
        print('%s & %d & %d & %d \\\\' % (index,
                                          it_row[0], it_row[1], it_row[2]))
    min_list = getMinRatio(it_list)
    max_list = getMaxRatio(it_list)
    avg_list = getAvgRatio(it_list)
    print_row3(min_list, 'Min')
    print_row3(max_list, 'Max')
    print_row3(avg_list, 'Average')


def test_alg3(df, eps, global_lip, repeat):
    # print('test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R')
    it_list = []
    time_list = []

    for test in df.itertuples():
        it_list_row = [None] * 8
        time_list_row = [0.] * 8
        points_db[test.Index] = {'bnb2_pslint_points_list': []}
        prob = uvpr.UniVarProblem(test.Index, test.objective, test.a, test.b, test.min_f, test.min_x,
                                  lambda x: log_point(x, points_db[test.Index]['bnb2_pslint_points_list']))
        T1 = time.perf_counter()
        for num in range(0, repeat):
            PC_N = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=1).nsteps
        T2 = time.perf_counter()
        time_PC_N = (T2 - T1) / repeat
        it_list_row[0] = PC_N
        time_list_row[0] = time_PC_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            PI_N = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=1).nsteps
        T2 = time.perf_counter()
        time_PI_N = (T2 - T1) / repeat
        it_list_row[1] = PI_N
        time_list_row[1] = time_PI_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QC_N = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=1).nsteps
        T2 = time.perf_counter()
        time_QC_N = (T2 - T1) / repeat
        it_list_row[2] = QC_N
        time_list_row[2] = time_QC_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QI_N = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=1).nsteps
        T2 = time.perf_counter()
        time_QI_N = (T2 - T1) / repeat
        it_list_row[3] = QI_N
        time_list_row[3] = time_QI_N

        T1 = time.perf_counter()
        for num in range(0, repeat):
            PC_R = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=2).nsteps
        T2 = time.perf_counter()
        time_PC_R = (T2 - T1) / repeat
        it_list_row[4] = PC_R
        time_list_row[4] = time_PC_R

        T1 = time.perf_counter()
        for num in range(0, repeat):
            PI_R = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=1,
                                    reduction=2).nsteps
        T2 = time.perf_counter()
        time_PI_R = (T2 - T1) / repeat
        it_list_row[5] = PI_R
        time_list_row[5] = time_PI_R

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QC_R = sfzcp.new_method(prob, symm=True, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=2).nsteps
        T2 = time.perf_counter()
        time_QC_R = (T2 - T1) / repeat
        it_list_row[6] = QC_R
        time_list_row[6] = time_QC_R

        T1 = time.perf_counter()
        for num in range(0, repeat):
            QI_R = sfzcp.new_method(prob, symm=False, epsilon=eps, global_lipschitz_interval=global_lip, estimator=2,
                                    reduction=2).nsteps
        T2 = time.perf_counter()
        time_QI_R = (T2 - T1) / repeat
        it_list_row[7] = QI_R
        time_list_row[7] = time_QI_R

        it_list.append(it_list_row)
        time_list.append(time_list_row)

        # print('%s & %d & %d & %d & %d & %d & %d & %d & %d \\\\' %(test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R))
        # print('%s & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f \\\\' % (
        # test.Index, time_PC_N, time_PI_N, time_QC_N, time_QI_N, time_PC_R, time_PI_R, time_QC_R, time_QI_R))
    index = 0
    print('test.Index,PC_2,PI_2,QC_2,QI_2,PC_3,PI_3,QC_3,QI_3')
    for time_row in time_list:
        index = index + 1
        print('%s & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f & %.5f \\\\' % (index,
                                                                                   time_row[0], time_row[1],
                                                                                   time_row[2], time_row[3],
                                                                                   time_row[4], time_row[5],
                                                                                   time_row[6], time_row[7]))
    min_list = getMinRatio(time_list)
    max_list = getMaxRatio(time_list)
    avg_list = getAvgRatio(time_list)
    print_row3(min_list, 'Min')
    print_row3(max_list, 'Max')
    print_row3(avg_list, 'Average')

    index = 0
    print('test.Index,PC_N,PI_N,QC_N,QI_N,PC_R,PI_R,QC_R,QI_R')
    for it_row in it_list:
        index = index + 1
        print('%s & %d & %d & %d & %d & %d & %d & %d & %d \\\\' % (index,
                                                                   it_row[0], it_row[1],
                                                                   it_row[2], it_row[3],
                                                                   it_row[4], it_row[5],
                                                                   it_row[6], it_row[7]))
    min_list = getMinRatio(it_list)
    max_list = getMaxRatio(it_list)
    avg_list = getAvgRatio(it_list)
    print_row3(min_list, 'Min')
    print_row3(max_list, 'Max')
    print_row3(avg_list, 'Average')
