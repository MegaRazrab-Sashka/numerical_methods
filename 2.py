'''2.	Решить СЛАУ с точностью  : А) методом прогонки,
Б) Методом Зейделя'''


import numpy as np

M = np.array([[7, -3, 0, 0, 12.4], [-2, 7, -3, 0, 6.7],
              [0, 3, 11, -5, 16.2], [0, 0, -2, 7, 10.5]])
k = len(M)  # для скорости программы, чтоб каждый раз его не вычислять
print(f'Исходная матрица: \n{M}')

# Метод прогонки:
A = [-M[0][1] / M[0][0]]
B = [M[0][-1] / M[0][0]]

for i in range(1, k):
    ai = M[i][i - 1]
    bi = M[i][i]
    ci = M[i][i + 1] if i < k - 1 else 0
    if i == k - 1:
        A.append(0)
    else:
        A.append(-ci / (bi + ai * A[-1]))
    B.append((M[i][-1] - ai * B[-1]) / (bi + ai * A[i - 1]))

x = [0] * k  # Список ответов
x[-1] = B[-1]
for i in range(k - 2, -1, -1):
    x[i] = A[i] * x[i + 1] + B[i]
x = [float(round(i, 2)) for i in x]
print(f'Метод прогонки: {x}')

# Метод Зейделя:
eps = 0.0001
max_it = 1000

b = [M[i][-1] for i in range(k)]
A_mat = [M[i][:-1] for i in range(k)]

x_new = [0] * k
x_old = [0] * k

for it in range(max_it):
    for i in range(k):
        sum1 = sum(A_mat[i][j] * x_new[j] for j in range(i))  # уже обновлённые
        sum2 = sum(A_mat[i][j] * x_old[j] for j in range(i + 1, k))  # старые
        x_new[i] = (b[i] - sum1 - sum2) / A_mat[i][i]
    err = max(abs(x_new[i] - x_old[i]) for i in range(k))
    if err < eps:
        break
    x_old = x_new.copy()

x_new = [float(i) for i in x_new]
print(f'Метод Зейделя: {x_new}')