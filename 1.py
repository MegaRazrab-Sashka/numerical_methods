'''1.	Решить СЛАУ с точностью: А) методом Гаусса,
 Б) методом простой итерации'''

import numpy as np

s1 = 6.42 - 1.13 + 1.05 + 6.15
s2 = 1.14 - 5.15 + 2.11 - 4.16
s3 = -0.71 + 0.81 - 3.02 - 0.17
A = np.array([
    [6.42, -1.13, 1.05, 6.15, s1],
    [1.14, -5.15, 2.11, -4.16, s2],
    [-0.71, 0.81, -3.02, -0.17, s3]
], dtype=float)
A1 = A.copy()
print(f'Матрица A: \n{A}\n')

print('Метод Гаусса:')
n, m = A.shape  # n - число строк (у меня 3), m - столбцов (5)
print('Прямой ход:')
for i in range(1, n):
    k = -(A[i][0] / A[0][0])
    A[i] = A[i] + A[0] * k
    print(A)
    print()
k = - (A[2][1] / A[1][1])
A[2] = A[2] + A[1] * k
print(A)
print()

print('Обратный ход:')
ans = []
x3 = A[2][3] / A[2][2]
print(f'x3 = {round(x3, 2)}')
ans.append(round(x3, 2))
x2 = (A[1][3] - A[1][2] * x3) / A[1][1]
print(f'x2 = {round(x2, 2)}')
ans.append(round(x2, 2))
x1 = (A[0][3] - A[0][2] * x3 - A[0][1] * x2) / A[0][0]
print(f'x1 = {round(x1, 2)}')
ans.append(round(x1, 2))

print('Метод простой итерации:')
A = np.array([
    [6.42, -1.13, 1.05],
    [1.14, -5.15, 2.11],
    [-0.71, 0.81, -3.02]
], dtype=float)
b = np.array([6.15, -4.16, -0.17], dtype=float) # столбец правых частей

C = np.zeros_like(A)  # матрица (то, что B в уравнении x = Bx + C)
d = np.zeros_like(b)  # новый столбец свободных членов (C в x = Bx + C)

for i in range(3):
    for j in range(3):
        if i != j:
            C[i,j] = -A[i,j] / A[i,i]
    d[i] = b[i] / A[i,i]

check = 1

# Проверка сходимости (А то не сойдётся, а мы тут пыхтим)
norm_inf = np.max(np.sum(np.abs(C), axis=1))   # ||C||_inf (по строкам)
norm_1 = np.max(np.sum(np.abs(C), axis=0))     # ||C||_1 (по столбцам)
if norm_inf < 1 or norm_1 < 1:
    print(f'Условие сходимости выполнено: ||C||_inf = {norm_inf:.4f}, ||C||_1 = {norm_1:.4f}')
else:
    # Проверим спектральный радиус на всякий случай
    spectral_radius = max(abs(np.linalg.eigvals(C)))
    if spectral_radius < 1:
        print(f'Нормы матрицы >= 1, но спектральный радиус {spectral_radius:.4f} < 1. Сходимость возможна.')
    else:
        print(f'Метод расходится! Спектральный радиус = {spectral_radius:.4f}. Прерывание.')
        check = 0

if check == 1:
    x = np.zeros(3)
    x_new = np.zeros(3)
    eps = 0.01

    for _ in range(100):
        x_new = C @ x + d
        if np.linalg.norm(x_new - x, np.inf) < eps:
            break
        x = x_new.copy()

    print(x_new)

    res = b - A @ x_new
    print(f'Невязка: {res}')
    print(f'Норма невязки: {np.linalg.norm(res)}')