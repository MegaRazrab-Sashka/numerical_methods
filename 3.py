"""3.	Методом вращения с точностью вычислить собственные значения и
собственные вектора симметричной матрицы"""

import numpy as np

A = np.array([[56, 8, 3], [8, 24, 8], [3, 8, 88]])
eps = 0.01

def method_vrascheniya(A, eps):
    n = len(A)
    V = np.eye(n)  # будущая матрица собственных векторов
    A_cur = A.copy().astype(float)
    max_iter = 1000

    for it in range(max_iter):
        # Поиск максимального внедиагонального элемента
        max_el = 0
        p, q = 0, 1
        for i in range(n):
            for j in range(i+1, n):
                if abs(A_cur[i][j]) > max_el:
                    max_el = abs(A_cur[i][j])
                    p, q = i, j

        # Проверка сходимости
        if max_el < eps:
            break

        # Вычисление угла поворота
        if abs(A_cur[p][p] - A_cur[q][q]) < 1e-10:
            theta = np.pi / 4
        else:
            theta = 0.5 * np.arctan2(2 * A_cur[p][q], A_cur[q][q] - A_cur[p][p])
        c = np.cos(theta)
        s = np.sin(theta)
        R = np.eye(n)
        R[p][p] = c
        R[q][q] = c
        R[p][q] = s
        R[q][p] = -s
        A_cur = R.T @ A_cur @ R
        V = V @ R
    lam = np.array([A_cur[i][i] for i in range(n)])
    return lam, V, it+1

lam, V, iter_count = method_vrascheniya(A, eps)

print(f'Симметричная матрица A:\n{A}\n')
print(f'Точность: {eps}\n')
print(f'Кол-во итераций: {iter_count}\n')
print('Собственные значения:')
for i, val in enumerate(lam):
    print(f'  λ{i+1} = {val:.6f}')

print('\nСобственные вектора (столбцы матрицы V):')
for i in range(len(A)):
    print(f'\nВектор v{i+1}:')
    for j in range(len(A)):
        print(f'  {V[j][i]:.6f}')
