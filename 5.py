'''5.	Выписать интерполяционные многочлены Лагранжа и Ньютона
для узловых значений {xi, yi}. Найти погрешность в точке x* = 1,3832'''


import numpy as np
import matplotlib.pyplot as plt

x_points = np.array([0.115, 0.120, 0.125, 0.130, 0.135, 0.140])
y_points = np.array([8.65729, 8.29329, 7.95829, 7.64893, 7.36235, 7.09613])

x_star = 0.1264  # точка для интерполяции


def lagrange_interpolation(x, x_points, y_points):
    """Вычисляет значение многочлена Лагранжа в точке x"""
    n = len(x_points)
    result = 0.0
    for i in range(n):
        # Базисный полином Li(x)
        li = 1.0
        for j in range(n):
            if j != i:
                li *= (x - x_points[j]) / (x_points[i] - x_points[j])
        result += y_points[i] * li
    return result


def newton_coefficients(x_points, y_points):
    #  Вычисляет коэффициенты многочлена Ньютона (разделённые разности)
    n = len(x_points)
    coef = y_points.copy()
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            coef[i] = (coef[i] - coef[i - 1]) / (x_points[i] - x_points[i - j])
    return coef


def newton_evaluation(x, x_points, coef):
    #  Вычисляет значение многочлена Ньютона в точке x по коэффициентам
    n = len(x_points)
    result = coef[0]
    product = 1.0
    for i in range(1, n):
        product *= (x - x_points[i - 1])
        result += coef[i] * product
    return result


# Лагранж
y_lagrange = lagrange_interpolation(x_star, x_points, y_points)

# Ньютон
newton_coef = newton_coefficients(x_points, y_points)
y_newton = newton_evaluation(x_star, x_points, newton_coef)


'''Погрешность (оценка через остаточный член)
Остаточный член: R_n(x) = f[x0,...,xn,x] * П(x - xi)
Для оценки используем разделённую разность на всех узлах + точка x*
Добавим x* во временный массив узлов для вычисления погрешности'''

def divided_difference_full(x_nodes, y_nodes):
    #  Вычисляет разделённые разности для массива узлов
    n = len(x_nodes)
    table = [y_nodes.copy()]
    for k in range(1, n):
        prev = table[-1]
        curr = []
        for i in range(n - k):
            denom = x_nodes[i + k] - x_nodes[i]
            if abs(denom) < 1e-15:
                curr.append(0.0)  # защита от деления на ноль
            else:
                curr.append((prev[i + 1] - prev[i]) / denom)
        table.append(curr)
    return table[-1][0] if len(table[-1]) > 0 else 0.0


# Оценим погрешность методом Рунге (через разность между интерполяциями разного порядка)
# Используем многочлены 5-й и 4-й степени
def error_runge(x, x_points, y_points):
    #  Оценка погрешности по Рунге: разность между интерполяциями с n и n-1 узлами
    # С 6 узлами (все)
    y_full = lagrange_interpolation(x, x_points, y_points)

    # С 5 узлами (без последнего)
    y_reduced = lagrange_interpolation(x, x_points[:-1], y_points[:-1])

    # Абсолютная погрешность примерно |y_full - y_reduced|
    return abs(y_full - y_reduced)


x_dense = np.linspace(0.114, 0.141, 200)
y_lagrange_dense = [lagrange_interpolation(x, x_points, y_points) for x in x_dense]
y_newton_dense = [newton_evaluation(x, x_points, newton_coef) for x in x_dense]

print("\nУзловые точки:")
for i in range(len(x_points)):
    print(f"  x_{i} = {x_points[i]:.3f}    y_{i} = {y_points[i]:.6f}")

print(f"\nТочка интерполяции x* = {x_star}")

print("\nРезультаты интерполяции")
print(f"Лагранж:  y = {y_lagrange:.10f}")
print(f"Ньютон:   y = {y_newton:.10f}")

# Проверка согласованности методов
print(f"\nСовпадение Лагранжа и Ньютона: {abs(y_lagrange - y_newton):.2e}")

# Оценка погрешности
error_estimate = error_runge(x_star, x_points, y_points)
print("\nОценка погрешности")
print(f"Метод Рунге (оценка): {error_estimate:.2e}")

# Дополнительная оценка через разделённую разность
# Построим разделённую разность 6-го порядка: f[x0,...,x5,x*]
x_extended = np.append(x_points, x_star)
# Значения в узлах + приближённое значение в x* (по интерполяции Ньютона)
y_extended = np.append(y_points, y_newton)
dd_6 = divided_difference_full(x_extended, y_extended)

# Произведение П(x*-xi)
product = np.prod(x_star - x_points)
error_dd = abs(dd_6 * product)
print(f"Через разделённую разность 6-го порядка: {error_dd:.2e}")

print(f"Приближённое значение функции в x* = {x_star}:")
print(f"y ≈ {y_newton:.8f}")
print(f"Оценка абсолютной погрешности: ~ {max(error_estimate, error_dd):.2e}")

plt.figure(figsize=(10, 6))
plt.plot(x_dense, y_lagrange_dense, 'b-', linewidth=2, label='Интерполяция Лагранжа/Ньютона')
plt.plot(x_points, y_points, 'ro', markersize=8, label='Узловые точки')
plt.plot(x_star, y_newton, 'gs', markersize=10, label=f'x* = {x_star}, y ≈ {y_newton:.6f}')
plt.title('Интерполяция многочленами Лагранжа и Ньютона')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# Вывод явного вида многочлена Ньютона (в символьной форме - приближённо)
print(f"\nКоэффициенты: {newton_coef}")
print("В виде: P(x) = c0 + c1*(x-x0) + c2*(x-x0)*(x-x1) + ...")
print(f"где c0 = {newton_coef[0]:.8f}")
print(f"c1 = {newton_coef[1]:.8f}")
print(f"c2 = {newton_coef[2]:.8f}")
print(f"c3 = {newton_coef[3]:.8f}")
print(f"c4 = {newton_coef[4]:.8f}")
print(f"c5 = {newton_coef[5]:.8f}")