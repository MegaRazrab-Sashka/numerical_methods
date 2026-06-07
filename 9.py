import numpy as np

def f(x):
    """Подынтегральная функция"""
    return 1.0 / np.sqrt(2.0 * x * x + 1.0)

def simpson(a, b, n):
    """
    Метод Симпсона для интеграла от a до b с n подотрезками (n должно быть чётным).
    """
    if n % 2 != 0:
        raise ValueError("n должно быть чётным")
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    # Сумма: f0 + fn + 4*(нечётные) + 2*(чётные)
    s = y[0] + y[-1]
    s += 4.0 * np.sum(y[1:-1:2])   # нечётные индексы (первый после 0)
    s += 2.0 * np.sum(y[2:-2:2])   # чётные индексы (второй после 0)
    return s * h / 3.0

# Параметры
a, b = 0.8, 1.6
eps = 1e-4

# Начальное число подотрезков (чётное)
n = 4
I_prev = simpson(a, b, n)
n *= 2
I_curr = simpson(a, b, n)

while abs(I_curr - I_prev) > eps:
    I_prev = I_curr
    n *= 2
    I_curr = simpson(a, b, n)

print(f"Число подотрезков n = {n}")
print(f"Приближённое значение интеграла: {I_curr:.8f}")
print(f"Оценка погрешности: {abs(I_curr - I_prev):.2e}")

# Для справки: точное значение (если нужно)
# Вычислим с очень малым шагом
I_exact = simpson(a, b, 10**5)
print(f"Справочное значение (очень точное): {I_exact:.8f}")