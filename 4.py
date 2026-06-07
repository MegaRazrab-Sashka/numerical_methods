'''4.	С точностью уточнить один из корней уравнения
А) методом простой итерации
Б) методом половинного деления (дихотомии)
В) Методом Ньютона'''


import math

# Исходная функция и производная для метода Ньютона
def f(x):
    return x * (2 ** x) - 1

def df(x):
    return 2 ** x + x * (2 ** x) * math.log(2)

# Функция для метода простой итерации: x = 1 / (2^x)
def phi(x):
    return 1 / (2 ** x)

eps = 1e-6
max_iter = 1000

def simple_iteration(x0, eps, max_iter):
    x_prev = x0
    for i in range(max_iter):
        x_next = phi(x_prev)
        if abs(x_next - x_prev) < eps:
            return x_next, i+1
        x_prev = x_next
    return x_prev, max_iter

def dikh(a, b, eps, max_iter):
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) и f(b) должны иметь разные знаки")
    for i in range(max_iter):
        c = (a + b) / 2
        if abs(f(c)) < eps or (b - a) / 2 < eps:
            return c, i+1
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2, max_iter

def newton(x0, eps, max_iter):
    x = x0
    for i in range(max_iter):
        fx = f(x)
        if abs(fx) < eps:
            return x, i+1
        dfx = df(x)
        if dfx == 0:
            raise ValueError("Производная равна нулю")
        x_new = x - fx / dfx
        if abs(x_new - x) < eps:
            return x_new, i+1
        x = x_new
    return x, max_iter

print("Решение уравнения x * 2^x = 1\n")

x00 = 0.5      # для простой итерации
a, b = 0, 1          # интервал для дихотомии
x0n = 0.6      # начальное приближение для Ньютона

root_si, iter_si = simple_iteration(x00, eps, max_iter)
print(f"А) Метод простой итерации: x = {root_si:.10f}, итераций: {iter_si}")

root_bs, iter_bs = dikh(a, b, eps, max_iter)
print(f"Б) Метод половинного деления: x = {root_bs:.10f}, итераций: {iter_bs}")

root_nw, iter_nw = newton(x0n, eps, max_iter)
print(f"В) Метод Ньютона: x = {root_nw:.10f}, итераций: {iter_nw}")