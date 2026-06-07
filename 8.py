# Исходные данные (из задания №5)
x = [0.115, 0.120, 0.125, 0.130, 0.135, 0.140]
y = [8.65729, 8.29329, 7.95829, 7.64893, 7.36235, 7.09613]

n = len(x) - 1                 # количество отрезков
h = [x[i+1] - x[i] for i in range(n)]   # шаги

# ---- 1. Построение системы для вторых производных (естественные границы) ----
# Система имеет вид: h[i-1]*c[i-1] + 2*(h[i-1]+h[i])*c[i] + h[i]*c[i+1] = 3*(dy[i]/h[i] - dy[i-1]/h[i-1])
# где c[i] = S''(x_i)/2, dy[i] = y[i+1]-y[i]

# Правая часть и диагонали для прогонки (метод Томаса)
alpha = [0.0] * (n+1)    # прогоночные коэффициенты
beta  = [0.0] * (n+1)
c = [0.0] * (n+1)        # искомые c[i] = S''(x_i)/2

# Естественные граничные условия: S''(x0)=0 => c[0]=0, S''(xn)=0 => c[n]=0
c[0] = 0.0
c[n] = 0.0

# Прямой ход прогонки для внутренних узлов (i = 1 .. n-1)
# Вычисляем правые части d[i] и прогоночные коэффициенты
d = [0.0] * (n+1)
for i in range(1, n):
    dy_left  = y[i] - y[i-1]
    dy_right = y[i+1] - y[i]
    d[i] = 3.0 * (dy_right / h[i] - dy_left / h[i-1])

# Настройка прогонки для трёхдиагональной матрицы:
# A[i]*c[i-1] + B[i]*c[i] + C[i]*c[i+1] = d[i]
# где A[i] = h[i-1], B[i] = 2*(h[i-1]+h[i]), C[i] = h[i]

A = [0.0] * (n+1)
B = [0.0] * (n+1)
C = [0.0] * (n+1)

for i in range(1, n):
    A[i] = h[i-1]
    B[i] = 2.0 * (h[i-1] + h[i])
    C[i] = h[i]

# Прямой ход (приведение к верхнетреугольному виду)
alpha[1] = -C[1] / B[1]
beta[1]  =  d[1] / B[1]
for i in range(2, n):
    denom = B[i] + A[i] * alpha[i-1]
    alpha[i] = -C[i] / denom
    beta[i]  = (d[i] - A[i] * beta[i-1]) / denom

# Обратный ход
c[n-1] = beta[n-1]
for i in range(n-2, 0, -1):
    c[i] = alpha[i] * c[i+1] + beta[i]

# Теперь c[i] – известны.

# ---- 2. Вычисление коэффициентов сплайна на каждом отрезке ----
# S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3, x in [x_i, x_{i+1}]
a_coef = y[:-1]           # a_i = y_i
b_coef = [0.0] * n
d_coef = [0.0] * n

for i in range(n):
    hi = h[i]
    b_coef[i] = (y[i+1] - y[i]) / hi - hi * (2*c[i] + c[i+1]) / 3.0
    d_coef[i] = (c[i+1] - c[i]) / (3.0 * hi)

# ---- 3. Вывод сплайнов (задание №6) ----
print("Кубические сплайны дефекта 1 (естественные граничные условия S''(x1)=S''(x6)=0)")
for i in range(n):
    print(f"\nОтрезок {i+1}: x ∈ [{x[i]:.3f}, {x[i+1]:.3f}]")
    print(f"S_{i+1}(x) = {a_coef[i]:.8f} + {b_coef[i]:.8f}*(x-{x[i]:.3f}) + {c[i]:.8f}*(x-{x[i]:.3f})^2 + {d_coef[i]:.8f}*(x-{x[i]:.3f})^3")
    # Раскрытый вид (для удобства)
    xi = x[i]
    coef0 = a_coef[i] - b_coef[i]*xi + c[i]*xi**2 - d_coef[i]*xi**3
    coef1 = b_coef[i] - 2*c[i]*xi + 3*d_coef[i]*xi**2
    coef2 = c[i] - 3*d_coef[i]*xi
    coef3 = d_coef[i]
    print(f"     = {coef0:.8f} + {coef1:.8f}·x + {coef2:.8f}·x^2 + {coef3:.8f}·x^3")

# ---- 4. Вычисление производных в заданных точках (задание №8) ----
def find_segment(xx):
    """Возвращает индекс отрезка, содержащего xx, или ближайший узел"""
    if xx <= x[0]:
        return 0
    if xx >= x[-1]:
        return n-1
    for i in range(n):
        if x[i] <= xx <= x[i+1]:
            return i
    return n-1

def spline_value(xx):
    """Значение сплайна в точке xx"""
    i = find_segment(xx)
    dx = xx - x[i]
    return a_coef[i] + b_coef[i]*dx + c[i]*dx**2 + d_coef[i]*dx**3

def spline_derivative1(xx):
    """Первая производная сплайна в точке xx"""
    i = find_segment(xx)
    dx = xx - x[i]
    return b_coef[i] + 2*c[i]*dx + 3*d_coef[i]*dx**2

def spline_derivative2(xx):
    """Вторая производная сплайна в точке xx"""
    i = find_segment(xx)
    dx = xx - x[i]
    return 2*c[i] + 6*d_coef[i]*dx

# Точки для вычисления
x_star = 0.1264
x2 = 0.120

print("Задание №8: значения 1-й и 2-й производных")

for label, xx in [("x* = 0.1264", x_star), ("x₂ = 0.120", x2)]:
    y1 = spline_derivative1(xx)
    y2 = spline_derivative2(xx)
    print(f"\n{label}")
    print(f"  y'  = {y1:.8f}")
    print(f"  y'' = {y2:.8f}")

# Проверка: в узлах сплайн точно проходит через табличные значения
print("\nПроверка в узлах (значения сплайна должны совпадать с y_i):")
for i in range(len(x)):
    val = spline_value(x[i])
    print(f"x={x[i]:.3f}  y_табл={y[i]:.8f}  y_сплайн={val:.8f}  разность={val-y[i]:.2e}")