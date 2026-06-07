import numpy as np
import matplotlib.pyplot as plt

# Исходные данные
x = np.array([0.00, 0.12, 0.19, 0.35, 0.40, 0.45, 0.62, 0.71, 0.84, 0.91, 1.00])
y = np.array([1.90, 1.80, 2.20, 2.10, 2.60, 2.10, 2.70, 2.50, 2.90, 2.70, 3.30])

#  Линейная аппроксимация (степень 1)
coef_lin = np.polyfit(x, y, 1)   # возвращает [a1, a0] для P(x)=a1*x + a0
a1, a0 = coef_lin
P1 = np.poly1d(coef_lin)          # полиномиальная функция

# Сумма квадратов отклонений (остаточная сумма квадратов)
res_lin = np.sum((y - P1(x))**2)

# Квадратичная аппроксимация (степень 2)
coef_quad = np.polyfit(x, y, 2)   # возвращает [b2, b1, b0] для P(x)=b2*x^2 + b1*x + b0
b2, b1, b0 = coef_quad
P2 = np.poly1d(coef_quad)

res_quad = np.sum((y - P2(x))**2)

print("МЕТОД НАИМЕНЬШИХ КВАДРАТОВ")

print("\n1. ЛИНЕЙНАЯ АППРОКСИМАЦИЯ")
print(f"   P₁(x) = {a0:.6f} + {a1:.6f}·x")
print(f"   Сумма квадратов невязок: {res_lin:.6f}")

print("\n2. КВАДРАТИЧНАЯ АППРОКСИМАЦИЯ")
print(f"   P₂(x) = {b0:.6f} + {b1:.6f}·x + {b2:.6f}·x²")
print(f"   Сумма квадратов невязок: {res_quad:.6f}")

# ---- Построение графиков (для наглядности) ----
x_plot = np.linspace(0, 1.05, 200)
y_lin_plot = P1(x_plot)
y_quad_plot = P2(x_plot)

plt.figure(figsize=(10, 6))
plt.scatter(x, y, color='red', label='Экспериментальные точки', zorder=5)
plt.plot(x_plot, y_lin_plot, 'b--', label=f'Линейная аппроксимация: {a0:.3f}+{a1:.3f}x')
plt.plot(x_plot, y_quad_plot, 'g-', label=f'Квадратичная аппроксимация: {b0:.3f}+{b1:.3f}x+{b2:.3f}x²')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Аппроксимация методом наименьших квадратов')
plt.legend()
plt.grid(True)
plt.show()