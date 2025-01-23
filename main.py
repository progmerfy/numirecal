import numpy as np
import math

# Функция для вычисления интегрального синуса с использованием ряда Тейлора
def si_taylor(x, epsilon):
    result = 0.0
    n = 0
    while True:
        term = ((-1)**n * x**(2*n + 1)) / ((2*n + 1) * math.factorial(2*n + 1))
        result += term
        if abs(term) < epsilon:
            break
        n += 1
    return result

# Функция для вычисления интерполяционного полинома Лагранжа
def lagrange_polynomial(x_values, func_values, x):
    result = 0.0
    n = len(x_values)
    for i in range(n):
        term = func_values[i]
        for j in range(n):
            if i != j:
                term *= (x - x_values[j]) / (x_values[i] - x_values[j])
        result += term
    return result

# Функция для вычисления узлов Чебышева
def chebyshev_nodes(a, b, n):
    nodes = []
    for i in range(n + 1):
        node = (a + b) / 2 + ((b - a) / 2) * math.cos((2 * i + 1) * math.pi / (2 * (n + 1)))
        nodes.append(node)
    return nodes

# Параметры задачи
a = 0.0
b = 2.0
h = 0.2  # Шаг табулирования
epsilon_si = 1e-6
points = [0.8, 2.0, 3.2]
# Табулирование Si(x) на интервале [0.0, 2.0] с шагом 0.2
x_values_tab = np.arange(a, b + h, h)
si_values_tab = [si_taylor(x, epsilon_si) for x in x_values_tab]

# Построение интерполяционного полинома Лагранжа по равномерным узлам
n = 6  # Количество узлов для интерполяции
h_uniform = (b - a) / n
x_values_uniform = np.arange(a, b + h_uniform, h_uniform)
si_values_uniform = [si_taylor(x, epsilon_si) for x in x_values_uniform]

# Вывод таблицы табулирования и интерполяции
print("Zi\tLn(Zi)\tS(Zi)\t|Ln(Zi)-S(Zi)|")
for i in range(len(x_values_tab)):
    zi = x_values_tab[i]
    sZi = si_values_tab[i]
    lnZi = lagrange_polynomial(x_values_uniform, si_values_uniform, zi)
    error = abs(lnZi - sZi)
    print(f"{zi:.1f}\t{lnZi:.5f}\t{sZi:.5f}\t{error:.7f}")

# Вывод таблицы с погрешностями для равномерных узлов и узлов Чебышева
print("\nn\tПогрешность (Лагранж)\tПогрешность (Чебышев)")
for n in range(6, 31):  # Увеличение узлов с 6 до 30
    # Равномерные узлы
    h_uniform = (b - a) / n
    x_values_uniform = np.arange(a, b + h_uniform, h_uniform)
    si_values_uniform = [si_taylor(x, epsilon_si) for x in x_values_uniform]
    
    # Узлы Чебышева
    x_values_cheb = chebyshev_nodes(a, b, n)
    si_values_cheb = [si_taylor(x, epsilon_si) for x in x_values_cheb]
    
    # Вычисление максимальной погрешности для равномерных узлов
    max_error_uniform = 0.0
    for i in range(len(x_values_tab)):
        zi = x_values_tab[i]
        sZi = si_values_tab[i]
        lnZi_uniform = lagrange_polynomial(x_values_uniform, si_values_uniform, zi)
        error_uniform = abs(lnZi_uniform - sZi)
        max_error_uniform = max(max_error_uniform, error_uniform)
    
    # Вычисление максимальной погрешности для узлов Чебышева
    max_error_cheb = 0.0
    for i in range(len(x_values_tab)):
        zi = x_values_tab[i]
        sZi = si_values_tab[i]
        lnZi_cheb = lagrange_polynomial(x_values_cheb, si_values_cheb, zi)
        error_cheb = abs(lnZi_cheb - sZi)
        max_error_cheb = max(max_error_cheb, error_cheb)
    
    # Вывод результатов для каждого n
    print(f"{n}\t{max_error_uniform:.10f}\t\t{max_error_cheb:.10f}")