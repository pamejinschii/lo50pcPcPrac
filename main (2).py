import numpy as np

# Данные
xi = np.array([0.1 * np.pi, 0.2 * np.pi, 0.3 * np.pi, 0.4 * np.pi])
yi = np.cos(xi)

# Интерполяция методом Ньютона
def newton_interpolation(xi, yi, x):
    n = len(xi)
    coef = np.zeros([n, n])
    coef[:,0] = yi
    
    for j in range(1, n):
        for i in range(n-j):
            coef[i,j] = (coef[i+1, j-1] - coef[i, j-1]) / (xi[i+j] - xi[i])
    
    x_term = 1
    y_interp = coef[0,0]
    
    for j in range(1, n):
        x_term *= (x - xi[j-1])
        y_interp += coef[0,j] * x_term
    
    return y_interp

# Вычисление значения функции в точке 0.25π
x_star = 0.25 * np.pi
y_star = newton_interpolation(xi, yi, x_star)
y_exact = np.cos(x_star)

# Вычисление погрешности
error = np.abs(y_star - y_exact)

print(f"Интерполяционное значение в точке {x_star:.5f}: {y_star:.5f}")
print(f"Точное значение в точке {x_star:.5f}: {y_exact:.5f}")
print(f"Погрешность: {error:.5e}")
