import numpy as np
import math
import matplotlib.pyplot as painter
# константы
N = 20

h = 1.0 / N

x = [h * i for i in range(N + 1)]


def real_function(x):
    return math.exp(x) + math.exp(-x) + 2.6 * (x ** 2) - 2.6 * x - 2


def f(x, y):
    return y - 2.6 * (x ** 2) + 2.6 * x + 7.2

# метод для разгона
def runge_kutta(y0, u0):
    y = [y0] * (N + 1)
    u = [u0] * (N + 1)

    for i in range(2):
        K1 = h * u[i]
        L1 = h * f(x[i], y[i])

        K2 = h * (u[i] + L1/2)
        L2 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * K1)

        K3 = h * (u[i] + L2/2)
        L3 = h * f(x[i] + 0.5 * h, y[i] + 0.5 * K2)

        K4 = h * (u[i] + L3)
        L4 = h * f(x[i] + h, y[i] + K3)

        y[i + 1] = y[i] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
        u[i + 1] = u[i] + (L1 + 2 * L2 + 2 * L3 + L4) / 6
    return y, u

# метод для решения задачи Коши
def Adams(y0, u0):
    y, u = runge_kutta(y0, u0) # разгоняем 2 значения

    for i in range(2, N):
        y[i + 1] = y[i] + h / 12 * (23 * u[i] - 16 * u[i-1] + 5 * u[i-2])
        u[i + 1] = u[i] + h / 12 * (23 * f(x[i], y[i]) - 16 * f(x[i - 1], y[i - 1])
                                    + 5 * f(x[i - 2], y[i - 2]))
    return y, u


def phi(mu):
    y0 = 0
    u0 = mu
    y, u = Adams(y0, u0)
    return y[-1] - math.exp(1) - math.exp(-1) + 2

# ищем mu
def method_hord(mu1, mu2):
    mu = mu1
    while np.abs(phi(mu)) > 0.00001:
        mu = mu - phi(mu) / (phi(mu) - phi(mu2)) * (mu - mu2)
    return mu

# ищем интервал, на котором потом будем искать mu
def find_start_mu(start=-100):
    while phi(start) < 0:
        start += 1
    return start - 1, start


def progonka():
    lambdas = [0.0] * (N + 1)
    mus = [0.0] * (N + 1)

    for i in range(1, N + 1):
        lambdas[i] = 1 / (h * h + 2 - lambdas[i-1]) # внесли минус в знаменатель
        mus[i] = (mus[i-1] - h * h * (2.6 * x[i] * (1 - x[i]) + 7.2)) / (h * h + 2 - lambdas[i-1]) # тут тоже

    lambdas[-1] = 0
    mus[-1] = math.e + 1 / math.e - 2

    y = [0.0] * (N + 1)
    y_temp = mus[-1]
    for i in range(N + 1, 0, -1):
        y_temp = lambdas[i - 1] * y_temp + mus[i - 1]
        y[i - 1] = y_temp
    return y


if __name__ == '__main__':
    mu1, mu2 = find_start_mu()
    mu = method_hord(mu1, mu2)
    y0 = 0
    u0 = mu
    y1, u = Adams(y0, u0)
    y2 = progonka()

    painter.title(f'N={N}')
    painter.plot(x, y2, label=f"Прогонка", color='blue', marker='o')
    painter.plot(x, y1, label=f"Стрельба, mu={np.round(mu, 6)}", color='green', marker='o')
    painter.plot(x, [real_function(x_i) for x_i in x], label="Точное решение", color='red', marker='o')
    painter.grid(True)
    painter.legend()
    painter.show()
