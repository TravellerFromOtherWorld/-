import math


def func(x):
    return 1/math.tan(x + 0.2) - x * x


def pr_func(x):
    return -1/(math.sin(x+0.2)**2) - 2 * x


def iter_fi(x):
    return (x * (1/math.tan(x + 0.2)))**(1/3)


def method_dihotomii(a, b, f, eps):
    n = 0
    x_n = 0
    a_n = a
    b_n = b
    while ((a+b)/2**n) > eps:
        x_n = (a_n + b_n)/2
        if (f(a_n) * f(x_n)) < 0:
            b_n = x_n
        else:
            a_n = x_n
        n += 1

    return x_n, n


def method_newton(x0, f, pr_f, m, eps):
    x_n = x0
    n = 0
    while (f(x_n)/pr_f(m)) > eps:
        n += 1
        x_n = x_n - (f(x_n)/pr_f(x_n))

    return x_n, n


def method_newtonm(x0, f, pr_f, m, eps):
    x_n = x0
    n = 0
    while (f(x_n)/pr_f(m)) > eps:
        n += 1
        x_n = x_n - (f(x_n)/pr_f(x0))

    return x_n, n


def method_hord(x0, x1, f, pr_f, m, eps):
    x_n = x1
    n = 0
    while abs((f(x_n)/pr_f(m))) > eps:
        n += 1
        x_n = x_n - (f(x_n)/(f(x_n) - f(x0))) * (x_n - x0)

    return x_n, n


def method_phord(x0, x1, f, pr_f, m, eps):
    x_n_1 = x0
    x_n = x1
    n = 0
    while abs((f(x_n)/pr_f(m))) > eps:
        n += 1
        x_cur = x_n
        x_n = x_n - (f(x_n)/(f(x_n) - f(x_n_1))) * (x_n - x_n_1)
        x_n_1 = x_cur

    return x_n, n


def method_pr_iter(x0, fi, q, eps):
    x_1 = fi(x0)
    n = 0
    x_n = x_1
    while (q**n/(1-q) * abs(x_1 - x0)) > eps:
        n += 1
        x_n = fi(x_n)

    return x_n, n


if __name__ == '__main__':
    ur = func
    pr_ur = pr_func
    start = 0
    end = math.pi/2 - 0.2
    e = 0.5 * 10**(-5)
    fi = iter_fi
    print('Метод дихотомии')
    print(method_dihotomii(start, end, ur, e))
    print('Метод Ньютона')
    print(method_newton(end, ur, pr_ur, end, e))
    print('Модифицированный метод Ньютона')
    print(method_newtonm(end, ur, pr_ur, end, e))
    print('Метод хорд')
    print(method_hord(end, start, ur, pr_ur, end, e))
    print('Метод подвижных хорд')
    print(method_phord(end, start, ur, pr_ur, end, e))
    print('Метод простой итерации')
    print(method_pr_iter(math.pi/12, fi, 1/3, e))
