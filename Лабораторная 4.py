import math


def expression(x):
    return math.e ** (math.sin(x))


def left_rectangle(h, a, b, f):
    k = 0
    m = (b - a) / h
    xi = a
    res = 0
    while k < m:
        res += f(xi)
        xi += h
        k += 1
    return h * res


def simpson(h, a, b, f):
    k = 0
    m = (b - a) / h
    xi = a
    res = -xi
    while k < m:
        res += 2 * f(xi)
        xi += h/2
        res += 4 * f(xi)
        xi += h / 2
        k += 1
    res += f(xi)
    return (h/6) * res


if __name__ == '__main__':
    print('Метод левых прямоугольников с h=0.1: ', left_rectangle(0.1, 2, 3, expression))
    print('Метод левых прямоугольников с h=0.05: ', left_rectangle(0.05, 2, 3, expression))
    print('Метод левых прямоугольников с h=0.025: ', left_rectangle(0.025, 2, 3, expression))

    print('Метод Симпсона с h=0.1: ', simpson(0.1, 2, 3, expression))
    print('Метод Симпсона с h=0.05: ', simpson(0.05, 2, 3, expression))
    print('Метод Симпсона с h=0.025: ', simpson(0.025, 2, 3, expression))

    print('Метод Гаусса: ', 0.444 * expression(2.5) + 0.278 * expression(2.887) + 0.278 * expression(2.113))
