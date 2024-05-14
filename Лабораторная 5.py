import matplotlib.pyplot as painter


def real_function(x):
    return 3 * x + 2


def f(xi, yi):
    return (yi + 1)/(xi + 1)


def Euler(xi, n):
    h = 1 / n
    results = [2]
    for x in xi:
        if x == 0:
            continue
        y = results[-1]
        results.append(y + h * f(x, y))
    return results


def Cauchy(xi, n):
    h = 1 / n
    results = [2]
    for x in xi:
        if x == 0:
            continue
        y = results[-1]
        y1 = y + h / 2 * f(x, y)
        results.append(y + h * f(x + h/2, y1))
    return results


def Adams(xi, n):
    h = 1 / n
    results = [2, 2 + h * f(0, 2) + h ** 3 / 3]
    for i in range(2, len(xi)):
        y = results[-1]
        y1 = results[-2]
        results.append((12*(xi[i]+1)) / (12*(xi[i]+1) - 5*h) * (y + 2*h/3 * f(xi[i-1], y) -
                                                                h/12*f(xi[i-2], y1) + 5*h/(12*(xi[i]+1))))
    return results


def list_of_xi(n):
    xi = [0]
    h = 1 / n
    for i in range(n):
        xi.append(xi[-1] + h)
    return xi


if __name__ == '__main__':
    xi1 = list_of_xi(10)
    xi2 = list_of_xi(20)
    xi3 = list_of_xi(30)

    painter.plot(xi1, Euler(xi1, 10), color='orange',
                 label='Метод Эйлера (n = 10)', marker='o')
    """
    painter.plot(xi2, Euler(xi2, 20), color='green',
                 label='Метод Эйлера (n = 20)', marker='o')
    
    painter.plot(xi3, Euler(xi3, 30), color='red',
                 label='Метод Эйлера (n = 30)', marker='o')
    """
    painter.plot(xi1, Cauchy(xi1, 10), color='red',
                 label='Метод Коши (n = 10)', marker='o')
    """
    painter.plot(xi2, Cauchy(xi2, 20), color='orange',
                 label='Метод Коши (n = 20)', marker='o')
    
    painter.plot(xi3, Cauchy(xi3, 30), color='pink',
                 label='Метод Коши (n = 30)', marker='o')
    """
    painter.plot(xi1, [real_function(x) for x in xi1], color='green',
                 label='Точное решение (n = 10)', marker='o')
    """
    painter.plot(xi2, [real_function(x) for x in xi2], color='blue',
                 label='Точное решение (n = 20)', marker='o')
    """
    painter.plot(xi1, Adams(xi1, 10), color='blue',
                 label='Метод Адамса-Мултона (n = 10)', marker='o')
    """
    painter.plot(xi2, Adams(xi2, 20), color='black',
                 label='Метод Адамса-Мултона (n = 20)', marker='o')

    painter.plot(xi3, Adams(xi3, 30), color='blue',
                 label='Метод Адамса-Мултона (n = 30)', marker='o')

    painter.plot(xi3, [real_function(x) for x in xi3], color='purple',
                 label='Точное решение (n = 30)', marker='o')
    """
    painter.xlabel('x')
    painter.ylabel('y')

    # painter.title('Метод Эйлера')
    # painter.title('Метод Коши')
    # painter.title('Метод Адамса-Мултона')
    painter.title('Сравнение')
    # painter.title('Точное решение')

    painter.grid(True)
    painter.legend()
    painter.show()
