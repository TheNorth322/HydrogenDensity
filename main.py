import numpy
import matplotlib.pyplot as plt
import scipy.special
import scipy.misc
import math
from mpl_toolkits import mplot3d

# Переход к сферическим координатам (r, theta, phi)
def r(x, y, z):
    return numpy.sqrt(x**2 + y**2 + z**2)

def theta(x, y, z):
    R = r(x,y,z)
    if R == 0:
        return 0
    elif z < -R or z > R:
        return numpy.nan
    else:
        return numpy.arccos(z / R)


def phi(x, y, z):
    if x == 0 and y == 0:
        return 0
    elif x == 0:
        if y > 0:
            return numpy.pi / 2
        else:
            return -numpy.pi / 2
    else:
        return numpy.arctan(y / x)


# Коэффициент нормализации
def normalization_factor(n, l, a0):
    return numpy.sqrt((2 / n / a0) ** 3 * math.factorial(n - l - 1) / (2 * n * math.factorial(n + l)))

# Радиальная функция
def R(r, n, l, a0):
   return normalization_factor(n, l, a0) * (2 * r / n / a0) ** l * numpy.exp(-r / n / a0) * scipy.special.genlaguerre(
       n - l - 1, 2 * l + 1)(2 * r / n / a0)

# Волновая функция
def wave_function(r, theta, phi, n, l, m, a0):
   return R(r, n, l, a0) * scipy.special.sph_harm(m, l, phi, theta)

# Плотность вероятности волновой функции
def absolute_wave_function(r, theta, phi, n, l, m, a0):
   return numpy.absolute(wave_function(r, theta, phi, n, l, m, a0)) ** 2

def w(x, y, z, n, l, m):
    a0 = 1
    return absolute_wave_function(r(x, y, z), theta(x, y, z), phi(x, y, z), n, l, m, a0)

# Визуализация данных
def visualize(elements, probability):
    probability = probability / sum(probability)
    coord = numpy.random.choice(elements, size=100000, replace=True, p=probability)
    elem_mat = [i.split(',') for i in coord]
    elem_mat = numpy.matrix(elem_mat)
    x_coords = [float(i.item()[1:]) for i in elem_mat[:, 0]]
    y_coords = [float(i.item()) for i in elem_mat[:, 1]]
    z_coords = [float(i.item()[0:-1]) for i in elem_mat[:, 2]]

    fig = plt.figure(figsize=(10, 10))
    fig.canvas.mpl_connect('close_event', on_close)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x_coords, y_coords, z_coords, alpha=0.05, s=2)
    ax.set_title("Электронная плотность атома водорода")
    plt.show()

# Метод валидации данных
def validate_data(n, l, m):
    if (n < -3 or n > 3):
        raise Exception("Неверное значение n!\n")
    elif (l < -3 or l >= n):
        raise Exception("Неверное значение l!\n")
    elif (m < -3 or m > l):
        raise Exception("Неверное значение m!\n")

# Метод расчета
def calculate(xs, ys, zs, n, l, m, elements, probability):
    for x in xs:
        for y in ys:
            for z in zs:
                elements.append(str((x, y, z)))
                probability.append(w(x, y, z, n, l, m))

def on_close(event):
    print("")

def main():
    # n - главное квантовое число,
    # l - квантовое число углового момента,
    # m - магнитное квантовое число
    while (True):
        n = input("Введите значение n: ")
        if (n.find("-") <= 1 and n.replace("-", "").isdigit()):
            n = int(n)
        else:
            continue

        l = input("Введите значение l: ")
        if (l.find("-") <= 1 and l.replace("-", "").isdigit()):
            l = int(l)
        else:
            continue

        m = input("Введите значение m: ")
        if (m.find("-") <= 1 and m.replace("-", "").isdigit()):
            m = int(m)
        else:
            continue

        step = 1
        xs = numpy.arange(-30, 30, step)
        ys = numpy.arange(-30, 30, step)
        zs = numpy.arange(-30, 30, step)

        elements = []
        probability = []

        try:
            validate_data(n, l, m)
        except Exception as e:
            print(e)
            continue

        print("Ожидайте построения графика электронной плотности атома водорода")

        calculate(xs, ys, zs, n, l, m, elements, probability)
        visualize(elements, probability)

main()