import matplotlib.pyplot as plt
import csv
import numpy
import scipy
from scipy import optimize, stats
from mpl_toolkits.mplot3d import axes3d, Axes3D

# Ten skrypt ma za zadanie znalezienia płaszczyzny, która najlepiej (w sensie błędu średniokwadratowego) przybliża dane
# opisujące ekspresję cząstek adhezyjnych. Przed uruchomieniem proszę pamiętać o zmianie ścieżek.

x = []
y = []
z = []
data = []

#Funkcje dla każdego typu cząstki.

def fit_PSelectin():

    with open('PSelectin_Fold1.csv', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    with open('PSelectin_Fold2.csv', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    guess0 = [0.8, 0.0, 0.8, 5.0, 2.8, 0.8, 22.0, 5.0, 0.2, 0.55, 3.5, 18.4]

    result = scipy.optimize.minimize(fun=residual_PSelectin, x0=guess0, args=data)
    coeff = result.x

    extend = 40.0
    x_1 = numpy.arange(min(x), max(x) + extend, 0.05)
    y_1 = numpy.arange(min(y), max(y), 0.05)
    X, Y = numpy.meshgrid(x_1, y_1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, fun_PSelectin(X, Y, coeff) + 1.0)

    new_z = [a + 1 for a in z]

    ax.scatter(x, y, new_z, marker='o', color='r')

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Dose (Gy)')
    ax.set_zlabel('Fold change')
    plt.gca().invert_xaxis()

    numpy.savetxt("PSelectincoeff.csv", coeff, delimiter=",")

    plt.title('P-Selectins expression change due to irradiation')
    plt.grid(which='minor', alpha=0.9)
    plt.grid(which='major', alpha=0.9)

    plt.savefig('PSelectin_fit.pdf')
    plt.show()

def fit_ESelectin():

    with open('ESelectin_Fold1.csv','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    with open('ESelectin_Fold2.csv', 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    guess0 = [0.6, -0.198, 3.38, 4.5, 14.4, 10.0, 0.5, 23.5, 1.05, 0.9, -0.3, 9.0, 55.0]

    result = scipy.optimize.minimize(fun=residual_ESelectin, x0=guess0, args=data)
    coeff = result.x

    extend = 40.0
    x_1 = numpy.arange(min(x), max(x) + extend, 0.05)
    y_1 = numpy.arange(min(y), max(y), 0.05)
    X, Y = numpy.meshgrid(x_1, y_1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, fun_ESelectin(X, Y, coeff) + 1.0)

    new_z = [a + 1 for a in z]

    ax.scatter(x, y, new_z, marker='o', color='r')

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Dose (Gy)')
    ax.set_zlabel('Fold change')
    plt.gca().invert_xaxis()

    numpy.savetxt("ESelectincoeff.csv", coeff, delimiter=",")

    plt.title('E-Selectins expression change due to irradiation')
    plt.grid(which='minor', alpha=0.9)
    plt.grid(which='major', alpha=0.9)

    plt.savefig('ESelectin_fit.pdf')
    plt.show()

def fit_VIntegrin():

    with open('VCAM_Fold1.csv','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    with open('VCAM_Fold2.csv','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    guess0 = [0.92, 0.035, 22.57, 205.4, 0.0, 450.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0]

    result = scipy.optimize.minimize(fun=residual_VIntegrin, x0=guess0, args=data)
    coeff = result.x

    extend = 40.0
    x_1 = numpy.arange(min(x), max(x) + extend, 0.05)
    y_1 = numpy.arange(min(y), max(y), 0.05)
    X, Y = numpy.meshgrid(x_1, y_1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, fun_VIntegrin(X, Y, coeff) + 1.0)

    new_z = [a + 1 for a in z]

    ax.scatter(x, y, new_z, marker='o', color='r')

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Dose (Gy)')
    ax.set_zlabel('Fold change')
    plt.gca().invert_xaxis()

    numpy.savetxt("VIntegrincoeff.csv", coeff, delimiter=",")

    plt.title('V-Integrins expression change due to irradiation')
    plt.grid(which='minor', alpha=0.9)
    plt.grid(which='major', alpha=0.9)

    plt.savefig('VIntegrin_fit.pdf')
    plt.show()

def fit_IIntegrin():

    with open('ICAM_Fold1.csv','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    with open('ICAM_Fold2.csv','r') as csvfile:
        plots = csv.reader(csvfile, delimiter=';')
        for row in plots:
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]) - 1.0)
            data.append((float(row[0]), float(row[1]), float(row[2]) - 1.0))

    guess0 = [1.0, -0.7, 18.5, 45.0, 0.0, 44.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0]

    result = scipy.optimize.minimize(fun=residual_IIntegrin, x0=guess0, args=data)
    coeff = result.x

    extend = 40.0
    x_1 = numpy.arange(min(x), max(x) + extend, 0.05)
    y_1 = numpy.arange(min(y), max(y), 0.05)
    X, Y = numpy.meshgrid(x_1, y_1)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot_surface(X, Y, fun_IIntegrin(X, Y, coeff) + 1.0)

    new_z = [a + 1 for a in z]

    ax.scatter(x, y, new_z, marker='o', color='r')

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('Dose (Gy)')
    ax.set_zlabel('Fold change')
    plt.gca().invert_xaxis()

    numpy.savetxt("IIntegrincoeff.csv", coeff, delimiter=",")

    plt.title('I-Integrins expression change due to irradiation')
    plt.grid(which='minor', alpha=0.9)
    plt.grid(which='major', alpha=0.9)

    plt.savefig('IIntegrin_fit.pdf')
    plt.show()

# Funkcje opsiujące szukaną płaszczyznę.

def fun_ESelectin(x, y, guess):

    a1 = guess[0]
    a2 = guess[1]
    a3 = guess[2]
    a4 = guess[3]
    a5 = guess[4]
    a6 = guess[5]
    a7 = guess[6]
    a8 = guess[7]
    a9 = guess[8]

    # Ustalone ze względu na wysokie wartości.
    a10 = 0.8

    a11 = guess[10]
    a12 = guess[11]
    a13 = guess[12]

    return (scipy.stats.lognorm.pdf(x+0.01,a1,a2*y+a3,a4) )  \
            * ( a5 * y + a6 * numpy.sqrt(y) ) \
            +(scipy.stats.lognorm.pdf(x+0.01, a7, a8, a9)) * (scipy.stats.lognorm.pdf(y,a10,a11, a12)) * a13

def fun_PSelectin(x, y, guess):

    #Parametry optymalizowane

    a1 = guess[0]
    a2 = guess[1]
    a3 = guess[2]
    a4 = guess[3]

    # Sztywne ustalenie pierwszego rozkładu z powodu zbyt wysokich wartości.
    a5 = 0.2
    a6 = 0.55
    a7 = 3.5

    a8 = guess[4]
    a9 = guess[5]
    a10 = guess[6]
    a11 = guess[7]
    a12 = guess[8]
    a13 = guess[9]
    a14 = guess[10]
    a15 = guess[11]

    return (scipy.stats.lognorm.pdf(x+0.01,a1,a2*y+a3,a4) )  \
            * ( scipy.stats.lognorm.pdf(y,a5,a6, a7)* a8 ) \
            +(scipy.stats.lognorm.pdf(x+0.01, a9, a10, a11)) * (scipy.stats.lognorm.pdf(y,a12,a13, a14)) * a15

def fun_IIntegrin(x, y, guess):

    a1 = guess[0]
    a2 = guess[1]
    a3 = guess[2]
    a4 = guess[3]
    a5 = guess[4]
    a6 = guess[5]
    a7 = guess[6]
    a8 = guess[7]
    a9 = guess[8]
    a10 = guess[9]
    a11 = guess[10]
    a12 = guess[11]
    a13 = guess[12]

    return (scipy.stats.lognorm.pdf(x+0.01,a1,a2*y+a3,a4) )  \
            * ( a5 * y + a6 * numpy.sqrt(y) ) \
            +(scipy.stats.lognorm.pdf(x+0.01, a7, a8, a9)) * (scipy.stats.lognorm.pdf(y,a10,a11, a12)) * a13

def fun_VIntegrin(x, y, guess):

    #Ustalone ze względu na zbyt wysokie wartości.
    a5 = 0.0
    a6 = 350.0

    a1 = guess[0]
    a2 = guess[1]
    a3 = guess[2]
    a4 = guess[3]
    a7 = guess[6]
    a8 = guess[7]
    a9 = guess[8]
    a10 = guess[9]
    a11 = guess[10]
    a12 = guess[11]
    a13 = guess[12]

    return (scipy.stats.lognorm.pdf(x+0.01,a1,a2*y+a3,a4) )  \
            * ( a5 * y + a6 * numpy.sqrt(y)) \
            +(scipy.stats.lognorm.pdf(x+0.01, a7, a8, a9)) * (scipy.stats.lognorm.pdf(y,a10,a11, a12)) * a13

#Funkcje opisujące residuum

def residual_PSelectin(params, points):
  residuals = [
    p[2] - fun_PSelectin(p[0], p[1], params) for p in points]

  return numpy.linalg.norm(residuals)

def residual_ESelectin(params, points):
  residuals = [
    p[2] - fun_ESelectin(p[0], p[1], params) for p in points]

  return numpy.linalg.norm(residuals)

def residual_VIntegrin(params, points):
  residuals = [
    p[2] - fun_VIntegrin(p[0], p[1], params) for p in points]

  return numpy.linalg.norm(residuals)

def residual_IIntegrin(params, points):
  residuals = [
    p[2] - fun_IIntegrin(p[0], p[1], params) for p in points]

  return numpy.linalg.norm(residuals)

def main():
    fit_ESelectin()
    fit_PSelectin()
    fit_IIntegrin()
    fit_VIntegrin()

if __name__ == "__main__":
    main()