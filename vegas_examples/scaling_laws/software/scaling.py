from matplotlib import pyplot
import numpy
from scipy.optimize import curve_fit


def function(x, A1, B1, C1, A2, B2, C2):
    return A1 * numpy.exp(-B1 * (x - C1)**2) + \
           A2 * numpy.exp(-B2 * (x - C2)**2)

def main():
    lengths = [10, 20, 30, 40, 50]
    
    V1_max_mean = []
    V2_max_mean = []
    for length in lengths:
        V1_max_arr = []
        V2_max_arr = []
        Chi_max_arr = []
        for i in range(1, 6):
            T, E, Cv, M, Chi, V1, V2 = numpy.loadtxt("length_%i_sim_%i_.mean" % (length, i), unpack=True)

            params, _ = curve_fit(function, T, V1,
                p0=(
                    1.0, 1.0, T[numpy.argmax(V1)],
                    1.0, 1.0, T[numpy.argmax(V1)],
                    ), maxfev=10000000)

            xnew = numpy.linspace(min(T), max(T), 1000)
            ynew = function(xnew, *params)
            V1_max = numpy.max(ynew)

            params, _ = curve_fit(function, T, V2,
                p0=(
                    1.0, 1.0, T[numpy.argmax(V2)],
                    1.0, 1.0, T[numpy.argmax(V2)],
                    ), maxfev=10000000)

            xnew = numpy.linspace(min(T), max(T), 1000)
            ynew = function(xnew, *params)
            V2_max = numpy.max(ynew)

            V1_max_arr.append(V1_max)
            V2_max_arr.append(V2_max)

        V1_max_mean.append(numpy.mean(V1_max_arr))
        V2_max_mean.append(numpy.mean(V2_max_arr))



    X = numpy.log(lengths)
    Y1 = numpy.log(V1_max_mean)
    Y2 = numpy.log(V2_max_mean)

    m1, B1 = numpy.polyfit(X, Y1, 1)
    nu1 = 1 / m1

    m2, B2 = numpy.polyfit(X, Y2, 1)
    nu2 = 1 / m2

    pyplot.figure()
    pyplot.plot(X, Y1, "or", label=r"$\ln V_{1}$")
    pyplot.plot(X, numpy.polyval((m1, B1), X), "--r")
    pyplot.plot(X, Y2, "sb", label=r"$\ln V_{2}$")
    pyplot.plot(X, numpy.polyval((m2, B2), X), "--b")
    pyplot.grid()
    pyplot.xlabel(r"$\ln L$", fontsize=20)
    pyplot.ylabel(r"$\ln V_1, \ \ln V_2$", fontsize=20)
    pyplot.legend(loc="best")
    pyplot.savefig("Vn.png")
    pyplot.close()


    nu = (nu1 + nu2) * 0.5

    print("nu1 = ", nu1)
    print("nu2 = ", nu2)
    print("nu = ", nu)



if __name__ == '__main__':
    main()