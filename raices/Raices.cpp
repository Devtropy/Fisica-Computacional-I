#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <tuple>

using namespace std;

tuple<double, double, int> Newton_Rhapson(function<double(double)> funcion,
                                          function<double(double)> derivada,
                                          double approx0, double tolerancia,
                                          int max_iteraciones) {
  double x = approx0;
  double e = 1;
  int i = 0;
  double aprox = 0;

  while (e > tolerancia && i <= max_iteraciones) {

    x -= funcion(x) / derivada(x);
    if (funcion(x) == 0) {
      break;
    }

    else {
      e = fabs(aprox - x) / fabs(aprox);
      aprox = x;
      i++;
    }
    if (i > max_iteraciones) {
      throw runtime_error("No se encontro en" + to_string(max_iteraciones) +
                          "iteraciones");
      break;
    }
  }
  return {x, e, i};
}

double f(double x) { return x * x * x + 4 * x * x - 10.0; }

tuple<double, double, int> biseccion(function<double(double)> funcion, double a,
                                     double b, double tolerancia,
                                     int max_iteraciones) {
  // Inicializacion de variables
  double m;       // raiz
  double m0;      // raiz anterior
  int i = 0;      // iteracion
  double e = 1.0; // error relativo
  // Algoritmo de biseccion con condicionde parada de error relativo
  while (e > tolerancia && i < max_iteraciones) {
    m0 = m;
    m = (b + a) / 2.0;
    e = fabs(m - m0) / fabs(m);
    // teorema de bolzano
    if (funcion(a) * funcion(m) < 0.0) {
      b = m;
    } else if (funcion(m) * funcion(b) < 0) {
      a = m;
    } else {
      break;
    }
    i++;
  }
  return {m, e, i};
}
tuple<double, double, int>
Metodo_de_la_Secante(function<double(double)> funcion, double x0, double x1,
                     double tolerancia, int max_iteraciones) {
  double e = 100;
  int i = 0;
  double x;

  while (e > tolerancia && i < max_iteraciones) {
    double f0 = funcion(x0);
    double f1 = funcion(x1);
    x = x1 - f1 * (x1 - x0) / (f1 - f0);
    x0 = x1;
    x1 = x;
    e = fabs(x1 - x0) / fabs(x1);
    i++;
  }
  return {x1, e, i};
}
