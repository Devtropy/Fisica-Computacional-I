#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>

// abreviacion de llamadas de libreria
using namespace std;

double f(double x) { return 8.0 - 4.5 * (x - sin(x)); }

double fx(double x) { return -4.5 * (1 - cos(x)); }

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
      break;
    }
  }
  return {x, e, i};
}

int main() {

  auto resultado = Newton_Rhapson(f, fx, 2.0, 1e-9, 100000);

  cout << "La raiz es:" << get<0>(resultado) << endl;
  cout << "con error relativo:" << get<1>(resultado) << endl;
  cout << "en la interacion:" << get<2>(resultado) << endl;
  return 0;
}
