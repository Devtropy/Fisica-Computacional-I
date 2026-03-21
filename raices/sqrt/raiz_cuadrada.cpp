#include <cmath>
#include <functional>
#include <iostream>
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

tuple<double, double, int> Raiz_cuadrada(double p) {
  if (p < 0) {
    throw runtime_error("El numero a calcular la raiz no puede ser negativo");
  } else {
    auto raiz = Newton_Rhapson([p](double x) { return x * x - p; },
                               [](double x) { return 2.0 * x; }, p, 1e-6, 20);
    return raiz;
  }
}

int main() {

  auto raiz1 = Raiz_cuadrada(729.0);
  auto raiz2 = Raiz_cuadrada(1500.0);
  cout << "--------------Raiz 1 ---------------" << endl;
  cout << "la raiz de 729 es:" << get<0>(raiz1) << endl;
  cout << "con error relativo:" << get<1>(raiz1) << endl;
  cout << "en la iteracion:" << get<2>(raiz1) << endl;
  cout << "--------------Raiz 2 ---------------" << endl;
  cout << "la raiz de 1500 es:" << get<0>(raiz2) << endl;
  cout << "con error relativo:" << get<1>(raiz2) << endl;
  cout << "en la iteracion:" << get<2>(raiz2) << endl;
  cout << "--------------Raiz 3 ---------------" << endl;
  auto raiz3 = Raiz_cuadrada(-72.0);

  return 0;
}
