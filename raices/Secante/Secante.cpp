#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>

// abreviacion de llamadas de libreria
using namespace std;

double f(double x) { return 8.0 - 4.5 * (x - sin(x)); }

double f2(double x) { return x * x * x - 3.0 * x + 2; }

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

int main() {

  auto cubica = Metodo_de_la_Secante(f2, -2.6, -2.4, 1e-9, 100000);

  cout << "La raiz es:" << get<0>(cubica) << endl;
  cout << "con error relativo:" << get<1>(cubica) << endl;
  cout << "en la interacion:" << get<2>(cubica) << endl;

  auto sinosuidal = Metodo_de_la_Secante(f, 2, 3, 1e-9, 100000);

  cout << "La raiz es:" << get<0>(sinosuidal) << endl;
  cout << "con error relativo:" << get<1>(sinosuidal) << endl;
  cout << "en la interacion:" << get<2>(sinosuidal) << endl;

  return 0;
}
