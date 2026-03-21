#include <cmath>
#include <functional>
#include <iostream>
#include <tuple>

// abreviacion de llamadas de libreria
using namespace std;

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

int main() {
  auto resultado = biseccion(f, 1, 2, 1e-5, 10000);
  cout << "La raiz es:" << get<0>(resultado) << endl;
  cout << "con error relativo:" << get<1>(resultado) << endl;
  cout << "en la interacion:" << get<2>(resultado) << endl;
  return 0;
}
