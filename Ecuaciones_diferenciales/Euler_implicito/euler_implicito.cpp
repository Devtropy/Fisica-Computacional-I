#include <cmath>
#include <fstream>
#include <functional>
#include <tuple>
#include <vector>

using namespace std;
template <typename Tipo>
tuple<Tipo, Tipo, int> Metodo_de_la_Secante(function<Tipo(Tipo)> funcion,
                                            Tipo x0, Tipo x1, Tipo tolerancia,
                                            int max_iteraciones) {
  Tipo e = 100;
  int i = 0;
  Tipo x;

  while (e > tolerancia && i < max_iteraciones) {
    Tipo f0 = funcion(x0);
    Tipo f1 = funcion(x1);
    if (fabs(f1 - f0) < 1e-14) {
      break;
    }
    x = x1 - f1 * (x1 - x0) / (f1 - f0);
    x0 = x1;
    x1 = x;
    if (fabs(x1) < 1e-14) {
      e = fabs(x1 - x0);
    } else {
      e = fabs(x1 - x0) / fabs(x1);
    }
    i++;
  }
  return {x1, e, i};
}
template <typename Tipo>
tuple<vector<Tipo>, vector<Tipo>>
Euler_implicito(const vector<Tipo> &r0, const function<Tipo(Tipo, Tipo)> f,
                Tipo h, vector<Tipo> &D) {
  int N = D[1] / h;
  vector<Tipo> x(N + 1), y(N + 1);
  x[0] = r0[0];
  y[0] = r0[1];
  for (int i = 0; i < N; i++) {
    x[i + 1] = x[i] + h;
    y[i + 1] = get<0>(
        Metodo_de_la_Secante(function<Tipo(Tipo)>([&](Tipo y1) {
                               return y1 - f(x[i + 1], y1) * h - y[i];
                             }),
                             y[i], y[i] + h * f(x[i + 1], y[i]), 1e-8, 1000));
  }
  return {x, y};
}

template <typename Tipo>
tuple<vector<Tipo>, vector<Tipo>>
Error(const vector<Tipo> x, const vector<Tipo> y,
      const function<Tipo(Tipo)> Solucion_exacta) {
  int N = x.size();
  vector<Tipo> y_exacta(N), error(N);
  for (int i = 0; i < N; i++) {
    y_exacta[i] = Solucion_exacta(x[i]);
    error[i] = abs((y_exacta[i] - y[i]));
  }
  return {y_exacta, error};
}

int main() {
  double pi = 355. / 113.;
  vector<double> r0 = {0, 1};
  vector<double> D = {0, 4 * pi};
  auto [x, y] = Euler_implicito<double>(
      r0, [](double x, double y) { return -1000 * (y - cos(x)) - sin(x); },
      0.2 * pi, D);
  auto [y_exacta, error] = Error<double>(x, y, [](double x) { return cos(x); });
  ofstream file("datos.csv");
  file << "x,y,y_exacta,error\n";
  for (int i = 0; i < x.size(); i++) {
    file << x[i] << "," << y[i] << "," << y_exacta[i] << "," << error[i]
         << "\n";
  }
  return 0;
}
