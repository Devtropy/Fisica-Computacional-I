#include <cmath>
#include <fstream>
#include <functional>
#include <tuple>
#include <vector>

using namespace std;

template <typename Tipo>
tuple<vector<Tipo>, vector<Tipo>> Euler(const vector<Tipo> r0,
                                        const function<Tipo(Tipo, Tipo)> f,
                                        Tipo h, vector<Tipo> D) {
  int N = D[1] / h;
  vector<Tipo> x(N + 1), y(N + 1);
  x[0] = r0[0];
  y[0] = r0[1];
  for (int i = 0; i < N; i++) {
    x[i + 1] = x[i] + h;
    y[i + 1] = y[i] + f(x[i], y[i]) * h;
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

template <typename Tipo> Tipo Solucion_exacta(const Tipo x) {
  return sqrt((Tipo(2) / Tipo(3) * x * x * x + Tipo(4)));
}
template <typename Tipo> Tipo Campo(const Tipo x, const Tipo y) {
  return (x * x) / y;
}

int main() {
  vector<double> r0 = {0, 2};
  vector<double> D = {0, 2.1};
  auto [x, y] = Euler<double>(r0, Campo<double>, 0.1, D);
  auto [y_exacta, error] = Error<double>(x, y, Solucion_exacta<double>);
  ofstream file("datos.csv");
  file << "x,y,y_exacta,error\n";
  for (int i = 0; i < x.size(); i++) {
    file << x[i] << "," << y[i] << "," << y_exacta[i] << "," << error[i]
         << "\n";
  }
  return 0;
}
