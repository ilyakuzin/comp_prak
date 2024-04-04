#include <iostream>
#include <cmath>
#include <vector>

const int M = 5; // разбиение по оси t
const int N = 10; // разбиение по оси x.
const double HALF_ONE = 0.5;
const double T = 1.0;
const double h = 1.0 / N;
const double tau = T / (M - 1);
const double hSquare = h * h;

const double A = tau / hSquare;
const double B =  2 * tau / hSquare + 1;
const double C = A;

double mu_1(double t) {
    return sin(t);
}

double mu_2(double t) {
    return cos(t) * sin(t);
}

double phi_1(double t) {
    return 0;
}

double phi_2(double x, double t) {
    return cos(x) * sin(t);
}

double f(double x, double t) {
    return cos(x) * sin(t) + cos(x) * cos(t);
}

double resultFunction(double x, double t) {
    return cos(x) * sin(t);
}

int main() {
    setlocale(LC_ALL, "Russian");
    std::vector<double> result_i(N); // значения функции-ответа на последнем слое (ответ).

    // Вычисляем y_i на самом нижнем слое. Заодно вычисляем результат на самом верхнем слое.
    std::vector<double> prevY_i(N); // y_(i-1).
    std::vector<double> y_i(N);     // y_i
    std::vector<double> nextY_i(N); // y_(i+1).
    double x = 0.;

    //результирующая функция + у на 0 слое + краевые условия
    for (int i = 0; i < N; i++) { // j = 0.
        result_i[i] = resultFunction(x, T);
        prevY_i[i] = phi_1(x); //вычисление на 0 слое
        x += h;
    }
    prevY_i[0] = mu_1(0); //краевые условия -- слева
    prevY_i[N - 1] = mu_2(T); //справа

    double t = tau; // t^j.
    std::vector<double> D_i(N);

    {
        std::vector<double> f_i(N);     // f^0_i.
        std::vector<double> f_i_jpp(N); // f^(j+1)_i.
        std::vector<double> phi_2_i(N); //краевые условия
        x = 0; // x_i.

        for (int i = 1; i < N; i++) { // Прямой ход метода прогонки: вычисление f_i.
            f_i[i] = f(x, 0);
            f_i_jpp[i] = f(x, t);
            phi_2_i[i] = phi_2(x, t);
            x += h;
        }

        x = 0;
        for (int i = 0; i < N - 1; i++) {
            y_i[i] = prevY_i[i] + tau * (phi_2_i[i] + f_i[i]); //вычислили значения на первом слое
            result_i[i] = resultFunction(x, tau);
            x += h;
        }

        for (int i = 0; i < N; i++) { // Обновление коэффициентов D_i по y_i и y_(i-1).
            D_i[i] = (2 * y_i[i] -  prevY_i[i] + tau * f_i_jpp[i]);
        }
    }

    for (int j = 2; j <= M; j++) { // Проход снизу вверх от слоя к слою.
        std::vector<double> f_i_jpp(N); // f^(j+1)_i.
        std::vector<double> alpha_i(N);
        std::vector<double> beta_i(N);
        alpha_i[0] = 0;
        beta_i[0] = mu_1(t);
        x = 0; // x_i.

        for (int i = 1; i < N - 1; i++) { // Прямой ход метода прогонки: вычисление f_i_jpp, alpha_i, beta_i.
            f_i_jpp[i] = f(x, t + tau);
            double denom = B - C * alpha_i[i - 1];
            alpha_i[i] = A / denom; //по формуле
            beta_i[i] = (D_i[i - 1] + C * beta_i[i - 1]) / denom; //по формуле стр.19
            x += h;
        }

        nextY_i[N - 1] = mu_2(t);
        for (int i = N - 2; i >= 0; i--) { // Обратный ход метода прогонки: вычисление y_i по y_(i+1), alpha_(i+1) и beta_(i+1).
            nextY_i[i] = alpha_i[i] * nextY_i[i + 1] + beta_i[i];
        }
        x = 0.;
        //проверка
        for (int i = 0; i < N; i++) { // Обновление коэффициентов D_i по y_i и y_(i-1).
            result_i[i] = resultFunction(x, T);
            x += h;
        }

        prevY_i = y_i;
        y_i = nextY_i;

        for (int i = 0; i < N; i++) { // Обновление коэффициентов D_i по y_i и y_(i-1).
            D_i[i] = (2 * y_i[i] -  prevY_i[i] + tau * f_i_jpp[i]);
        }
        t += tau;
    }

    double z = 0.;
    for (int i = 0; i < N; i++) {
        z = std::max(z, std::abs(y_i[i] - result_i[i]));
    }

    double result_norma = 0.;
    for (int i = 0; i < N; i++) {
        result_norma = std::max(result_norma, std::abs(result_i[i]));
        std::cout << "|" << nextY_i[i] << " - " << result_i[i] << "| = " << std::abs(nextY_i[i] - result_i[i]) << std::endl;
    }

    std::cout << "Кубическая норма разности полученного решения и правильного, z: " << std::endl;
    std::cout << z << std::endl;

    std::cout << "Отношение z к норме правильного решения: " << std::endl;
    std::cout << z / result_norma << std::endl;

    return 0;
}
