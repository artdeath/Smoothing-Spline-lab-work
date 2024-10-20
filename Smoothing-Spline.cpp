#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>



class Spline {

	double p;

	std::vector<double> x;

	std::vector<double> y;

	std::vector<double> w;

	std::vector<double> alpha;

	std::vector<double> main_diag;

	std::vector<double> p_diag;

	std::vector<double> right;

public:


	Spline() {

		p = 0.;

	}

	void set_p(double Input) {

		p = Input;

	}

	double get_p() {

		return p;

	}

	void push_x(double Input) {

		x.push_back(Input);

	}

	void set_x(int Index, double Input) {

		if ((Index >= x.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		x.at(Index) = Input;

	}

	double get_x(int Index) {

		if ((Index >= x.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		return x.at(Index);

	}

	size_t size_x() {

		return x.size();

	}

	void push_y(double Input) {

		y.push_back(Input);

	}

	void set_y(int Index, double Input) {

		if ((Index >= y.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		y.at(Index) = Input;

	}

	double get_y(int Index) {

		if ((Index >= y.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		return y.at(Index);

	}

	void push_w(double Input) {

		w.push_back(Input);

	}

	void set_w(int Index, double Input) {

		if ((Index >= w.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		w.at(Index) = Input;

	}

	double get_w(int Index) {

		if ((Index >= w.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		return w.at(Index);

	}


	double get_alpha(int Index) {

		if ((Index >= alpha.size()) or (Index < 0)) {

			throw std::runtime_error("Out of range\n");

		}

		return alpha.at(Index);

	}

	double First_Phi(double ksi) {

		return (1. - ksi) / 2.;

	}

	double Second_Phi(double ksi) {

		return (1. + ksi) / 2.;

	}

	double Ksi(double x, double x_first, double x_second) {

		return (2. * (x - x_first) / (x_second - x_first)) - 1.;

	}

	double md_first(int i) { //главная диагональ

		return (w.at(i) * First_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) * First_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) + w.at(i + 1) * First_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))) * First_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))));

	}

	double md_second(int i) { //главная диагональ

		return (w.at(i) * Second_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) * Second_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) + w.at(i + 1) * Second_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))) * Second_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))));

	}

	double pd(int i) { //побочная

		return (w.at(i) * First_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) * Second_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) + w.at(i + 1)* First_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1)))* Second_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))));

	}

	double right_first(int i) {

		return (w.at(i) * y.at(i) * First_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1)))+ w.at(i + 1)* y.at(i + 1)* First_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))));

	}

	double right_second(int i) {

		return (w.at(i) * y.at(i) * Second_Phi(Ksi(x.at(i), x.at(i), x.at(i + 1))) + w.at(i + 1) * y.at(i + 1) * Second_Phi(Ksi(x.at(i + 1), x.at(i), x.at(i + 1))));

	}


	void gMake() {

		double Math = 0.;

		for (int i = 0; i < y.size(); i++) {

			Math += y.at(i);

		}

		Math /= x.size();


		main_diag.push_back(0.);

		right.push_back(0.);

		for (int i = 0; i < x.size() - 1; i++) { 

			main_diag.push_back(0.);

			right.push_back(0.);

			main_diag.at(i) += (1. - p) * md_first(i);

			main_diag.at(i + 1) += (1. - p) * md_second(i);

			right.at(i) += (1. - p) * right_first(i);

			right.at(i + 1) += (1. - p) * right_second(i);


		}

		right.at(0) += p * Math;

		for (int i = 1; i < x.size() - 1; i++) {

			right.at(i) += p * 2. * Math;

		}

		right.at(right.size() - 1) += p * Math;

		for (int i = 0; i < x.size() - 1; i++) {

			main_diag.at(i) += p / (x.at(i + 1) - x.at(i));

			main_diag.at(i + 1) += p / (x.at(i + 1) - x.at(i));

		}

		for (int i = 0; i < x.size() - 1; i++) {

			p_diag.push_back((1. - p) * pd(i));

			p_diag.back() -= p / (2 * x.at(i + 1) - x.at(i));

		}

		std::vector<double> v; std::vector<double> u;

		double z; //знаменатель метода прогонки

		z = main_diag.at(0); v.push_back(-p_diag.at(0) / z); u.push_back(right.at(0) / z);

		

		for (int i = 1; i < main_diag.size() - 1; i++) {

			z = main_diag.at(i) + p_diag.at(i - 1) * v.at(i - 1);

			v.push_back(-p_diag.at(i - 1) / z);

			u.push_back((right.at(i) - p_diag.at(i - 1) * u.at(i - 1)) / z);

		}

		z = main_diag.at(main_diag.size() - 1) + p_diag.at(p_diag.size() - 1) * v.at(v.size() - 1); v.push_back(0.); u.push_back((right.at(right.size() - 1) - p_diag.at(p_diag.size() - 1) * u.at(u.size() - 1)) / z);


		alpha.push_back(u.at(u.size() - 1));


		for (int i = u.size() - 2; i >= 0; i--) {

			alpha.push_back(v.at(i) * alpha.back() + u.at(i));

		}

		std::reverse(alpha.begin(), alpha.end()); //~~записывает значение в обратном порядке


	}

	double d_First_Phi(int i) {

		return (-1 / (x.at(i + 1) - x.at(i)));

	}

	double d_Second_Phi(int i) {

		return (1 / (x.at(i + 1) - x.at(i)));

	}
};



int main() {

	Spline a;

	std::ifstream ifs("Input.txt");

	if (ifs.is_open()) {

		double temp = 0.;

		int i = 1;

		while (ifs >> temp) {

			a.push_x(i);

			a.push_y(temp);

			a.push_w(1.);

			if (std::fabs(a.get_y(i - 1)) > 7) {

				a.set_w(i - 1, 0.1);

			}

			i++;

		}


	}

	a.set_p(0.5);

	a.gMake();

	ifs.close();
	ifs.clear();


	std::ofstream ofs("Output.txt");

	if (ofs.is_open()) {

		for (int i = 0; i < a.size_x(); i++) {

			int k = 1;

			while (!(a.get_x(i) < a.get_x(k))  and (std::fabs((a.get_x(i) - a.get_x(k))) > 1e-7)) {

				k++;

			}

			ofs << (a.get_alpha(k - 1) * a.First_Phi(a.Ksi(a.get_x(i), a.get_x(k - 1), a.get_x(k))) + a.get_alpha(k) * a.Second_Phi(a.Ksi(a.get_x(i), a.get_x(k - 1), a.get_x(k)))) << std::endl;

		}


	}

	return 0;
}


//0.532396358
//-3.508040518
//1.19580284
//3.190684949
//3.138327355