/*******************************
UI3-N500-05   Jason Ma
Computation of Frechet Distance
(Baseline)
7/7/2017 Friday
V1.0.0
********************************/

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<ctime>
#include<iomanip>

using namespace std;

int read_queryset(string fname, string * & query_no, double * & range, int pos_det) {
	ifstream fin(fname.c_str());
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();

	int count = 0;
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	//cout << pos_last << endl;
	while (pos_nl != string::npos) {
		count++;
		pos_nl = lines.find('\n', pos_safe + 1);
		if (pos_safe == pos_nl - 1) {
			count = count - 1;
			break;
		}
		if (pos_last == pos_nl + 1) {
			break;
		}
		pos_safe = pos_nl;
		//cout << pos_nl << endl;
	}

	//cout << count << endl;

	query_no = new string[count];
	range = new double[count];
	pos_nl = 0;
	int pos_space = 0;
	int pos_int = 0;
	for (int i = 0; i < count; i++) {
		pos_space = lines.find(' ', pos_nl);
		query_no[i] = lines.substr(pos_nl, pos_space - pos_nl).c_str();
		if (pos_det != string::npos) {
			pos_nl = lines.find('\r', pos_space);
			while (((int)lines[pos_space]) > 58 || ((int)lines[pos_space]) < 48) {
				pos_space++;
			}
			range[i] = atof(lines.substr(pos_space, pos_nl - pos_space).c_str());
			pos_nl = pos_nl + 2;
		}
		else {
			pos_nl = lines.find('\n', pos_space);
			while (((int)lines[pos_space]) > 58 || ((int)lines[pos_space]) < 48) {
				pos_space++;
			}
			range[i] = atof(lines.substr(pos_space, pos_nl - pos_space).c_str());
			pos_nl = pos_nl + 1;
		}
	}

	//for (int i = 0; i < count; i++) {
		//cout << range[i] << endl;
	//}

	return count;
}


int read_dataset(string fname, string * & data_no, int & pos_det) {
	ifstream fin(fname.c_str());
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();

	int count = 0;
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	//cout << pos_last<<endl;
	while (pos_nl != string::npos) {
		count++;
		pos_nl = lines.find('\n', pos_safe + 1);
		if (pos_safe == pos_nl - 1) {
			count = count - 1;
			break;
		}
		if (pos_last == pos_nl + 1) {
			break;
		}
		pos_safe = pos_nl;
		//cout << pos_nl << endl;
	}

	//cout <<count<< endl;
	
	pos_det = lines.find('\r');
	data_no = new string[count];
	pos_nl = 0;
	int pos_zero = 0;
	for (int i = 0; i < count; i++) {
		if (pos_det != string::npos) {
			pos_zero = lines.find('\r', pos_nl);
			data_no[i] = lines.substr(pos_nl, pos_zero - pos_nl).c_str();
			pos_nl = pos_zero + 2;
		}
		else {
			pos_zero = lines.find('\n', pos_nl);
			data_no[i] = lines.substr(pos_nl, pos_zero - pos_nl).c_str();
			pos_nl = pos_zero + 1;
		}
	}

	//for (int i = 0; i < count; i++) {
		//cout << data_no[i] << endl;
	//}

	return count;
}

int read_trj(string fname, double * & x_value, double * & y_value) {
	ifstream fin(fname.c_str());
	//cout<<fin.is_open()<<endl;
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();
	
	//cout<<lines<<endl;

	int count = -1;
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	while (pos_nl != string::npos) {
		count++;
		pos_nl = lines.find('\n', pos_safe + 1);
		if (pos_safe == pos_nl - 1) {
			count = count - 1;
			break;
		}
		if (pos_last == pos_nl + 1) {
			break;
		}
		pos_safe = pos_nl;
	}
	
	//cout<<count<<endl;

	x_value = new double[count];
	y_value = new double[count];
	pos_nl = 0;
	int pos_xs = 0;
	int pos_ys = 0;
	for (int i = 0; i < count; i++) {
		pos_nl = lines.find('\n', pos_ys);
		pos_xs = lines.find(' ', pos_nl);
		x_value[i] = atof(lines.substr(pos_nl + 1, pos_xs - pos_nl - 1).c_str());
		pos_ys = lines.find(' ', pos_xs + 1);
		y_value[i] = atof(lines.substr(pos_xs + 1, pos_ys - pos_xs - 1).c_str());
	}

	return count;
}

double cal_distance(double * & x_value_one, double * & y_value_one, double * & x_value_two, double * & y_value_two, int i, int j) {
	double distance = sqrt(pow((x_value_one[i] - x_value_two[j]), 2) + pow((y_value_one[i] - y_value_two[j]), 2));
	return distance;
}



void com_cell(double * & cell, bool * & label, double * & x_value_one, double * & y_value_one, double * & x_value_two, double * & y_value_two, int a, int d, int b, int c, int p, double range) {
	/*Compute the left empty space of a cell*/
	double x_one = x_value_one[a];
	double y_one = y_value_one[a];
	double x_two = x_value_two[b];
	double y_two = y_value_two[b];
	double x_three = x_value_two[c];
	double y_three = y_value_two[c];
	double m;
	double k;
	double A;
	double B;
	double C;
	double det;
	if (x_three != x_two) {
		if (x_value_two[b] > x_value_two[c]) {
			x_two = x_value_two[c];
			y_two = y_value_two[c];
			x_three = x_value_two[b];
			y_three = y_value_two[b];
		}
		m = (y_three - y_two) / (x_three - x_two);
		k = y_two - m*x_two;
		A = 1 + pow(m, 2);
		B = 2 * m*(k - y_one) - 2 * x_one;
		C = pow(x_one, 2) + pow(k, 2) - 2 * k*y_one + pow(y_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;

		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= x_three&&solution >= x_two) {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = (x_three - solution) / (x_three - x_two);
					cell[4 * (b*p + a) + 3] = (x_three - solution) / (x_three - x_two);
				}
				else {
					cell[4 * (b*p + a) + 2] = (solution - x_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 3] = (solution - x_two) / (x_three - x_two);
				}
			}
			else {
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= x_two && solution_two >= x_three) {
				cell[4 * (b*p + a) + 2] = 0;
				cell[4 * (b*p + a) + 3] = 1;
			}
			else if ((solution_one < x_two && solution_two < x_two) || (solution_one > x_three && solution_two > x_three)) {
			}
			else if (solution_one >= x_two&&solution_one <= x_three&&solution_two >= x_three) {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 0;
					cell[4 * (b*p + a) + 3] = (x_three - solution_one) / (x_three - x_two);
				}
				else {
					cell[4 * (b*p + a) + 2] = (solution_one - x_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 3] = 1;
				}
			}
			else if (solution_two >= x_two&&solution_two <= x_three&&solution_one <= x_two) {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = (x_three - solution_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 3] = 1;
				}
				else {
					cell[4 * (b*p + a) + 2] = 0;
					cell[4 * (b*p + a) + 3] = (solution_two - x_two) / (x_three - x_two);
				}
			}
			else {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = (x_three - solution_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 3] = (x_three - solution_one) / (x_three - x_two);
				}
				else {
					cell[4 * (b*p + a) + 2] = (solution_one - x_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 3] = (solution_two - x_two) / (x_three - x_two);
				}
			}
		}
	}
	else {
		if (y_value_two[b] > y_value_two[c]) {
			x_two = x_value_two[c];
			y_two = y_value_two[c];
			x_three = x_value_two[b];
			y_three = y_value_two[b];
		}
		A = 1;
		B = -2 * y_one;
		C = pow(y_one, 2) + pow(x_two, 2) - 2 * x_one*x_two + pow(x_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;
		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= y_three&&solution >= y_two) {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = (y_three - solution) / (y_three - y_two);
					cell[4 * (b*p + a) + 3] = (y_three - solution) / (y_three - y_two);
				}
				else {
					cell[4 * (b*p + a) + 2] = (solution - y_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 3] = (solution - y_two) / (y_three - y_two);
				}
			}
			else {
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= y_two && solution_two >= y_three) {
				cell[4 * (b*p + a) + 2] = 0;
				cell[4 * (b*p + a) + 3] = 1;
			}
			else if ((solution_one < y_two && solution_two < y_two) || (solution_one > y_three && solution_two > y_three)) {
			}
			else if (solution_one >= y_two&&solution_one <= y_three&&solution_two >= y_three) {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 0;
					cell[4 * (b*p + a) + 3] = (y_three - solution_one) / (y_three - y_two);
				}
				else {
					cell[4 * (b*p + a) + 2] = (solution_one - y_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 3] = 1;
				}
			}
			else if (solution_two >= y_two&&solution_two <= y_three&&solution_one <= y_two) {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = (y_three - solution_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 3] = 1;
				}
				else {
					cell[4 * (b*p + a) + 2] = 0;
					cell[4 * (b*p + a) + 3] = (solution_two - y_two) / (y_three - y_two);
				}
			}
			else {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = (y_three - solution_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 3] = (y_three - solution_one) / (y_three - y_two);
				}
				else {
					cell[4 * (b*p + a) + 2] = (solution_one - y_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 3] = (solution_two - y_two) / (y_three - y_two);
				}
			}
		}
	}

	//cout << m << endl << k << endl << A << endl << B << endl << C << endl << det;



	/*Compute the bottom empty space of a cell*/
	x_one = x_value_two[b];
	y_one = y_value_two[b];
	x_two = x_value_one[a];
	y_two = y_value_one[a];
	x_three = x_value_one[d];
	y_three = y_value_one[d];
	if (x_three != x_two) {
		if (x_value_one[a] > x_value_one[d]) {
			x_two = x_value_one[d];
			y_two = y_value_one[d];
			x_three = x_value_one[a];
			y_three = y_value_one[a];
		}
		m = (y_three - y_two) / (x_three - x_two);
		k = y_two - m*x_two;
		A = 1 + pow(m, 2);
		B = 2 * m*(k - y_one) - 2 * x_one;
		C = pow(x_one, 2) + pow(k, 2) - 2 * k*y_one + pow(y_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;

		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= x_three&&solution >= x_two) {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = (x_three - solution) / (x_three - x_two);
					cell[4 * (b*p + a) + 1] = (x_three - solution) / (x_three - x_two);
				}
				else {
					cell[4 * (b*p + a)] = (solution - x_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 1] = (solution - x_two) / (x_three - x_two);
				}
			}
			else {
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= x_two && solution_two >= x_three) {
				cell[4 * (b*p + a)] = 0;
				cell[4 * (b*p + a) + 1] = 1;
			}
			else if ((solution_one < x_two && solution_two < x_two) || (solution_one > x_three && solution_two > x_three)) {
			}
			else if (solution_one >= x_two&&solution_one <= x_three&&solution_two >= x_three) {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = 0;
					cell[4 * (b*p + a) + 1] = (x_three - solution_one) / (x_three - x_two);
				}
				else {
					cell[4 * (b*p + a)] = (solution_one - x_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 1] = 1;
				}
			}
			else if (solution_two >= x_two&&solution_two <= x_three&&solution_one <= x_two) {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = (x_three - solution_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 1] = 1;
				}
				else {
					cell[4 * (b*p + a)] = 0;
					cell[4 * (b*p + a) + 1] = (solution_two - x_two) / (x_three - x_two);
				}
			}
			else {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = (x_three - solution_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 1] = (x_three - solution_one) / (x_three - x_two);
				}
				else {
					cell[4 * (b*p + a)] = (solution_one - x_two) / (x_three - x_two);
					cell[4 * (b*p + a) + 1] = (solution_two - x_two) / (x_three - x_two);
				}
			}
		}
	}
	else {
		if (y_value_one[a] > y_value_one[d]) {
			x_two = x_value_one[d];
			y_two = y_value_one[d];
			x_three = x_value_one[a];
			y_three = y_value_one[a];
		}
		A = 1;
		B = -2 * y_one;
		C = pow(y_one, 2) + pow(x_two, 2) - 2 * x_one*x_two + pow(x_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;
		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= y_three&&solution >= y_two) {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = (y_three - solution) / (y_three - y_two);
					cell[4 * (b*p + a) + 1] = (y_three - solution) / (y_three - y_two);
				}
				else {
					cell[4 * (b*p + a)] = (solution - y_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 1] = (solution - y_two) / (y_three - y_two);
				}
			}
			else {
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= y_two && solution_two >= y_three) {
				cell[4 * (b*p + a)] = 0;
				cell[4 * (b*p + a) + 1] = 1;
			}
			else if ((solution_one < y_two && solution_two < y_two) || (solution_one > y_three && solution_two > y_three)) {
			}
			else if (solution_one >= y_two&&solution_one <= y_three&&solution_two >= y_three) {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = 0;
					cell[4 * (b*p + a) + 1] = (y_three - solution_one) / (y_three - y_two);
				}
				else {
					cell[4 * (b*p + a)] = (solution_one - y_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 1] = 1;
				}
			}
			else if (solution_two >= y_two&&solution_two <= y_three&&solution_one <= y_two) {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = (y_three - solution_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 1] = 1;
				}
				else {
					cell[4 * (b*p + a)] = 0;
					cell[4 * (b*p + a) + 1] = (solution_two - y_two) / (y_three - y_two);
				}
			}
			else {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = (y_three - solution_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 1] = (y_three - solution_one) / (y_three - y_two);
				}
				else {
					cell[4 * (b*p + a)] = (solution_one - y_two) / (y_three - y_two);
					cell[4 * (b*p + a) + 1] = (solution_two - y_two) / (y_three - y_two);
				}
			}
		}
	}

	label[b*p + a] = 1;
}

void dec_frechet(bool * & dp_label, double * & cell, bool * & label, int i, int j, int p, int q, double * & x_value_one, double * & y_value_one, double * & x_value_two, double * & y_value_two, double range, int & deter, int last) {
	if (deter == 1) {
		return;
	}

	if ((i == (p - 1) && j == (q - 1))) {
		deter = 1;
	}
	else if (i == (p - 1)) {
		if (label[(j + 1)*p + i] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i, i + 1, j + 1, j + 2, p, range);
		}

		if (cell[4 * ((j + 1)*p + i)] != -1) {
			if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 1) {
				dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);

				dp_label[2 * (j*p + i) + 1] = 0;
			}

			else if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 1) {
				if (cell[4 * (j*p + i) + 1] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				else if (cell[4 * (j*p + i) + 1] >= cell[4 * ((j + 1)*p + i)] && cell[4 * (j*p + i)] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}

				dp_label[2 * (j*p + i)] = 0;
			}
		}
	}
	else if (j == (q - 1)) {
		if (label[j*p + i + 1] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i + 1, i + 2, j, j + 1, p, range);
		}

		if (cell[4 * (j*p + i + 1) + 2] != -1) {
			if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 1) {
				dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				dp_label[2 * (j*p + i)] = 0;
			}
			else if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 1) {
				if (cell[4 * (j*p + i) + 3] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
				else if (cell[4 * (j*p + i) + 3] >= cell[4 * (j*p + i + 1) + 2] && cell[4 * (j*p + i) + 2] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
				dp_label[2 * (j*p + i) + 1] = 0;
			}
		}
	}
	else {
		if (label[j*p + i + 1] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i + 1, i + 2, j, j + 1, p, range);
		}

		if (cell[4 * (j*p + i + 1) + 2] != -1) {
			if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 1) {
				dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
			}
			else if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 1) {
				if (cell[4 * (j*p + i) + 3] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
				else if (cell[4 * (j*p + i) + 3] >= cell[4 * (j*p + i + 1) + 2] && cell[4 * (j*p + i) + 2] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
			}
		}

		if (label[(j + 1)*p + i] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i, i + 1, j + 1, j + 2, p, range);
		}

		if (cell[4 * ((j + 1)*p + i)] != -1) {
			if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 1) {
				dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				dp_label[2 * (j*p + i) + 1] = 0;
			}
			else if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 1) {
				if (cell[4 * (j*p + i) + 1] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				else if (cell[4 * (j*p + i) + 1] >= cell[4 * ((j + 1)*p + i)] && cell[4 * (j*p + i)] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				dp_label[2 * (j*p + i)] = 0;
			}
		}
	}
}

int main(void) {

	ofstream outResult;
	outResult.open("jason/baseline_result.txt");

	clock_t start_time=clock();

	int query_count;
	int dataset_count;
	string * data_no;
	string * query_no;
	double * range;
	int pos_det;
	
	//cout<<"checkpoint1"<<endl;
	
	dataset_count = read_dataset("jason/files/dataset.txt", data_no, pos_det);
	
	//cout<<"checkpoint2"<<endl;
	
	query_count = read_queryset("jason/files/query.txt", query_no, range, pos_det);
	
	//cout<<"checkpoint3"<<endl;

	clock_t query_close_time = 0;
	clock_t data_close_time=0;
	clock_t query_time=0;
	clock_t Frechet_time_end = 0;

	for (int i = 0; i < query_count; i++) {
		
		clock_t loop_start = clock();

		ofstream outFile;
		string answer_name;
		if (i >= 0 && i <= 9) {
			answer_name = "jason/files/baseline_result-000" + to_string(i)+".txt";
		}
		else if (i >= 10 && i <= 99) {
			answer_name = "jason/files/baseline_result-00" + to_string(i) + ".txt";
		}
		else if (i >= 100 && i <= 999) {
			answer_name = "jason/files/baseline_result-0" + to_string(i) + ".txt";
		}
		else {
			answer_name = "jason/files/baseline_result-" + to_string(i) + ".txt";
		}
		outFile.open(answer_name);

		clock_t query_open_time = clock();

		double * x_value_one;
		double * y_value_one;
		
		//cout<<"checkpoint4"<<endl;
		
		//cout<<query_no[i]<<endl;
		
		int count_one = read_trj(query_no[i], x_value_one, y_value_one);
		
		//cout<<count_one<<endl;
		
		//cout<<"checkpoint5"<<endl;

		query_close_time = query_close_time + clock() - query_open_time;

		for (int j = 0; j < dataset_count; j++) {

			clock_t data_open_time = clock(); //time to open the file

			double * x_value_two;
			double * y_value_two;
			
			//cout<<"checkpoint6"<<endl;
			
			//cout<<data_no[j]<<endl;
			
			int count_two = read_trj(data_no[j], x_value_two, y_value_two);
			
			//for(int m=0;m<count_two;m++){
				//cout<<setprecision(10)<<setiosflags(ios::showpoint)<<x_value_two[m]<<" "<<y_value_two[m]<<endl;
			//}
			
			//cout<<count_two<<endl;
			
			//cout<<"checkpoint7"<<endl;
			
			data_close_time = data_close_time + clock() - data_open_time;

			clock_t Frechet_time_start = clock();
			
			//cout<<"checkpoint8"<<endl;
			
			bool * label = new bool[(count_one - 1)*(count_two - 1)];
			
			//cout<<"checkpoint8.5\n";
			
			double * cell = new double[(count_one - 1)*(count_two - 1) * 4];

			for (int k = 0; k < (count_one - 1)*(count_two - 1) * 4; k++) {
				cell[k] = -1;
			}

			for (int k = 0; k < (count_one - 1)*(count_two - 1); k++) {
				label[k] = 0;
			}
			
			//cout<<"checkpoint9"<<endl;

			//com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, 19, 20, 0, 1, count_one, range);

			//for (int i = 0; i < 80; i++) {
			//cout << cell[i] << endl;
			//}

			int deter = 0;
			
			//cout<<"checkpoint10"<<endl;
			
			bool * dp_label = new bool[2*(count_one - 1)*(count_two - 1)];

			for (int i = 0; i < 2*(count_one - 1)*(count_two - 1); i++) {
				dp_label[i] = 1;
			}
			
			//cout<<"checkpoint11"<<endl;

			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, 0, 1, 0, 1, count_one - 1, range[i]);
			double end_point = cal_distance(x_value_one, y_value_one, x_value_two, y_value_two, count_one - 1, count_two - 1);

			if (end_point > range[i] || cell[0] == -1) {
			}
			else {
				dec_frechet(dp_label, cell, label, 0, 0, count_one - 1, count_two - 1, x_value_one, y_value_one, x_value_two, y_value_two, range[i], deter,1);
				if (deter == 1) {
					outFile << data_no[j] << endl;
				}
			}

			//for (int i = 0; i < (count_one - 1)*(count_two - 1); i++) {
				//cout << i << "    " << label[i] << endl;
			//}

			delete[] x_value_two;
			delete[] y_value_two;
			delete[] cell;
			delete[] label;
			delete[] dp_label;

			Frechet_time_end = Frechet_time_end + clock() - Frechet_time_start;
		}
		outFile.close();

		delete[] x_value_one;
		delete[] y_value_one;

		query_time = query_time + clock() - loop_start;
	}

	delete[] data_no;
	delete[] range;
	delete[] query_no;

	clock_t end_time = clock() - start_time; //overall running time

	outResult << "Average reading query time:" << " " << (float)query_close_time / (CLOCKS_PER_SEC*query_count)<<endl
		<<"Average reading data time:"<<" "<<(float)data_close_time/ (CLOCKS_PER_SEC*dataset_count) << endl
		<<"Average Frechet Distance Computation time"<<"  "<< (float)Frechet_time_end / (CLOCKS_PER_SEC*dataset_count*query_count) << endl
		<<"Average query time:" << " " << (float)query_close_time / (CLOCKS_PER_SEC*query_count) << endl
		<<"Overall running time:"<<" "<<(float)end_time/ CLOCKS_PER_SEC;

	outResult.close();

	return 0;
}

