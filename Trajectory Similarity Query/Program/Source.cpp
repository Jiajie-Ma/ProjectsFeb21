/*******************************
UI3-N500-05   Jason Ma
Computation of Frechet Distance
(Trajectory MMB R-tree Indexing)
(Parallelization)
7/26/2017 Wednesday
********************************/

//Header files for Computations
#include<boost/geometry.hpp>
#include<cmath>

//Header files for I/O
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

//Header files for parallelization
#include<omp.h>

//Boost Library (V1.6.2) namespace neat definition
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//Geometry Type Definition for Range Query
typedef bg::model::point<double, 2, bg::cs::cartesian> point; //point geometry
typedef bg::model::multi_point<point> pointset; //A set of points geometry
typedef bg::model::box<point> traj_box; //box geometry (aligned with the defined axis)
typedef std::pair<traj_box, int> value; //value stored in r-tree with pair feature
typedef bgi::rtree<value, bgi::rstar<20, 6, 0>> rtree; //r*-tree with maximum entry of 20, minimum entry of 6, and disabled re-insertion

//A robust function that reads all points coodinates (x,y respectively) from a trajectory into point geometry, 
//record four extreme points of the trajectory to construct MMB, working for both Windows and Linux System
//Return: # of points of a trajectory
int read_trj(std::string fname, point & up_point, point & down_point, point & left_point, point & right_point, int pos_det) {

	//Open the file and read the entire file into a string
	std::ifstream fin(fname.c_str());
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();

	//A dope module that computes # of lines in the file with string::find(), no matter it is:
	//1. end with several dead lines
	//2. end with no dead lines
	int count = -1;//# of lines starting with -1 for the first friendly line
	//A set of variables help compute # of lines
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	//For Windows system (without '\r')
	if (pos_det != std::string::npos) {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 2) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}
	//For Linux system (with'\r')
	else {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 1) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}

	//Read required information with string::find() and store them into corresponding point geometry
	//A set of variables that help find required information
	pos_nl = 0;
	int pos_xs = 0;
	int pos_ys = 0;
	for (int i = 0; i < count; i++) {
		//Read x coordinate
		pos_nl = lines.find('\n', pos_ys);
		pos_xs = lines.find(' ', pos_nl);
		double x_value = atof(lines.substr(pos_nl + 1, pos_xs - pos_nl - 1).c_str());
		//Read y coordinate
		pos_ys = lines.find(' ', pos_xs + 1);
		double y_value = atof(lines.substr(pos_xs + 1, pos_ys - pos_xs - 1).c_str());
		if (y_value > up_point.get<1>()) {
			up_point.set<0>(x_value);
			up_point.set<1>(y_value);
		}
		if (y_value < down_point.get<1>()) {
			down_point.set<0>(x_value);
			down_point.set<1>(y_value);
		}
		if (x_value < left_point.get<0>()) {
			left_point.set<0>(x_value);
			left_point.set<1>(y_value);
		}
		if (x_value > right_point.get<0>()) {
			right_point.set<0>(x_value);
			right_point.set<1>(y_value);
		}
	}

	return count;
}

//A robust function that reads all trajectory file names from "dataset.txt" ,
//construct MBBs for all trajectories and inserts them into R-tree by calling read_trj(),
//working for both Windows and Linux System
//Return: # of trajecotries in the whole dataset
int read_dataset(std::string dataset_name, std::string * & data_no, int * & traj_count, rtree & trj_rtree,  int & pos_det) {

	//Open the file and read the entire file into a string
	std::ifstream fin(dataset_name.c_str());
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();

	pos_det = lines.find('\r');//System Determinant: 1. Windows: "\n" 2. Linux: "\r\n"

	//A dope module that computes # of lines in the file with string::find(), no matter it is:
	//1. end with several dead lines
	//2. end with no dead lines
	int count = 0;//# of lines
	//A set of variables help compute # of lines
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	//For Windows system (without '\r')
	if (pos_det != std::string::npos) {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 2) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}
	//For Linux system (with'\r')
	else {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 1) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}

	//Read required information with string::find() and store them into a string array
	data_no = new std::string[count]; //String array storing all trajectory file names
	//A set of variables that help find required information
	pos_nl = 0;
	int pos_zero = 0;
	for (int i = 0; i < count; i++) {
		//For Linux System
		if (pos_det != std::string::npos) {
			pos_zero = lines.find('\r', pos_nl);
			data_no[i] = lines.substr(pos_nl, pos_zero - pos_nl).c_str();
			pos_nl = pos_zero + 2;
		}
		//For Windows System
		else {
			pos_zero = lines.find('\n', pos_nl);
			data_no[i] = lines.substr(pos_nl, pos_zero - pos_nl).c_str();
			pos_nl = pos_zero + 1;
		}
	}


	traj_count = new int[count];//int array that stores # of points of each trajectory


	//Read points from each trajectory data files, 
	//construct MBBs based on their extreme points
	//and insert them into R-tree
	for (int i = 0; i < count; ++i) {
		//Initialize the extreme points for record
		point up_point = { 0,-INFINITY };
		point down_point = { 0,INFINITY };
		point left_point = { INFINITY,0 };
		point right_point = { -INFINITY,0 };
		//Read points from trajectory data files and record their extreme points
		traj_count[i] = read_trj(data_no[i], up_point, down_point, left_point, right_point, pos_det);
		pointset ps; //pointset geometry used to store the extreme points for construction of MBB
		boost::geometry::append(ps, up_point);
		boost::geometry::append(ps, down_point);
		boost::geometry::append(ps, left_point);
		boost::geometry::append(ps, right_point);
		//Construct MBB
		traj_box traj_MMB; //result MBB
		boost::geometry::envelope(ps, traj_MMB);
		//Insert the MMB into R-tree with id of the trajectory
		trj_rtree.insert(std::make_pair(traj_MMB, i));
	}

	return count;
}

//A robust function that reads all query file names and query ranges from 'query.txt' into a string array, working for both Windows and Linux System
//Return: # of query trajectories
int read_queryset(std::string fname, std::string * & query_no, double * & range, int pos_det) {

	//Open the file and read the entire file into a string
	std::ifstream fin(fname.c_str());
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();

	//A dope module that computes # of lines in the file with string::find(), no matter it is:
	//1. end with several dead lines
	//2. end with no dead lines
	int count = 0; //# of lines in the file
	//A set of variables help compute # of lines
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	//For Windows system (without '\r')
	if (pos_det != std::string::npos) {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 2) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}
	//For Linux system (with'\r')
	else {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 1) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}

	//Read required information with string::find() and store them into corresponding arrays, working for both Windows and Linux system
	query_no = new std::string[count]; //array storing query file names
	range = new double[count]; //array storing query ranges
	//A set of variables that help find required information
	pos_nl = 0;
	int pos_space = 0;
	int pos_int = 0;
	for (int i = 0; i < count; i++) {
		//Read query file names
		pos_space = lines.find(' ', pos_nl);
		query_no[i] = lines.substr(pos_nl, pos_space - pos_nl).c_str();
		//Read corresponding query ranges
		//For Linux System
		if (pos_det != std::string::npos) {
			pos_nl = lines.find('\r', pos_space);
			while (((int)lines[pos_space]) > 58 || ((int)lines[pos_space]) < 48) {
				pos_space++;
			}
			range[i] = atof(lines.substr(pos_space, pos_nl - pos_space).c_str());
			pos_nl = pos_nl + 2;
		}
		//For Windows System
		else {
			pos_nl = lines.find('\n', pos_space);
			while (((int)lines[pos_space]) > 58 || ((int)lines[pos_space]) < 48) {
				pos_space++;
			}
			range[i] = atof(lines.substr(pos_space, pos_nl - pos_space).c_str());
			pos_nl = pos_nl + 1;
		}
	}

	return count;
}

//A robust function that reads all points of the query trajectory into two double arrays, 
//recording the four extreme points to make bounding box for tree qeury, working for both Windows and Linux System
//Return: # of points of the query trajectory
int read_query(std::string fname, double * & x_value, double * & y_value, point & up_point, point & down_point, point & left_point, point & right_point, int pos_det) {

	//Open the file and read the entire file into a string
	std::ifstream fin(fname.c_str());
	std::string lines((std::istreambuf_iterator<char>(fin)),
		std::istreambuf_iterator<char>());
	fin.close();

	//A dope module that computes # of lines in the file with string::find(), no matter it is:
	//1. end with several dead lines
	//2. end with no dead lines
	int count = -1;//# of lines starting with -1 for the first friendly line
	//A set of variables help compute # of lines
	int pos_safe = 0;
	int pos_nl = 0;
	int pos_last = lines.size();
	//For Windows system (without '\r')
	if (pos_det != std::string::npos) {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 2) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}
	//For Linux system (with'\r')
	else {
		while (pos_nl != std::string::npos) {
			count++;
			pos_nl = lines.find('\n', pos_safe + 1);
			//Check for several dead lines condition
			if (pos_safe == pos_nl - 1) {
				count = count - 1;
				break;
			}
			//Check for no dead lines condition
			if (pos_last == pos_nl + 1) {
				break;
			}
			pos_safe = pos_nl;
		}
	}

	//Read points (x,y respectively) and record four extreme ones
	x_value = new double[count]; //Double array storing the x coordinate 
	y_value = new double[count]; //Double array storing the y coordinate
	//A set of variables that help find required information
	pos_nl = 0;
	int pos_xs = 0;
	int pos_ys = 0;
	for (int i = 0; i < count; i++) {
		//Read the x coordinate
		pos_nl = lines.find('\n', pos_ys);
		pos_xs = lines.find(' ', pos_nl);
		x_value[i] = atof(lines.substr(pos_nl + 1, pos_xs - pos_nl - 1).c_str());
		//Read the y coordinate
		pos_ys = lines.find(' ', pos_xs + 1);
		y_value[i] = atof(lines.substr(pos_xs + 1, pos_ys - pos_xs - 1).c_str());
		//Record the uppest point
		if (y_value[i] > up_point.get<1>()) {
			up_point.set<0>(x_value[i]);
			up_point.set<1>(y_value[i]);
		}
		//Record the downest point
		if (y_value[i] < down_point.get<1>()) {
			down_point.set<0>(x_value[i]);
			down_point.set<1>(y_value[i]);
		}
		//Record the leftest point
		if (x_value[i] < left_point.get<0>()) {
			left_point.set<0>(x_value[i]);
			left_point.set<1>(y_value[i]);
		}
		//Record the rightest point
		if (x_value[i] > right_point.get<0>()) {
			right_point.set<0>(x_value[i]);
			right_point.set<1>(y_value[i]);
		}
	}

	return count;
}

//A simple function that computes the Euclidean distance between two points
double cal_distance(double x_value_one, double y_value_one, double x_value_two, double y_value_two) {
	double distance = sqrt(pow((x_value_one - x_value_two), 2) + pow((y_value_one - y_value_two), 2));
	return distance;
}

//A cheesy function that implements free-space diagram cell computation in geometrical way:
//Consider a cell for two segments of two trajectories, we want to compute its left empty space and bottom empty space.
//To achieve this, we draw a circle with *radius* of the range query, *centre* of either segement's beginning vertices and analyze the intersection condition.
//The cell array will store the terminal points of the free space on the left & bottom edge of a cell (Valid: 1<=cell[i]<=2)
void com_cell(double * & cell, bool * & label, double * & x_value_one, double * & y_value_one, double * & x_value_two, double * & y_value_two, int a, int d, int b, int c, int p, double range) {

	//Compute the left empty space of a cell, we draw circle on segment one and make segment two as our studied segment
	//For the sake of convenience
	double x_one = x_value_one[a];
	double y_one = y_value_one[a];
	double x_two = x_value_two[b];
	double y_two = y_value_two[b];
	double x_three = x_value_two[c];
	double y_three = y_value_two[c];
	//Variables used for algebra computation
	double m;
	double k;
	double A;
	double B;
	double C;
	double det;
	//When the studied segment is not a vertical line
	if (x_three != x_two) {
		//For the sake of generalization, make sure x_three > x_two 
		if (x_value_two[b] > x_value_two[c]) {
			x_two = x_value_two[c];
			y_two = y_value_two[c];
			x_three = x_value_two[b];
			y_three = y_value_two[b];
		}
		//Compute the determinant for intersection analysis
		m = (y_three - y_two) / (x_three - x_two);
		k = y_two - m*x_two;
		A = 1 + pow(m, 2);
		B = 2 * m*(k - y_one) - 2 * x_one;
		C = pow(x_one, 2) + pow(k, 2) - 2 * k*y_one + pow(y_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;

		//When the circle doesn't intersect the studied segment line
		if (det < 0) {
		}
		//When the circle intersects the studied segment line with only one point (tangent)
		else if (det == 0) {
			//Compute the x coordinate of the intersection
			double solution = (-B + sqrt(det)) / (2 * A);
			//When the intersection is on the studied segment. Otherwise, the circle doesn't have true intersection with the segment
			if (solution <= x_three&&solution >= x_two) {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 1 + ((x_three - solution) / (x_three - x_two));
					cell[4 * (b*p + a) + 3] = 1 + ((x_three - solution) / (x_three - x_two));
				}
				else {
					cell[4 * (b*p + a) + 2] = 1 + ((solution - x_two) / (x_three - x_two));
					cell[4 * (b*p + a) + 3] = 1 + ((solution - x_two) / (x_three - x_two));
				}
			}
		}
		//When the circle intersects the studied segment line with two points
		else {
			//Compute the x coordinate of the intersections
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			//When the intersections are both not on the studied segment, but are on the "right" positions
			if (solution_one <= x_two && solution_two >= x_three) {
				cell[4 * (b*p + a) + 2] = 1;
				cell[4 * (b*p + a) + 3] = 2;
			}
			//When the intersections are both not on the studied segment, but aren't on the "right" positions
			else if ((solution_one < x_two && solution_two < x_two) || (solution_one > x_three && solution_two > x_three)) {
			}
			//When the smaller intersection is on the studied segment
			else if (solution_one >= x_two&&solution_one <= x_three&&solution_two >= x_three) {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 1;
					cell[4 * (b*p + a) + 3] = 1 + ((x_three - solution_one) / (x_three - x_two));
				}
				else {
					cell[4 * (b*p + a) + 2] = 1 + ((solution_one - x_two) / (x_three - x_two));
					cell[4 * (b*p + a) + 3] = 2;
				}
			}
			//When the larger intersection is on the studied segment
			else if (solution_two >= x_two&&solution_two <= x_three&&solution_one <= x_two) {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 1 + ((x_three - solution_two) / (x_three - x_two));
					cell[4 * (b*p + a) + 3] = 2;
				}
				else {
					cell[4 * (b*p + a) + 2] = 1;
					cell[4 * (b*p + a) + 3] = ((solution_two - x_two) / (x_three - x_two)) + 1;
				}
			}
			//When both intersections are on the studied segement
			else {
				if (x_value_two[b] > x_value_two[c]) {
					cell[4 * (b*p + a) + 2] = ((x_three - solution_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 3] = ((x_three - solution_one) / (x_three - x_two)) + 1;
				}
				else {
					cell[4 * (b*p + a) + 2] = ((solution_one - x_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 3] = ((solution_two - x_two) / (x_three - x_two)) + 1;
				}
			}
		}
	}
	//When the studied segment is a vertical line
	else {
		//For the sake of generalization, make sure y_three > y_two 
		if (y_value_two[b] > y_value_two[c]) {
			x_two = x_value_two[c];
			y_two = y_value_two[c];
			x_three = x_value_two[b];
			y_three = y_value_two[b];
		}
		//Compute the determinant for intersection analysis
		A = 1;
		B = -2 * y_one;
		C = pow(y_one, 2) + pow(x_two, 2) - 2 * x_one*x_two + pow(x_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;

		//Same discussion as above
		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= y_three&&solution >= y_two) {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = ((y_three - solution) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 3] = ((y_three - solution) / (y_three - y_two)) + 1;
				}
				else {
					cell[4 * (b*p + a) + 2] = ((solution - y_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 3] = ((solution - y_two) / (y_three - y_two)) + 1;
				}
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= y_two && solution_two >= y_three) {
				cell[4 * (b*p + a) + 2] = 1;
				cell[4 * (b*p + a) + 3] = 2;
			}
			else if ((solution_one < y_two && solution_two < y_two) || (solution_one > y_three && solution_two > y_three)) {
			}
			else if (solution_one >= y_two&&solution_one <= y_three&&solution_two >= y_three) {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 1;
					cell[4 * (b*p + a) + 3] = 1 + ((y_three - solution_one) / (y_three - y_two));
				}
				else {
					cell[4 * (b*p + a) + 2] = 1 + ((solution_one - y_two) / (y_three - y_two));
					cell[4 * (b*p + a) + 3] = 2;
				}
			}
			else if (solution_two >= y_two&&solution_two <= y_three&&solution_one <= y_two) {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = 1 + ((y_three - solution_two) / (y_three - y_two));
					cell[4 * (b*p + a) + 3] = 2;
				}
				else {
					cell[4 * (b*p + a) + 2] = 1;
					cell[4 * (b*p + a) + 3] = ((solution_two - y_two) / (y_three - y_two)) + 1;
				}
			}
			else {
				if (y_value_two[b] > y_value_two[c]) {
					cell[4 * (b*p + a) + 2] = ((y_three - solution_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 3] = ((y_three - solution_one) / (y_three - y_two)) + 1;
				}
				else {
					cell[4 * (b*p + a) + 2] = ((solution_one - y_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 3] = ((solution_two - y_two) / (y_three - y_two)) + 1;
				}
			}
		}
	}


	//Compute the bottom empty space of a cell, we draw circle on segment two and make segment one as our studied segment
	//For the sake of convenience
	x_one = x_value_two[b];
	y_one = y_value_two[b];
	x_two = x_value_one[a];
	y_two = y_value_one[a];
	x_three = x_value_one[d];
	y_three = y_value_one[d];
	//When the studied segment is not a vertical line
	if (x_three != x_two) {
		//For the sake of generalization, make sure y_three > y_two 
		if (x_value_one[a] > x_value_one[d]) {
			x_two = x_value_one[d];
			y_two = y_value_one[d];
			x_three = x_value_one[a];
			y_three = y_value_one[a];
		}
		//Compute the determinant for intersection analysis
		m = (y_three - y_two) / (x_three - x_two);
		k = y_two - m*x_two;
		A = 1 + pow(m, 2);
		B = 2 * m*(k - y_one) - 2 * x_one;
		C = pow(x_one, 2) + pow(k, 2) - 2 * k*y_one + pow(y_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;

		//Same discussion as above
		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= x_three&&solution >= x_two) {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = ((x_three - solution) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((x_three - solution) / (x_three - x_two)) + 1;
				}
				else {
					cell[4 * (b*p + a)] = ((solution - x_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((solution - x_two) / (x_three - x_two)) + 1;
				}
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= x_two && solution_two >= x_three) {
				cell[4 * (b*p + a)] = 1;
				cell[4 * (b*p + a) + 1] = 2;
			}
			else if ((solution_one < x_two && solution_two < x_two) || (solution_one > x_three && solution_two > x_three)) {
			}
			else if (solution_one >= x_two&&solution_one <= x_three&&solution_two >= x_three) {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = 1;
					cell[4 * (b*p + a) + 1] = ((x_three - solution_one) / (x_three - x_two)) + 1;
				}
				else {
					cell[4 * (b*p + a)] = ((solution_one - x_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 1] = 2;
				}
			}
			else if (solution_two >= x_two&&solution_two <= x_three&&solution_one <= x_two) {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = ((x_three - solution_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 1] = 2;
				}
				else {
					cell[4 * (b*p + a)] = 1;
					cell[4 * (b*p + a) + 1] = ((solution_two - x_two) / (x_three - x_two)) + 1;
				}
			}
			else {
				if (x_value_one[a] > x_value_one[d]) {
					cell[4 * (b*p + a)] = ((x_three - solution_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((x_three - solution_one) / (x_three - x_two)) + 1;
				}
				else {
					cell[4 * (b*p + a)] = ((solution_one - x_two) / (x_three - x_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((solution_two - x_two) / (x_three - x_two)) + 1;
				}
			}
		}
	}
	//When the studied segment is a vertical line
	else {
		//For the sake of generalization, make sure y_three > y_two 
		if (y_value_one[a] > y_value_one[d]) {
			x_two = x_value_one[d];
			y_two = y_value_one[d];
			x_three = x_value_one[a];
			y_three = y_value_one[a];
		}
		//Compute the determinant for intersection analysis
		A = 1;
		B = -2 * y_one;
		C = pow(y_one, 2) + pow(x_two, 2) - 2 * x_one*x_two + pow(x_one, 2) - pow(range, 2);
		det = pow(B, 2) - 4 * A*C;

		//Same discussion as above
		if (det < 0) {
		}
		else if (det == 0) {
			double solution = (-B + sqrt(det)) / (2 * A);
			if (solution <= y_three&&solution >= y_two) {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = ((y_three - solution) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((y_three - solution) / (y_three - y_two)) + 1;
				}
				else {
					cell[4 * (b*p + a)] = ((solution - y_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((solution - y_two) / (y_three - y_two)) + 1;
				}
			}
		}
		else {
			double solution_one = (-B - sqrt(det)) / (2 * A);
			double solution_two = (-B + sqrt(det)) / (2 * A);
			if (solution_one <= y_two && solution_two >= y_three) {
				cell[4 * (b*p + a)] = 1;
				cell[4 * (b*p + a) + 1] = 2;
			}
			else if ((solution_one < y_two && solution_two < y_two) || (solution_one > y_three && solution_two > y_three)) {
			}
			else if (solution_one >= y_two&&solution_one <= y_three&&solution_two >= y_three) {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = 1;
					cell[4 * (b*p + a) + 1] = ((y_three - solution_one) / (y_three - y_two)) + 1;
				}
				else {
					cell[4 * (b*p + a)] = ((solution_one - y_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 1] = 2;
				}
			}
			else if (solution_two >= y_two&&solution_two <= y_three&&solution_one <= y_two) {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = ((y_three - solution_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 1] = 2;
				}
				else {
					cell[4 * (b*p + a)] = 1;
					cell[4 * (b*p + a) + 1] = ((solution_two - y_two) / (y_three - y_two)) + 1;
				}
			}
			else {
				if (y_value_one[a] > y_value_one[d]) {
					cell[4 * (b*p + a)] = ((y_three - solution_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((y_three - solution_one) / (y_three - y_two)) + 1;
				}
				else {
					cell[4 * (b*p + a)] = ((solution_one - y_two) / (y_three - y_two)) + 1;
					cell[4 * (b*p + a) + 1] = ((solution_two - y_two) / (y_three - y_two)) + 1;
				}
			}
		}
	}

	//A determinant array (bool) to make sure we don't compute the same cell twice or moar!
	label[b*p + a] = 1;
}

//A recursive function that implements the depth-first search algorithm for the decision problem
void dec_frechet(bool * & dp_label, double * & cell, bool * & label, int i, int j, int p, int q, double * & x_value_one, double * & y_value_one, double * & x_value_two, double * & y_value_two, double range, int & deter, int last) {

	//Decision determinant: If it turns to 1, stop the recursion
	if (deter == 1) {
		return;
	}

	//When the search reaches the right-top cell of the diagram, stop the recursion. i.e. Make Decision determinant to 1
	if ((i == (p - 1) && j == (q - 1))) {
		deter = 1;
	}
	//When the search reaches the right edge of the diagram, only search upwards
	else if (i == (p - 1)) {
		//If the next cell is not computed, compute it!
		if (label[(j + 1)*p + i] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i, i + 1, j + 1, j + 2, p, range);
		}
		//Make sure there's some empty space at the bottom of the next cell
		if (cell[4 * ((j + 1)*p + i)] != 0) {
			//When the search comes from the left edge of the current cell:
			//If the dp_label is 1, stop
			if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 0) {
				dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				//Make dp_label 1 as we know this route is blocked so that we won't compute it again
				dp_label[2 * (j*p + i) + 1] = 1;
			}
			//When the search comes from the left edge of the current cell:
			//If the dp_label is 1, stop
			else if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 0) {
				//When the free space at the bottom of the current cell is gernarlly lefter than the one of the next cell
				if (cell[4 * (j*p + i) + 1] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				//When there's a "spatial intersection" between the free space at the bottom of the current cell and the one of the next cell 
				else if (cell[4 * (j*p + i) + 1] >= cell[4 * ((j + 1)*p + i)] && cell[4 * (j*p + i)] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				//Make dp_label 1 as we know this route is blocked so that we won't compute it again
				dp_label[2 * (j*p + i)] = 1;
			}
		}
	}
	//When the search reaches the top edge of the diagram, only search rightwards
	//Same discussion
	else if (j == (q - 1)) {

		if (label[j*p + i + 1] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i + 1, i + 2, j, j + 1, p, range);
		}

		if (cell[4 * (j*p + i + 1) + 2] != 0) {
			if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 0) {
				dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				dp_label[2 * (j*p + i)] = 1;
			}
			else if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 0) {
				if (cell[4 * (j*p + i) + 3] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
				else if (cell[4 * (j*p + i) + 3] >= cell[4 * (j*p + i + 1) + 2] && cell[4 * (j*p + i) + 2] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
				dp_label[2 * (j*p + i) + 1] = 1;
			}
		}
	}
	//When the search is at any other parts of the diagram, search both upwards and rightwards
	//Here we will do rightwards search first
	//Same discussion
	else {

		//Search rightwards
		if (label[j*p + i + 1] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i + 1, i + 2, j, j + 1, p, range);
		}

		if (cell[4 * (j*p + i + 1) + 2] != 0) {
			if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 0) {
				dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
			}
			else if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 0) {
				if (cell[4 * (j*p + i) + 3] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
				else if (cell[4 * (j*p + i) + 3] >= cell[4 * (j*p + i + 1) + 2] && cell[4 * (j*p + i) + 2] <= cell[4 * (j*p + i + 1) + 2]) {
					dec_frechet(dp_label, cell, label, i + 1, j, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 2);
				}
			}
		}

		//Search upwards
		if (label[(j + 1)*p + i] != 1) {
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, i, i + 1, j + 1, j + 2, p, range);
		}

		if (cell[4 * ((j + 1)*p + i)] != 0) {
			if ((last == 2 || last == 1) && dp_label[2 * (j*p + i) + 1] == 0) {
				dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				dp_label[2 * (j*p + i) + 1] = 1;
			}
			else if ((last == 0 || last == 1) && dp_label[2 * (j*p + i)] == 0) {
				if (cell[4 * (j*p + i) + 1] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				else if (cell[4 * (j*p + i) + 1] >= cell[4 * ((j + 1)*p + i)] && cell[4 * (j*p + i)] <= cell[4 * ((j + 1)*p + i)]) {
					dec_frechet(dp_label, cell, label, i, j + 1, p, q, x_value_one, y_value_one, x_value_two, y_value_two, range, deter, 0);
				}
				dp_label[2 * (j*p + i)] = 1;
			}
		}
	}
}

//The magic begins here
int main(void) {

	omp_set_num_threads(omp_get_num_procs()); //assign parallelization CPUs

	//Preparation Step (Tree Construction)
	int query_count; //# of query trajectories
	int dataset_count; //# of trajectories in the whole dataset
	int pos_det; //Position determinant that recognizes the machine system
	int * traj_count; //int elements storing # of points of each trajectory
	double * range; //double elements storing range queries of corresponding trajectories
	std::string * data_no; //string elements storing the names of the trajectory files in the whole dataset
	std::string * query_no; //string elements storing the names of the query trajectory files
	rtree trj_rtree; //An r*-stree for trajectory indexing and retrieval
	clock_t rtree_time_end;
	//Read "dataset.txt", store all trajectories' names in the dataset and insert all the points into our r-tree
	dataset_count = read_dataset("dataset.txt", data_no, traj_count, trj_rtree, pos_det);
	//Read "query.txt", store all query trajectories' names and their corresponding Frechet Distance range queries
	query_count = read_queryset("query.txt", query_no, range, pos_det);


	//This is where a query starts
    #pragma omp parallel for
	for (int i = 0; i < query_count; i++) {
		
		//Filter step
		double * x_value_one; //x coordinate of the query trajectory
		double * y_value_one; //y coordinate of the query trajectory
							  //Initialize the extreme points of the query trajecotry to make the Bounding box
		point up_point = { 0,-INFINITY };
		point down_point = { 0,INFINITY };
		point left_point = { INFINITY,0 };
		point right_point = { -INFINITY,0 };
		//Get # of points in the query trajectory
		//Read the coordinates and find its extreme points
		int count_one = read_query(query_no[i], x_value_one, y_value_one, up_point, down_point, left_point, right_point, pos_det);

		//Extend the positions of extreme points by the query range in their representative directions 
		//so that the box bounds the buffer of the trajectory
		up_point.set<1>(up_point.get<1>() + range[i]);
		down_point.set<1>(down_point.get<1>() - range[i]);
		left_point.set<0>(left_point.get<0>() - range[i]);
		right_point.set<0>(right_point.get<0>() + range[i]);

		//Put extreme points into a pointset to make the bounding box
		pointset ps;
		boost::geometry::append(ps, up_point);
		boost::geometry::append(ps, down_point);
		boost::geometry::append(ps, left_point);
		boost::geometry::append(ps, right_point);

		//Construct the bounding box by envelope()
		traj_box range_box;
		boost::geometry::envelope(ps, range_box);

		//Tree Query starts (Filter)
		std::vector<value> index_result; //Vector used to store the candiate trajectories
		trj_rtree.query(bgi::within(range_box), std::back_inserter(index_result)); //Tree query step

		//Refine step starts here!
		std::string answers; //string used to store answer trajectory file names
		int true_no = index_result.size(); //# of the candidate trajectories

		for (int k = 0; k < true_no; ++k) {

			int m = std::get<1>(index_result[k]); //record id of the candidate trajectory
			double * x_value_two; //x-coordinate of the candidate trajectory
			double * y_value_two; //y-coordinate of the candidate trajectory

			//Read the candidate trajectory
			std::ifstream fin(data_no[m]);
			std::string lines((std::istreambuf_iterator<char>(fin)),
				std::istreambuf_iterator<char>());
			fin.close();

			int count_two = traj_count[m];
			x_value_two = new double[traj_count[m]];
			y_value_two = new double[traj_count[m]];
			int pos_nl = 0;
			int pos_xs = 0;
			int pos_ys = 0;
			for (int i = 0; i < traj_count[m]; i++) {
				pos_nl = lines.find('\n', pos_ys);
				pos_xs = lines.find(' ', pos_nl);
				x_value_two[i] = atof(lines.substr(pos_nl + 1, pos_xs - pos_nl - 1).c_str());
				pos_ys = lines.find(' ', pos_xs + 1);
				y_value_two[i] = atof(lines.substr(pos_xs + 1, pos_ys - pos_xs - 1).c_str());
			}

			//Preparation & initialization
			int deter = 0; //Decision determinant used to stop the recursion when our search reaches the right-top cell of the diagram
			double * cell = new double[(count_one - 1)*(count_two - 1) * 4](); //Used to store the information of cells
			bool * label = new bool[(count_one - 1)*(count_two - 1)](); //Cell determinant used to make sure we don't compute the same cell twice or moar!
			bool * dp_label = new bool[(count_one - 1)*(count_two - 1) * 2](); //Search determinant used to make sure we don't search the same route twice

			//A little trick to filter moar
			//Compute the first cell ((0,0))
			com_cell(cell, label, x_value_one, y_value_one, x_value_two, y_value_two, 0, 1, 0, 1, count_one - 1, range[i]);
			//Compute the distance between the last points of two trajectories
			double end_point = cal_distance(x_value_one[count_one - 1], y_value_one[count_one - 1], x_value_two[count_two - 1], y_value_two[count_two - 1]);

			//If the distance between the last points of two trajectories is below the query range and,
			//   the first cell ((0,0)) is not filled
			//Begin depth-first search
			if (end_point <= range[i] && cell[0] != 0) {
				dec_frechet(dp_label, cell, label, 0, 0, count_one - 1, count_two - 1, x_value_one, y_value_one, x_value_two, y_value_two, range[i], deter, 1);
				//Save answer trajectories into the string
				if (deter == 1) {
					answers.append(data_no[m]);
					answers.append("\n");
				}
			}
			//Clean up!
			delete[] cell;
			delete[] label;
			delete[] dp_label;
			delete[] x_value_two;
			delete[] y_value_two;
		}

		//Open the output file to write answers in
		std::ofstream outFile;
		std::string answer_name;
		if (i >= 0 && i <= 9) {
			answer_name = "result-000" + std::to_string(i) + ".txt";
		}
		else if (i >= 10 && i <= 99) {
			answer_name = "result-00" + std::to_string(i) + ".txt";
		}
		else if (i >= 100 && i <= 999) {
			answer_name = "result-0" + std::to_string(i) + ".txt";
		}
		else {
			answer_name = "result-" + std::to_string(i) + ".txt";
		}
		outFile.open(answer_name);

		outFile << answers; //Write answers into result file

		outFile.close(); //close the file

		//Clean up!
		delete[] x_value_one;
		delete[] y_value_one;
	}
	//Clean up!
	delete[] data_no;
	delete[] query_no;
	delete[] range;
	delete[] traj_count;

	return 0;
}