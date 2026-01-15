#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>

#include "BasicTSPheuristics.hpp"

enum class EdgeWeightType {
	EUC_2D,
	ATT,
	GEO,
	CEIL_2D
};

double g_calc_distance(double x1, double y1, double x2, double y2, EdgeWeightType type) {
	double dx = x1 - x2;
	double dy = y1 - y2;

	switch (type) {
		case EdgeWeightType::EUC_2D:
			return sqrt(dx * dx + dy * dy);

		case EdgeWeightType::ATT: {
			double rij = sqrt((dx * dx + dy * dy) / 10.0);
			int tij = static_cast<int>(round(rij));
			return (tij < rij) ? tij + 1 : tij;
		}

		case EdgeWeightType::CEIL_2D:
			return ceil(sqrt(dx * dx + dy * dy));

		case EdgeWeightType::GEO: {
			double PI = 3.141592;
			auto to_geo = [PI](double x) {
				int deg = static_cast<int>(x);
				double min = x - deg;
				return PI * (deg + 5.0 * min / 3.0) / 180.0;
			};
			double lat1 = to_geo(x1);
			double lon1 = to_geo(y1);
			double lat2 = to_geo(x2);
			double lon2 = to_geo(y2);
			double RRR = 6378.388;
			double q1 = cos(lon1 - lon2);
			double q2 = cos(lat1 - lat2);
			double q3 = cos(lat1 + lat2);
			return static_cast<int>(RRR * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
		}

		default:
			return sqrt(dx * dx + dy * dy);
	}
}

int get_tsplib(const string &as_name, vector<CityPair> &av_dat, vector<vector<double>> &ad_dist, vector<vector<double>> &ad_xy) {
	string s_filename = "tsplib/" + as_name + ".tsp";
	ifstream file(s_filename);

	if (!file.is_open()) {
		throw runtime_error("Could not open file: " + s_filename);
	}

	string line;
	int n = 0;
	EdgeWeightType edge_type = EdgeWeightType::EUC_2D;

	// Parse header
	while (getline(file, line)) {
		// Remove leading/trailing whitespace
		size_t start = line.find_first_not_of(" \t");
		if (start == string::npos) continue;
		line = line.substr(start);

		if (line.find("DIMENSION") != string::npos) {
			size_t pos = line.find(':');
			if (pos == string::npos) pos = line.find_first_of(" \t", 9);
			if (pos != string::npos) {
				n = stoi(line.substr(pos + 1));
			}
		}
		else if (line.find("EDGE_WEIGHT_TYPE") != string::npos) {
			if (line.find("EUC_2D") != string::npos) {
				edge_type = EdgeWeightType::EUC_2D;
			}
			else if (line.find("ATT") != string::npos) {
				edge_type = EdgeWeightType::ATT;
			}
			else if (line.find("GEO") != string::npos) {
				edge_type = EdgeWeightType::GEO;
			}
			else if (line.find("CEIL_2D") != string::npos) {
				edge_type = EdgeWeightType::CEIL_2D;
			}
			else if (line.find("EXPLICIT") != string::npos) {
				throw runtime_error("EXPLICIT edge weight type not yet supported: " + as_name);
			}
		}
		else if (line.find("NODE_COORD_SECTION") != string::npos) {
			break;
		}
	}

	if (n == 0) {
		throw runtime_error("Could not find DIMENSION in file: " + s_filename);
	}

	// Resize data structures
	vector<double> d_x(n);
	vector<double> d_y(n);

	ad_dist.resize(n);
	ad_xy.resize(n);
	for (int i = 0; i < n; ++i) {
		ad_dist[i].resize(n);
		ad_xy[i].resize(2);
	}

	// Parse coordinates
	int nodes_read = 0;
	while (getline(file, line) && nodes_read < n) {
		size_t start = line.find_first_not_of(" \t");
		if (start == string::npos) continue;
		line = line.substr(start);

		if (line.find("EOF") != string::npos) break;

		istringstream iss(line);
		int id;
		double x, y;
		if (iss >> id >> x >> y) {
			int idx = id - 1;  // TSPLIB uses 1-based indexing
			d_x[idx] = x;
			d_y[idx] = y;
			ad_xy[idx][0] = x;
			ad_xy[idx][1] = y;
			++nodes_read;
		}
	}

	file.close();

	if (nodes_read != n) {
		throw runtime_error("Expected " + to_string(n) + " nodes but read " + to_string(nodes_read));
	}

	// Calculate distances and populate av_dat
	av_dat.resize(n * n);
	int k = 0;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			av_dat[k].i = i;
			av_dat[k].j = j;
			av_dat[k].dij = g_calc_distance(d_x[i], d_y[i], d_x[j], d_y[j], edge_type);
			ad_dist[i][j] = av_dat[k].dij;
			++k;
		}
	}

	return n;
}

// Overload that reads instance name from config file
// CSV format: param,value,comment (third column ignored)
// Expected row: instance,<name>,<optional comment>
//
int get_tsplib(vector<CityPair> &av_dat, vector<vector<double>> &ad_dist, vector<vector<double>> &ad_xy, string config_fname="") {
	string s_fname = config_fname.empty() ? "tsplib_option_default.csv" : config_fname;
	ifstream cfg(s_fname);
	if (!cfg.is_open()) {
		throw runtime_error("Could not open config file: " + s_fname);
	}

	string line, s_name;
	getline(cfg, line); // skip header
	while (getline(cfg, line)) {
		istringstream iss(line);
		string param, value;
		if (getline(iss, param, ',') && getline(iss, value, ',')) {
			if (param == "instance") {
				s_name = value;
				break;
			}
		}
	}
	cfg.close();

	if (s_name.empty()) {
		throw runtime_error("Could not find 'instance' in tsplib_option.csv");
	}

	return get_tsplib(s_name, av_dat, ad_dist, ad_xy);
}

