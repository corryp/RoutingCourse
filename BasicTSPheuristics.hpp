#pragma once

#include <cstdio>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <cfloat>
using namespace std;

const double gd_tol = 0.000000005;
int gi_itctr;		//counter to count number of iterations

class CityPair {
public:
	int i;
	int j;
	double dij;

	CityPair(int ai_i, int ai_j, double ad_dij) : i(ai_i), j(ai_j), dij(ad_dij) {}
	CityPair() {}
};//CityPair

void g_gen_tsp_instance(const int ai_n, vector<CityPair> &av_dat, vector<vector<double>> &ad_dist, vector<vector<double>> &ad_xy) {
	vector<double> d_x(ai_n);
	vector<double> d_y(ai_n);

	ad_dist.resize(ai_n);
	ad_xy.resize(ai_n);
	for (int i = 0; i < ai_n; ++i) {
		ad_dist[i].resize(ai_n);
		ad_xy[i].resize(2);
	}//i

	for (int i = 0; i < ai_n; ++i) {
		d_x[i] = rand()*100.0 / RAND_MAX;
		d_y[i] = rand()*100.0 / RAND_MAX;
		ad_xy[i][0] = d_x[i];
		ad_xy[i][1] = d_y[i];
	}//i

	av_dat.resize(ai_n * ai_n);
	int k = 0;
	for (int i = 0; i < ai_n; ++i) {
		for (int j = 0; j < ai_n; ++j) {
			av_dat[k].i = i;
			av_dat[k].j = j;
			av_dat[k].dij = sqrt((d_x[i] - d_x[j])*(d_x[i] - d_x[j]) + (d_y[i] - d_y[j])*(d_y[i] - d_y[j]));
			ad_dist[i][j] = av_dat[k].dij;
			++k;
		}//j
	}//i
}//g_gen_tsp_instance

void g_output_tour_csv(const vector<int> &ai_tour, const vector<vector<double>> &ad_xy) {
	ofstream ofs("tsp.csv");
	ofs << "order,city,x,y" << endl;
	for (int k = 0; k < ai_tour.size(); ++k) {
		int city = ai_tour[k];
		ofs << k << "," << city << "," << ad_xy[city][0] << "," << ad_xy[city][1] << endl;
	}
	ofs.close();
}//g_output_tour_csv

double gd_calc_z(int ai_n, const vector<vector<double>> &ad_dist, const vector<int> &ai_tour) {
	double z = ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];

	for (int k = 1; k < ai_n; ++k)
		z += ad_dist[ai_tour[k]][ai_tour[k - 1]];

	return z;
}//gd_calc_z

//generates a random integer between 0 and ai_n - 1
//
int gi_rand(int ai_n) {
	if (ai_n == 0) return 0;
	int i = ai_n * (1.0 * rand() / RAND_MAX);
	i = (i < ai_n ? i : ai_n - 1);	//for the extremely rare occurence of sampling a value of exactly 1
	return i;
}//gi_rand

double gd_swap(int ai_n, const vector<vector<double>> &ad_dist, double ad_z, vector<int> &ai_tour, int ai_sw1, int ai_sw2) {
	double z;

	assert(ai_sw1 < ai_sw2);

	z = ad_z;

	if (ai_sw1 + 1 == ai_sw2) {
		z -= ai_sw1 > 0 ? ad_dist[ai_tour[ai_sw1 - 1]][ai_tour[ai_sw1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];
		z -= ad_dist[ai_tour[ai_sw1]][ai_tour[ai_sw2]];
		z -= ai_sw2 < ai_n - 1 ? ad_dist[ai_tour[ai_sw2]][ai_tour[ai_sw2 + 1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];

		z += ai_sw1 > 0 ? ad_dist[ai_tour[ai_sw1 - 1]][ai_tour[ai_sw2]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[ai_sw2]];
		z += ad_dist[ai_tour[ai_sw2]][ai_tour[ai_sw1]];
		z += ai_sw2 < ai_n - 1 ? ad_dist[ai_tour[ai_sw1]][ai_tour[ai_sw2 + 1]] : ad_dist[ai_tour[ai_sw1]][ai_tour[0]];
	}
	else if (ai_sw1 + ai_n == ai_sw2 + 1) {
		z -= ad_dist[ai_tour[ai_sw2 - 1]][ai_tour[ai_sw2]];
		z -= ad_dist[ai_tour[ai_sw2]][ai_tour[ai_sw1]];
		z -= ad_dist[ai_tour[ai_sw1]][ai_tour[ai_sw1 + 1]];

		z += ad_dist[ai_tour[ai_sw2 - 1]][ai_tour[ai_sw1]];
		z += ad_dist[ai_tour[ai_sw1]][ai_tour[ai_sw2]];
		z += ad_dist[ai_tour[ai_sw2]][ai_tour[ai_sw1 + 1]];
	}//if
	else {
		z -= ai_sw1 > 0 ? ad_dist[ai_tour[ai_sw1 - 1]][ai_tour[ai_sw1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];
		z -= ad_dist[ai_tour[ai_sw1]][ai_tour[ai_sw1 + 1]];
		z -= ad_dist[ai_tour[ai_sw2 - 1]][ai_tour[ai_sw2]];
		z -= ai_sw2 < ai_n - 1 ? ad_dist[ai_tour[ai_sw2]][ai_tour[ai_sw2 + 1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];

		z += ad_dist[ai_tour[ai_sw2 - 1]][ai_tour[ai_sw1]];
		z += ai_sw2 < ai_n - 1 ? ad_dist[ai_tour[ai_sw1]][ai_tour[ai_sw2 + 1]] : ad_dist[ai_tour[ai_sw1]][ai_tour[0]];
		z += ai_sw1 > 0 ? ad_dist[ai_tour[ai_sw1 - 1]][ai_tour[ai_sw2]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[ai_sw2]];
		z += ad_dist[ai_tour[ai_sw2]][ai_tour[ai_sw1 + 1]];
	}//else

	int i_tmp = ai_tour[ai_sw1];
	ai_tour[ai_sw1] = ai_tour[ai_sw2];
	ai_tour[ai_sw2] = i_tmp;

	return z;
}//gd_swap

//NOTE: z calc assume symmetric distance matrix
//
double gd_two_opt(int ai_n, const vector<vector<double>> &ad_dist, double ad_z, const vector<int> &ai_tour, int ai_hd1, int ai_hd2, vector<int> &ai_nbr) {
	double z;
	int k = ai_hd1;
	int h = ai_hd2 - 1;
	int i_tmp;

	z = ad_z;
	z -= ai_hd1 > 0 ? ad_dist[ai_tour[ai_hd1 - 1]][ai_tour[ai_hd1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];
	z -= ai_hd2 > 0 ? ad_dist[ai_tour[ai_hd2 - 1]][ai_tour[ai_hd2]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];
	z += ad_dist[ai_tour[ai_hd1]][ai_tour[ai_hd2]];
	z += ad_dist[ai_tour[ai_hd1 > 0 ? ai_hd1 - 1 : ai_n - 1]][ai_tour[ai_hd2 > 0 ? ai_hd2 - 1 : ai_n - 1]];

	ai_nbr = ai_tour;
	if (h < 0)
		h = ai_n - 1;
	while (k != h) {
		i_tmp = ai_nbr[k];
		ai_nbr[k] = ai_nbr[h];
		ai_nbr[h] = i_tmp;

		++k;
		if (k == ai_n)
			k = 0;
		if (k != h) {
			--h;
			if (h == -1)
				h = ai_n - 1;
		}//if
	}//while

	return z;
}//gd_two_opt

//ai_tour is the current tour, ai_nbr should be pre-sized and will be overwritten
//
double gd_insertion(int ai_n, const vector<vector<double>> &ad_dist, double ad_z, const vector<int> &ai_tour, int ai_rmv, int ai_ins, vector<int> &ai_nbr) {
	double z;

	int i_rmv = ai_rmv;
	int i_ins = ai_ins;

	int k, j;
	if (i_ins < i_rmv) {
		j = ai_tour[i_rmv];
		for (k = 0; k <= i_ins; ++k)
			ai_nbr[k] = ai_tour[k];
		for (k = i_rmv; k > i_ins + 1; --k)
			ai_nbr[k] = ai_tour[k - 1];
		ai_nbr[i_ins + 1] = j;
		for (k = i_rmv + 1; k < ai_n; ++k)
			ai_nbr[k] = ai_tour[k];
	}
	else {
		j = ai_tour[i_rmv];
		for (k = 0; k < i_rmv; ++k)
			ai_nbr[k] = ai_tour[k];
		for (k = i_rmv; k < i_ins; ++k)
			ai_nbr[k] = ai_tour[k + 1];
		ai_nbr[i_ins] = j;
		for (k = i_ins + 1; k < ai_n; ++k)
			ai_nbr[k] = ai_tour[k];
	}//else

	z = ad_z;
	z -= i_rmv > 0 ? ad_dist[ai_tour[i_rmv - 1]][ai_tour[i_rmv]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];				//-d_R-1,R
	z -= i_rmv < ai_n - 1 ? ad_dist[ai_tour[i_rmv]][ai_tour[i_rmv + 1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];		//-d_R,R+1
	z -= i_ins < ai_n - 1 ? ad_dist[ai_tour[i_ins]][ai_tour[i_ins + 1]] : ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];		//-d_I,I+1
	z += i_rmv == 0 ? ad_dist[ai_tour[ai_n - 1]][ai_tour[1]] : 
		(i_rmv == ai_n - 1 ? ad_dist[ai_tour[ai_n - 2]][ai_tour[0]] : ad_dist[ai_tour[i_rmv - 1]][ai_tour[i_rmv + 1]]);	//+d_R-1,R+1
	z += ad_dist[ai_tour[i_ins]][ai_tour[i_rmv]];																		//+d_I,R
	z += i_ins < ai_n - 1 ? ad_dist[ai_tour[i_rmv]][ai_tour[i_ins + 1]] : ad_dist[ai_tour[i_rmv]][ai_tour[0]];			//+d_R,I+1

	return z;
}//gd_insertion

 //ai_tour should contain a valid initial solution for the local search
 //
double gd_local_swap_search(int ai_n, vector<vector<double>> &ad_dist, vector<int> &ai_tour, bool ab_full_nbrhd) {
	double z, d_znbr, d_znbr_best;
	int i_sw1, i_sw2, i_sw1_best, i_sw2_best;
	int i_itctr = 0;

	d_znbr_best = gd_calc_z(ai_n, ad_dist, ai_tour);
	do {
		z = d_znbr_best;
		d_znbr_best = z;
		++i_itctr;
		for (i_sw1 = 0; i_sw1 < ai_n - 1; ++i_sw1) {
			for (i_sw2 = i_sw1 + 1; i_sw2 < ai_n; ++i_sw2) {
				d_znbr = gd_swap(ai_n, ad_dist, z, ai_tour, i_sw1, i_sw2);
				if (d_znbr < d_znbr_best) {
					d_znbr_best = d_znbr;
					if (!ab_full_nbrhd)
						break;	//accept first improvement found
					else {
						i_sw1_best = i_sw1;
						i_sw2_best = i_sw2;
					}//else
				}//if
				gd_swap(ai_n, ad_dist, d_znbr, ai_tour, i_sw1, i_sw2);	//undo the previous swap before trying the next one
			}//i_sw2
			if (!ab_full_nbrhd && d_znbr_best < z)
				break;
		}//i_sw1
		if (ab_full_nbrhd && d_znbr_best < z)
			d_znbr = gd_swap(ai_n, ad_dist, z, ai_tour, i_sw1_best, i_sw2_best);
	} while (d_znbr_best + gd_tol < z);

	gi_itctr = i_itctr;
	return z;
}//gd_local_swap_search

//ai_tour should contain a valid initial solution for the local search
//
double gd_local_insertion_search(int ai_n, vector<vector<double>> &ad_dist, vector<int> &ai_tour, bool ab_full_nbrhd) {	
	double z, d_znbr, d_znbr_best;
	vector<int> i_nbr(ai_n);
	int i_rmv, i_ins, i_rmv_best, i_ins_best;
	int i_itctr = 0;
	
	d_znbr_best = gd_calc_z(ai_n, ad_dist, ai_tour);
	do {
		++i_itctr;
		z = d_znbr_best;
		d_znbr_best = z;
		for (i_rmv = 0; i_rmv < ai_n; ++i_rmv) {
			for (i_ins = 0; i_ins < ai_n; ++i_ins) {
				if (i_rmv != i_ins && i_rmv != i_ins + 1 && i_rmv != i_ins - ai_n + 1) {
					d_znbr = gd_insertion(ai_n, ad_dist, z, ai_tour, i_rmv, i_ins, i_nbr);
					if (d_znbr < d_znbr_best) {
						d_znbr_best = d_znbr;
						if (!ab_full_nbrhd) {
							ai_tour = i_nbr;
							break;	//accept first improvement found
						}//if
						else {
							i_rmv_best = i_rmv;
							i_ins_best = i_ins;
						}//else
					}//if
				}//if
			}//i_ins
			if (!ab_full_nbrhd && d_znbr_best < z)
				break;
		}//i_rmv
		if (ab_full_nbrhd && d_znbr_best < z) {
			d_znbr = gd_insertion(ai_n, ad_dist, z, ai_tour, i_rmv_best, i_ins_best, i_nbr);
			ai_tour = i_nbr;
		}//if
	} while (d_znbr_best + gd_tol < z);

	gi_itctr = i_itctr;
	return z;
}//gd_local_insertion_search

 //ai_tour should contain a valid initial solution for the local search
 //
double gd_local_2opt_search(int ai_n, vector<vector<double>> &ad_dist, vector<int> &ai_tour, bool ab_full_nbrhd) {
	double z, d_znbr, d_znbr_best;
	vector<int> i_nbr(ai_n);
	int i_hd1, i_hd2, i_hd1_best, i_hd2_best;
	int i_itctr = 0;

	d_znbr_best = gd_calc_z(ai_n, ad_dist, ai_tour);
	do {
		++i_itctr;
		z = d_znbr_best;
		for (i_hd1 = 0; i_hd1 < ai_n; ++i_hd1) {
//			for (i_hd2 = 0; i_hd2 < ai_n; ++i_hd2) {
			for (i_hd2 = i_hd1 + 2; i_hd2 < ai_n; ++i_hd2) {
//			if (abs(i_hd1 - i_hd2) > 1 && abs(i_hd1 - i_hd2) < ai_n) {
					d_znbr = gd_two_opt(ai_n, ad_dist, z, ai_tour, i_hd1, i_hd2, i_nbr);
					if (d_znbr < d_znbr_best) {
						d_znbr_best = d_znbr;
						if (!ab_full_nbrhd) {
							ai_tour = i_nbr;
							break;	//accept first improvement found
						}//if
						else {
							i_hd1_best = i_hd1;
							i_hd2_best = i_hd2;
						}//else
					}//if
//				}//if
			}//i_hd2
			if (!ab_full_nbrhd && d_znbr_best < z)
				break;
		}//i_hd1
		if (ab_full_nbrhd && d_znbr_best < z) {
			d_znbr = gd_two_opt(ai_n, ad_dist, z, ai_tour, i_hd1_best, i_hd2_best, i_nbr);
			ai_tour = i_nbr;
		}//if
	} while (d_znbr_best + gd_tol < z);

	gi_itctr = i_itctr;
	return z;
}//gd_local_2opt_search

//ai_tour will be cleared/resized
//
double gd_nearest_neighbour(int ai_n, vector<vector<double>> &ad_dist, int i_i0, vector<int> &ai_tour) {
	int i = i_i0;
	int j;
	vector<bool> b_visited(ai_n);
	double d_min;
	int j_min;
	double z = 0;

	for (int k = 0; k < ai_n; ++k)
		b_visited[k] = false;

	ai_tour.resize(ai_n);
	ai_tour[0] = i;
	for (int k = 1; k < ai_n; ++k) {
		b_visited[i] = true;
		d_min = DBL_MAX;
		for (j = 0; j < ai_n; ++j) {
			if (!b_visited[j]) {
				if (ad_dist[i][j] < d_min) {
					d_min = ad_dist[i][j];
					j_min = j;
				}//if
			}//if
		}//j
		ai_tour[k] = j_min;
		z += d_min;
		i = j_min;
	}//k
	z += ad_dist[ai_tour[ai_n - 1]][ai_tour[0]];
	return z;
}//gd_nearest_neighbour

 //ai_tour will be cleared/resized
 //
double gd_random_tour(int ai_n, const vector<vector<double>> &ad_dist, vector<int> &ai_tour) {
	double z;

	vector<bool> b_used(ai_n);
	for (int i = 0; i < ai_n; ++i)
		b_used[i] = false;

	ai_tour.resize(ai_n);

	int n = ai_n - 1;
	int j, k, i_ctr;
	for (int i = 0; i < ai_n; ++i) {
		k = gi_rand(n);
		i_ctr = 0;
		for (j = 0; j < ai_n; ++j) {
			if (!b_used[j]) {
				if (i_ctr == k)
					break;
				++i_ctr;
			}//if
		}//j
		b_used[j] = true;
		ai_tour[i] = j;
		--n;
	}//i

	z = gd_calc_z(ai_n, ad_dist, ai_tour);

	return z;
}//gd_random_tour


///////////////////////////TEST STUBS///////////////////////////////

void g_test_rand_tour() {
	vector<CityPair> v_dat;
	vector<vector<double>> d_dist, d_xy;
	vector<int> i_tour;
	int n = 8;
	double z;

	g_gen_tsp_instance(n, v_dat, d_dist, d_xy);
	for (int i = 0; i < 200; ++i) {
		z = gd_random_tour(n, d_dist, i_tour);
		cout << z;
		for (int j = 0; j < n; ++j)
			cout << " " << i_tour[j];
		cout << endl;
	}//i
	cin >> n;
}//g_test_rand_tour

void g_test_heuristics() {
	vector<CityPair> v_dat;
	vector<vector<double>> d_dist, d_xy;
	vector<int> i_tour, i_nbr;
	int n = 8;
	i_nbr.resize(n);

	g_gen_tsp_instance(n, v_dat, d_dist, d_xy);
	//double z = gd_nearest_neighbour(n, d_dist, 0, i_tour);
	i_tour.resize(n);
	for (int i = 0; i < n; ++i)
		i_tour[i] = i;
	double z = gd_calc_z(n, d_dist, i_tour);

	double d_znbr, d_zcheck, d_zcheck2;
	d_zcheck = gd_calc_z(n, d_dist, i_tour);
	cout << z << "==" << d_zcheck << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (abs(i - j) > 1 && abs(i - j) < n - 1) {
				//d_znbr = gd_insertion(n, d_dist, z, i_tour, i, j, i_nbr);
				//d_znbr = gd_swap(n, d_dist, z, i_tour, i, j);
				d_znbr = gd_two_opt(n, d_dist, z, i_tour, i, j, i_nbr);
				d_zcheck = gd_calc_z(n, d_dist, i_nbr);
				//d_zcheck2 = gd_swap(n, d_dist, d_znbr, i_tour, i, j);	//undo the swap
				//cout << i << " " << j << " " << d_znbr << "==" << d_zcheck <<  "..." << d_zcheck2 << "==" << z << endl;
				cout << i << " " << j << " " << d_znbr << "==" << d_zcheck << endl;
				if (abs(d_znbr - d_zcheck) > 0.000005)
					i = i;
			}//if
		}//j
	}//i

	cin >> n;

}//g_test_heuristics


/*
int main() {
	g_test_heuristics();
	return 0;
}
	*/
