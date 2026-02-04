#pragma once

#include "BasicTSPheuristics.hpp"
#include "metaheuristics.hpp"

class TSPsoln : public Soln {
public:
	vector<int> mi_x;

	TSPsoln(int ai_n) : mi_x(ai_n), Soln(ai_n) {}
	TSPsoln(const vector<int> &ai_tour, double ad_z) : mi_x(ai_tour), Soln(ai_tour.size()) {
		md_z = ad_z;
	}
	
	Soln& operator=(const Soln &a_rhs) {	//instantiation of base class assignment operator
		const TSPsoln *p_rhs = dynamic_cast<const TSPsoln*>(&a_rhs);
		*this = *p_rhs;		//this calls the TSPsoln assignment oeprator below
		return *this;
	}//=

	TSPsoln& operator=(const TSPsoln &a_rhs) {
		if (mi_n != a_rhs.mi_n)
			mi_x.resize(mi_n);
		for (int i = 0; i < mi_n; ++i)
			mi_x[i] = a_rhs.mi_x[i];

		mi_n = a_rhs.mi_n;
		md_z = a_rhs.md_z;

		return *this;
	}//=
};//TSPsoln

class TSPnbrOp : public NbrOp {
public:
	const int mi_n;
	const vector<vector<double>> &md_dist;

	int mi_i;				//current cities/or tour positions in full neighbourhood enumeration
	int mi_j;

	int mi_i_prev;			//most recent pair of cities evaluated
	int mi_j_prev;
	double md_z_prev;

	int mi_i_held;			//pair of cities held for later build of neighbour
	int mi_j_held;
	double md_z_held;

	TSPnbrOp(const int ai_n, const vector<vector<double>> &ad_dist) : mi_n(ai_n), md_dist(ad_dist) {}

	void m_hold_nbr() {
		mi_i_held = mi_i_prev;
		mi_j_held = mi_j_prev;
		md_z_held = md_z_prev;
	}//m_hold_nbr

	virtual double md_nbr_z(const TSPsoln &a_x) = 0;		//calculate z-value for current mi_i and mi_j

	Soln* mp_new_soln() {
		return new TSPsoln(mi_n);
	}//mp_new_soln

	void m_randomise(Soln &a_x) {
		TSPsoln &x = dynamic_cast<TSPsoln&>(a_x);
		x.md_z = gd_random_tour(mi_n, md_dist, x.mi_x);
	}//m_randomise
};//TSPnbrOp

class Op2opt : public TSPnbrOp {
public:
	void m_reset_nhd() {
		mi_i = 0;
		mi_j = 2;
	}//m_reset_nhd

	double md_rand_nbr(const Soln &a_x, Soln* ap_nbr = 0) {
		int i = gi_rand(mi_n);
		int j = gi_rand(mi_n);
		while (abs(i - j) < 2 || abs(i - j) == mi_n - 1)
			j = gi_rand(mi_n);

		mi_i = i;
		mi_j = j;
		mi_i_prev = i;
		mi_j_prev = j;

		if (ap_nbr != 0) {
			m_hold_nbr();
			m_build_held_nbr(a_x, ap_nbr);
		}//ap_nbr

		return md_nbr_z(dynamic_cast<const TSPsoln&>(a_x));
	}//md_rand_nbr;

	double md_next_nbr(const Soln &a_x, bool &ab_stop, Soln* ap_nbr = 0) {
		double z = md_nbr_z(dynamic_cast<const TSPsoln&>(a_x));
		mi_i_prev = mi_i;
		mi_j_prev = mi_j;

		if (mi_j == mi_n - 1) {		//index to next neighbour
			if (mi_i == mi_n - 3) {
				m_reset_nhd();
				ab_stop = true;
				return z;
			}//if
			++mi_i;
			mi_j = mi_i + 2;
		}
		else
			++mi_j;
		ab_stop = false;
		return z;
	}//md_next_nbr

	void m_build_held_nbr(const Soln &a_x, Soln* ap_nbr) {
		const TSPsoln& x = dynamic_cast<const TSPsoln&>(a_x);
		TSPsoln *p_nbr = dynamic_cast<TSPsoln*>(ap_nbr);

		p_nbr->md_z = gd_two_opt(mi_n, md_dist, x.md_z, x.mi_x, mi_i_held, mi_j_held, p_nbr->mi_x);
	}//m_build_held_nbr

	double md_nbr_z(const TSPsoln &a_x) {
		double z = a_x.md_z;
		z -= mi_i > 0 ? md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_i]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];
		z -= mi_j > 0 ? md_dist[a_x.mi_x[mi_j - 1]][a_x.mi_x[mi_j]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];
		z += md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_j]];
		z += md_dist[a_x.mi_x[mi_i > 0 ? mi_i - 1 : mi_n - 1]][a_x.mi_x[mi_j > 0 ? mi_j - 1 : mi_n - 1]];

		return z;
	}//md_nbr_z

	Op2opt(const int ai_n, const vector<vector<double>> &ad_dist) : TSPnbrOp(ai_n, ad_dist) {}
};//Op2opt

class OpSwap : public TSPnbrOp {
public:
	double md_rand_nbr(const Soln &a_x, Soln* ap_nbr = 0) {
		mi_i = gi_rand(mi_n - 1);
		mi_j = mi_i + 1 + gi_rand(mi_n - 1 - (mi_i + 1));
		mi_i_prev = mi_i;
		mi_j_prev = mi_j;

		if (ap_nbr != 0) {
			m_hold_nbr();
			m_build_held_nbr(a_x, ap_nbr);
		}//ap_nbr

		return md_nbr_z(dynamic_cast<const TSPsoln&>(a_x));
	}//g_rand_nbr

	void m_reset_nhd() {
		mi_i = 0;
		mi_j = 1;
	}//m_reset_nhd

	double md_next_nbr(const Soln &a_x, bool &ab_stop, Soln* ap_nbr = 0) {
		double z = md_nbr_z(dynamic_cast<const TSPsoln&>(a_x));
		mi_i_prev = mi_i;
		mi_j_prev = mi_j;

		if (mi_j == mi_n - 1) {		//index to next neighbour
			if (mi_i == mi_n - 2) {
				m_reset_nhd();
				ab_stop = true;
				return z;
			}//if
			++mi_i;
			mi_j = mi_i + 1;
		}
		else
			++mi_j;
		ab_stop = false;
		return z;
	}//g_next_nbr

	void m_build_held_nbr(const Soln &a_x, Soln* ap_nbr) {
		const TSPsoln& x = dynamic_cast<const TSPsoln&>(a_x);
		TSPsoln *p_nbr = dynamic_cast<TSPsoln*>(ap_nbr);
		*p_nbr = a_x;
		p_nbr->md_z = gd_swap(mi_n, md_dist, p_nbr->md_z, p_nbr->mi_x, mi_i_held, mi_j_held);
	}//m_build_held_nbr

	double md_nbr_z(const TSPsoln &a_x) {
		double z;

		z =  a_x.md_z;

		if (mi_i + 1 == mi_j) {
			z -= mi_i > 0 ? md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_i]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];
			z -= md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_j]];
			z -= mi_j < mi_n - 1 ? md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_j + 1]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];

			z += mi_i > 0 ? md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_j]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[mi_j]];
			z += md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_i]];
			z += mi_j < mi_n - 1 ? md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_j + 1]] : md_dist[a_x.mi_x[mi_i]][a_x.mi_x[0]];
		}
		else if (mi_i + mi_n == mi_j + 1) {
			z -= md_dist[a_x.mi_x[mi_j - 1]][a_x.mi_x[mi_j]];
			z -= md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_i]];
			z -= md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_i + 1]];

			z += md_dist[a_x.mi_x[mi_j - 1]][a_x.mi_x[mi_i]];
			z += md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_j]];
			z += md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_i + 1]];
		}//if
		else {
			z -= mi_i > 0 ? md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_i]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];
			z -= md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_i + 1]];
			z -= md_dist[a_x.mi_x[mi_j - 1]][a_x.mi_x[mi_j]];
			z -= mi_j < mi_n - 1 ? md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_j + 1]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];

			z += md_dist[a_x.mi_x[mi_j - 1]][a_x.mi_x[mi_i]];
			z += mi_j < mi_n - 1 ? md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_j + 1]] : md_dist[a_x.mi_x[mi_i]][a_x.mi_x[0]];
			z += mi_i > 0 ? md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_j]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[mi_j]];
			z += md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_i + 1]];
		}//else

		return z;
	}//md_nbr_z

	OpSwap(const int ai_n, const vector<vector<double>> &ad_dist) : TSPnbrOp(ai_n, ad_dist) {}
};//OpSwap

class OpIns : public TSPnbrOp {
public:
	double md_rand_nbr(const Soln &a_x, Soln* ap_nbr = 0) {
		mi_i = gi_rand(mi_n - 1);
		mi_j = gi_rand(mi_n - 1);
		while (mi_i == mi_j || mi_i == mi_j + 1 || mi_i == mi_j - mi_n + 1)
			mi_j = gi_rand(mi_n - 1);
		mi_i_prev = mi_i;
		mi_j_prev = mi_j;

		if (ap_nbr != 0) {
			m_hold_nbr();
			m_build_held_nbr(a_x, ap_nbr);
		}//ap_nbr

		return md_nbr_z(dynamic_cast<const TSPsoln&>(a_x));
	}//g_rand_nbr

	void m_reset_nhd() {
		mi_i = 0;
		mi_j = 2;
	}//m_reset_nhd

	double md_next_nbr(const Soln &a_x, bool &ab_stop, Soln* ap_nbr = 0) {
		double z = md_nbr_z(dynamic_cast<const TSPsoln&>(a_x));
		mi_i_prev = mi_i;
		mi_j_prev = mi_j;

		++mi_j;
		while (mi_i == mi_j || mi_i == mi_j + 1 || mi_i == mi_j - mi_n + 1)
			++mi_j;
		if (mi_j >= mi_n) {
			if (mi_i == mi_n - 1) {
				m_reset_nhd();
				ab_stop = true;
				return z;
			}//if
			++mi_i;
			mi_j = 0;
			while (mi_i == mi_j || mi_i == mi_j + 1 || mi_i == mi_j - mi_n + 1)
				++mi_j;
		}//if
		ab_stop = false;
		return z;
	}//md_next_nbr

	void m_build_held_nbr(const Soln &a_x, Soln* ap_nbr) {
		const TSPsoln& x = dynamic_cast<const TSPsoln&>(a_x);
		TSPsoln *p_nbr = dynamic_cast<TSPsoln*>(ap_nbr);
		p_nbr->md_z = gd_insertion(mi_n, md_dist, x.md_z, x.mi_x, mi_i_held, mi_j_held, p_nbr->mi_x);
	}//m_build_held_nbr

	double md_nbr_z(const TSPsoln &a_x) {
		double z = a_x.md_z;
		z -= mi_i > 0 ? md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_i]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];				//-d_R-1,R
		z -= mi_i < mi_n - 1 ? md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_i + 1]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];		//-d_R,R+1
		z -= mi_j < mi_n - 1 ? md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_j + 1]] : md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[0]];		//-d_I,I+1
		z += mi_i == 0 ? md_dist[a_x.mi_x[mi_n - 1]][a_x.mi_x[1]] :
			(mi_i == mi_n - 1 ? md_dist[a_x.mi_x[mi_n - 2]][a_x.mi_x[0]] : md_dist[a_x.mi_x[mi_i - 1]][a_x.mi_x[mi_i + 1]]);	//+d_R-1,R+1
		z += md_dist[a_x.mi_x[mi_j]][a_x.mi_x[mi_i]];																		//+d_I,R
		z += mi_j < mi_n - 1 ? md_dist[a_x.mi_x[mi_i]][a_x.mi_x[mi_j + 1]] : md_dist[a_x.mi_x[mi_i]][a_x.mi_x[0]];			//+d_R,I+1
		return z;
	}//md_nbr_z

	OpIns(const int ai_n, const vector<vector<double>> &ad_dist) : TSPnbrOp(ai_n, ad_dist) {}
};//OpIns


class OpMulti : public TSPnbrOp {
public:
	vector<NbrOp*> mv_nhd;
	vector<double> md_pnhd;

	int mi_nhd;
	int mi_prev_nhd;
	int mi_held_nhd;

	void m_add_nhd(NbrOp &a_nhd, double ad_p) {
		mv_nhd.push_back(&a_nhd);
		md_pnhd.push_back(ad_p);
	}//m_add_nhd

	double md_rand_nbr(const Soln &a_x, Soln* ap_nbr = 0) {
		double d_rnd = md_rand();
		double d_sum = 0;
		for (int k = 0; k < mv_nhd.size() - 1; ++k) {
			d_sum += md_pnhd[k];
			if (d_rnd <= d_sum) {
				mi_prev_nhd = k;
				return mv_nhd[k]->md_rand_nbr(a_x);
			}//if
		}//k
		mi_prev_nhd = mv_nhd.size() - 1;

		if (ap_nbr != 0) {
			m_hold_nbr();
			m_build_held_nbr(a_x, ap_nbr);
		}//ap_nbr

		return mv_nhd.back()->md_rand_nbr(a_x);
	}//md_rand_nbr

	double md_next_nbr(const Soln &a_x, bool &ab_stop, Soln* ap_nbr = 0) {
		bool b_stop;
		double z = mv_nhd[mi_nhd]->md_next_nbr(a_x, b_stop);
		mi_prev_nhd = mi_nhd;

		if (b_stop)
			++mi_nhd;
		ab_stop = (mi_nhd >= mv_nhd.size());

		return z;
	}//md_next_nbr

	void m_reset_nhd() {
		mi_nhd = 0;
		for (int i = 0; i < mv_nhd.size(); ++i)
			mv_nhd[i]->m_reset_nhd();
	}//m_reset_nhd

	void m_build_held_nbr(const Soln &a_x, Soln* ap_nbr) {
		mv_nhd[mi_held_nhd]->m_build_held_nbr(a_x, ap_nbr);
	}//m_build_held_nbr

	void m_hold_nbr() {
		mi_held_nhd = mi_prev_nhd;
		mv_nhd[mi_held_nhd]->m_hold_nbr();
	}//m_hold_nbr

	double md_nbr_z(const TSPsoln &a_x) {
		return 0;
	}//md_nbr_z

	OpMulti(const int ai_n, const vector<vector<double>> &ad_dist) : TSPnbrOp(ai_n, ad_dist) {}
};//OpMulti
