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

#include "LogFile.hpp"

const double gd_mhtol = 0.000000005;

class Soln {
public:
	int mi_n;										//number of degrees of freedom in the solution
	double md_z;									//objective value

	Soln(int ai_n) : mi_n(ai_n) {}
	virtual Soln& operator=(const Soln &a_rhs) = 0;	//Note: dynamic cast will be required in override of this
};//Soln

class NbrOp {
public:
	virtual double md_rand_nbr(const Soln &a_x, Soln* ap_nbr = 0) = 0;	//generates a random neighbour (ensure disposal of ap_nbr)
	virtual double md_next_nbr(const Soln &a_x, bool &ab_stop, Soln* ap_nbr = 0) = 0;	//ab_stop set to true when end of neighbourhood is reached
	virtual void m_hold_nbr() {}			//if required, this procedure stores data required to generate the most recent neighbour
	virtual void m_build_held_nbr(const Soln &a_x, Soln* ap_nbr) {}	//if required, this procedure builds the neighbor which was held

	virtual Soln* mp_new_soln() { return 0; }	//provides mechanism to create solution of a derived class
	virtual void m_reset_nhd() {};									//used to force a reset of a neighbourhood evaluation (works with m_next_nbr)
	virtual void m_randomise(Soln &a_x) {};			//used by GA to generate an initial random population

	//helper functions
	//
	virtual double md_rand() { return 1.0 * rand() / RAND_MAX; }	//creates [0,1] uniform random variable
	virtual int mi_rand(int ai_n) {	//creates integer random variable between 0 and ai_n - 1
		if (ai_n == 0) return 0;
		int i = ai_n * md_rand();
		i = (i < ai_n ? i : ai_n - 1);	//for the extremely rare occurence of sampling a value of exactly 1
		return i;
	}//mi_rand

	NbrOp() { m_reset_nhd(); }
};//NbrOp


////////////////////////////////////////////LOCAL SEARCH...
struct LSsummary {
	double md_z0;
	double md_zbest;
	int mi_nhd_ctr;
	int mi_soln_ctr;
};//LSsummary

//TODO: add logging
//
LSsummary gx_local_search(bool ab_fullnhd, NbrOp &a_nbrhd, Soln &a_x, Soln &a_xbest, Logger *ap_log = 0) {
	LSsummary stats;
	stats.md_z0 = a_x.md_z;
	stats.mi_nhd_ctr = 0;
	stats.mi_soln_ctr = 0;

	Soln* p_nbr = a_nbrhd.mp_new_soln();
	a_xbest = a_x;
	bool b_stop, b_best;
	double z, d_zbest;
	do {
		++stats.mi_nhd_ctr;
		b_stop = false;
		b_best = false;
		a_nbrhd.m_reset_nhd();
		d_zbest = a_xbest.md_z;
		while (!b_stop) {
			++stats.mi_soln_ctr;
			z = a_nbrhd.md_next_nbr(a_xbest, b_stop);
			if (z + gd_mhtol < d_zbest) {
				b_best = true;
				d_zbest = z;
				a_nbrhd.m_hold_nbr();
				if (!ab_fullnhd)
					break;
			}//if
		}//wend
		if (b_best) {
			a_nbrhd.m_build_held_nbr(a_xbest, p_nbr);
			a_xbest = *p_nbr;
		}//if
	} while (b_best);
	stats.md_zbest = a_xbest.md_z;
	delete p_nbr;
	return stats;
}//gx_local_search


////////////////////////////////////////////SIMULATED ANNEALING...

struct SActrl {
	double md_init_accept_rate;	//this is used to calculate the initial T0
	int mi_accept_limit_coef;	//next temp step if mi_accept_limit_coef * N nbrs accepted
	int mi_nbr_limit_coef;		//next temp step if mi_nbr_limit_coef * N nbrs attempted
	double md_tfact;
	int mi_stale_step_lim;		//algorithm will stop after this many steps without improvment;
	int mi_max_iter;			//maximum number of iterations (0 = no limit)

	//logfile settings
	bool mb_every_it;			//true => every iteration logged : false => log end of temp step only
	int mi_log_gap;				//log every Nth record to file and screen (1 = every record)
};//SActrl

struct SAsummary {
	double md_z0;
	double md_zbest;
	double md_t0;
	double md_tfinal;
	int mi_tstep_ctr;
	int mi_nbr_ctr;
	int mi_mv_ctr;
};//SAsummary

struct SAtrace {
	int mi_totit;
	int mi_tstep;
	double md_t;
	int mi_stepit;
	double md_zbest;
	double md_z;
	double md_znbr;
	int mi_accpt;
	double md_paccpt;
};//SAtrace

void g_init_sa_log(Logger *ap_log) {
	if (ap_log == 0)
		return;

	ap_log->m_f("totit", "totit");
	ap_log->m_f("tstep", "tstep");
	ap_log->m_f("temp", "temp");
	ap_log->m_f("stepit", "stepit");
	ap_log->m_f("zbest", "zbest");
	ap_log->m_f("z", "z");
	ap_log->m_f("znbr", "znbr");
	ap_log->m_f("accpt", "accpt");
	ap_log->m_f("paccpt", "paccpt");
	ap_log->m_write_header();
}//g_init_sa_log

void g_write_sa_log(Logger *ap_log, const SActrl &a_ctrl, SAtrace a_dat) {
	if (ap_log == 0)
		return;

	ap_log->m_log("totit", a_dat.mi_totit);
	ap_log->m_log("tstep", a_dat.mi_tstep);
	ap_log->m_log("temp", a_dat.md_t);
	ap_log->m_log("stepit", a_dat.mi_stepit);
	ap_log->m_log("zbest", a_dat.md_zbest);
	ap_log->m_log("z", a_dat.md_z);
	ap_log->m_log("znbr", a_dat.md_znbr);
	ap_log->m_log("accpt", a_dat.mi_accpt);
	ap_log->m_log("paccpt", a_dat.md_paccpt);
	ap_log->m_write_record();
}//g_write_sa_log

//a_x should contain the initial solution, and will be overwritten for use as incumbent solution
//
SAsummary gx_simulated_annealing(const SActrl &a_ctrl, NbrOp &a_nbrhd, Soln &a_x, Soln &a_xbest, Logger *ap_log = 0) {
	g_init_sa_log(ap_log);
	//initialisation: calculation of T0
	//
	double d_z;
	double d_t0;
	double d_deltaE_avg = 0;
	for (int k = 0; k < 100; ++k) {
		d_z = a_nbrhd.md_rand_nbr(a_x);
		d_deltaE_avg += 0.01 * (d_z - a_x.md_z);
	}//k
	d_t0 = -abs(d_deltaE_avg) / log(a_ctrl.md_init_accept_rate);	//NB: abs in case most moves improving from poor starting solution

	SAsummary stats;
	stats = { a_x.md_z, a_x.md_z, d_t0, 0, 0, 0 };

	Soln* p_nbr = a_nbrhd.mp_new_soln();
	int i_nbr_ctr, i_accpt_ctr, i_accpt, i_totit = 0, i_stale_ctr = 0;
	int n = a_x.mi_n;
	double p_accpt;
	double d_temp = d_t0;
	a_xbest = a_x;
	while (i_stale_ctr < a_ctrl.mi_stale_step_lim && (a_ctrl.mi_max_iter == 0 || i_totit < a_ctrl.mi_max_iter)) {
		i_accpt_ctr = 0;
		i_nbr_ctr = 0;
		d_z = a_x.md_z;
		while (i_accpt_ctr < a_ctrl.mi_accept_limit_coef * n && i_nbr_ctr < a_ctrl.mi_nbr_limit_coef * n) {
			d_z = a_nbrhd.md_rand_nbr(a_x);

			if (d_z < a_x.md_z)
				p_accpt = 1.0;
			else
				p_accpt = exp(-(d_z - a_x.md_z) / d_temp);
			if (p_accpt == 1.0 || a_nbrhd.md_rand() < p_accpt) {
				a_nbrhd.m_hold_nbr();
				a_nbrhd.m_build_held_nbr(a_x, p_nbr);
				a_x = *p_nbr;
				if (a_x.md_z < a_xbest.md_z)
					a_xbest = a_x;
				++i_accpt_ctr;
				i_accpt = 1;
				++stats.mi_mv_ctr;
			}//if
			else
				i_accpt = 0;
			++i_nbr_ctr;
			++stats.mi_nbr_ctr;

			++i_totit;
			if (a_ctrl.mb_every_it || i_totit == 1)
				g_write_sa_log(ap_log, a_ctrl, {i_totit, stats.mi_tstep_ctr, d_temp, i_nbr_ctr, a_xbest.md_z, d_z, d_z, i_accpt, (1.0 * i_accpt_ctr) / i_nbr_ctr });
		}//wend
		++stats.mi_tstep_ctr;
		if (i_accpt_ctr == 0)
			++i_stale_ctr;
		else
			i_stale_ctr = 0;
		if (!a_ctrl.mb_every_it)
			g_write_sa_log(ap_log, a_ctrl, { i_totit, stats.mi_tstep_ctr, d_temp, i_nbr_ctr, a_xbest.md_z, d_z, d_z, i_accpt, (1.0 * i_accpt_ctr) / i_nbr_ctr });
		d_temp *= a_ctrl.md_tfact;
	}//wend

	delete p_nbr;
	stats.md_zbest = a_xbest.md_z;
	stats.md_tfinal = d_temp;
	return stats;
}//gx_simulated_annealing

///////////////////Tabu search...

class TabuListBase {
public:
	virtual void m_add_tabu(int ai_iter, const Soln &a_x, const NbrOp &a_nhd) = 0;
	virtual bool mb_is_tabu(int ai_iter, double ad_z, const NbrOp &a_nhd) { return false; }
	virtual void m_reset() {}
};//TabuListBase

class TabuListZhash : public TabuListBase {
public:
	int mi_hlen;					//length of hash table
	int mi_tl;						//number of iterations a solution remains tabu

	vector<int> mi_end_tabu;		//first iteration where solutin is no-longer tabu

	TabuListZhash(int ai_hlen, int ai_tl) : mi_hlen(ai_hlen), mi_end_tabu(ai_hlen), mi_tl(ai_tl) {
		m_reset();
	}//TabuListZhash

	void m_add_tabu(int ai_iter, const Soln &a_x, const NbrOp &a_nhd) {
		int k = mi_hash(a_x.md_z);
		mi_end_tabu[k] = ai_iter + mi_tl;
	}//m_add_tabu

	bool mb_is_tabu(int ai_iter, double ad_z, const NbrOp &a_nhd) {
		int k = mi_hash(ad_z);
		return mi_end_tabu[k] >= ai_iter;
	}//mb_is_tabu

	void m_reset() {
		for (int i = 0; i < mi_hlen; ++i)
			mi_end_tabu[i] = 0;
	}//m_reset

private:
	int mi_hash(double ad_z) {
		return (int)(mi_hlen * ad_z) % mi_hlen;
	}//md_hash
};//TabuListZhash

class LongMemBase {
public:
	virtual double md_freq_penalty(int ai_iter, double ad_z, const NbrOp &a_nhd) { return 0; };
	virtual void m_update_freq(int ai_iter, double ad_z, const NbrOp &a_nhd) {}
	virtual void m_reset() {}
};//LongMemBase

struct TSctrl {
	int mi_itmax;		//max number of iterations
	int mi_log_gap;		//log every Nth record to file and screen (1 = every record)
};//TSctrl

struct TSsummary {
	double md_z0;
	double md_zbest;
	int mi_nhd_ctr;
};//TSsummary

struct TStrace {
	int mi_nhd_ctr;
	int mi_cnddt_ctr;
	int mi_tabu_ctr;
	double md_zbest;
	double md_z;
};//TStrace

void g_init_ts_log(Logger *ap_log) {
	if (ap_log == 0)
		return;

	ap_log->m_f("iter", "iter");
	ap_log->m_f("nbrn", "candidates");
	ap_log->m_f("tabu", "tabu count");
	ap_log->m_f("zbest", "zbest");
	ap_log->m_f("z", "z");
	ap_log->m_write_header();
}//g_init_ts_log

void g_write_ts_log(Logger *ap_log, TStrace a_dat) {
	if (ap_log == 0)
		return;

	ap_log->m_log("iter", a_dat.mi_nhd_ctr);
	ap_log->m_log("nbrn", a_dat.mi_cnddt_ctr);
	ap_log->m_log("tabu", a_dat.mi_tabu_ctr);
	ap_log->m_log("zbest", a_dat.md_zbest);
	ap_log->m_log("z", a_dat.md_z);
	ap_log->m_write_record();
}//g_write_ts_log

TSsummary gx_tabu_search(const TSctrl &a_ctrl, TabuListBase &a_tl, LongMemBase *ap_mem, NbrOp &a_nbrhd, Soln &a_x, Soln &a_xbest, Logger *ap_log = 0) {
	g_init_ts_log(ap_log);
	TSsummary stats = { a_x.md_z, a_x.md_z, 0 };

	Soln* p_nbr = a_nbrhd.mp_new_soln();
	a_nbrhd.m_reset_nhd();
	a_tl.m_reset();
	LongMemBase *p_mem = (ap_mem == 0 ? new LongMemBase : ap_mem);	//just populates with a do-nothing empty shell
	p_mem->m_reset();
	a_xbest = a_x;

	bool b_stop, b_move;
	double d_z, d_zpen, d_zbest;
	TStrace trc = {0, 0, 0, a_x.md_z, a_x.md_z};
	g_write_ts_log(ap_log, trc);
	a_xbest = a_x;
	for (trc.mi_nhd_ctr = 1; trc.mi_nhd_ctr <= a_ctrl.mi_itmax; ++trc.mi_nhd_ctr) {	//main loop
		trc.mi_cnddt_ctr = 0;
		trc.mi_tabu_ctr = 0;
		++stats.mi_nhd_ctr;
		b_stop = false;
		b_move = false;
		d_zbest = DBL_MAX;
		while (!b_stop) {	//loop to evaluate neighbourhood
			d_z = a_nbrhd.md_next_nbr(a_x, b_stop);
			if (!a_tl.mb_is_tabu(trc.mi_nhd_ctr, d_z, a_nbrhd)) {
				++trc.mi_cnddt_ctr;
				d_zpen = d_z + p_mem->md_freq_penalty(trc.mi_nhd_ctr, d_z, a_nbrhd);
				if (d_zpen + gd_mhtol < d_zbest) {		//test for best neighbour found
					a_nbrhd.m_hold_nbr();
					d_zbest = d_zpen;
					b_move = true;
				}//if
			}//if
			else
				++trc.mi_tabu_ctr;
		}//wend
		if (b_move) {	//build solution for the best neighbour found
			a_nbrhd.m_build_held_nbr(a_x, p_nbr);
			a_x = *p_nbr;
			a_tl.m_add_tabu(trc.mi_nhd_ctr, a_x, a_nbrhd);	//update tabu list
			p_mem->m_update_freq(trc.mi_nhd_ctr, a_x.md_z, a_nbrhd);
			if (a_x.md_z + gd_mhtol < a_xbest.md_z)	//check for best known solution
				a_xbest = a_x;
		}//if

		trc.md_z = a_x.md_z;
		trc.md_zbest = a_xbest.md_z;
		g_write_ts_log(ap_log, trc);
	}//trc.mi_nhd_ctr

	delete p_nbr;
	if (ap_mem == 0)
		delete p_mem;

	stats.md_zbest = a_xbest.md_z;
	return stats;
}//gx_tabu_search

////////////////////////////////Genetic Algorithm....

struct GAtrace {
	int mi_gen_ctr;
	double md_zbest;
	double md_zmin;
	double md_zavg;
	double md_zmax;
};//GAtrace

struct GActrl {
	int mi_max_gens;
	int mi_pop_size;
	int mi_tourn_size;
	int mi_n_xover;		//generally assumed 2 offspring per xover
	double md_pmutate;
	int mi_elitist;		//1 => elitist strategy used, => 0 otherwise
	int mi_log_gap;		//log every Nth record to file and screen (1 = every record)
};//GActrl

class XoverBase {
public:
	virtual void m_xover(const Soln& a_p1, const Soln& a_p2, vector<Soln*> ap_o) = 0;
};//XoverBase

void g_init_ga_log(Logger *ap_log) {
	if (ap_log == 0)
		return;

	ap_log->m_f("iter", "iter");
	ap_log->m_f("zbest", "zbest");
	ap_log->m_f("zmin", "zmin");
	ap_log->m_f("zavg", "zavg");
	ap_log->m_f("zmax", "zmax");
	ap_log->m_write_header();
}//g_init_ga_log

void g_write_ga_log(Logger *ap_log, GAtrace a_dat) {
	if (ap_log == 0)
		return;

	ap_log->m_log("iter", a_dat.mi_gen_ctr);
	ap_log->m_log("zbest", a_dat.md_zbest);
	ap_log->m_log("zmin", a_dat.md_zmin);
	ap_log->m_log("zavg", a_dat.md_zavg);
	ap_log->m_log("zmax", a_dat.md_zmax);
	ap_log->m_write_record();
}//g_write_ga_log

GAtrace gx_genetic_alg(const GActrl &a_ctrl, NbrOp &a_mutt, XoverBase &a_xover, Soln &a_xbest, Logger *ap_log = 0) {
	GAtrace trc = { 0, DBL_MAX, DBL_MAX, 0, -DBL_MAX };
	Soln* p_nbr = a_mutt.mp_new_soln();
	g_init_ga_log(ap_log);
	vector<Soln*> popln(a_ctrl.mi_pop_size);

	//create initial population
	//
	int i_kbest;
	for (int k = 0; k < a_ctrl.mi_pop_size; ++k) {
		popln[k] = a_mutt.mp_new_soln();
		a_mutt.m_randomise(*popln[k]);
		if (popln[k]->md_z < trc.md_zbest)
			i_kbest = k;
		trc.md_zavg += popln[k]->md_z / a_ctrl.mi_pop_size;
		trc.md_zmin = popln[k]->md_z < trc.md_zmin ? popln[k]->md_z : trc.md_zmin;
		trc.md_zbest = popln[k]->md_z < trc.md_zbest ? popln[k]->md_z : trc.md_zbest;
		trc.md_zmax = popln[k]->md_z > trc.md_zmax ? popln[k]->md_z : trc.md_zmax;
	}//k
	a_xbest = *(popln[i_kbest]);

	vector<int> i_tourn(a_ctrl.mi_tourn_size);
	vector<int> i_parents(2);	//winner of the two tournaments to decide the two parents
	vector<Soln*> p_offspring(2);
	p_offspring[0] = a_mutt.mp_new_soln();
	p_offspring[1] = a_mutt.mp_new_soln();

	bool b_notok;
	int i_soln, i_winner;
	double md_zwinner;
	vector<double> d_zunfit(2);
	vector<int> i_unfit(2);
	for (trc.mi_gen_ctr = 0; trc.mi_gen_ctr < a_ctrl.mi_max_gens; ++trc.mi_gen_ctr) {
		for (int k = 0; k < a_ctrl.mi_n_xover; ++k) {
			//determine parents for the tournament
			//
			i_parents[0] = -1;
			for (int p = 0; p < 2; ++p) {
				//build the next tournament
				//
				md_zwinner = DBL_MAX;
				for (int i = 0; i < a_ctrl.mi_tourn_size; ++i) {
					do {	//randomly pick tournament member, making sure it's not already selected
						i_soln = a_mutt.mi_rand(a_ctrl.mi_pop_size);
						b_notok = i_soln == i_parents[0];
						for (int j = 0; j < i; ++j)
							b_notok = b_notok || i_tourn[j] == i_soln;
					} while (b_notok);
					i_tourn[i] = i_soln;
					if (popln[i_soln]->md_z < md_zwinner) {
						md_zwinner = popln[i_soln]->md_z;
						i_winner = i_soln;
					}//if
				}//i
				i_parents[p] = i_winner;
			}//p (parent selection loop)
			a_xover.m_xover(*(popln[i_parents[0]]), *(popln[i_parents[1]]), p_offspring);

			d_zunfit[0] = -DBL_MAX;
			d_zunfit[1] = -DBL_MAX;
			for (int i = 0; i < a_ctrl.mi_pop_size; ++i) {	//determine two most unfit individuals to replace
				for (int p = 0; p < 2; ++p) {
					if (popln[i]->md_z - gd_mhtol > d_zunfit[p]) {
						if (p == 0) {
							i_unfit[1] = i_unfit[0];
							d_zunfit[1] = d_zunfit[0];
						}//if
						d_zunfit[p] = popln[i]->md_z;
						i_unfit[p] = i;
						break;
					}//if
				}//p
			}//i

			for (int p = 0; p < 2; ++p) {	//replace worst individuals with the new offspring
				if (p_offspring[p]->md_z < trc.md_zbest) {
					trc.md_zbest = p_offspring[p]->md_z;
					a_xbest = *(p_offspring[p]);
				}//if
				*(popln[i_unfit[p]]) = *(p_offspring[p]);
			}//p
		}//k (xovers loop)

		//randomly apply mutations to population
		//
		trc.md_zavg = 0;
		trc.md_zmin = DBL_MAX;
		trc.md_zmax = -DBL_MAX;
		d_zunfit[0] = -DBL_MAX;
		d_zunfit[1] = -DBL_MAX;
		for (int i = 0; i < a_ctrl.mi_pop_size; ++i) {
			if (a_mutt.md_rand() < a_ctrl.md_pmutate) {
				a_mutt.md_rand_nbr(*(popln[i]), p_nbr);
				*(popln[i]) = *p_nbr;
				if (p_nbr->md_z < trc.md_zbest) {
					trc.md_zbest = p_nbr->md_z;
					a_xbest = *p_nbr;
				}//if
			}//if
			trc.md_zavg += popln[i]->md_z / a_ctrl.mi_pop_size;
			trc.md_zmin = popln[i]->md_z < trc.md_zmin ? popln[i]->md_z : trc.md_zmin;
			trc.md_zmax = popln[i]->md_z > trc.md_zmax ? popln[i]->md_z : trc.md_zmax;

			if (a_ctrl.mi_elitist) {
				for (int p = 0; p < 2; ++p) {	//TODO: put this code in a function
					if (popln[i]->md_z - gd_mhtol > d_zunfit[p]) {
						if (p == 0) {
							i_unfit[1] = i_unfit[0];
							d_zunfit[1] = d_zunfit[0];
						}//if
						d_zunfit[p] = popln[i]->md_z;
						i_unfit[p] = i;
						break;
					}//if
				}//p
			}//if
		}//i
		if (a_ctrl.mi_elitist) {
			*(popln[i_unfit[0]]) = a_xbest;		//replace worst with best known solution
			trc.md_zmax = d_zunfit[1];			//update zmax
		}//if

		g_write_ga_log(ap_log, trc);
	}//trc.mi_gen_ctr

	for (int k = 0; k < a_ctrl.mi_pop_size; ++k)
		delete popln[k];
	delete p_offspring[0];
	delete p_offspring[1];
	delete p_nbr;
	return trc;
}//gx_genetic_alg

////////////////////////////////Variable Neighbourhood Search....

struct VNSctrl {
	int mi_kmax;			//max number of neighbourhoods to cycle through
	int mi_max_iter;		//maximum number of VNS iterations (0 = no limit)
	int mi_shake_steps;		//number of random moves per shake step

	//logfile settings
	bool mb_every_it;		//true => every k-step logged : false => log end of iteration only
	int mi_log_gap;			//log every Nth record to file and screen (1 = every record)
};//VNSctrl

struct VNSsummary {
	double md_z0;
	double md_zbest;
	int mi_iter_ctr;
	int mi_ls_ctr;
};//VNSsummary

struct VNStrace {
	int mi_iter;
	int mi_k;
	double md_zbest;
	double md_z;
	double md_zshake;
	double md_zls;
};//VNStrace

void g_init_vns_log(Logger *ap_log) {
	if (ap_log == 0)
		return;

	ap_log->m_f("iter", "iter");
	ap_log->m_f("k", "k");
	ap_log->m_f("zbest", "zbest");
	ap_log->m_f("z", "z");
	ap_log->m_f("zshake", "zshake");
	ap_log->m_f("zls", "zls");
	ap_log->m_write_header();
}//g_init_vns_log

void g_write_vns_log(Logger *ap_log, VNStrace a_dat) {
	if (ap_log == 0)
		return;

	ap_log->m_log("iter", a_dat.mi_iter);
	ap_log->m_log("k", a_dat.mi_k);
	ap_log->m_log("zbest", a_dat.md_zbest);
	ap_log->m_log("z", a_dat.md_z);
	ap_log->m_log("zshake", a_dat.md_zshake);
	ap_log->m_log("zls", a_dat.md_zls);
	ap_log->m_write_record();
}//g_write_vns_log

//a_x should contain the initial solution, and will be overwritten for use as incumbent solution
//a_nhds is an ordered vector of neighbourhood operators (determines VNS cycling order)
//
VNSsummary gx_basic_vns(const VNSctrl &a_ctrl, vector<NbrOp*> &a_nhds, Soln &a_x, Soln &a_xbest, Logger *ap_log = 0) {
	g_init_vns_log(ap_log);
	VNSsummary stats = { a_x.md_z, a_x.md_z, 0, 0 };

	int i_kmax = a_ctrl.mi_kmax < (int)a_nhds.size() ? a_ctrl.mi_kmax : (int)a_nhds.size();
	Soln* p_xshake = a_nhds[0]->mp_new_soln();
	Soln* p_xls = a_nhds[0]->mp_new_soln();
	Soln* p_nbr = a_nhds[0]->mp_new_soln();
	a_xbest = a_x;

	VNStrace trc = { 0, 0, a_x.md_z, a_x.md_z, 0, 0 };
	g_write_vns_log(ap_log, trc);

	for (trc.mi_iter = 1; a_ctrl.mi_max_iter == 0 || trc.mi_iter <= a_ctrl.mi_max_iter; ++trc.mi_iter) {
		int k = 0;
		while (k < i_kmax) {
			//shake: apply shake_steps random moves using neighbourhood k
			//
			*p_xshake = a_x;
			for (int s = 0; s < a_ctrl.mi_shake_steps; ++s) {
				a_nhds[k]->md_rand_nbr(*p_xshake, p_nbr);
				*p_xshake = *p_nbr;
			}//s
			trc.md_zshake = p_xshake->md_z;

			//local search: first improvement using neighbourhood k
			//
			gx_local_search(false, *a_nhds[k], *p_xshake, *p_xls);
			++stats.mi_ls_ctr;
			trc.md_zls = p_xls->md_z;

			//neighbourhood change
			//
			trc.mi_k = k;
			if (p_xls->md_z + gd_mhtol < a_x.md_z) {
				a_x = *p_xls;
				k = 0;		//improvement found: reset to first neighbourhood
				if (a_x.md_z + gd_mhtol < a_xbest.md_z)
					a_xbest = a_x;
			}//if
			else
				++k;		//no improvement: try next neighbourhood

			trc.md_z = a_x.md_z;
			trc.md_zbest = a_xbest.md_z;
			if (a_ctrl.mb_every_it)
				g_write_vns_log(ap_log, trc);
		}//wend (k loop)
		++stats.mi_iter_ctr;
		if (!a_ctrl.mb_every_it)
			g_write_vns_log(ap_log, trc);
	}//trc.mi_iter

	delete p_xshake;
	delete p_xls;
	delete p_nbr;

	stats.md_zbest = a_xbest.md_z;
	return stats;
}//gx_basic_vns