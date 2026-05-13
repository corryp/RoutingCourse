#include "tsp_nhoods.hpp"

class TSPnhdMem : public LongMemBase {
public:
	int mi_n;
	TSPnbrOp &m_nhd;

	double md_wt;	//weighting of the penalty relative to average of two most recent z values

	vector<vector<vector<int>>> mi_usectr;	//[nhd_id][i][j] count of how many times an (i,j) move was used per neighbourhood type
	double md_z1, md_z2;					//tracks two most recent z values

	double md_freq_penalty(int ai_iter, double ad_z, const NbrOp &a_nhd) {
		double d_zavg = 0.5 * (md_z1 + md_z2);
		double d_freq = mi_usectr[m_nhd.mi_k_prev][m_nhd.mi_i_prev][m_nhd.mi_j_prev] / (1.0 * ai_iter);
		return md_wt * d_zavg * d_freq;
	}//md_freq_penalty

	void m_update_freq(int ai_iter, double ad_z, const NbrOp &a_nhd) {
		++mi_usectr[m_nhd.mi_k_held][m_nhd.mi_i_held][m_nhd.mi_j_held];
		md_z2 = md_z1;
		md_z1 = ad_z;
	}//m_update_freq

	void m_reset() {
		for (int k = 0; k < NHD_COUNT; ++k) {
			for (int i = 0; i < mi_n; ++i) {
				for (int j = 0; j < mi_n; ++j)
					mi_usectr[k][i][j] = 0;
			}//i
		}//k
		md_z1 = 0;
		md_z2 = 0;
	}//m_reset

	TSPnhdMem(TSPnbrOp &a_nhd, double ad_wt) : mi_usectr(NHD_COUNT), mi_n(a_nhd.mi_n), m_nhd(a_nhd), md_wt(ad_wt) {
		for (int k = 0; k < NHD_COUNT; ++k) {
			mi_usectr[k].resize(mi_n);
			for (int i = 0; i < mi_n; ++i)
				mi_usectr[k][i].resize(mi_n);
		}//k
		m_reset();
	}//ctor
};//TSPnhdMem