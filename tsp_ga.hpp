#include "tsp_nhoods.hpp"

class TSPxover : public XoverBase {
public:
	int mi_n;
	vector<vector<double>> &md_dist;

	vector<const TSPsoln*> m_p;
	vector<TSPsoln*> m_ch;
	vector<int> mi_template;
	vector<vector<bool>> mb_xother;
	vector<int> mi_idx;

	void m_xover(const Soln& a_p1, const Soln& a_p2, vector<Soln*> ap_o) {
		m_p[0] = dynamic_cast<const TSPsoln*>(&a_p1);
		m_p[1] = dynamic_cast<const TSPsoln*>(&a_p2);
		m_ch[0] = dynamic_cast<TSPsoln*> (ap_o[0]);
		m_ch[1] = dynamic_cast<TSPsoln*> (ap_o[1]);

		for (int i = 0; i < mi_n; ++i) {
			mb_xother[0][i] = false;
			mb_xother[1][i] = false;
		}//i

		for (int i = 0; i < mi_n; ++i)
			mi_template[i] = gi_rand(2);

		int i_c, i_p, i_xp;
		for (int i = 0; i < mi_n; ++i) {
			i_c = mi_template[i];
			i_xp = i_c == 0 ? 1 : 0;
			m_ch[i_c]->mi_x[i] = m_p[i_c]->mi_x[i];
			mb_xother[i_xp][m_p[i_xp]->mi_x[i]] = true;
		}//i

		int c;
		mi_idx[0] = 0;
		mi_idx[1] = 0;
		for (int i = 0; i < mi_n; ++i) {
			i_xp = mi_template[i];
			i_c = i_xp == 0 ? 1 : 0;
			while (!mb_xother[i_c][m_p[i_xp]->mi_x[mi_idx[i_c]]])
				++mi_idx[i_c];
			m_ch[i_c]->mi_x[i] = m_p[i_xp]->mi_x[mi_idx[i_c]];
			++mi_idx[i_c];
		}//i

		for (int p = 0; p < 2; ++p) {
			m_ch[p]->md_z = gd_calc_z(mi_n, md_dist, m_ch[p]->mi_x);
		}//p

	}//m_xover

	TSPxover(int ai_n, vector<vector<double>> &ad_dist) : 
		mi_n(ai_n), m_p(2), mi_template(ai_n), m_ch(2), mi_idx(2), mb_xother(2), md_dist(ad_dist) {
		
		for (int c = 0; c < 2; ++c) {
			mb_xother[c].resize(mi_n);
		}//c
	}//ctor
};//TSPxover