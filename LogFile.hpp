#pragma once

#include<iostream>
#include<fstream>
#include<string>
#include<vector>

using namespace std;

class LogField {
public:
	string mc_id;
	string mc_title;
	string mc_val;
	bool mb_to_screen;
	LogField(string ac_id, string ac_title, bool ab_to_screen) : mc_id(ac_id), mc_title(ac_title), mb_to_screen(ab_to_screen) {}
};//LogField

class Logger {
private:
	ofstream &mf;
	vector<LogField> mv_cols;

	int mi_rec;
	vector<LogField>::iterator m_it;
public:
	bool mb_to_screen;
	bool mb_to_file;
	int mi_screen_gap;

	Logger(ofstream &af, bool ab_to_screen = false, int ai_screen_gap = 1) : mf(af), mb_to_file(true), mb_to_screen(ab_to_screen), mi_screen_gap(ai_screen_gap) {}

	void m_set_options(bool ab_to_screen, bool ab_to_file, int ai_screen_gap) {
		mb_to_screen = ab_to_screen;
		mb_to_file = ab_to_file;
		mi_screen_gap = ai_screen_gap;
	}//m_set_options

	void m_f(string ac_id, string ac_title) {
		mv_cols.emplace_back(ac_id, ac_title, true);
	}//m_add_field

	void m_f(string ac_id, string ac_title, bool ab_to_screen) {
		mv_cols.emplace_back(ac_id, ac_title, ab_to_screen);
	}//m_add_field


	void m_write_header() {
		for (m_it = mv_cols.begin(); m_it != mv_cols.end(); m_it++) {
			if (m_it != mv_cols.begin()) mf << ",";
			mf << m_it->mc_title;
		}//m_it
		mf << endl;
		mi_rec = 0;
		if (mb_to_screen) m_write_screen_header();
	}//m_write_header

	void m_write_screen_header() {
		for (m_it = mv_cols.begin(); m_it != mv_cols.end(); m_it++) {
			cout << " " << m_it->mc_title;
		}//m_it
		cout << endl;
	}//m_write_screen_header

	void m_log(string ac_id, string ac_val) {
		m_load_field(ac_id);
		m_it->mc_val = ac_val;
	}//m_log

	void m_log(string ac_id, double ad_val) {
		float d_tmp = (float)ad_val;
		m_load_field(ac_id);
		m_it->mc_val = to_string(d_tmp);
	}//m_log

	void m_log(string ac_id, int ai_val) {
		m_load_field(ac_id);
		m_it->mc_val = to_string(ai_val);
	}//m_log

	void m_log(string ac_id, bool ab_val) {
		m_load_field(ac_id);
		m_it->mc_val = ab_val ? "true" : "false";
	}//m_log

	void m_write_record() {
		++mi_rec;
		if (!mb_to_file && !mb_to_screen)
			return;
		for (m_it = mv_cols.begin(); m_it != mv_cols.end(); m_it++) {
			if (mb_to_file) {
				if (m_it != mv_cols.begin()) mf << ",";
				mf << m_it->mc_val;
			}//if
			if (mb_to_screen && m_it->mb_to_screen && (mi_screen_gap == 1 || mi_rec % mi_screen_gap == 1))
				cout << " " << m_it->mc_val;
		}//m_it
		if (mb_to_file) mf << endl;
		if (mb_to_screen && (mi_screen_gap == 1 || mi_rec % mi_screen_gap == 1))
			cout << endl;
	}//m_write_record

private:
	void m_load_field(string ac_id) {
		for (m_it = mv_cols.begin(); m_it != mv_cols.end(); m_it++) {
			if (m_it->mc_id == ac_id) return;
		}//m_it
	}//m_load_field

};//Logger