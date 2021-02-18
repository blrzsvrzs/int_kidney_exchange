#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/adaptors.h>
#include <lemon/core.h>
#include <time.h>
#include <glpk.h>
#include <iomanip>

#include <math.h>
#include <stdio.h>
#include "windows.h"
#include "psapi.h"

using namespace lemon;
using namespace std;

double cpuTime();
bool is_next_char_digit(string& line, unsigned int l);
unsigned int char2uint(char& p);
void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, vector<ListGraph::Node>& c, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, ListGraph& g_ideal, vector<unsigned short int>& node_arrives, unsigned short int& no_of_nodes);
void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& dispy, unsigned int& S, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving);
void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& y, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& target, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes);
void shapley(vector<double>& shapl, vector<unsigned short int>& v, unsigned short int& n, unsigned int& s);
void de2bi_card(unsigned int& k, vector<bool>& a, unsigned short int& n, unsigned short int& card);
void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n);
void insertion_sort(vector<unsigned short int>& w, vector<double>& y, unsigned short int& N);

void nucl(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& x, vector<unsigned short int>& v, double& prec, double& min_satisfaction);
void excess_init(vector<double>& exc, vector<bool>& unsettled, vector<vector<bool>>& A, vector<double>& x, vector<unsigned short int>& v, unsigned int& s, unsigned short int& n);
void excess_init_mem(vector<double>& exc, vector<bool>& unsettled, vector<bool>& a, vector<double>& x, vector<unsigned short int>& v, unsigned int& s, unsigned short int& n);
void A_mx(vector<vector<bool>>& A, unsigned short int& n, unsigned int& s);
void nucl_comp(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& excess, double& prec, vector<bool>& unsettled, unsigned short int& iter, unsigned int& piv, unsigned int& sr, double& t, vector<double>& x, vector<vector<bool>>& A, double& t1, vector<double>& singleton_bounds, bool& nlsu, double& min_satisfaction);
void nucl_comp_mem(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& excess, double& prec, vector<bool>& unsettled, unsigned short int& iter, unsigned int& piv, unsigned int& sr, double& t, vector<double>& x, vector<bool>& a, double& t1, vector<double>& singleton_bounds, bool& nlsu, double& min_satisfaction);
void pivot(double& epsi, unsigned int& s, vector<double>& excess, double& prec, unsigned short int& n, vector<vector<bool>>& A, vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& unsettled, unsigned short int& rank, vector<double>& d, vector<double>& x, bool& disp, vector<vector<bool>>& Asettled, unsigned int& piv, unsigned int& sr_count, unsigned short int& iter, vector<bool>& unsettled_p, vector<double>& singleton_bounds, double& epsi_old, bool& nlsu, double& min_satisfaction);
void pivot_mem(double& epsi, unsigned int& s, vector<double>& excess, double& prec, unsigned short int& n, vector<bool>& a, vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& unsettled, unsigned short int& rank, vector<double>& d, vector<double>& x, bool& disp, vector<vector<bool>>& Asettled, unsigned int& piv, unsigned int& sr_count, unsigned short int& iter, vector<bool>& unsettled_p, vector<double>& singleton_bounds, double& epsi_old, bool& nlsu, double& min_satisfaction);
void subroutine(vector<bool>& U, vector<bool>& U2, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<vector<double>>& Arref, vector<bool>& J, double& prec, unsigned short int& n, unsigned int& tight_size, unsigned short int& tight2_size, unsigned short int& rank, bool& disp, vector<vector<bool>>& Asettled, unsigned int& sr_count, bool& u, unsigned int& s, vector<unsigned int>& T_coord, vector<unsigned int>& T2_coord, vector<bool>& unsettled, double& epsi_old, double& epsi, vector<bool>& unsettled_p, bool& settled, bool& nlsu);
void subr_upd(vector<vector<double>>& Arref, vector<bool>& J, unsigned int& i, unsigned short int& n, double& prec, vector<bool>& U, vector<bool>& U2, unsigned int& sumt, unsigned int& sumt2, vector<bool>& t, vector<bool>& t2, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, unsigned int& tight_size, unsigned short int& tight2_size, unsigned short int& rank, vector<bool>& unsettled, vector<vector<bool>>& Asettled, bool& disp, unsigned int& s, vector<unsigned int>& T_coord, vector<unsigned int>& T2_coord, double& epsi_old, double& epsi, vector<bool>& unsettled_p, bool& settled, glp_prob*& lp, vector<bool>& ar0pos);
void step(vector<bool>& T, vector<bool>& T2, vector<bool>& unsettled, vector<bool>& unsettled_p, unsigned int& s, vector<vector<bool>>& A, double& epsi, vector<double>& excess, vector<double>& d, unsigned short int& n, vector<double>& x, vector<double>& singleton_bounds, bool& disp, double& prec);
void step_mem(vector<bool>& T, vector<bool>& T2, vector<bool>& unsettled, vector<bool>& unsettled_p, unsigned int& s, vector<bool>& a, double& epsi, vector<double>& excess, vector<double>& d, unsigned short int& n, vector<double>& x, vector<double>& singleton_bounds, bool& disp, double& prec);
void imprdir(vector<double>& d, unsigned short int& n, unsigned int& t_size, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2, unsigned short int& rank, vector<vector<bool>>& Asettled, bool& disp);
void vec_min_uns(double& m, vector<double>& x, vector<bool>& unsettled, unsigned int& s);
void tight_coal(vector<bool>& T, vector<double>& excess, double& epsi, double& prec, unsigned int& s, vector<unsigned int>& T_coord, vector<bool>& unsettled, unsigned int& t_size);
void tight_coal2(vector<bool>& T2, vector<double>& x, vector<double>& singleton_bounds, double& prec, unsigned short int& n, vector<unsigned int>& T2_coord, vector<bool>& unsettled_p, unsigned short int& t2_size);
bool nonz_vec(vector<double>& x, double& prec);
void sum_vecb(unsigned int& s, vector<bool>& x);
bool binrank(vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& b, unsigned short int& n);
void sc_vec_prod(vector<double>& y, double a, vector<double>& x);
void vec_subtract(vector<double>& z, vector<double>& x, vector<double>& y);
void rowechform_piv2(vector<vector<double>>& rref, unsigned int& i, unsigned short int& n);
void rowechform_piv(vector<vector<double>>& rref, vector<int>& nonz, unsigned int& i, unsigned short int& j, unsigned int& k, unsigned short int& n);
void rowechform_loop(vector<vector<double>>& rref, vector<bool>& J, unsigned int& i, unsigned short int& j, unsigned short int& rank, double& prec, unsigned short int& n);
void rowechform(vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& B, unsigned short int& n, unsigned short int& rank);
void swap_ith_and_firstone(vector<vector<double>>& rref, vector<int>& ones, vector<int>& nonz, unsigned int& i);
void swap_ith_and_firstnz(vector<vector<double>>& rref, vector<int>& nonz, unsigned int& i);

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods);
void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<unsigned short int>& s);
void no_lex_min(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, vector<unsigned short int>& s_NLM, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_nucl, vector<double>& target, vector<unsigned short int>& v, unsigned int& S, double& prec, double& min_satisfaction, vector<double>& target_NLM, bool& empty_core, vector<unsigned short int>& s_noncoop);
void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, vector<unsigned short int>& s_LMWC, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, bool& target_nucl, double& min_satisfaction, vector<double>& target, vector<double>& target_LMWC, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, bool& empty_core);
void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_nucl, bool& disp, vector<double>& target, double& prec, double& min_satisfaction, bool& empty_core, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, unsigned short int& ideal_sum, vector<unsigned short int>& s);

int main() {
	cout << "I solemnly swear that I am up to no good." << endl;
	double t0 = cpuTime();
	// input parameters and data
	unsigned short int N = 4;
	unsigned short int graph_size = 1000;
	bool target_nucl = true;
	bool disp = false;
	unsigned short int years = 3;
	unsigned short int periods_per_year = 4;

	string line;
	ifstream inp;
	inp.open("genxml-0.xml");
	getline(inp, line);
	inp.close();

	// building the (undirected) compatibility and the (directed) 'matching' graphs
	unsigned short int Vp = 4 * (unsigned short int)((graph_size / 4) / N);
	unsigned short int no_of_nodes = N * Vp;
	vector<unsigned int> arc_out(0, 0);
	vector<unsigned int> arc_in(0, 0);
	unsigned int m = 0;
	unsigned short int k = 0;
	vector<unsigned short int> node_labels(no_of_nodes, 0);
	vector<unsigned short int> label_positions(graph_size, graph_size + 1);
	unsigned int S = pow(2, N) - 2;
	ListGraph g;
	vector<ListGraph::Node> c(no_of_nodes);
	xml_parser(line, node_labels, label_positions, c, k, g, arc_in, arc_out, m, no_of_nodes);
	double t1 = cpuTime();
	double read_time = t1 - t0;

	// determining starting pairs and arrival times of others
	unsigned short int periods = years * periods_per_year;
	vector<unsigned short int> no_of_active_nodes(N, Vp / 4);
	ListGraph::NodeMap<bool> active_nodes(g);
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = 0; j < Vp; j++) {
			active_nodes[c[i * Vp + j]] = false;
		}
	}

	unsigned int seed = GetTickCount64();
	//inp.open("seed.txt");
	//inp >> seed;
	//inp.close();
	srand(seed);
	ofstream res;
	initial_pairs(Vp, N, active_nodes, c);
	vector<unsigned short int> node_arrives(no_of_nodes, 0);
	arrival_times(node_arrives, Vp, N, active_nodes, c, periods);


	t0 = t1;
	ListGraph::EdgeMap<unsigned short int> edge_card_weight(g, 0);
	ListGraph g_ideal;
	GraphCopy<ListGraph, ListGraph> gr_cop(g, g_ideal);
	gr_cop.run();
	undi_lemon(m, arc_in, arc_out, label_positions, g, c, edge_card_weight, g_ideal, node_arrives, no_of_nodes);
	t1 = cpuTime();
	double graph_time = t1 - t0;
	t0 = t1;

	vector<unsigned short int> v(S + 1, 0);
	vector<unsigned short int> s(N, 0);
	vector<unsigned short int> s_cumm(N, 0);
	double prec = pow(10, -7);
	vector<double> target(N, 0);
	vector<double> target_cumm(N, 0);
	double min_satisfaction;
	double game_time = 0;
	double target_time = 0;
	double matching_time = 0;
	vector<double> credit(N, 0);
	vector<double> y(N, 0);
	vector<bool> pos(N, false);
	vector<unsigned short int> w(N, 0);
	unsigned short int p;
	vector<unsigned short int> lb(N, 0);
	vector<unsigned short int> ub(N, 0);
	double opt = 0;
	vector<bool> leaving(no_of_nodes, false);
	unsigned short int Q = 0;
	unsigned short int sum_s_cumm = 0;
	unsigned short int sum_s_NLM = 0;
	unsigned short int sum_s_LMWC = 0;
	unsigned short int sum_s_NC = 0;
	MaxMatching<ListGraph> ms(g);
	ms.run();
	unsigned short int max_sum = 2 * ms.matchingSize();
	bool empty_core = false;

	t0 = cpuTime();
	unsigned short int ideal_sum;
	vector<unsigned short int> s_ideal(N, 0);
	vector<double> target_ideal(N, 0);
	ideal_matching(S, g_ideal, N, Vp, c, target_nucl, disp, target_ideal, prec, min_satisfaction, empty_core, y, pos, w, p, lb, ub, opt, no_of_nodes, ideal_sum, s_ideal);
	t1 = cpuTime();
	double ideal_time = t1 - t0;
	t0 = t1;

	while (Q < periods) {
		if (disp) {
			cout << endl << "--== PERIOD " << Q + 1 << " ==--" << endl << endl;
		}
		if (disp) {
			cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << no_of_active_nodes[i] << " ";
			cout << endl;
		}
		// cooperative game and target
		coop_game(g, v, s, c, disp, S, Vp, N, active_nodes, leaving);
		t1 = cpuTime();
		game_time += t1 - t0;
		t0 = t1;
		if (target_nucl) {
			nucl(disp, N, S, target, v, prec, min_satisfaction);
			if (disp) {
				if (min_satisfaction < -prec) {
					cout << " --== The core of the matching game is EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
					empty_core = true;
				}
				else {
					cout << " --== The core of the matching game is NON-EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
				}
			}
		}
		else {
			shapley(target, v, N, S);
		}
		for (unsigned short int i = 0; i < N; i++)
			target_cumm[i] += target[i];
		t1 = cpuTime();
		target_time += t1 - t0;
		t0 = t1;
		if (disp) {
			if (target_nucl) {
				cout << "Nucleolus: ";
			}
			else {
				cout << "Shapley: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				cout << target[i] << " ";
			}
			cout << endl;
		}
		lex_min_matching(N, v[S], S, target, s, no_of_active_nodes, g, edge_card_weight, c, prec, disp, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, Vp);
		t1 = cpuTime();
		matching_time += t1 - t0;
		t0 = t1;
		for (unsigned short int i = 0; i < N; i++)
			credit[i] += target[i] - s[i];
		if (disp) {
			cout << "Credits: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << credit[i] << " ";
			cout << endl;
		}
		for (unsigned short int i = 0; i < N; i++)
			s_cumm[i] += s[i];
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		if (disp)
			cin.get();
	}
	t0 = cpuTime();
	cout << "lexmin+credits done... ";
	vector<unsigned short int> s_NLM(N, 0);
	vector<double> target_NLM(N, 0);
	vector<unsigned short int> s_NC(N, 0);
	no_lex_min(node_arrives, g, leaving, active_nodes, c, s_NLM, no_of_active_nodes, N, Vp, periods, disp, s, target_nucl, target, v, S, prec, min_satisfaction, target_NLM, empty_core, s_NC);
	t1 = cpuTime();
	cout << "random matching done... ";
	double no_lex_time = t1 - t0;
	t0 = t1;
	vector<unsigned short int> s_LMWC(N, 0);
	vector<double> target_LMWC(N, 0);
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, s_LMWC, no_of_active_nodes, N, Vp, periods, disp, s, Q, v, S, target_nucl, min_satisfaction, target, target_LMWC, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, empty_core);
	t1 = cpuTime();
	cout << "lexmin done!" << endl;
	double without_credits_time = t1 - t0;
	double sum_err = 0;
	double sum_err_LMWC = 0;
	double sum_err_NLM = 0;
	double max_err = 0;
	double max_err_LMWC = 0;
	double max_err_NLM = 0;
	for (unsigned short int i = 0; i < N; i++) {
		if (abs(s_cumm[i] - target_cumm[i]) > max_err)
			max_err = abs(s_cumm[i] - target_cumm[i]);
		if (abs(s_LMWC[i] - target_LMWC[i]) > max_err_LMWC)
			max_err_LMWC = abs(s_LMWC[i] - target_LMWC[i]);
		if (abs(s_NLM[i] - target_NLM[i]) > max_err_NLM)
			max_err_NLM = abs(s_NLM[i] - target_NLM[i]);
		sum_err += abs(s_cumm[i] - target_cumm[i]);
		sum_err_LMWC += abs(s_LMWC[i] - target_LMWC[i]);
		sum_err_NLM += abs(s_NLM[i] - target_NLM[i]);
		sum_s_cumm += s_cumm[i];
		sum_s_LMWC += s_LMWC[i];
		sum_s_NLM += s_NLM[i];
		sum_s_NC += s_NC[i];
	}
	res.open("results.txt", ofstream::out | ofstream::trunc);
	res << seed << endl << endl;
	res << sum_err << endl << sum_err_LMWC << endl << sum_err_NLM << endl << max_err << endl << max_err_LMWC << endl << max_err_NLM << endl << endl;
	res << max_sum << endl << ideal_sum << endl << sum_s_cumm << endl << sum_s_NLM << endl << sum_s_LMWC << endl << sum_s_NC << endl << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(17) << target_cumm[i] << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(0) << s_cumm[i] << endl;
	res << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(17) << target_LMWC[i] << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(0) << s_LMWC[i] << endl;
	res << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(17) << target_NLM[i] << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(0) << s_NLM[i] << endl;
	res << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << s_NC[i] << endl;
	res << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(17) << target_ideal[i] << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(0) << s_ideal[i] << endl;
	res << endl << fixed << setprecision(17) << read_time << endl << graph_time << endl << game_time << endl << target_time << endl << matching_time << endl;
	res << no_lex_time << endl << without_credits_time << endl << ideal_time;
	res << endl << endl << fixed << setprecision(0) << empty_core;
	res.close();
	if (disp) {
		cout << "read time: " << read_time << endl;
		cout << "graph time: " << graph_time << endl;
		cout << "game time: " << game_time << endl;
		if (target_nucl) {
			cout << "Nucleolus time: ";
		}
		else {
			cout << "Shapley time: ";
		}
		cout << target_time << endl;
		cout << "matching time: " << matching_time << endl;
		cout << "no lex min matching time: " << no_lex_time << endl;
		cout << "without credits time: " << without_credits_time << endl;
		cout << "ideal time: " << ideal_time << endl;
		cout << "Mischief managed! BUT CAN YOU PROVE IT CABRÓN?!" << endl;
		cout << endl << "sum of errors: " << sum_err << " (max error=" << max_err << ")" << endl;
		cout << "  LMWC errors: " << sum_err_LMWC << " (max error=" << max_err_LMWC << ")" << endl;
		cout << "   NLM errors: " << sum_err_NLM << " (max error=" << max_err_NLM << ")" << endl << endl;
		cout << " max sum: " << max_sum << endl;
		cout << "   ideal: " << ideal_sum << endl;
		cout << "     sum: " << sum_s_cumm << endl;
		cout << "sum LMWC: " << sum_s_LMWC << endl;
		cout << " sum NLM: " << sum_s_NLM << endl;
		cout << "  sum NC: " << sum_s_NC << endl << endl;
		cout << "  Target: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(4) << target_cumm[i] << " ";
		cout << endl;
		cout << "Matching: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(0) << s_cumm[i] << "      ";
		cout << endl << endl;
		cout << "  Target LMWC: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(4) << target_LMWC[i] << " ";
		cout << endl;
		cout << "Matching LMWC: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(0) << s_LMWC[i] << "      ";
		cout << endl << endl;
		cout << "   Target NLM: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(4) << target_NLM[i] << " ";
		cout << endl;
		cout << " Matching NLM: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(0) << s_NLM[i] << "      ";
		cout << endl << endl;
		cout << "  Matching NC: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << s_NC[i] << "      ";
		cout << endl << endl;
		cout << " Ideal Target: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << fixed << setprecision(4) << target_ideal[i] << " ";
		cout << endl;
		cout << "  Ideal match: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << s_ideal[i] << "      ";
		cout << endl << endl;
		cout << "TIMES" << endl;
		cout << "Init: " << fixed << setprecision(4) << read_time + graph_time << endl;
		cout << "IKEP: " << game_time + target_time + matching_time << endl;
		cout << "LMWC: " << without_credits_time << endl;
		cout << " NLM: " << no_lex_time << endl;
		cout << endl << "empty core: " << empty_core << endl;
	}
	return 0;
}

void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_nucl, bool& disp, vector<double>& target, double& prec, double& min_satisfaction, bool& empty_core, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, unsigned short int& ideal_sum, vector<unsigned short int>& s) {
	vector<unsigned short int> v(S + 1, 0);
	vector<bool> a(N, false);
	for (unsigned int i = 0; i < S; i++) {
		ListGraph::NodeMap<bool> coal(g, false);
		de2bi(i, a, N);
		for (unsigned short int j = 0; j < N; j++) {
			if (a[j]) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++) {
					coal[c[k]] = true;
				}
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m(FilterNodes<ListGraph>(g, coal));
		coal_m.run();
		v[i] = 2 * coal_m.matchingSize();
	}
	MaxMatching<ListGraph> grand_coal(g);
	grand_coal.run();
	v[S] = 2 * grand_coal.matchingSize();
	ideal_sum = v[S];
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = i * Vp; j < (i + 1) * Vp; j++) {
			if (!(grand_coal.mate(c[j]) == INVALID)) {
				s[i]++;
			}
		}
	}
	if (target_nucl) {
		nucl(disp, N, S, target, v, prec, min_satisfaction);
		if (disp) {
			if (min_satisfaction < -prec) {
				cout << " --== The core of the matching game is EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
				empty_core = true;
			}
			else {
				cout << " --== The core of the matching game is NON-EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
			}
		}
	}
	else {
		shapley(target, v, N, S);
	}
	vector<unsigned short int> no_of_active_nodes(N, Vp);
	ListGraph::NodeMap<bool> active_nodes(g, true);
	ListGraph::EdgeMap<unsigned short int> edge_card_weight(g, 1);
	vector<double> credit(N, 0);
	vector<bool> leaving(no_of_nodes, false);
	lex_min_matching(N, v[S], S, target, s, no_of_active_nodes, g, edge_card_weight, c, prec, disp, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, Vp);
	return;
}

void no_lex_min(vector<unsigned short int> &node_arrives, ListGraph &g, vector<bool> &leaving, ListGraph::NodeMap<bool> &active_nodes, vector<ListGraph::Node> &c, vector<unsigned short int> &s_NLM, vector<unsigned short int> &no_of_active_nodes, unsigned short int &N, unsigned short int &Vp, unsigned short int &periods, bool &dispy, vector<unsigned short int> &s, bool &target_nucl, vector<double> &target, vector<unsigned short int> &v, unsigned int &S, double &prec, double &min_satisfaction, vector<double> &target_NLM, bool &empty_core, vector<unsigned short int> &s_noncoop) {
	bool disp = false;
	if (dispy)
		cout << " --== Without lex min matching == -- " << endl;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		no_of_active_nodes[i] = Vp / 4;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (node_arrives[i*Vp + j] == 0) {
				active_nodes[c[i*Vp + j]] = true;
			}
			else {
				active_nodes[c[i*Vp + j]] = false;
			}
		}
	}
	for (unsigned short int Q = 0; Q < periods; Q++) {
		if (Q > 0)
			changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		coop_game(g, v, s, c, disp, S, Vp, N, active_nodes, leaving);
		if (target_nucl) {
			nucl(disp, N, S, target, v, prec, min_satisfaction);
			if (dispy) {
				if (min_satisfaction < -prec) {
					cout << " --== The core of the matching game is EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
					empty_core = true;
				}
				else {
					cout << " --== The core of the matching game is NON-EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
				}
			}
		}
		else {
			shapley(target, v, N, S);
		}
		for (unsigned short int i = 0; i < N; i++) {
			s_noncoop[i] += v[pow(2, i) - 1];
			target_NLM[i] += target[i];
		}
		if (dispy) {
			if (target_nucl) {
				cout << "Nucleolus: ";
			}
			else {
				cout << "Shapley: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				cout << target[i] << " ";
			}
			cout << endl;
		}
		if (dispy)
			cout << "s[" << Q << "]=";
		for (unsigned int i = 0; i < N; i++) {
			s_NLM[i] += s[i];
			if (dispy)
				cout << s[i] << " ";
		}
		if (dispy)
			cout << endl;
	}
	return;
}

void coop_game(ListGraph &g, vector<unsigned short int> &v, vector<unsigned short int> &s, vector<ListGraph::Node> &c, bool &dispy, unsigned int &S, unsigned short int &Vp, unsigned short int &N, ListGraph::NodeMap<bool> &active_nodes, vector<bool> &leaving) {
	vector<bool> a(N, false);
	for (unsigned int i = 0; i < S; i++) {
		ListGraph::NodeMap<bool> coal(g, false);
		de2bi(i, a, N);
		for (unsigned short int j = 0; j < N; j++) {
			if (a[j]) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++) {
					if (active_nodes[c[k]]) {
						coal[c[k]] = true;
					}
				}
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m(FilterNodes<ListGraph>(g, coal));
		coal_m.run();
		v[i] = 2 * coal_m.matchingSize();
	}
	FilterNodes<ListGraph> sg(g, active_nodes);
	MaxMatching<FilterNodes<ListGraph>> grand_coal(sg);
	grand_coal.run();
	v[S] = 2 * grand_coal.matchingSize();
	if (dispy)
		cout << "grand coal: " << v[S] << endl;
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = i * Vp; j < (i + 1) * Vp; j++) {
			if (active_nodes[c[j]]) {
				if (!(grand_coal.mate(c[j]) == INVALID)) {
					s[i]++;
					leaving[j] = true;
				}
				else {
					leaving[j] = false;
				}
			}
			else {
				leaving[j] = false;
			}
		}
	}
	if (dispy) {
		cout << "s: ";
		for (unsigned short int i = 0; i < N; i++) {
			cout << s[i] << " ";
		}
		cout << endl;
	}
	return;
}

void lex_min_without_credits(vector<unsigned short int> &node_arrives, ListGraph &g, vector<bool> &leaving, ListGraph::NodeMap<bool> &active_nodes, vector<ListGraph::Node> &c, vector<unsigned short int> &s_LMWC, vector<unsigned short int> &no_of_active_nodes, unsigned short int &N, unsigned short int &Vp, unsigned short int &periods, bool &dispy, vector<unsigned short int> &s, unsigned short int &Q, vector<unsigned short int> &v, unsigned int &S, bool &target_nucl, double &min_satisfaction, vector<double> &target, vector<double> &target_LMWC, vector<double> &credit, double &prec, vector<double> &y, vector<bool> &pos, vector<unsigned short int> &w, unsigned short int &p, double &opt, vector<unsigned short int> &lb, vector<unsigned short int> &ub, ListGraph::EdgeMap<unsigned short int> &edge_card_weight, bool &empty_core) {
	Q = 0;
	if (dispy)
		cout << " --== Without lex min matching == -- " << endl;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		credit[i] = 0;
		no_of_active_nodes[i] = Vp / 4;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (node_arrives[i*Vp + j] == 0) {
				active_nodes[c[i*Vp + j]] = true;
			}
			else {
				active_nodes[c[i*Vp + j]] = false;
			}
		}
	}
	while (Q < periods) {
		if (dispy) {
			cout << endl << "--== PERIOD " << Q + 1 << " ==--" << endl << endl;
		}
		if (dispy) {
			cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << no_of_active_nodes[i] << " ";
			cout << endl;
		}
		// cooperative game and target
		coop_game(g, v, s, c, dispy, S, Vp, N, active_nodes, leaving);
		if (target_nucl) {
			nucl(dispy, N, S, target, v, prec, min_satisfaction);
			if (dispy) {
				if (min_satisfaction < -prec) {
					cout << " --== The core of the matching game is EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
					empty_core = true;
				}
				else {
					cout << " --== The core of the matching game is NON-EMPTY! (min satisfaction = " << min_satisfaction << " )==-- " << endl;
				}
			}
		}
		else {
			shapley(target, v, N, S);
		}
		for (unsigned short int i = 0; i < N; i++)
			target_LMWC[i] += target[i];
		if (dispy) {
			if (target_nucl) {
				cout << "Nucleolus: ";
			}
			else {
				cout << "Shapley: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				cout << target[i] << " ";
			}
			cout << endl;
		}
		lex_min_matching(N, v[S], S, target, s, no_of_active_nodes, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, Vp);
		for (unsigned short int i = 0; i < N; i++)
			s_LMWC[i] += s[i];
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		if (dispy)
			cin.get();
	}
	return;
}

void lex_min_matching(unsigned short int &N, unsigned short int &grandcoal, unsigned int &S, vector<double> &target, vector<unsigned short int> &s, vector<unsigned short int> &no_of_active_nodes, ListGraph &g, ListGraph::EdgeMap<unsigned short int> &edge_card_weight, vector<ListGraph::Node> &c, double &prec, bool &dispy, vector<double> &y, vector<bool> &pos, vector<unsigned short int> &w, unsigned short int &p, vector<unsigned short int> &lb, vector<unsigned short int> &ub, double &opt, ListGraph::NodeMap<bool> &active_nodes, vector<bool> &leaving, vector<double> &credit, unsigned short int &Vp) {
	p = 0;
	for (unsigned short int i = 0; i < N; i++) {
		if (target[i] + credit[i] - s[i] >= -prec) {
			y[i] = target[i] + credit[i] - s[i];
			pos[i] = true;
		}
		else {
			y[i] = s[i] - target[i] - credit[i];
			pos[i] = false;
		}
	}
	insertion_sort(w, y, N);
	if (dispy) {
		cout << "y: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << y[i] << " ";
		cout << endl << "largest difference: " << y[w[0]] << endl;
		cout << "order: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << w[i] << " ";
		cout << endl;
	}
	if (y[w[0]] <= 0.5 + prec)
		p = N;
	while (p < N) {
		for (unsigned short int i = p; i < N; i++) {
			if (i == p || abs(y[w[i]] - y[w[p]]) < prec) {
				if (pos[w[i]]) {
					lb[w[i]] = s[w[i]];
					if (target[w[i]] + credit[w[i]] + y[w[p]] > (double)no_of_active_nodes[w[i]] - prec) {
						ub[w[i]] = no_of_active_nodes[w[i]];
					}
					else {
						ub[w[i]] = (short int)(target[w[i]] + credit[w[i]] + y[w[p]] + prec);
					}
				}
				else {
					ub[w[i]] = s[w[i]];
					if (target[w[i]] + credit[w[i]] - y[w[p]] < -prec) {
						lb[w[i]] = 0;
					}
					else {
						if (abs((short int)(target[w[i]] + credit[w[i]] - y[w[p]] + prec) - target[w[i]] - credit[w[i]] + y[w[p]]) < prec) {
							lb[w[i]] = (short int)(target[w[i]] + credit[w[i]] - y[w[p]] + prec);
						}
						else {
							lb[w[i]] = (short int)(target[w[i]] + credit[w[i]] - y[w[p]]) + 1;
						}
					}
				}
			}
			else {
				if (target[w[i]] + credit[w[i]] - y[w[p]] < -prec) {
					lb[w[i]] = 0;
				}
				else {
					lb[w[i]] = (short int)(target[w[i]] + credit[w[i]] - y[w[p]] + prec) + 1;
				}
				if (target[w[i]] + credit[w[i]] + y[w[p]] > (double)no_of_active_nodes[w[i]] - prec) {
					ub[w[i]] = no_of_active_nodes[w[i]];
				}
				else {
					if (abs((short int)(target[w[i]] + credit[w[i]] + y[w[p]] + prec) - target[w[i]] - credit[w[i]] - y[w[p]]) < prec) {
						ub[w[i]] = (short int)(target[w[i]] + credit[w[i]] + y[w[p]] + prec) - 1;
					}
					else {
						ub[w[i]] = (short int)(target[w[i]] + credit[w[i]] + y[w[p]]);
					}
				}
			}
		}
		new_matching(lb, ub, g, edge_card_weight, w, p, y, c, opt, s, grandcoal, prec, dispy, target, no_of_active_nodes, N, active_nodes, pos, leaving, credit, Vp);
	}
	return;
}

void new_matching(vector<unsigned short int> &lb, vector<unsigned short int> &ub, ListGraph &g, ListGraph::EdgeMap<unsigned short int> &edge_card_weight, vector<unsigned short int> &w, unsigned short int &p, vector<double> &y, vector<ListGraph::Node> &c, double &opt, vector<unsigned short int> &s, unsigned short int &max_match, double &prec, bool &dispy, vector<double> &target, vector<unsigned short int> &no_of_active_nodes, unsigned short int &N, ListGraph::NodeMap<bool> &active_nodes, vector<bool> &pos, vector<bool> &leaving, vector<double> &credit, unsigned short int &Vp) {
	unsigned short int lb_old = lb[w[p]];
	unsigned short int ub_old = ub[w[p]];
	if (y[w[p]] >= 1 - prec) {
		if ((lb[w[p]] >= target[w[p]] + credit[w[p]] - y[w[p]] - prec) && (lb[w[p]] <= target[w[p]] + credit[w[p]] - y[w[p]] + prec) && (lb[w[p]] < no_of_active_nodes[w[p]])) {
			lb[w[p]] += 1;
		}
		if ((ub[w[p]] >= target[w[p]] + credit[w[p]] + y[w[p]] - prec) && (ub[w[p]] <= target[w[p]] + credit[w[p]] + y[w[p]] + prec) && (ub[w[p]] > prec)) {
			ub[w[p]] -= 1;
		}
	}
	else {
		if ((target[w[p]] + credit[w[p]] - s[w[p]] >= y[w[p]] - prec) && (target[w[p]] + credit[w[p]] - s[w[p]] <= y[w[p]] + prec) && (s[w[p]] < no_of_active_nodes[w[p]])) {
			lb[w[p]] = s[w[p]] + 1;
		}
		if ((s[w[p]] - target[w[p]] - credit[w[p]] >= y[w[p]] - prec) && (s[w[p]] - target[w[p]] - credit[w[p]] <= y[w[p]] + prec) && (s[w[p]] > prec)) {
			ub[w[p]] = s[w[p]] - 1;
		}
	}
	if (dispy)
		cout << "COUNTRY " << w[p] << " lb: " << lb[w[p]] << " ; ub: " << ub[w[p]] << endl;
	if ((lb[w[p]] == ub[w[p]]) && (s[w[p]] == ub[w[p]]) && (lb[w[p]] == s[w[p]])) {
		p++;
		if (dispy) {
			cout << "p=" << p << endl;
			if (p < N)
				cout << "largest (unfinished) difference: " << y[w[p]] << endl;
		}
	}
	else {
		unsigned short int Apmax = ub[w[0]] - lb[w[0]];
		unsigned short int Bpmax = no_of_active_nodes[w[0]] - ub[w[0]];
		unsigned short int count = 0;
		for (unsigned short int i = 0; i < p + 1; i++) {
			count += no_of_active_nodes[w[i]];
			if (dispy)
				cout << "(COUNTRY " << w[i] << " lb: " << lb[w[i]] << " ; ub: " << ub[w[i]] << ")" << endl;
			if (lb[w[i]] > ub[w[i]])
				cout << endl << " *** BOUND ERROR! *** " << endl << endl;
			if (no_of_active_nodes[w[i]] - ub[w[i]] > Bpmax)
				Bpmax = no_of_active_nodes[w[i]] - ub[w[i]];
			if (ub[w[i]] - lb[w[i]] > Apmax)
				Apmax = ub[w[i]] - lb[w[i]];
		}
		for (unsigned short int i = p + 1; i < N; i++) {
			count += no_of_active_nodes[w[i]];
			if (dispy)
				cout << "COUNTRY " << w[i] << " lb: " << lb[w[i]] << " ; ub: " << ub[w[i]] << endl;
			if (lb[w[i]] > ub[w[i]])
				cout << endl << " *** BOUND ERROR! *** " << endl << endl;
			if (ub[w[i]] - lb[w[i]] > Apmax)
				Apmax = ub[w[i]] - lb[w[i]];
			if (no_of_active_nodes[w[i]] - ub[w[i]] > Bpmax)
				Bpmax = no_of_active_nodes[w[i]] - ub[w[i]];
		}
		vector<vector<ListGraph::Node>> Ap(N, vector<ListGraph::Node>(Apmax));
		vector<vector<ListGraph::Node>> Bp(N, vector<ListGraph::Node>(Bpmax));
		for (unsigned short int i = 0; i < N; i++) {
			for (unsigned short int j = 0; j < no_of_active_nodes[w[i]] - ub[w[i]]; j++) {
				Bp[w[i]][j] = g.addNode();
				active_nodes[Bp[w[i]][j]] = true;
				for (unsigned short int k = 0; k < Vp; k++) {
					if (active_nodes[c[k + w[i] * Vp]]) {
						ListGraph::Edge e = g.addEdge(Bp[w[i]][j], c[k + w[i] * Vp]);
						edge_card_weight[e] = 0;
					}
				}
			}
			count += no_of_active_nodes[w[i]] - ub[w[i]];
			for (unsigned short int j = 0; j < ub[w[i]] - lb[w[i]]; j++) {
				Ap[w[i]][j] = g.addNode();
				active_nodes[Ap[w[i]][j]] = true;
				for (unsigned short int k = 0; k < Vp; k++) {
					if (active_nodes[c[k + w[i] * Vp]]) {
						ListGraph::Edge e = g.addEdge(Ap[w[i]][j], c[k + w[i] * Vp]);
						edge_card_weight[e] = 0;
					}
				}
				for (unsigned short int k = 0; k < i; k++) {
					for (unsigned short int l = 0; l < ub[w[k]] - lb[w[k]]; l++) {
						ListGraph::Edge e = g.addEdge(Ap[w[i]][j], Ap[w[k]][l]);
						edge_card_weight[e] = 0;
					}
				}
			}
			for (unsigned short int j = 0; j < ub[w[i]] - lb[w[i]]; j++) {
				for (unsigned short int l = j + 1; l < ub[w[i]] - lb[w[i]]; l++) {
					ListGraph::Edge e = g.addEdge(Ap[w[i]][j], Ap[w[i]][l]);
					edge_card_weight[e] = 0;
				}
			}
			count += ub[w[i]] - lb[w[i]];
		}
		ListGraph::Node ap;
		if (!(count % 2 == 0)) {
			ap = g.addNode();
			active_nodes[ap] = true;
			for (unsigned short int k = 0; k < N; k++) {
				for (unsigned short int l = 0; l < ub[w[k]] - lb[w[k]]; l++) {
					ListGraph::Edge e = g.addEdge(ap, Ap[w[k]][l]);
					edge_card_weight[e] = 0;
				}
			}
		}
		MaxWeightedPerfectMatching<FilterNodes<ListGraph>, ListGraph::EdgeMap<unsigned short int>> max_weight(FilterNodes<ListGraph>(g, active_nodes), edge_card_weight);
		bool feas = max_weight.run();
		if (dispy)
			cout << "Matching feasibility: " << feas << endl;
		if (feas) {
			opt = max_weight.matchingWeight();
			opt *= 2;
		}
		else {
			opt = -DBL_MAX;
			if (dispy) {
				unsigned short int county = 0;
				for (ListGraph::NodeIt v(g); v != INVALID; ++v) {
					if (active_nodes[v])
						county++;
				}
				cout << "number of nodes: " << county << endl;
			}
		}
		if (dispy)
			cout << "opt: " << opt << "; max_match: " << max_match << endl;
		if (abs(opt - max_match) < prec) {
			vector<unsigned short int> t(N, 0);
			for (unsigned short int i = 0; i < N*Vp; i++) {
				if (active_nodes[c[i]]) {
					if (edge_card_weight[max_weight.matching(c[i])] == 1) {
						for (unsigned short int j = 0; j < N; j++) {
							if (j*Vp <= i && i < (j + 1)*Vp) {
								t[j]++;
								break;
							}
						}
					}
				}
			}
			for (unsigned short int i = 0; i < N; i++)
				s[i] = t[i];
			if (dispy) {
				cout << "s: ";
				for (unsigned short int i = 0; i < N; i++) {
					cout << s[i] << " ";
				}
				cout << endl;
			}
			for (unsigned short int i = 0; i < N*Vp; i++) {
				if (active_nodes[c[i]]) {
					if (edge_card_weight[max_weight.matching(c[i])] == 1) {
						leaving[i] = true;
					}
					else {
						leaving[i] = false;
					}
				}
				else {
					leaving[i] = false;
				}
			}
			for (unsigned short int i = 0; i < N; i++) {
				for (unsigned short int j = 0; j < no_of_active_nodes[w[i]] - ub[w[i]]; j++)
					g.erase(Bp[w[i]][j]);
				for (unsigned short int j = 0; j < ub[w[i]] - lb[w[i]]; j++)
					g.erase(Ap[w[i]][j]);
			}
			for (unsigned short int i = 0; i < N; i++) {
				if (target[i] + credit[i] - s[i] >= 0) {
					y[i] = target[i] + credit[i] - s[i];
					pos[i] = true;
				}
				else {
					y[i] = s[i] - target[i] - credit[i];
					pos[i] = false;
				}
			}
			insertion_sort(w, y, N);
			if (dispy) {
				cout << "y: ";
				for (unsigned short int i = 0; i < N; i++)
					cout << y[i] << " ";
				cout << endl << "largest (unfinished) difference: " << y[w[p]] << endl;
				cout << "order: ";
				for (unsigned short int i = 0; i < N; i++)
					cout << w[i] << " ";
				cout << endl;
			}
		}
		else {
			for (unsigned short int i = 0; i < N; i++) {
				for (unsigned short int j = 0; j < no_of_active_nodes[w[i]] - ub[w[i]]; j++)
					g.erase(Bp[w[i]][j]);
				for (unsigned short int j = 0; j < ub[w[i]] - lb[w[i]]; j++)
					g.erase(Ap[w[i]][j]);
			}
			if (s[w[p]] < lb[w[p]]) {
				if (abs(target[w[p]] + credit[w[p]] - lb_old) > y[w[p]]) {
					lb[w[p]] = s[w[p]];
				}
				else {
					lb[w[p]] = lb_old;
				}
			}
			if (s[w[p]] > ub[w[p]]) {
				if (abs(target[w[p]] + credit[w[p]] - ub_old) > y[w[p]]) {
					ub[w[p]] = s[w[p]];
				}
				else {
					ub[w[p]] = ub_old;
				}
			}
			p++;
			if (dispy) {
				cout << "p=" << p << endl;
				if (p < N)
					cout << "largest (unfinished) difference: " << y[w[p]] << endl;
			}
		}
		if (!(count % 2 == 0))
			g.erase(ap);
	}
	if (p < N) {
		if (y[w[p]] <= 0.5 + prec)
			p = N;
	}
	return;
}

void changing_nodes(ListGraph::NodeMap<bool> &active_nodes, vector<bool> &leaving, vector<unsigned short int> &no_of_active_nodes, unsigned short int &N, unsigned short int &Vp, vector<unsigned short int> &node_arrives, unsigned short int &Q, vector<ListGraph::Node> &c, vector<unsigned short int> &s) {
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (leaving[i*Vp + j]) {
				active_nodes[c[i*Vp + j]] = false;
				no_of_active_nodes[i]--;
				leaving[i*Vp + j] = false;
			}
			else {
				if (active_nodes[c[i*Vp + j]] && node_arrives[i*Vp + j] == Q - 4) {
					active_nodes[c[i*Vp + j]] = false;
					no_of_active_nodes[i]--;
				}
			}
			if (node_arrives[i*Vp + j] == Q) {
				active_nodes[c[i*Vp + j]] = true;
				no_of_active_nodes[i]++;
			}
		}
	}
	return;
}

void initial_pairs(unsigned short int &Vp, unsigned short int &N, ListGraph::NodeMap<bool> &active_nodes, vector<ListGraph::Node> &c) {
	unsigned short int coal = rand() % Vp;
	unsigned short int count = 0;
	for (unsigned short int i = 0; i < N; i++) {
		while (count < Vp / 4) {
			if (active_nodes[c[i*Vp + coal]]) {
				coal = rand() % Vp;
			}
			else {
				active_nodes[c[i*Vp + coal]] = true;
				count++;
				coal = rand() % Vp;
			}
		}
		count = 0;
	}
	return;
}

void arrival_times(vector<unsigned short int> &node_arrives, unsigned short int &Vp, unsigned short int &N, ListGraph::NodeMap<bool> &active_nodes, vector<ListGraph::Node> &c, unsigned short int &periods) {
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = 0; j < Vp; j++) {
			if (!(active_nodes[c[i*Vp + j]])) {
				node_arrives[i*Vp + j] = rand() % (periods - 1) + 1;
			}
		}
	}
	return;
}

void insertion_sort(vector<unsigned short int> &w, vector<double> &y, unsigned short int &N) {
	w[0] = 0;
	for (unsigned short int i = 1; i < N; i++) {
		if (y[i] <= y[w[i - 1]]) {
			w[i] = i;
		}
		else {
			w[i] = w[i - 1];
			if (i == 1) {
				w[0] = i;
			}
			else {
				for (unsigned short int j = i - 2; j >= 0; j--) {
					if (y[i] <= y[w[j]]) {
						w[j + 1] = i;
						break;
					}
					else {
						w[j + 1] = w[j];
						if (j == 0) {
							w[0] = i;
							break;
						}
					}
				}
			}
		}
	}
	return;
}

void undi_lemon(unsigned int &m, vector<unsigned int> &arc_in, vector<unsigned int> &arc_out, vector<unsigned short int> &label_positions, ListGraph &g, vector<ListGraph::Node> &c, ListGraph::EdgeMap<unsigned short int> &edge_card_weight, ListGraph &g_ideal, vector<unsigned short int> &node_arrives, unsigned short int &no_of_nodes) {
	bool halt = false;
	for (unsigned int i = 0; i < m - 1; i++) {
		if (label_positions[arc_in[i]] < no_of_nodes) {
			if (arc_out[i] < arc_in[i]) {
				for (unsigned int j = i + 1; arc_out[j] < arc_in[i] + 1; j++) {
					if (arc_out[j] == arc_in[i]) {
						for (unsigned int k = j; arc_out[k] == arc_out[j]; k++) {
							if (arc_out[i] == arc_in[k]) {
								ListGraph::Edge e = g.addEdge(c[label_positions[arc_out[i]]], c[label_positions[arc_in[i]]]);
								edge_card_weight[e] = 1;
								if (!(abs(node_arrives[label_positions[arc_out[i]]] - node_arrives[label_positions[arc_in[i]]]) > 3)) {
									g_ideal.addEdge(c[label_positions[arc_out[i]]], c[label_positions[arc_in[i]]]);
								}
								halt = true;
							}
							if ((halt) || (k == m - 1))
								break;
						}
					}
					if ((halt) || (j == m - 1)) {
						halt = false;
						break;
					}
				}
			}
		}
	}
	return;
}

void de2bi(unsigned int &k, vector<bool>&a, unsigned short int &n) {
	vector<bool> zero(n, false);
	a = zero;
	unsigned int i = 2;
	for (unsigned short int c = 0; c < n - 2; c++)
		i += i;
	unsigned int j = k + 1;
	unsigned short int l = n - 1;
	while (j > 0) {
		if (j >= i) {
			a[l] = true;
			j -= i;
		}
		i /= 2;
		l--;
	}
	return;
}

void shapley(vector<double> &shapl, vector<unsigned short int> &v, unsigned short int &n, unsigned int &s) {
	vector<double> w(n, 0);
	w[0] = 1 / (double)n;
	vector<double> expo(n, 1);
	for (unsigned short int j = 1; j < n; j++) {
		w[j] = w[j - 1] * j / (n - j);
		expo[j] = pow(2, j);
	}
	vector<bool> a(n, false);
	unsigned short int k = 0;
	for (unsigned short int i = 0; i < n - 1; i++)
		shapl[i] = (double)(v[expo[i] - 1]) / n;
	for (unsigned int i = 0; i < s; i++) {
		de2bi_card(i, a, n, k);
		for (unsigned short int j = 0; j < n - 1; j++) {
			if (!a[j])
				shapl[j] += w[k] * (v[i + expo[j]] - v[i]);
		}
	}
	shapl[n - 1] = v[s];
	for (unsigned short int j = 0; j < n - 1; j++)
		shapl[n - 1] -= shapl[j];
	return;
}

void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes) {
	unsigned int l = 6;
	unsigned short int n = 0;
	while (l < line.size() - 7) {
		if (line[l] == '<' && line[l + 1] == 'e') {
			l = l + 17;
			n++;
			if (!is_next_char_digit(line, l)) {
				node_labels[n - 1] = char2uint(line[l]);
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					node_labels[n - 1] = 10 * char2uint(line[l]) + char2uint(line[l + 1]);
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						node_labels[n - 1] = 100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]);
						l = l + 2;
					}
					else {
						node_labels[n - 1] = 1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]);
						l = l + 3;
					}
				}
			}
			if (n + k - 1 == node_labels[n - 1]) {
				label_positions[n + k - 1] = n - 1;
			}
			else {
				while (n + k - 1 < node_labels[n - 1]) {
					label_positions[n + k - 1] = 65535;
					label_positions.push_back(0);
					k++;
				}
				label_positions[n + k - 1] = n - 1;
			}

			c[n - 1] = g.addNode();
			l = l + 9;
			if (!is_next_char_digit(line, l)) {
				////donor_ages.push_back(char2uint(line[l]));
				//donor_ages[n - 1] = char2uint(line[l]);
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					////donor_ages.push_back(10*char2uint(line[l])+char2uint(line[l+1]));
					//donor_ages[n - 1] = 10*char2uint(line[l])+char2uint(line[l+1]);
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						////donor_ages.push_back(100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]));
						//donor_ages[n - 1] = 100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]);
						l = l + 2;
					}
					else {
						////if (!is_next_char_digit(line, l + 3)){
						////donor_ages.push_back(1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]));
						//donor_ages[n - 1] = 1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]);
						l = l + 3;
						////}
					}
				}
			}
			l = l + 25;
			if (!is_next_char_digit(line, l)) {
				if (node_labels[n - 1] != char2uint(line[l]))
					cout << "ID ERROR!" << endl;
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					if (node_labels[n - 1] != 10 * char2uint(line[l]) + char2uint(line[l + 1]))
						cout << "ID ERROR!" << endl;
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						if (node_labels[n - 1] != 100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]))
							cout << "ID ERROR!" << endl;
						l = l + 2;
					}
					else {
						//if (!is_next_char_digit(line, l + 3)){
						if (node_labels[n - 1] != 1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]))
							cout << "ID ERROR!" << endl;
						l = l + 3;
						//}
					}
				}
			}
			if (line[l + 21] == 'm')
				l = l + 29;
			else
				l = l + 28;
		}
		while (line[l] == '<' && line[l + 1] == 'm' && line[l + 6] == '>') {
			m++;
			l = l + 18;
			arc_out.push_back(node_labels[n - 1]);
			if (!is_next_char_digit(line, l)) {
				arc_in.push_back(char2uint(line[l]));
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					arc_in.push_back(10 * char2uint(line[l]) + char2uint(line[l + 1]));
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						arc_in.push_back(100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]));
						l = l + 2;
					}
					else {
						//if (!is_next_char_digit(line, l + 3)){
						arc_in.push_back(1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]));
						l = l + 3;
						//}
					}
				}
			}
			l = l + 20;
			if (!is_next_char_digit(line, l)) {
				//arc_weight.push_back(char2uint(line[l]));
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					//arc_weight.push_back(10*char2uint(line[l])+char2uint(line[l+1]));
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						//arc_weight.push_back(100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]));
						l = l + 2;
					}
					else {
						////if (!is_next_char_digit(line, l + 3)){
						//arc_weight.push_back(1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]));
						l = l + 3;
						////}
					}
				}
			}
			l = l + 17;
		}
		if (!(line[l] == '<' && line[l + 1] == 'e')) {
			l = l + 18;
		}
		if (n == no_of_nodes)
			break;
	}
	return;
}

bool is_next_char_digit(string &line, unsigned int l) {
	if (line[l + 1] == '0' || line[l + 1] == '1' || line[l + 1] == '2' || line[l + 1] == '3' || line[l + 1] == '4' || line[l + 1] == '5' || line[l + 1] == '6' || line[l + 1] == '7' || line[l + 1] == '8' || line[l + 1] == '9')
		return true;
	return false;
}

unsigned int char2uint(char &p) {
	if (p == '1')
		return 1;
	else
		if (p == '2')
			return 2;
		else
			if (p == '3')
				return 3;
			else
				if (p == '4')
					return 4;
				else
					if (p == '5')
						return 5;
					else
						if (p == '6')
							return 6;
						else
							if (p == '7')
								return 7;
							else
								if (p == '8')
									return 8;
								else
									if (p == '9')
										return 9;
									else
										return 0;
}

void de2bi_card(unsigned int &k, vector<bool>&a, unsigned short int &n, unsigned short int &card) {
	vector<bool> zero(n, false);
	card = 0;
	a = zero;
	unsigned int i = 2;
	for (unsigned short int c = 0; c < n - 2; c++)
		i += i;
	unsigned int j = k + 1;
	unsigned short int l = n - 1;
	while (j > 0) {
		if (j >= i) {
			a[l] = true;
			card++;
			j -= i;
		}
		i /= 2;
		l--;
	}
	return;
}

double cpuTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}

void nucl(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &x, vector<unsigned short int> &v, double &prec, double &min_satisfaction) {
	bool memo = false;
	bool nlsu = false;
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	vector<double> excess(s, 0);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	unsigned short int iter = 0;
	unsigned int piv = 0;
	unsigned int sr = 0;
	double t = 0;
	double t1 = cpuTime();
	for (unsigned short int i = 0; i < n; i++) {
		singleton_bounds[i] = v[pow(2, i) - 1];
		impu += singleton_bounds[i];
	}
	x = singleton_bounds;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += (v[s] - impu) / n;
	if (memo) {
		vector<bool> a(n, false);
		excess_init_mem(excess, unsettled, a, x, v, s, n);
		nucl_comp_mem(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, a, t1, singleton_bounds, nlsu, min_satisfaction);
	}
	else {
		vector<vector<bool>> A(s + 1, vector<bool>(n, false));
		A_mx(A, n, s);
		excess_init(excess, unsettled, A, x, v, s, n);
		nucl_comp(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu, min_satisfaction);
	}
	return;
}

void nucl_comp(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec, vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, vector<vector<bool>> &A, double &t1, vector<double> &singleton_bounds, bool &nlsu, double &min_satisfaction) {
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<vector<bool>> Asettled(n, vector<bool>(n, 0));
	Asettled[0] = vector<bool>(n, true);
	if (disp) {
		cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
		cout << "Starting point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	vector<double> d(n, 0);
	double epsi = 0;
	double epsi_old = -DBL_MAX;
	while (rank < n)
		pivot(epsi, s, excess, prec, n, A, Arref, J, unsettled, rank, d, x, disp, Asettled, piv, sr, iter, unsettled_p, singleton_bounds, epsi_old, nlsu, min_satisfaction);
	t = cpuTime() - t1;
	if (disp) {
		cout << "The nucleolus solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Time needed: " << t << " seconds" << endl;
		cout << "Iterations needed: " << iter << endl;
		cout << "Pivots needed: " << piv << endl;
		cout << "Subroutine solves needed: " << sr << endl;
	}
	return;
}

void pivot(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<vector<bool>>&A, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv, unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds, double &epsi_old, bool &nlsu, double &min_satisfaction) {
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned int t_size = 0;
	tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size);
	unsigned short int t2_size = 0;
	tight_coal2(T2, x, singleton_bounds, prec, n, T2_coord, unsettled_p, t2_size);
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t_size; i++)
		Atight[i] = A[T_coord[i]];
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		Atight2[i] = A[T2_coord[i]];
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	bool u = true;
	bool settled = false;
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sr_count, u, s, T_coord, T2_coord, unsettled, epsi_old, epsi, unsettled_p, settled, nlsu);
	if (disp)
		cout << endl << "  --==  subroutine finished  ==--  " << endl << endl;
	if (settled) {
		iter++;
		if (iter == 1)
			min_satisfaction = epsi;
	}
	if (disp) {
		cout << "T:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "U:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "T0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (!U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
		cout << "U0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
	}
	if (u) {
		piv++;
		if (disp)
			cout << endl << "  --==  solving improving direction LP  ==--  " << endl << endl;
		imprdir(d, n, t_size, t2_size, Atight, Atight2, U, U2, rank, Asettled, disp);
		if (disp)
			cout << endl << "  --==  improving direction obtained  ==--  " << endl << endl;
		if (disp) {
			cout << "Improving direction:" << endl;
			for (unsigned short int i = 0; i < n; i++) {
				cout << d[i] << "    ";
			}
			cout << endl;
		}
		if (disp)
			cout << endl << "  --==  computing step size  ==--  " << endl << endl;
		step(T, T2, unsettled, unsettled_p, s, A, epsi, excess, d, n, x, singleton_bounds, disp, prec);
	}
	else {
		if (disp)
			cout << endl << "  --==  minimal tight set found!  ==--  " << endl << endl;
		if (rank == n)
			return;
		if (!nlsu) {
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					if (!(binrank(Arref, J, A[i], n))) {
						unsettled[i] = false;
						unsettled[s - 1 - i] = false;
					}
				}
			}
		}
		for (unsigned short int i = 0; i < n; i++)
			if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
				unsettled_p[i] = false;
	}
	if (disp && settled)
		cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
	return;
}

void subroutine(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &tight_size, unsigned short int &tight2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sr_count, bool &u, unsigned int &s, vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, vector<bool> &unsettled, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled, bool &nlsu) {
	unsigned int sumt = 0;
	vector<bool> t(tight_size, false);
	unsigned int sumt2 = 0;
	vector<bool> t2(tight2_size, false);
	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_rows(lp, n + 1);
	glp_add_cols(lp, tight_size + tight2_size + rank);
	for (unsigned short int i = 1; i < n + 1; i++)
		glp_set_row_bnds(lp, i, GLP_FX, 0, 0);
	glp_set_row_bnds(lp, n + 1, GLP_FX, 1, 1);
	for (unsigned int i = 1; i < tight_size + tight2_size + 1; i++) {
		glp_set_col_bnds(lp, i, GLP_LO, 0, DBL_MAX);
		glp_set_obj_coef(lp, i, 1);
	}
	for (unsigned int i = tight_size + tight2_size + 1; i < tight_size + tight2_size + rank + 1; i++) {
		glp_set_col_bnds(lp, i, GLP_FR, -DBL_MAX, DBL_MAX);
		glp_set_obj_coef(lp, i, 0);
	}
	vector<int> ia((n + 1)*(tight_size + tight2_size + rank) + 1, 0);
	vector<int> ja((n + 1)*(tight_size + tight2_size + rank) + 1, 0);
	vector<double> ar((n + 1)*(tight_size + tight2_size + rank) + 1, 0);
	vector<bool> ar0pos(tight_size, false);
	unsigned int count = 0;
	for (unsigned int j = 1; j < tight_size + 1; j++) {
		for (unsigned short int i = 1; i < n + 1; i++) {
			count++;
			ia[count] = i;
			ja[count] = j;
			if (Atight[j - 1][i - 1]) {
				ar[count] = 1;
			}
			else {
				ar[count] = 0;
			}
		}
	}
	for (unsigned short int j = 1; j < tight2_size + 1; j++) {
		for (unsigned short int i = 1; i < n + 1; i++) {
			count++;
			ia[count] = i;
			ja[count] = j + tight_size;
			if (Atight2[j - 1][i - 1]) {
				ar[count] = 1;
			}
			else {
				ar[count] = 0;
			}
		}
	}
	for (unsigned short int j = 1; j < rank + 1; j++) {
		for (unsigned short int i = 1; i < n + 1; i++) {
			count++;
			ia[count] = i;
			ja[count] = j + tight_size + tight2_size;
			if (Asettled[j - 1][i - 1]) {
				ar[count] = 1;
			}
			else {
				ar[count] = 0;
			}
		}
	}
	for (unsigned int j = 1; j < tight_size + 1; j++) {
		count++;
		ia[count] = n + 1;
		ja[count] = j;
		ar[count] = 1;
	}
	for (unsigned int j = tight_size + 1; j < tight_size + tight2_size + rank + 1; j++) {
		count++;
		ia[count] = n + 1;
		ja[count] = j;
		ar[count] = 0;
	}
	int* ia_arr = ia.data();
	int* ja_arr = ja.data();
	double* ar_arr = ar.data();
	glp_load_matrix(lp, count, ia_arr, ja_arr, ar_arr);
	if (disp)
		cout << endl << "  --==  solving subroutine LP  ==--  " << endl << endl;
	glp_smcp parm;
	glp_init_smcp(&parm);
	if (!disp)
		parm.msg_lev = GLP_MSG_OFF;
	parm.presolve = GLP_ON;
	glp_simplex(lp, &parm);
	bool feas = false;
	if (glp_get_prim_stat(lp) == 2)
		feas = true;
	if (disp)
		cout << "subroutine feasibility: " << feas << endl;
	if (feas && nlsu)
		settled = true;
	sr_count++;
	unsigned int i;
	unsigned short int rank_old = rank;
	while (feas) {
		subr_upd(Arref, J, i, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, tight_size, tight2_size, rank, unsettled, Asettled, disp, s, T_coord, T2_coord, epsi_old, epsi, unsettled_p, settled, lp, ar0pos);
		if (rank == n) {
			u = false;
			glp_delete_prob(lp);
			glp_free_env();
			return;
		}
		for (unsigned int i = 0; i < tight_size; i++) {
			if (ar0pos[i]) {
				ar[count - tight_size - tight2_size - rank_old + i + 1] = 0;
				ar0pos[i] = false;
			}
		}
		if (sumt < tight_size) {
			i = 0;
			while (i < tight_size) {
				if (t[i] == false) {
					if (!(binrank(Arref, J, Atight[i], n))) {
						U[i] = false;
						t[i] = true;
						ar[count - tight_size - tight2_size - rank_old + i + 1] = 0;
						glp_set_obj_coef(lp, i + 1, 0);
						sumt++;
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
						if (sumt == tight_size && sumt2 == tight2_size) {
							u = false;
							glp_delete_prob(lp);
							glp_free_env();
							return;
						}
					}
				}
				i++;
			}
			i = 0;
			while (i < tight2_size) {
				if (t2[i] == false) {
					if (!(binrank(Arref, J, Atight2[i], n))) {
						U2[i] = false;
						t2[i] = true;
						glp_set_obj_coef(lp, tight_size + i + 1, 0);
						sumt2++;
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						if (disp)
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
						if (sumt == tight_size && sumt2 == tight2_size) {
							u = false;
							glp_delete_prob(lp);
							glp_free_env();
							return;
						}
					}
				}
				i++;
			}
			for (unsigned short int i = 0; i < n; i++)
				if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
					unsettled_p[i] = false;
			glp_load_matrix(lp, count, ia_arr, ja_arr, ar_arr);
			if (disp)
				cout << endl << "  --==  solving subroutine LP again  ==--  " << endl << endl;
			glp_simplex(lp, &parm);
			feas = false;
			if (glp_get_prim_stat(lp) == 2)
				feas = true;
			if (disp)
				cout << "subroutine feasibility: " << feas << endl;
			sr_count++;
		}
		else {
			u = false;
			glp_delete_prob(lp);
			glp_free_env();
			return;
		}
	}
	glp_delete_prob(lp);
	glp_free_env();
	return;
}

void subr_upd(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &tight_size, unsigned short int &tight2_size, unsigned short int &rank, vector<bool> &unsettled, vector<vector<bool>> &Asettled, bool &disp, unsigned int &s, vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled, glp_prob* &lp, vector<bool> &ar0pos) {
	i = 0;
	vector<double> lambdi(tight_size + tight2_size, 0);
	for (unsigned int j = 0; j < tight_size; j++) {
		if (t[j] == false)
			lambdi[j] = glp_get_col_prim(lp, j + 1);
	}
	for (unsigned short int j = 0; j < tight2_size; j++) {
		if (t2[j] == false)
			lambdi[j + tight_size] = glp_get_col_prim(lp, tight_size + j + 1);
	}
	while (i < tight_size && sumt < tight_size) {
		if (lambdi[i] > prec) {
			U[i] = false;
			t[i] = true;
			ar0pos[i] = true;
			glp_set_obj_coef(lp, i + 1, 0);
			sumt++;
			unsettled[T_coord[i]] = false;
			unsettled[s - 1 - T_coord[i]] = false;
			if (binrank(Arref, J, Atight[i], n)) {
				rank++;
				if (epsi > epsi_old) {
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T_coord[i] << " settled as well)" << endl;
				rowechform(Arref, J, Atight[i], n, rank);
				Asettled[rank - 1] = Atight[i];
				if (rank == n) {
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				glp_set_col_bnds(lp, i + 1, GLP_FR, -DBL_MAX, DBL_MAX);
			}
			else {
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, got settled (with " << s - T_coord[i] << ") without rank increase" << endl;
			}
		}
		i++;
	}
	i = 0;
	while (i < tight2_size && sumt2 < tight2_size) {
		if (lambdi[i + tight_size] > prec) {
			U2[i] = false;
			t2[i] = true;
			sumt2++;
			glp_set_obj_coef(lp, tight_size + i + 1, 0);
			unsettled[T2_coord[i]] = false;
			unsettled[s - 1 - T2_coord[i]] = false;
			if (binrank(Arref, J, Atight2[i], n)) {
				rank++;
				if (epsi > epsi_old) {
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
					cout << "lambda_" << T2_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T2_coord[i] << " settled as well)" << endl;
				rowechform(Arref, J, Atight2[i], n, rank);
				Asettled[rank - 1] = Atight2[i];
				if (rank == n) {
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				glp_set_col_bnds(lp, tight_size + i + 1, GLP_FR, -DBL_MAX, DBL_MAX);
			}
			else {
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, got settled (with " << s - T_coord[i] << ") without rank increase" << endl;
			}
		}
		i++;
	}
	for (unsigned short int i = 0; i < n; i++)
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			unsettled_p[i] = false;
	return;
}

void step(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s, vector<vector<bool>> &A, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec) {
	double alpha = DBL_MAX;
	double Ad;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha)
				alpha = (epsi - excess[i]) / (Ad - 1);
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (!T2[i] && unsettled_p[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha)
				alpha = (singleton_bounds[i] - x[i]) / d[i];
		}
	}
	if (disp)
		cout << "Step size: " << alpha << endl;
	if (disp)
		cout << endl << "  --==  step size obtained  ==--  " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New point: " << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			excess[i] += alpha*Ad;
		}
	}
	return;
}

void imprdir(vector<double>&d, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<bool> &U, vector<bool> &U2, unsigned short int &rank, vector<vector<bool>>&Asettled, bool &disp) {
	glp_prob *dir_lp;
	dir_lp = glp_create_prob();
	glp_set_obj_dir(dir_lp, GLP_MIN);
	glp_add_cols(dir_lp, n);
	glp_add_rows(dir_lp, t_size + t2_size + rank);
	for (unsigned int i = 1; i < t_size + 1; i++) {
		if (U[i - 1]) {
			glp_set_row_bnds(dir_lp, i, GLP_LO, 1, DBL_MAX);
		}
		else {
			glp_set_row_bnds(dir_lp, i, GLP_FX, 0, 0);
		}
	}
	for (unsigned short int i = 1; i < t2_size + 1; i++)
		glp_set_row_bnds(dir_lp, t_size + i, GLP_LO, 0, DBL_MAX);
	for (unsigned short int i = 1; i < rank + 1; i++)
		glp_set_row_bnds(dir_lp, t_size + t2_size + i, GLP_FX, 0, 0);
	for (unsigned short int i = 1; i < n + 1; i++)
		glp_set_col_bnds(dir_lp, i, GLP_FR, -DBL_MAX, DBL_MAX);
	vector<unsigned int> sumd(n, 0);
	for (unsigned int i = 0; i < t_size; i++) {
		if (U[i]) {
			for (unsigned short int j = 0; j < n; j++) {
				if (Atight[i][j])
					sumd[j]++;
			}
		}
	}
	for (unsigned short int i = 1; i < n + 1; i++)
		glp_set_obj_coef(dir_lp, i, sumd[i - 1]);
	vector<int> ia(n*(t_size + t2_size + rank) + 1, 0);
	vector<int> ja(n*(t_size + t2_size + rank) + 1, 0);
	vector<double> ar(n*(t_size + t2_size + rank) + 1, 0);
	unsigned int count = 0;
	for (unsigned int i = 1; i < t_size + 1; i++) {
		for (unsigned short int j = 1; j < n + 1; j++) {
			count++;
			ia[count] = i;
			ja[count] = j;
			if (Atight[i - 1][j - 1])
				ar[count] = 1;
			else
				ar[count] = 0;
		}
	}
	for (unsigned int i = 1; i < t2_size + 1; i++) {
		for (unsigned short int j = 1; j < n + 1; j++) {
			count++;
			ia[count] = i + t_size;
			ja[count] = j;
			if (Atight2[i - 1][j - 1])
				ar[count] = 1;
			else
				ar[count] = 0;
		}
	}
	for (unsigned int i = 1; i < rank + 1; i++) {
		for (unsigned short int j = 1; j < n + 1; j++) {
			count++;
			ia[count] = i + t_size + t2_size;
			ja[count] = j;
			if (Asettled[i - 1][j - 1])
				ar[count] = 1;
			else
				ar[count] = 0;
		}
	}
	int* ia_arr = ia.data();
	int* ja_arr = ja.data();
	double* ar_arr = ar.data();
	glp_load_matrix(dir_lp, count, ia_arr, ja_arr, ar_arr);
	glp_smcp parm;
	glp_init_smcp(&parm);
	if (!disp)
		parm.msg_lev = GLP_MSG_OFF;
	glp_simplex(dir_lp, &parm);
	for (unsigned short int j = 1; j < n + 1; j++)
		d[j - 1] = glp_get_col_prim(dir_lp, j);
	glp_delete_prob(dir_lp);
	glp_free_env();
	return;
}

void nucl_comp_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec, vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, vector<bool> &a, double &t1, vector<double> &singleton_bounds, bool &nlsu, double &min_satisfaction) {
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<vector<bool>> Asettled(n, vector<bool>(n, 0));
	Asettled[0] = vector<bool>(n, true);
	if (disp) {
		cout << "Starting point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	vector<double> d(n, 0);
	double epsi = 0;
	double epsi_old = -DBL_MAX;
	while (rank < n)
		pivot_mem(epsi, s, excess, prec, n, a, Arref, J, unsettled, rank, d, x, disp, Asettled, piv, sr, iter, unsettled_p, singleton_bounds, epsi_old, nlsu, min_satisfaction);
	t = cpuTime() - t1;
	if (disp) {
		cout << "The nucleolus solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Time needed: " << t << " seconds" << endl;
		cout << "Iterations needed: " << iter << endl;
		cout << "Pivots needed: " << piv << endl;
		cout << "Subroutine solves needed: " << sr << endl;
	}
	return;
}

void pivot_mem(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<bool>&a, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv, unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds, double &epsi_old, bool &nlsu, double &min_satisfaction) {
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned int t_size = 0;
	tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size);
	unsigned short int t2_size = 0;
	tight_coal2(T2, x, singleton_bounds, prec, n, T2_coord, unsettled_p, t2_size);
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t_size; i++)
		de2bi(T_coord[i], Atight[i], n);
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		de2bi(T2_coord[i], Atight2[i], n);
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	bool u = true;
	bool settled = false;
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sr_count, u, s, T_coord, T2_coord, unsettled, epsi_old, epsi, unsettled_p, settled, nlsu);
	if (disp)
		cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
	if (settled) {
		iter++;
		if (iter == 1)
			min_satisfaction = epsi;
	}
	if (disp) {
		cout << "T:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "U:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "T0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (!U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
		cout << "U0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
	}
	if (u) {
		piv++;
		if (disp)
			cout << endl << "   ---===   SOLVING IMPROVING DIRECTION LP   ===---   " << endl << endl;
		imprdir(d, n, t_size, t2_size, Atight, Atight2, U, U2, rank, Asettled, disp);
		if (disp)
			cout << endl << "   ---===   IMPROVING DIRECTION OBTAINED   ===---   " << endl << endl;
		if (disp) {
			cout << "Improving direction:" << endl;
			for (unsigned short int i = 0; i < n; i++) {
				cout << d[i] << "    ";
			}
			cout << endl;
		}
		if (disp)
			cout << endl << "   ---===   COMPUTING STEP SIZE   ===---   " << endl << endl;
		step_mem(T, T2, unsettled, unsettled_p, s, a, epsi, excess, d, n, x, singleton_bounds, disp, prec);
	}
	else {
		if (disp)
			cout << "Min tight set found! Rank increased to: " << rank << endl;
		if (rank == n)
			return;
		if (!nlsu) {
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					de2bi(i, a, n);
					if (!(binrank(Arref, J, a, n))) {
						unsettled[i] = false;
						unsettled[s - 1 - i] = false;
					}
				}
			}
		}
		for (unsigned short int i = 0; i < n; i++)
			if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
				unsettled_p[i] = false;
	}
	return;
}

void step_mem(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s, vector<bool> &a, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec) {
	double alpha = DBL_MAX;
	double Ad;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha)
				alpha = (epsi - excess[i]) / (Ad - 1);
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (!T2[i] && unsettled_p[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha)
				alpha = (singleton_bounds[i] - x[i]) / d[i];
		}
	}
	if (disp)
		cout << "Step size: " << alpha << endl;
	if (disp)
		cout << endl << "   ---===   STEP SIZE OBTAINED   ===---   " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New x point: " << endl;
		for (unsigned short int i = 0; i < n; i++) {
			cout << x[i] << endl;
		}
	}
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			Ad = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					Ad += d[j];
			}
			excess[i] += alpha*Ad;
		}
	}
	return;
}

void excess_init(vector<double> &exc, vector<bool> &unsettled, vector<vector<bool>> &A, vector<double> &x, vector<unsigned short int> &v, unsigned int &s, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			exc[i] = -(short int)v[i];
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					exc[i] += x[j];
			}
		}
		else
			exc[i] = DBL_MAX;
	}
	return;
}

void excess_init_mem(vector<double> &exc, vector<bool> &unsettled, vector<bool> &a, vector<double> &x, vector<unsigned short int> &v, unsigned int &s, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			exc[i] = -(short int)v[i];
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					exc[i] += x[j];
			}
		}
		else
			exc[i] = DBL_MAX;
	}
	return;
}

void tight_coal(vector<bool>&T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (abs(excess[i] - epsi) < prec) {
				t_size++;
				T[i] = true;
				T_coord.push_back(i);
			}
		}
	}
	return;
}

void tight_coal2(vector<bool>&T2, vector<double> &x, vector<double> &singleton_bounds, double &prec, unsigned short int &n, vector<unsigned int> &T2_coord, vector<bool> &unsettled_p, unsigned short int &t2_size) {
	for (unsigned int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (abs(x[i] - singleton_bounds[i]) < prec) {
				t2_size++;
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
			}
		}
	}
	return;
}

void rowechform(vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &B, unsigned short int &n, unsigned short int &rank) {
	double prec = pow(10, -10);
	vector<vector<double>> rref(rank + 1, vector<double>(n, 0));
	rref[0] = Arref[0];// first row done
	for (unsigned int i = 1; i < rank + 1; i++) {
		for (unsigned short int j = 0; j < n; j++) {
			if (i < rank) {
				if (Arref[i][j] > prec || Arref[i][j] < -prec)
					rref[i][j] = Arref[i][j];
			}
			else {
				if (B[j])
					rref[i][j] = 1;
			}
		}
	}
	for (unsigned int i = 1; i < rank + 1; i++) {
		if (rref[i][0] > prec || rref[i][0] < -prec) {
			if (rref[i][0] < 1 + prec && rref[i][0] > 1 - prec)
				vec_subtract(rref[i], rref[0], rref[i]);
			else {
				rowechform_piv2(rref, i, n);
			}
		}
	}//first column done
	unsigned int l = 1;
	unsigned short int j = 1;
	while (l < rank + 1 && j < n) {
		rowechform_loop(rref, J, l, j, rank, prec, n);
	}
	if (rank + 1 < n) {
		for (unsigned short int l = 1; l < rank + 1; l++) {
			Arref[l] = rref[l];
		}
	}
	else {
		for (unsigned short int l = 1; l < n; l++) {
			Arref[l] = rref[l];
		}
	}
}

void rowechform_loop(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned short int &rank, double &prec, unsigned short int &n) {
	vector<int> nonz(0, 0);
	vector<int> ones(0, 0);
	for (unsigned int k = i; k < rank + 1; k++) {
		if (rref[k][j] > prec || rref[k][j] < -prec) {
			if (rref[k][j] < 1 + prec && rref[k][j] > 1 - prec) {
				ones.push_back(k);
			}
			else {
				nonz.push_back(k);
			}
		}
	}
	if (ones.size() == 0 && nonz.size() == 0) {
		j++; // all zero column, we can move on
	}
	else {
		if (ones.size() == 0) { // there are no 1s, but have non-zeros
			if (nonz[0] != i) { // if the first non-zero is not in the i-th row => Arref[i][j]=0
				swap_ith_and_firstnz(rref, nonz, i);// swap i-th and first non-zero row
			}
			sc_vec_prod(rref[i], 1 / rref[i][j], rref[i]);
		}
		else { // there are 1s => if first 1 is in the 1-th row, there's nothing to do
			if (ones[0] != i) { // if it's not in the i-th row, we swap the i-th and the first 1
				swap_ith_and_firstone(rref, ones, nonz, i);
			}
		}
		// eliminate all the pos with i-th row, then i++, j++ and move on
		if (ones.size() > 0) {
			if (ones[0] == i) {
				for (unsigned int k = 1; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
			else {
				for (unsigned int k = 0; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
		}
		if (nonz.size() > 0) {
			if (nonz[0] == i) {
				for (unsigned int k = 1; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
			else {
				for (unsigned int k = 0; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
		}
		i++;
		J[j] = false;
		j++;
	}
}

void rowechform_piv(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i, unsigned short int &j, unsigned int &k, unsigned short int &n) {
	vector<double> aux(n, 0);
	sc_vec_prod(aux, rref[nonz[k]][j], rref[i]);
	vec_subtract(rref[nonz[k]], rref[nonz[k]], aux);
}

void rowechform_piv2(vector<vector<double>> &rref, unsigned int &i, unsigned short int &n) {
	vector<double> aux(n, 0);
	sc_vec_prod(aux, rref[i][0], rref[0]);
	vec_subtract(rref[i], aux, rref[i]);
}

void swap_ith_and_firstnz(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i) {
	vector<double> aux = rref[nonz[0]]; // swap i-th and first non-zero row
	rref[nonz[0]] = rref[i];
	rref[i] = aux;
	nonz[0] = i;
}

void swap_ith_and_firstone(vector<vector<double>> &rref, vector<int> &ones, vector<int> &nonz, unsigned int &i) {
	vector<double> aux = rref[ones[0]];
	rref[ones[0]] = rref[i];
	rref[i] = aux;
	if (nonz.size() > 0) {
		if (nonz[0] == i) {
			nonz[0] = ones[0];
		}
	}
	ones[0] = i;
}

bool binrank(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n) {
	double prec = pow(10, -10);
	vector<double> B(n, 0);
	for (unsigned short int i = 0; i < n; i++) {
		if (b[i] == true)
			B[i] = 1;
	}
	unsigned int m = 0;
	bool size = true;
	while (size) {
		if (nonz_vec(Arref[m], prec))
			m++;
		else
			size = false;
	}
	if (m >= n)
		return false;
	else {
		vector<bool> pivot_col(n, false);
		for (unsigned short int i = 0; i < n; i++) {
			if (J[i] == false)
				pivot_col[i] = true;
		}
		unsigned short int j = 0;
		vector<bool> piv(n, false);
		vector<double> aux(n, 0);
		unsigned short int k = 0;
		unsigned short int I = 0;
		unsigned int s = 0;
		unsigned short int ind = 0;
		unsigned short int count = 0;
		while (j < n) {
			for (unsigned short i = 0; i < n; i++) {
				if (B[i] > prec || B[i] < -prec)
					piv[i] = true;
			}
			sum_vecb(s, piv);
			if (s == 0)
				return false;
			else {
				while (k == 0) {
					if (piv[I] == true)
						k = I + 1;
					I++;
				}
				k--;
				I = 0;
				if (J[k] == true)
					return true;
				else {
					while (count < k + 1) {
						if (pivot_col[count])
							ind++;
						count++;
					}
					ind--;
					count = 0;
					sc_vec_prod(aux, B[k] / Arref[ind][k], Arref[ind]);
					vec_subtract(B, B, aux);
					j++;
				}
			}
			for (unsigned short int l = 0; l < n; l++)
				piv[l] = false;
			k = 0;
			ind = 0;
		}
		return false;
	}
}

void A_mx(vector<vector<bool>>&A, unsigned short int &n, unsigned int &s) {
	// creates boolean matrix A containing all the possible n-length boolean vectors (except for full zeros)
	for (unsigned int k = 0; k != s + 1; k++) {
		unsigned int i = 2;
		for (unsigned short int c = 0; c < n - 2; c++)
			i += i;
		unsigned int j = k + 1;
		unsigned short int l = n - 1;
		while (j > 0) {
			if (j >= i) {
				A[k][l] = true;
				j -= i;
			}
			i /= 2;
			l--;
		}
	}
	return;
}

void sum_vecb(unsigned int&s, vector<bool>&x) {
	// sums up the values of boolean x
	s = 0;
	for (unsigned int i = 0; i < x.size(); i++)
		s += x[i];
	return;
}

bool nonz_vec(vector<double>&x, double &prec) {
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] > prec || x[i] < -prec)
			return true;
	}
	return false;
}

void vec_min_uns(double&m, vector<double>&x, vector<bool> &unsettled, unsigned int &s) {
	m = DBL_MAX;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i] && x[i] < m)
			m = x[i];
	}
	return;
}

void vec_subtract(vector<double>&z, vector<double>&x, vector<double>&y) {
	// subtracts vector (double) y from vector (double) x
	for (unsigned int i = 0; i != x.size(); i++)
		z[i] = x[i] - y[i];
	return;
}

void sc_vec_prod(vector<double>& y, double a, vector<double> &x) {
	for (unsigned int i = 0; i < x.size(); i++)
		y[i] = a*x[i];
	return;
}
