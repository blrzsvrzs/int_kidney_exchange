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
void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& v_impu, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& dispy, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving);
void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& y, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& target, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes);
void insertion_sort(vector<unsigned short int>& w, vector<double>& y, unsigned short int& N);

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods);
void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<unsigned short int>& s);
void no_lex_min(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, vector<unsigned short int>& s_NLM, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_benefit, vector<double>& target, vector<unsigned short int>& v, vector<unsigned short int>& v_impu, double& prec, double& min_satisfaction, vector<double>& target_NLM, vector<unsigned short int>& s_noncoop);
void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, vector<unsigned short int>& s_LMWC, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, vector<unsigned short int>& v_impu, bool& target_benefit, double& min_satisfaction, vector<double>& target, vector<double>& target_LMWC, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight);
void ideal_matching(ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_benefit, bool& disp, vector<double>& target, double& prec, double& min_satisfaction, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, unsigned short int& ideal_sum, vector<unsigned short int>& s);

int main() {
	cout << "I solemnly swear that I am up to no good." << endl;
	double t0 = cpuTime();
	// input parameters and data
	unsigned short int N = 4;
	unsigned short int graph_size = 1000;
	bool target_benefit = true;
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

	vector<unsigned short int> v_impu(N, 0);
	vector<unsigned short int> v(N + 1, 0);
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

	t0 = cpuTime();
	unsigned short int ideal_sum;
	vector<unsigned short int> s_ideal(N, 0);
	vector<double> target_ideal(N, 0);
	ideal_matching(g_ideal, N, Vp, c, target_benefit, disp, target_ideal, prec, min_satisfaction, y, pos, w, p, lb, ub, opt, no_of_nodes, ideal_sum, s_ideal);
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
		coop_game(g, v, v_impu, s, c, disp, Vp, N, active_nodes, leaving);
		t1 = cpuTime();
		game_time += t1 - t0;
		t0 = t1;
		double suma = 0;
		double sumimpu = 0;
		if (target_benefit) {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i] - v_impu[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
		}
		else {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
		}
		for (unsigned short int i = 0; i < N; i++)
			target_cumm[i] += target[i];
		t1 = cpuTime();
		target_time += t1 - t0;
		t0 = t1;
		if (disp) {
			if (target_benefit) {
				cout << "Benefit: ";
			}
			else {
				cout << "Contribution: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				cout << target[i] << " ";
			}
			cout << endl;
		}
		lex_min_matching(N, v[N], target, s, no_of_active_nodes, g, edge_card_weight, c, prec, disp, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, Vp);
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
	no_lex_min(node_arrives, g, leaving, active_nodes, c, s_NLM, no_of_active_nodes, N, Vp, periods, disp, s, target_benefit, target, v, v_impu, prec, min_satisfaction, target_NLM, s_NC);
	t1 = cpuTime();
	cout << "random matching done... ";
	double no_lex_time = t1 - t0;
	t0 = t1;
	vector<unsigned short int> s_LMWC(N, 0);
	vector<double> target_LMWC(N, 0);
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, s_LMWC, no_of_active_nodes, N, Vp, periods, disp, s, Q, v, v_impu, target_benefit, min_satisfaction, target, target_LMWC, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight);
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
	res.close();
	cout << "Mischief managed! BUT CAN YOU PROVE IT CABRÓN?!" << endl;
	if (disp) {
		cout << "read time: " << read_time << endl;
		cout << "graph time: " << graph_time << endl;
		cout << "game time: " << game_time << endl;
		if (target_benefit) {
			cout << "Benefit time: ";
		}
		else {
			cout << "Contribution time: ";
		}
		cout << target_time << endl;
		cout << "matching time: " << matching_time << endl;
		cout << "no lex min matching time: " << no_lex_time << endl;
		cout << "without credits time: " << without_credits_time << endl;
		cout << "ideal time: " << ideal_time << endl;
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
	}
	return 0;
}

void ideal_matching(ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_benefit, bool& disp, vector<double>& target, double& prec, double& min_satisfaction, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, unsigned short int& ideal_sum, vector<unsigned short int>& s) {
	vector<unsigned short int> v_impu(N, 0);
	vector<unsigned short int> v(N + 1, 0);
	vector<bool> a(N, false);
	for (unsigned int i = 0; i < N; i++) {
		ListGraph::NodeMap<bool> coal1(g, false);
		a[i] = true;
		if (i > 0) {
			a[i - 1] = false;
		}
		for (unsigned short int k = i * Vp; k < (i + 1) * Vp; k++) {
			coal1[c[k]] = true;
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m1(FilterNodes<ListGraph>(g, coal1));
		coal_m1.run();
		v_impu[i] = 2 * coal_m1.matchingSize();
		ListGraph::NodeMap<bool> coal2(g, false);
		for (unsigned int j = 0; j < N; j++) {
			if (i != j) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++) {
					coal2[c[k]] = true;
				}
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m2(FilterNodes<ListGraph>(g, coal2));
		coal_m2.run();
		v[i] = 2 * coal_m2.matchingSize();
	}
	MaxMatching<ListGraph> grand_coal(g);
	grand_coal.run();
	v[N] = 2 * grand_coal.matchingSize();
	ideal_sum = v[N];
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = i * Vp; j < (i + 1) * Vp; j++) {
			if (!(grand_coal.mate(c[j]) == INVALID)) {
				s[i]++;
			}
		}
	}
	double suma = 0;
	double sumimpu = 0;
	if (target_benefit) {
		for (unsigned short int i = 0; i < N; i++) {
			suma += v[N] - v[i] - v_impu[i];
			sumimpu += v_impu[i];
		}
		for (unsigned short int i = 0; i < N; i++)
			target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
	}
	else {
		for (unsigned short int i = 0; i < N; i++) {
			suma += v[N] - v[i];
			sumimpu += v_impu[i];
		}
		for (unsigned short int i = 0; i < N; i++)
			target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
	}
	vector<unsigned short int> no_of_active_nodes(N, Vp);
	ListGraph::NodeMap<bool> active_nodes(g, true);
	ListGraph::EdgeMap<unsigned short int> edge_card_weight(g, 1);
	vector<double> credit(N, 0);
	vector<bool> leaving(no_of_nodes, false);
	lex_min_matching(N, v[N], target, s, no_of_active_nodes, g, edge_card_weight, c, prec, disp, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, Vp);
	return;
}

void no_lex_min(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, vector<unsigned short int>& s_NLM, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_benefit, vector<double>& target, vector<unsigned short int>& v, vector<unsigned short int>& v_impu, double& prec, double& min_satisfaction, vector<double>& target_NLM, vector<unsigned short int>& s_noncoop) {
	bool disp = false;
	if (dispy)
		cout << " --== Without lex min matching == -- " << endl;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		no_of_active_nodes[i] = Vp / 4;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (node_arrives[i * Vp + j] == 0) {
				active_nodes[c[i * Vp + j]] = true;
			}
			else {
				active_nodes[c[i * Vp + j]] = false;
			}
		}
	}
	for (unsigned short int Q = 0; Q < periods; Q++) {
		if (Q > 0)
			changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		coop_game(g, v, v_impu, s, c, disp, Vp, N, active_nodes, leaving);
		double suma = 0;
		double sumimpu = 0;
		if (target_benefit) {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i] - v_impu[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
		}
		else {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
		}
		for (unsigned short int i = 0; i < N; i++) {
			s_noncoop[i] += v_impu[i];
			target_NLM[i] += target[i];
		}
		if (dispy) {
			if (target_benefit) {
				cout << "Benefit: ";
			}
			else {
				cout << "Contribution: ";
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

void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& v_impu, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& dispy, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving) {
	vector<bool> a(N, false);
	for (unsigned int i = 0; i < N; i++) {
		ListGraph::NodeMap<bool> coal1(g, false);
		a[i] = true;
		if (i > 0) {
			a[i - 1] = false;
		}
		for (unsigned short int k = i * Vp; k < (i + 1) * Vp; k++) {
			if (active_nodes[c[k]]) {
				coal1[c[k]] = true;
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m1(FilterNodes<ListGraph>(g, coal1));
		coal_m1.run();
		v_impu[i] = 2 * coal_m1.matchingSize();
		ListGraph::NodeMap<bool> coal2(g, false);
		for (unsigned int j = 0; j < N; j++) {
			if (i != j) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++) {
					if (active_nodes[c[k]]) {
						coal2[c[k]] = true;
					}
				}
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m2(FilterNodes<ListGraph>(g, coal2));
		coal_m2.run();
		v[i] = 2 * coal_m2.matchingSize();
	}
	FilterNodes<ListGraph> sg(g, active_nodes);
	MaxMatching<FilterNodes<ListGraph>> grand_coal(sg);
	grand_coal.run();
	v[N] = 2 * grand_coal.matchingSize();
	if (dispy)
		cout << "grand coal: " << v[N] << endl;
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

void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, vector<unsigned short int>& s_LMWC, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, vector<unsigned short int>& v_impu, bool& target_benefit, double& min_satisfaction, vector<double>& target, vector<double>& target_LMWC, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight) {
	Q = 0;
	if (dispy)
		cout << " --== Without lex min matching == -- " << endl;
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		credit[i] = 0;
		no_of_active_nodes[i] = Vp / 4;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (node_arrives[i * Vp + j] == 0) {
				active_nodes[c[i * Vp + j]] = true;
			}
			else {
				active_nodes[c[i * Vp + j]] = false;
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
		coop_game(g, v, v_impu, s, c, dispy, Vp, N, active_nodes, leaving);
		double suma = 0;
		double sumimpu = 0;
		if (target_benefit) {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i] - v_impu[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
		}
		else {
			for (unsigned short int i = 0; i < N; i++) {
				suma += v[N] - v[i];
				sumimpu += v_impu[i];
			}
			for (unsigned short int i = 0; i < N; i++)
				target[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
		}
		for (unsigned short int i = 0; i < N; i++)
			target_LMWC[i] += target[i];
		if (dispy) {
			if (target_benefit) {
				cout << "Benefit: ";
			}
			else {
				cout << "Contribution: ";
			}
			for (unsigned short int i = 0; i < N; i++) {
				cout << target[i] << " ";
			}
			cout << endl;
		}
		lex_min_matching(N, v[N], target, s, no_of_active_nodes, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, Vp);
		for (unsigned short int i = 0; i < N; i++)
			s_LMWC[i] += s[i];
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		if (dispy)
			cin.get();
	}
	return;
}

void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp) {
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

void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& y, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& target, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp) {
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
			for (unsigned short int i = 0; i < N * Vp; i++) {
				if (active_nodes[c[i]]) {
					if (edge_card_weight[max_weight.matching(c[i])] == 1) {
						for (unsigned short int j = 0; j < N; j++) {
							if (j * Vp <= i && i < (j + 1) * Vp) {
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
			for (unsigned short int i = 0; i < N * Vp; i++) {
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

void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<unsigned short int>& s) {
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (leaving[i * Vp + j]) {
				active_nodes[c[i * Vp + j]] = false;
				no_of_active_nodes[i]--;
				leaving[i * Vp + j] = false;
			}
			else {
				if (active_nodes[c[i * Vp + j]] && node_arrives[i * Vp + j] == Q - 4) {
					active_nodes[c[i * Vp + j]] = false;
					no_of_active_nodes[i]--;
				}
			}
			if (node_arrives[i * Vp + j] == Q) {
				active_nodes[c[i * Vp + j]] = true;
				no_of_active_nodes[i]++;
			}
		}
	}
	return;
}

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c) {
	unsigned short int coal = rand() % Vp;
	unsigned short int count = 0;
	for (unsigned short int i = 0; i < N; i++) {
		while (count < Vp / 4) {
			if (active_nodes[c[i * Vp + coal]]) {
				coal = rand() % Vp;
			}
			else {
				active_nodes[c[i * Vp + coal]] = true;
				count++;
				coal = rand() % Vp;
			}
		}
		count = 0;
	}
	return;
}

void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods) {
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = 0; j < Vp; j++) {
			if (!(active_nodes[c[i * Vp + j]])) {
				node_arrives[i * Vp + j] = rand() % (periods - 1) + 1;
			}
		}
	}
	return;
}

void insertion_sort(vector<unsigned short int>& w, vector<double>& y, unsigned short int& N) {
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

void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, vector<ListGraph::Node>& c, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, ListGraph& g_ideal, vector<unsigned short int>& node_arrives, unsigned short int& no_of_nodes) {
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

bool is_next_char_digit(string& line, unsigned int l) {
	if (line[l + 1] == '0' || line[l + 1] == '1' || line[l + 1] == '2' || line[l + 1] == '3' || line[l + 1] == '4' || line[l + 1] == '5' || line[l + 1] == '6' || line[l + 1] == '7' || line[l + 1] == '8' || line[l + 1] == '9')
		return true;
	return false;
}

unsigned int char2uint(char& p) {
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

double cpuTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}
