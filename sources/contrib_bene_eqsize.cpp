/*
*    int_kidney_exchange
*    contrib_bene_eqsize.cpp
*    Purpose: computational study for
*             Benedek et al. (2021) - Computing Balanced Solutions for
*             Large Kidney Exchange Schemes
*             https://arxiv.org/abs/2109.06788
*             benefit and contribution value allocations with equal country sizes
*
*    @author Marton Benedek
*    @version 1.0 19/10/2021
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program. If not, see:
*    <https://github.com/blrzsvrzs/int_kidney_exchange>.
*/

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
void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& y, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& target, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp, bool& disp, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, unsigned int& m, unsigned short int& no_of_nodes);
void insertion_sort(vector<unsigned short int>& w, vector<double>& y, unsigned short int& N);

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods);
void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, bool& disp, unsigned short int& Vp, bool& d1);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, vector<ListGraph::Node>& c, vector<unsigned short int>& s);
void arbitrary_matching(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_omega, vector<double>& target, vector<unsigned short int>& v, unsigned int& S, double& prec, double& t0, vector<vector<unsigned short int>>& actual_alloc_NC, vector<vector<unsigned short int>>& actual_alloc, vector<vector<double>>& initial_alloc, vector<unsigned short int>& v_impu);
void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, bool& target_omega, vector<double>& target, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, double& t0, bool& d1, unsigned short int& I, vector<vector<double>>& init_alloc, vector<vector<unsigned short int>>& actual_alloc, vector<unsigned short int>& v_impu);
void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_omega, bool& disp, bool& dispy, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, vector<unsigned short int>& s_ideal, double& t0, bool& d1, vector<double>& target_ideal, vector<unsigned short int>& s_ideal_d1, double& ideal_time, double& ideal_d1_time, vector<unsigned short int>& v_impu);
void LMC(unsigned short int& Q, unsigned short int& periods, bool& dispy, unsigned short int& N, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, unsigned int& S, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, unsigned short int& Vp, double& t1, double& game_time, bool& target_omega, double& t0, vector<double>& target, double& prec, double& target_time, vector<vector<double>>& init_alloc, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, bool& d1, vector<double>& credit, double& matching_time, vector<vector<unsigned short int>>& actual_alloc, vector<unsigned short int>& node_arrives, vector<unsigned short int>& v_impu);
void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n);

int main() {
	cout << "I solemnly swear that I am up to no good." << endl;
	double t0 = cpuTime();
	// input parameters and data
	unsigned short int N = 4;
	unsigned short int inst = 0; // instance number, integer between 0 and 99
	bool target_omega = false; // true: first to evaluate is the benefit value
	bool second = true; // false: only the first of contribution/benefit value will be evaluated
	bool dispy = false; // true: information in terminal while running
	bool disp = false; // true: extremely detailed information while running, avoid with large graphs

	bool d1 = false;
	unsigned short int years = 6;
	unsigned short int periods_per_year = 4;
	string line;
	ifstream inp;
	unsigned short int graph_size = 2000;
	inp.open("genxml-" + to_string(inst) + ".xml"); // 1 out of the 100 instances generated by William Pettersson's web tool: https://wpettersson.github.io/kidney-webapp/#/
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
	t0 = t1;

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
	seed = 18397106; // for instance genxml-0.xml it is 18397106 for N=4, 20469263 for N=5, 22805501 for N=6, 25083567 for N=7, 27432197 for N=8, 30095162 for N=9, 33411331 for N=10, 6368187 for N=11, 13109406 for N=12, 23969593 for N=13, 43358281 for N=14, 79289906 for N=15
	//inp.open("n" + to_string(N) + "inst" + to_string(inst) + ".txt");
	//inp >> seed;
	//inp.close();
	srand(seed);
	ofstream res;
	initial_pairs(Vp, N, active_nodes, c);
	vector<unsigned short int> node_arrives(no_of_nodes, 0);
	arrival_times(node_arrives, Vp, N, active_nodes, c, periods);
	t1 = cpuTime();
	double rand_time = t1 - t0;


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
	vector<vector<unsigned short int>> actual_alloc_LMC(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> actual_alloc_LM(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> actual_alloc_d1C(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> actual_alloc_d1(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> actual_alloc_rand(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> actual_alloc_NC(periods, vector<unsigned short int>(N, 0));
	double prec = pow(10, -7);
	vector<double> target(N, 0);
	vector<vector<double>> init_alloc_LMC(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_LM(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_d1C(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_d1(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_rand(periods, vector<double>(N, 0));
	double game_time = 0;
	double target_time = 0;
	double matching_time = 0;
	unsigned short int I1 = 0;
	unsigned short int I11 = 0;
	unsigned short int I2 = 0;
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
	t1 = cpuTime();
	double init_time = t1 - t0;
	t0 = t1;

	vector<unsigned short int> s_ideal(N, 0);
	vector<double> init_alloc_ideal(N, 0);
	vector<unsigned short int> s_ideal_d1(N, 0);
	double ideal_time = 0;
	double ideal_d1_time = 0;
	t0 = cpuTime();
	ideal_matching(S, g_ideal, N, Vp, c, target_omega, disp, dispy, prec, y, pos, w, p, lb, ub, opt, no_of_nodes, s_ideal, t0, d1, init_alloc_ideal, s_ideal_d1, ideal_time, ideal_d1_time, v_impu);
	ideal_time += init_time + read_time + rand_time + graph_time;
	ideal_d1_time += init_time + read_time + rand_time + graph_time;
	cout << "ideal done... ";
	t0 = cpuTime();

	LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, game_time, target_omega, t0, target, prec, target_time, init_alloc_LMC, I1, I11, I2, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, matching_time, actual_alloc_LMC, node_arrives, v_impu);
	double LMC_time = init_time + game_time + target_time + matching_time + read_time + rand_time + graph_time;
	cout << "lexmin+c done... ";

	t0 = cpuTime();
	Q = 0;
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
	double tmp = 0;
	unsigned short int I = 0;
	d1 = true;
	LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, tmp, target_omega, tmp, target, prec, tmp, init_alloc_d1C, I, I, I, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, tmp, actual_alloc_d1C, node_arrives, v_impu);
	t1 = cpuTime();
	double d1C_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	d1 = false;
	cout << "d1+c done... ";

	t0 = cpuTime();
	arbitrary_matching(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, target_omega, target, v, S, prec, t0, actual_alloc_NC, actual_alloc_rand, init_alloc_rand, v_impu);
	t1 = cpuTime();
	double arbitrary_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	cout << "arbitrary matching done... ";

	t0 = cpuTime();
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_omega, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_LM, actual_alloc_LM, v_impu);
	t1 = cpuTime();
	double LM_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	cout << "lexmin done...";

	t0 = cpuTime();
	Q = 0;
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
	d1 = true;
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_omega, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_d1, actual_alloc_d1, v_impu);
	t1 = cpuTime();
	double d1_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	d1 = false;
	cout << "d1 done!" << endl;

	res.open("results.txt", ofstream::out | ofstream::trunc);
	res << seed << endl << endl;
	res << read_time << endl << graph_time << endl << rand_time << endl << init_time << endl << game_time << endl << target_time << endl << matching_time << endl << endl;
	res << LMC_time << endl << d1C_time << endl << LM_time << endl << d1_time << endl << arbitrary_time << endl << ideal_time << endl << ideal_d1_time << endl << endl;
	res << I1 << endl << I11 << endl << I2 << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_LMC[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << actual_alloc_LMC[Q][i] << endl;
		res << endl;
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_d1C[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << actual_alloc_d1C[Q][i] << endl;
		res << endl;
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_LM[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << actual_alloc_LM[Q][i] << endl;
		res << endl;
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_d1[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << actual_alloc_d1[Q][i] << endl;
		res << endl;
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_rand[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << actual_alloc_rand[Q][i] << endl;
		res << endl;
	}
	res << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(17) << init_alloc_ideal[i] << endl;
	res << endl;
	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(0) << s_ideal[i] << endl;
	res << endl;

	for (unsigned short int i = 0; i < N; i++)
		res << fixed << setprecision(0) << s_ideal_d1[i] << endl;
	res << endl;

	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << actual_alloc_NC[Q][i] << endl;
		res << endl;
	}
	res.close();

	if (second) {
		target_omega = true;
		I1 = 0;
		I11 = 0;
		I2 = 0;
		Q = 0;

		t0 = cpuTime();
		ideal_matching(S, g_ideal, N, Vp, c, target_omega, disp, dispy, prec, y, pos, w, p, lb, ub, opt, no_of_nodes, s_ideal, t0, d1, init_alloc_ideal, s_ideal_d1, ideal_time, ideal_d1_time, v_impu);
		ideal_time += init_time + read_time + rand_time + graph_time;
		ideal_d1_time += init_time + read_time + rand_time + graph_time;
		cout << "ideal done... ";
		t0 = cpuTime();

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

		LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, game_time, target_omega, t0, target, prec, target_time, init_alloc_LMC, I1, I11, I2, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, matching_time, actual_alloc_LMC, node_arrives, v_impu);
		LMC_time = init_time + game_time + target_time + matching_time + read_time + rand_time + graph_time;
		cout << "lexmin+c done... ";

		t0 = cpuTime();
		Q = 0;
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
		tmp = 0;
		I = 0;
		d1 = true;
		LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, tmp, target_omega, tmp, target, prec, tmp, init_alloc_d1C, I, I, I, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, tmp, actual_alloc_d1C, node_arrives, v_impu);
		t1 = cpuTime();
		d1C_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		d1 = false;
		cout << "d1+c done... ";

		t0 = cpuTime();
		arbitrary_matching(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, target_omega, target, v, S, prec, t0, actual_alloc_NC, actual_alloc_rand, init_alloc_rand, v_impu);
		t1 = cpuTime();
		arbitrary_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		cout << "arbitrary matching done... ";

		t0 = cpuTime();
		lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_omega, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_LM, actual_alloc_LM, v_impu);
		t1 = cpuTime();
		LM_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		cout << "lexmin done...";

		t0 = cpuTime();
		Q = 0;
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
		d1 = true;
		lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_omega, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_d1, actual_alloc_d1, v_impu);
		t1 = cpuTime();
		d1_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		d1 = false;
		cout << "d1 done!" << endl;

		res.open("results.txt", ofstream::out | ofstream::trunc);
		res << seed << endl << endl;
		res << read_time << endl << graph_time << endl << rand_time << endl << init_time << endl << game_time << endl << target_time << endl << matching_time << endl << endl;
		res << LMC_time << endl << d1C_time << endl << LM_time << endl << d1_time << endl << arbitrary_time << endl << ideal_time << endl << ideal_d1_time << endl << endl;
		res << I1 << endl << I11 << endl << I2 << endl;
		for (unsigned short int Q = 0; Q < periods; Q++) {
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(17) << init_alloc_LMC[Q][i] << endl;
			res << endl;
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(0) << actual_alloc_LMC[Q][i] << endl;
			res << endl;
		}
		res << endl;
		for (unsigned short int Q = 0; Q < periods; Q++) {
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(17) << init_alloc_d1C[Q][i] << endl;
			res << endl;
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(0) << actual_alloc_d1C[Q][i] << endl;
			res << endl;
		}
		res << endl;
		for (unsigned short int Q = 0; Q < periods; Q++) {
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(17) << init_alloc_LM[Q][i] << endl;
			res << endl;
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(0) << actual_alloc_LM[Q][i] << endl;
			res << endl;
		}
		res << endl;
		for (unsigned short int Q = 0; Q < periods; Q++) {
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(17) << init_alloc_d1[Q][i] << endl;
			res << endl;
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(0) << actual_alloc_d1[Q][i] << endl;
			res << endl;
		}
		res << endl;
		for (unsigned short int Q = 0; Q < periods; Q++) {
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(17) << init_alloc_rand[Q][i] << endl;
			res << endl;
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(0) << actual_alloc_rand[Q][i] << endl;
			res << endl;
		}
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_ideal[i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << s_ideal[i] << endl;
		res << endl;

		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << s_ideal_d1[i] << endl;
		res << endl;

		for (unsigned short int Q = 0; Q < periods; Q++) {
			for (unsigned short int i = 0; i < N; i++)
				res << fixed << setprecision(0) << actual_alloc_NC[Q][i] << endl;
			res << endl;
		}
		res.close();
	}
	cout << "Mischief managed!" << endl;
	return 0;
}

void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n) {
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

void LMC(unsigned short int& Q, unsigned short int& periods, bool& dispy, unsigned short int& N, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, unsigned int& S, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, unsigned short int& Vp, double& t1, double& game_time, bool& target_omega, double& t0, vector<double>& target, double& prec, double& target_time, vector<vector<double>>& init_alloc, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, bool& d1, vector<double>& credit, double& matching_time, vector<vector<unsigned short int>>& actual_alloc, vector<unsigned short int>& node_arrives, vector<unsigned short int>&v_impu) {
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
		t0 = cpuTime();
		double suma = 0;
		double sumimpu = 0;
		if (target_omega) {
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
		if (disp) {
			if (target_omega) {
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
		t1 = cpuTime();
		target_time += t1 - t0;
		init_alloc[Q] = target;
		t0 = cpuTime();
		lex_min_matching(N, v[N], S, target, s, no_of_active_nodes, I1, I11, I2, g, edge_card_weight, c, prec, disp, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
		t1 = cpuTime();
		matching_time += t1 - t0;
		actual_alloc[Q] = s;
		t0 = cpuTime();
		for (unsigned short int i = 0; i < N; i++)
			credit[i] += target[i] - s[i];
		if (disp) {
			cout << "Credits: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << credit[i] << " ";
			cout << endl;
		}
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		if (disp)
			cin.get();
	}
	return;
}

void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_omega, bool& disp, bool& dispy, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, vector<unsigned short int>& s_ideal, double& t0, bool& d1,  vector<double>& target_ideal, vector<unsigned short int>& s_ideal_d1, double& ideal_time, double& ideal_d1_time, vector<unsigned short int> &v_impu) {
	vector<bool> a(N, false);
	vector<unsigned short int> v(N + 1, 0);
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
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = i * Vp; j < (i + 1) * Vp; j++) {
			if (!(grand_coal.mate(c[j]) == INVALID)) {
				s_ideal[i]++;
			}
		}
	}
	s_ideal_d1 = s_ideal;
	vector<unsigned short int> no_of_active_nodes(N, Vp);
	ListGraph::NodeMap<bool> active_nodes(g, true);
	unsigned short int I = 0;
	ListGraph::EdgeMap<unsigned short int> edge_card_weight(g, 1);
	vector<double> credit(N, 0);
	vector<bool> leaving(no_of_nodes, false);
	double suma = 0;
	double sumimpu = 0;
	if (target_omega) {
		for (unsigned short int i = 0; i < N; i++) {
			suma += v[N] - v[i] - v_impu[i];
			sumimpu += v_impu[i];
		}
		for (unsigned short int i = 0; i < N; i++)
			target_ideal[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i] - v_impu[i]) / suma);
	}
	else {
		for (unsigned short int i = 0; i < N; i++) {
			suma += v[N] - v[i];
			sumimpu += v_impu[i];
		}
		for (unsigned short int i = 0; i < N; i++)
			target_ideal[i] = v_impu[i] + (v[N] - sumimpu) * ((v[N] - v[i]) / suma);
	}
	ideal_time = cpuTime() - t0;
	ideal_d1_time = ideal_time;
	t0 = cpuTime();
	lex_min_matching(N, v[N], S, target_ideal, s_ideal, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
	ideal_time += cpuTime() - t0;
	t0 = cpuTime();
	d1 = true;
	lex_min_matching(N, v[N], S, target_ideal, s_ideal_d1, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
	ideal_d1_time += cpuTime() - t0;
	d1 = false;
	return;
}

void arbitrary_matching(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_omega, vector<double>& target, vector<unsigned short int>& v, unsigned int& S, double& prec, double& t0, vector<vector<unsigned short int>>& actual_alloc_NC, vector<vector<unsigned short int>>& actual_alloc, vector<vector<double>>& initial_alloc, vector<unsigned short int>& v_impu) {
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
		if (target_omega) {
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
		initial_alloc[Q] = target;
		if (dispy) {
			if (target_omega) {
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
		actual_alloc[Q] = s;
		if (dispy)
			cout << "s[" << Q << "]=";
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

void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, bool& target_omega, vector<double>& target, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, double& t0, bool& d1, unsigned short int& I, vector<vector<double>>& init_alloc, vector<vector<unsigned short int>>& actual_alloc, vector<unsigned short int>& v_impu) {
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
		if (target_omega) {
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
		init_alloc[Q] = target;
		if (dispy) {
			if (target_omega) {
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
		lex_min_matching(N, v[N], S, target, s, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
		actual_alloc[Q] = s;
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, c, s);
		if (dispy)
			cin.get();
	}
	return;
}

void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, bool& disp, unsigned short int& Vp, bool& d1) {
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
		new_matching(lb, ub, g, edge_card_weight, w, p, y, c, opt, s, grandcoal, prec, dispy, target, no_of_active_nodes, N, active_nodes, pos, leaving, credit, Vp, disp, I1, I11, I2);
		if (p > 0 && d1)
			break;
	}
	return;
}

void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& y, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& target, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp, bool& disp, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2) {
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
		if (feas)
			opt = 2*max_weight.matchingWeight();
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
			I1++;
			if (p == 0)
				I11++;
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
			I2++;
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
