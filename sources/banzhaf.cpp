/*
*    int_kidney_exchange
*	 banzhaf.cpp
*    Purpose: computational study for
*             M. Benedek, P. Biró, W. Kern and D. Paulusma (2022)
*			  Computing Balanced Solutions for Large Kidney Exchange Schemes
*             https://dl.acm.org/doi/10.5555/3535850.3535861
*			  and for
*			  M. Benedek, P. Biró, D. Paulusma and X. Ye (2023)
*			  Computing Balanced Solutions for Large Kidney Exchange Schemes
*			  https://arxiv.org/abs/2109.06788
* 
*             Banzhaf value allocation with equal country sizes
*
*    @author Márton Benedek
*    @version 1.0 21/03/2023
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
#include <lemon/list_graph.h> // needed for ListDigraph
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

void LMC(unsigned short int& Q, unsigned short int& periods, bool& dispy, unsigned short int& N, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, unsigned int& S, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, unsigned short int& Vp, double& t1, double& game_time, double& t0, vector<double>& init_alloc, double& prec, double& init_alloc_time, vector<vector<double>>& init_allocations, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, bool& d1, vector<double>& credit, double& matching_time, vector<vector<unsigned short int>>& solution, vector<unsigned short int>& node_arrives, vector<bool>& unique_impu_LMC, unsigned short int& unique_impus_LMC, vector<unsigned short int>& v_accum, vector<double>& init_alloc_accum, vector<double>& s_accum, double& y_core_dist, double& s_core_dist, bool& core_dist, bool& unique_imputation);
void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, vector<double>& init_alloc, vector<double>& credit, double& prec, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, double& t0, bool& d1, unsigned short int& I, vector<vector<double>>& init_allocations, vector<vector<unsigned short int>>& solution, vector<bool>& unique_impu_LM, unsigned short int& unique_impus_LM, vector<unsigned short int>& v_accum, vector<double>& init_alloc_accum, vector<double>& s_accum, double& y_core_dist, double& s_core_dist, bool& core_dist, bool& unique_imputation);
void arbitrary_matching(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, vector<double>& init_alloc, vector<unsigned short int>& v, unsigned int& S, double& prec, double& t0, vector<vector<unsigned short int>>& solution_NC, vector<vector<unsigned short int>>& solution, vector<vector<double>>& initial_allocations, vector<bool>& unique_impu_rand, unsigned short int& unique_impus_rand, vector<unsigned short int>& v_accum, vector<double>& init_alloc_accum, vector<double>& s_accum, double& y_core_dist, double& s_core_dist, bool& core_dist, bool& unique_imputation);
void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& disp, bool& dispy, double& prec, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, vector<unsigned short int>& s_ideal, double& t0, bool& d1, vector<unsigned short int>& v_ideal, vector<double>& init_alloc_ideal, vector<unsigned short int>& s_ideal_d1, double& ideal_time, double& ideal_d1_time, bool& ideal_unique_impu, double y_core_dist, double& s_core_dist, double& s_core_dist_d1, bool& core_dist, bool& unique_imputation);
void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& init_alloc, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, bool& disp, unsigned short int& Vp, bool& d1);
void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& deviation, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& init_alloc, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp, bool& disp, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2);
void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, bool& dispy, unsigned int& S, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, bool& unique_impu, double& prec, bool& unique_imputation);
void norm_banzhaf(vector<double>& banz, vector<unsigned short int>& v, unsigned short int& n, unsigned int& s);
void insertion_sort(vector<unsigned short int>& w, vector<double>& deviation, unsigned short int& N);
void de2bi_card(unsigned int& k, vector<bool>& a, unsigned short int& n, unsigned short int& card);
void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n);
double core_distance(vector<double>& x, vector<unsigned short int>& v, unsigned short int& n, unsigned int& S);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, bool& disp, vector<ListGraph::Node>& c, vector<unsigned short int>& s);
void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods, bool& disp);
void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, vector<ListGraph::Node>& c, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, ListGraph& g_ideal, vector<unsigned short int>& node_arrives, unsigned short int& no_of_nodes);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, bool& disp, unsigned int& m, unsigned short int& no_of_nodes);
bool is_next_char_digit(string& line, unsigned int l);
unsigned int char2uint(char& p);
double cpuTime();

int main() {
	cout << "I solemnly swear that I am up to no good." << endl;
	double t0 = cpuTime();
	// input parameters and data
	unsigned short int N = 4; // number of countries/players
	unsigned short int inst = 0; // instance number, integer between 0 and 99
	bool dispy = false; // true: information in terminal while running
	bool disp = false; // true: extremely detailed information while running, avoid with large graphs
	bool unique_imputation = false;
	bool core_dist = false;
	bool ideal_scenario = false;

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
	ListGraph g;
	vector<ListGraph::Node> c(no_of_nodes);
	xml_parser(line, node_labels, label_positions, c, k, g, arc_in, arc_out, disp, m, no_of_nodes);
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
	initial_pairs(Vp, N, active_nodes, c, disp);
	vector<unsigned short int> node_arrives(no_of_nodes, 0);
	arrival_times(node_arrives, Vp, N, active_nodes, c, periods, disp);
	if (disp) {
		cout << endl;
		vector<unsigned short int> bla(periods, 0);
		for (unsigned short int i = 0; i < no_of_nodes; i++)
			bla[node_arrives[i]]++;
		cout << "no of arrivals: ";
		for (unsigned short int i = 0; i < periods; i++)
			cout << bla[i] << " ";
		cout << endl << endl;
	}
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

	unsigned int S = pow(2, N) - 2;
	vector<unsigned short int> v(S + 1, 0);
	vector<unsigned short int> v_accum(S + 1, 0);
	vector<unsigned short int> s(N, 0);
	vector<double> s_accum(N, 0);
	vector<vector<unsigned short int>> solution_LMC(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> solution_LM(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> solution_d1C(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> solution_d1(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> solution_rand(periods, vector<unsigned short int>(N, 0));
	vector<vector<unsigned short int>> solution_NC(periods, vector<unsigned short int>(N, 0));
	double prec = pow(10, -7);
	vector<double> init_alloc(N, 0);
	vector<double> init_alloc_accum(N, 0);
	vector<vector<double>> init_alloc_LMC(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_LM(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_d1C(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_d1(periods, vector<double>(N, 0));
	vector<vector<double>> init_alloc_rand(periods, vector<double>(N, 0));
	vector<bool> unique_impu_LMC(periods, false);
	vector<bool> unique_impu_LM(periods, false);
	vector<bool> unique_impu_d1C(periods, false);
	vector<bool> unique_impu_d1(periods, false);
	vector<bool> unique_impu_rand(periods, false);
	unsigned short int unique_impus_LMC = 0;
	unsigned short int unique_impus_LM = 0;
	unsigned short int unique_impus_d1C = 0;
	unsigned short int unique_impus_d1 = 0;
	unsigned short int unique_impus_rand = 0;
	double y_core_dist_LMC = 0;
	double y_core_dist_LM = 0;
	double y_core_dist_d1C = 0;
	double y_core_dist_d1 = 0;
	double y_core_dist_rand = 0;
	double s_core_dist_LMC = 0;
	double s_core_dist_LM = 0;
	double s_core_dist_d1C = 0;
	double s_core_dist_d1 = 0;
	double s_core_dist_rand = 0;
	double game_time = 0;
	double init_alloc_time = 0;
	double matching_time = 0;
	unsigned short int I1 = 0;
	unsigned short int I11 = 0;
	unsigned short int I2 = 0;
	vector<double> credit(N, 0);
	vector<double> deviation(N, 0);
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

	// IDEAL scenario:
	vector<unsigned short int> s_ideal(N, 0);
	vector<double> init_alloc_ideal(N, 0);
	vector<unsigned short int> s_ideal_d1(N, 0);
	bool ideal_unique_impu = false;
	double ideal_time = 0;
	double ideal_d1_time = 0;
	double y_core_dist_ideal = 0;
	double s_core_dist_ideal = 0;
	double s_core_dist_ideal_d1 = 0;
	if (ideal_scenario) {
		t0 = cpuTime();
		ideal_matching(S, g_ideal, N, Vp, c, disp, dispy, prec, deviation, pos, w, p, lb, ub, opt, no_of_nodes, s_ideal, t0, d1, v, init_alloc_ideal, s_ideal_d1, ideal_time, ideal_d1_time, ideal_unique_impu, y_core_dist_ideal, s_core_dist_ideal, s_core_dist_ideal_d1, core_dist, unique_imputation);
		ideal_time += init_time + read_time + rand_time + graph_time;
		ideal_d1_time += init_time + read_time + rand_time + graph_time;
		cout << "ideal done... ";
		t0 = cpuTime();
	}

	LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, game_time, t0, init_alloc, prec, init_alloc_time, init_alloc_LMC, I1, I11, I2, edge_card_weight, deviation, pos, w, p, lb, ub, opt, d1, credit, matching_time, solution_LMC, node_arrives, unique_impu_LMC, unique_impus_LMC, v_accum, init_alloc_accum, s_accum, y_core_dist_LMC, s_core_dist_LMC, core_dist, unique_imputation);
	double LMC_time = init_time + game_time + init_alloc_time + matching_time + read_time + rand_time + graph_time;
	cout << "lexmin+c done... ";

	t0 = cpuTime();
	Q = 0;
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
	double tmp = 0;
	unsigned short int I = 0;
	d1 = true;
	LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, tmp, tmp, init_alloc, prec, tmp, init_alloc_d1C, I, I, I, edge_card_weight, deviation, pos, w, p, lb, ub, opt, d1, credit, tmp, solution_d1C, node_arrives, unique_impu_d1C, unique_impus_d1C, v_accum, init_alloc_accum, s_accum, y_core_dist_d1C, s_core_dist_d1C, core_dist, unique_imputation);
	t1 = cpuTime();
	double d1C_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	d1 = false;
	cout << "d1+c done... ";

	t0 = cpuTime();
	arbitrary_matching(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, init_alloc, v, S, prec, t0, solution_NC, solution_rand, init_alloc_rand, unique_impu_rand, unique_impus_rand, v_accum, init_alloc_accum, s_accum, y_core_dist_rand, s_core_dist_rand, core_dist, unique_imputation);
	t1 = cpuTime();
	double arbitrary_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	cout << "arbitrary matching done... ";

	t0 = cpuTime();
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, init_alloc, credit, prec, deviation, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_LM, solution_LM, unique_impu_LM, unique_impus_LM, v_accum, init_alloc_accum, s_accum, y_core_dist_LM, s_core_dist_LM, core_dist, unique_imputation);
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
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, init_alloc, credit, prec, deviation, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_d1, solution_d1, unique_impu_d1, unique_impus_d1, v_accum, init_alloc_accum, s_accum, y_core_dist_d1, s_core_dist_d1, core_dist, unique_imputation);
	t1 = cpuTime();
	double d1_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	d1 = false;
	cout << "d1 done!" << endl;

	res.open("results.txt");
	//res.open("results.txt", ofstream::out | ofstream::trunc);
	res << fixed << setprecision(0) << seed << endl << endl;
	res << fixed << setprecision(17) << read_time << endl << graph_time << endl << rand_time << endl << init_time << endl << game_time << endl << init_alloc_time << endl << matching_time << endl << endl;
	res << fixed << setprecision(17) << LMC_time << endl << d1C_time << endl << LM_time << endl << d1_time << endl << arbitrary_time << endl;
	if (ideal_scenario)
		res << ideal_time << endl << ideal_d1_time << endl;
	res << endl;
	if (unique_imputation)
		res << fixed << setprecision(0) << unique_impus_LMC << endl << unique_impus_d1C << endl << unique_impus_LM << endl << unique_impus_d1 << endl << unique_impus_rand << endl << endl;
	if (core_dist) {
		res << fixed << setprecision(17) << y_core_dist_LMC << endl << y_core_dist_d1C << endl << y_core_dist_LM << endl << y_core_dist_d1 << endl << y_core_dist_rand << endl;
		if (ideal_scenario)
			res << y_core_dist_ideal << endl;
		res << endl;
		res << fixed << setprecision(17) << s_core_dist_LMC << endl << s_core_dist_d1C << endl << s_core_dist_LM << endl << s_core_dist_d1 << endl << s_core_dist_rand << endl;
		if (ideal_scenario)
			res << s_core_dist_ideal << endl << s_core_dist_ideal_d1 << endl;
		res << endl;
	}
	res << fixed << setprecision(0) << I1 << endl << I11 << endl << I2 << endl << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_LMC[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << solution_LMC[Q][i] << endl;
		res << endl;
		if (unique_imputation) {
			res << fixed << setprecision(0) << unique_impu_LMC[Q] << endl;
			res << endl;
		}
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_d1C[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << solution_d1C[Q][i] << endl;
		res << endl;
		if (unique_imputation) {
			res << fixed << setprecision(0) << unique_impu_d1C[Q] << endl;
			res << endl;
		}
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_LM[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << solution_LM[Q][i] << endl;
		res << endl;
		if (unique_imputation) {
			res << fixed << setprecision(0) << unique_impu_LM[Q] << endl;
			res << endl;
		}
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_d1[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << solution_d1[Q][i] << endl;
		res << endl;
		if (unique_imputation) {
			res << fixed << setprecision(0) << unique_impu_d1[Q] << endl;
			res << endl;
		}
	}
	res << endl;
	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_rand[Q][i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << solution_rand[Q][i] << endl;
		res << endl;
		if (unique_imputation) {
			res << fixed << setprecision(0) << unique_impu_rand[Q] << endl;
			res << endl;
		}
	}
	res << endl;
	if (ideal_scenario) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(17) << init_alloc_ideal[i] << endl;
		res << endl;
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << s_ideal[i] << endl;
		res << endl;

		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << s_ideal_d1[i] << endl;
		res << endl;
		if (unique_imputation) {
			res << fixed << setprecision(0) << ideal_unique_impu << endl;
			res << endl;
		}
	}

	for (unsigned short int Q = 0; Q < periods; Q++) {
		for (unsigned short int i = 0; i < N; i++)
			res << fixed << setprecision(0) << solution_NC[Q][i] << endl;
		res << endl;
	}
	res.close();

	cout << "Mischief managed!" << endl;
	return 0;
}

void LMC(unsigned short int& Q, unsigned short int& periods, bool& dispy, unsigned short int& N, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, unsigned int& S, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, unsigned short int& Vp, double& t1, double& game_time, double& t0, vector<double>& init_alloc, double& prec, double& init_alloc_time, vector<vector<double>>& init_allocations, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, bool& d1, vector<double>& credit, double& matching_time, vector<vector<unsigned short int>>& solution, vector<unsigned short int>& node_arrives, vector<bool>& unique_impu_LMC, unsigned short int& unique_impus_LMC, vector<unsigned short int>& v_accum, vector<double>& init_alloc_accum, vector<double>& s_accum, double& y_core_dist, double& s_core_dist, bool &core_dist, bool& unique_imputation) {
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
		bool unique_impu = false;
		coop_game(g, v, s, c, disp, dispy, S, Vp, N, active_nodes, leaving, unique_impu, prec, unique_imputation);
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < S + 1; i++)
					v_accum[i] = v[i];
			}
			else {
				for (unsigned short int i = 0; i < S + 1; i++)
					v_accum[i] += v[i];
			}
		}
		if (unique_imputation) {
			unique_impu_LMC[Q] = unique_impu;
			unique_impus_LMC += unique_impu;
		}
		t1 = cpuTime();
		game_time += t1 - t0;
		t0 = cpuTime();
		norm_banzhaf(init_alloc, v, N, S);
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < N; i++)
					init_alloc_accum[i] = init_alloc[i];
			}
			else {
				for (unsigned short int i = 0; i < N; i++)
					init_alloc_accum[i] += init_alloc[i];
			}
		}
		if (dispy) {
			cout << "Banzhaf: ";
			for (unsigned short int i = 0; i < N; i++) {
				cout << init_alloc[i] << " ";
			}
			cout << endl;
		}
		t1 = cpuTime();
		init_alloc_time += t1 - t0;
		init_allocations[Q] = init_alloc;
		t0 = cpuTime();
		lex_min_matching(N, v[S], S, init_alloc, s, no_of_active_nodes, I1, I11, I2, g, edge_card_weight, c, prec, dispy, deviation, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
		t1 = cpuTime();
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < N; i++)
					s_accum[i] = s[i];
			}
			else {
				for (unsigned short int i = 0; i < N; i++)
					s_accum[i] += s[i];
			}
		}
		matching_time += t1 - t0;
		solution[Q] = s;
		t0 = cpuTime();
		for (unsigned short int i = 0; i < N; i++)
			credit[i] += init_alloc[i] - s[i];
		if (dispy) {
			cout << "Credits: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << credit[i] << " ";
			cout << endl;
		}
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, disp, c, s);
		if (dispy) {
			cout << endl << "Press enter to continue." << endl;
			cin.get();
		}
	}
	if (core_dist) {
		y_core_dist = core_distance(init_alloc_accum, v_accum, N, S);
		s_core_dist = core_distance(s_accum, v_accum, N, S);
	}
	return;
}

void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, vector<double>& init_alloc, vector<double>& credit, double& prec, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, double& t0, bool& d1, unsigned short int& I, vector<vector<double>>& init_allocations, vector<vector<unsigned short int>>& solution, vector<bool>& unique_impu_LM, unsigned short int& unique_impus_LM, vector<unsigned short int>& v_accum, vector<double>& init_alloc_accum, vector<double>& s_accum, double& y_core_dist, double& s_core_dist, bool &core_dist, bool &unique_imputation) {
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
		bool unique_impu = false;
		coop_game(g, v, s, c, disp, dispy, S, Vp, N, active_nodes, leaving, unique_impu, prec, unique_imputation);
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < S + 1; i++)
					v_accum[i] = v[i];
			}
			else {
				for (unsigned short int i = 0; i < S + 1; i++)
					v_accum[i] += v[i];
			}
		}
		if (unique_imputation) {
			unique_impu_LM[Q] = unique_impu;
			unique_impus_LM += unique_impu;
		}
		norm_banzhaf(init_alloc, v, N, S);
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < N; i++)
					init_alloc_accum[i] = init_alloc[i];
			}
			else {
				for (unsigned short int i = 0; i < N; i++)
					init_alloc_accum[i] += init_alloc[i];
			}
		}
		init_allocations[Q] = init_alloc;
		if (dispy) {
			cout << "Banzhaf: ";
			for (unsigned short int i = 0; i < N; i++) {
				cout << init_alloc[i] << " ";
			}
			cout << endl;
		}
		lex_min_matching(N, v[S], S, init_alloc, s, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, deviation, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
		solution[Q] = s;
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < N; i++)
					s_accum[i] = s[i];
			}
			else {
				for (unsigned short int i = 0; i < N; i++)
					s_accum[i] += s[i];
			}
		}
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, disp, c, s);
		if (dispy) {
			cout << endl << "Press enter to continue." << endl;
			cin.get();
		}
	}
	if (core_dist) {
		y_core_dist = core_distance(init_alloc_accum, v_accum, N, S);
		s_core_dist = core_distance(s_accum, v_accum, N, S);
	}
	return;
}

void arbitrary_matching(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, vector<double>& init_alloc, vector<unsigned short int>& v, unsigned int& S, double& prec, double& t0, vector<vector<unsigned short int>>& solution_NC, vector<vector<unsigned short int>>& solution, vector<vector<double>>& initial_allocations, vector<bool>& unique_impu_rand, unsigned short int& unique_impus_rand, vector<unsigned short int>& v_accum, vector<double>& init_alloc_accum, vector<double>& s_accum, double& y_core_dist, double& s_core_dist, bool &core_dist, bool &unique_imputation) {
	if (dispy)
		cout << " --== Arbitrary matching == -- " << endl;
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
			changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, disp, c, s);
		if (disp) {
			cout << "Number of active nodes: ";
			for (unsigned short int i = 0; i < N; i++)
				cout << no_of_active_nodes[i] << " ";
			cout << endl;
		}
		bool unique_impu = false;
		coop_game(g, v, s, c, disp, disp, S, Vp, N, active_nodes, leaving, unique_impu, prec, unique_imputation);
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < S + 1; i++)
					v_accum[i] = v[i];
			}
			else {
				for (unsigned short int i = 0; i < S + 1; i++)
					v_accum[i] += v[i];
			}
		}
		if (unique_imputation) {
			unique_impu_rand[Q] = unique_impu;
			unique_impus_rand += unique_impu;
		}
		norm_banzhaf(init_alloc, v, N, S);
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < N; i++)
					init_alloc_accum[i] = init_alloc[i];
			}
			else {
				for (unsigned short int i = 0; i < N; i++)
					init_alloc_accum[i] += init_alloc[i];
			}
		}
		initial_allocations[Q] = init_alloc;
		for (unsigned short int i = 0; i < N; i++)
			solution_NC[Q][i] = v[pow(2, i) - 1];
		if (dispy) {
			cout << "Banzhaf: ";
			for (unsigned short int i = 0; i < N; i++) {
				cout << init_alloc[i] << " ";
			}
			cout << endl;
		}
		solution[Q] = s;
		if (dispy) {
			cout << "s[" << Q << "]=";
			for (unsigned int i = 0; i < N; i++)
				cout << s[i] << " ";
			cout << endl;
		}
		if (core_dist) {
			if (Q == 0) {
				for (unsigned short int i = 0; i < N; i++)
					s_accum[i] = s[i];
			}
			else {
				for (unsigned short int i = 0; i < N; i++)
					s_accum[i] += s[i];
			}
		}
	}
	if (core_dist) {
		y_core_dist = core_distance(init_alloc_accum, v_accum, N, S);
		s_core_dist = core_distance(s_accum, v_accum, N, S);
	}
	return;
}

void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& disp, bool& dispy, double& prec, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, vector<unsigned short int>& s_ideal, double& t0, bool& d1, vector<unsigned short int>& v_ideal, vector<double>& init_alloc_ideal, vector<unsigned short int>& s_ideal_d1, double& ideal_time, double& ideal_d1_time, bool& ideal_unique_impu, double y_core_dist, double& s_core_dist, double& s_core_dist_d1, bool &core_dist, bool& unique_imputation) {
	vector<bool> a(N, false);
	for (unsigned int i = 0; i < S; i++) {
		ListGraph::NodeMap<bool> coal(g, false);
		de2bi(i, a, N);
		for (unsigned short int j = 0; j < N; j++) {
			if (a[j]) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++)
					coal[c[k]] = true;
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m(FilterNodes<ListGraph>(g, coal));
		coal_m.run();
		v_ideal[i] = 2 * coal_m.matchingSize();
	}
	MaxMatching<ListGraph> grand_coal(g);
	grand_coal.run();
	v_ideal[S] = 2 * grand_coal.matchingSize();
	if (unique_imputation) {
		double surplus = v_ideal[S];
		for (unsigned short int i = 0; i < N; i++)
			surplus -= v_ideal[pow(2, i) - 1];
		if (surplus > prec)
			ideal_unique_impu = false;
		else
			ideal_unique_impu = true;
	}
	for (unsigned short int i = 0; i < N; i++) {
		for (unsigned short int j = i * Vp; j < (i + 1) * Vp; j++) {
			if (!(grand_coal.mate(c[j]) == INVALID))
				s_ideal[i]++;
		}
	}
	s_ideal_d1 = s_ideal;
	vector<unsigned short int> no_of_active_nodes(N, Vp);
	ListGraph::NodeMap<bool> active_nodes(g, true);
	unsigned short int I = 0;
	ListGraph::EdgeMap<unsigned short int> edge_card_weight(g, 1);
	vector<double> credit(N, 0);
	vector<bool> leaving(no_of_nodes, false);
	norm_banzhaf(init_alloc_ideal, v_ideal, N, S);
	if (core_dist)
		y_core_dist = core_distance(init_alloc_ideal, v_ideal, N, S);
	ideal_time = cpuTime() - t0;
	ideal_d1_time = ideal_time;
	t0 = cpuTime();
	lex_min_matching(N, v_ideal[S], S, init_alloc_ideal, s_ideal, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, deviation, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
	vector<double> s_ideal_d(N, s_ideal[0]);
	for (unsigned short int i = 1; i < N; i++)
		s_ideal_d[i] = s_ideal[i];
	if (core_dist)
		s_core_dist = core_distance(s_ideal_d, v_ideal, N, S);
	ideal_time += cpuTime() - t0;
	t0 = cpuTime();
	d1 = true;
	lex_min_matching(N, v_ideal[S], S, init_alloc_ideal, s_ideal_d1, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, deviation, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
	for (unsigned short int i = 0; i < N; i++)
		s_ideal_d[i] = s_ideal_d1[i];
	if (core_dist)
		s_core_dist_d1 = core_distance(s_ideal_d, v_ideal, N, S);
	ideal_d1_time += cpuTime() - t0;
	d1 = false;
	return;
}

void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& init_alloc, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& deviation, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, bool& disp, unsigned short int& Vp, bool& d1) {
	p = 0;
	if (disp)
		cout << "w*: " << grandcoal << endl;
	for (unsigned short int i = 0; i < N; i++) {
		if (init_alloc[i] + credit[i] - s[i] >= -prec) {
			deviation[i] = init_alloc[i] + credit[i] - s[i];
			pos[i] = true;
		}
		else {
			deviation[i] = s[i] - init_alloc[i] - credit[i];
			pos[i] = false;
		}
	}
	insertion_sort(w, deviation, N);
	if (dispy) {
		cout << "deviations: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << deviation[i] << " ";
		cout << endl << "largest difference: " << deviation[w[0]] << endl;
		cout << "order: ";
		for (unsigned short int i = 0; i < N; i++)
			cout << w[i] << " ";
		cout << endl;
	}
	if (deviation[w[0]] <= 0.5 + prec)
		p = N;
	while (p < N) {
		for (unsigned short int i = p; i < N; i++) {
			if (i == p || abs(deviation[w[i]] - deviation[w[p]]) < prec) {
				if (pos[w[i]]) {
					lb[w[i]] = s[w[i]];
					if (init_alloc[w[i]] + credit[w[i]] + deviation[w[p]] > (double)no_of_active_nodes[w[i]] - prec) {
						ub[w[i]] = no_of_active_nodes[w[i]];
					}
					else {
						ub[w[i]] = (short int)(init_alloc[w[i]] + credit[w[i]] + deviation[w[p]] + prec);
					}
				}
				else {
					ub[w[i]] = s[w[i]];
					if (init_alloc[w[i]] + credit[w[i]] - deviation[w[p]] < -prec) {
						lb[w[i]] = 0;
					}
					else {
						if (abs((short int)(init_alloc[w[i]] + credit[w[i]] - deviation[w[p]] + prec) - init_alloc[w[i]] - credit[w[i]] + deviation[w[p]]) < prec) {
							lb[w[i]] = (short int)(init_alloc[w[i]] + credit[w[i]] - deviation[w[p]] + prec);
						}
						else {
							lb[w[i]] = (short int)(init_alloc[w[i]] + credit[w[i]] - deviation[w[p]]) + 1;
						}
					}
				}
			}
			else {
				if (init_alloc[w[i]] + credit[w[i]] - deviation[w[p]] < -prec) {
					lb[w[i]] = 0;
				}
				else {
					lb[w[i]] = (short int)(init_alloc[w[i]] + credit[w[i]] - deviation[w[p]] + prec) + 1;
				}
				if (init_alloc[w[i]] + credit[w[i]] + deviation[w[p]] > (double)no_of_active_nodes[w[i]] - prec) {
					ub[w[i]] = no_of_active_nodes[w[i]];
				}
				else {
					if (abs((short int)(init_alloc[w[i]] + credit[w[i]] + deviation[w[p]] + prec) - init_alloc[w[i]] - credit[w[i]] - deviation[w[p]]) < prec) {
						ub[w[i]] = (short int)(init_alloc[w[i]] + credit[w[i]] + deviation[w[p]] + prec) - 1;
					}
					else {
						ub[w[i]] = (short int)(init_alloc[w[i]] + credit[w[i]] + deviation[w[p]]);
					}
				}
			}
			if (disp)
				cout << "COUNTRY " << w[i] << " lb: " << lb[w[i]] << " ; ub: " << ub[w[i]] << endl;
		}
		new_matching(lb, ub, g, edge_card_weight, w, p, deviation, c, opt, s, grandcoal, prec, dispy, init_alloc, no_of_active_nodes, N, active_nodes, pos, leaving, credit, Vp, disp, I1, I11, I2);
		if (p > 0 && d1)
			break;
	}
	return;
}

void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& deviation, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& init_alloc, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp, bool& disp, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2) {
	unsigned short int lb_old = lb[w[p]];
	unsigned short int ub_old = ub[w[p]];
	if (deviation[w[p]] >= 1 - prec) {
		if ((lb[w[p]] >= init_alloc[w[p]] + credit[w[p]] - deviation[w[p]] - prec) && (lb[w[p]] <= init_alloc[w[p]] + credit[w[p]] - deviation[w[p]] + prec) && (lb[w[p]] < no_of_active_nodes[w[p]])) {
			lb[w[p]] += 1;
		}
		if ((ub[w[p]] >= init_alloc[w[p]] + credit[w[p]] + deviation[w[p]] - prec) && (ub[w[p]] <= init_alloc[w[p]] + credit[w[p]] + deviation[w[p]] + prec) && (ub[w[p]] > prec)) {
			ub[w[p]] -= 1;
		}
	}
	else {
		if ((init_alloc[w[p]] + credit[w[p]] - s[w[p]] >= deviation[w[p]] - prec) && (init_alloc[w[p]] + credit[w[p]] - s[w[p]] <= deviation[w[p]] + prec) && (s[w[p]] < no_of_active_nodes[w[p]])) {
			lb[w[p]] = s[w[p]] + 1;
		}
		if ((s[w[p]] - init_alloc[w[p]] - credit[w[p]] >= deviation[w[p]] - prec) && (s[w[p]] - init_alloc[w[p]] - credit[w[p]] <= deviation[w[p]] + prec) && (s[w[p]] > prec)) {
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
				cout << "largest (unfinished) difference: " << deviation[w[p]] << endl;
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
			opt = 2 * max_weight.matchingWeight();
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
				if (init_alloc[i] + credit[i] - s[i] >= 0) {
					deviation[i] = init_alloc[i] + credit[i] - s[i];
					pos[i] = true;
				}
				else {
					deviation[i] = s[i] - init_alloc[i] - credit[i];
					pos[i] = false;
				}
			}
			insertion_sort(w, deviation, N);
			if (dispy) {
				cout << "deviation: ";
				for (unsigned short int i = 0; i < N; i++)
					cout << deviation[i] << " ";
				cout << endl << "largest (unfinished) difference: " << deviation[w[p]] << endl;
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
				if (abs(init_alloc[w[p]] + credit[w[p]] - lb_old) > deviation[w[p]]) {
					lb[w[p]] = s[w[p]];
				}
				else {
					lb[w[p]] = lb_old;
				}
			}
			if (s[w[p]] > ub[w[p]]) {
				if (abs(init_alloc[w[p]] + credit[w[p]] - ub_old) > deviation[w[p]]) {
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
					cout << "largest (unfinished) difference: " << deviation[w[p]] << endl;
			}
		}
		if (!(count % 2 == 0))
			g.erase(ap);
	}
	if (p < N) {
		if (deviation[w[p]] <= 0.5 + prec)
			p = N;
	}
	return;
}

void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, bool& dispy, unsigned int& S, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, bool& unique_impu, double& prec, bool &unique_imputation) {
	vector<bool> a(N, false);
	for (unsigned int i = 0; i < S; i++) {
		ListGraph::NodeMap<bool> coal(g, false);
		de2bi(i, a, N);
		for (unsigned short int j = 0; j < N; j++) {
			if (a[j]) {
				for (unsigned short int k = j * Vp; k < (j + 1) * Vp; k++) {
					if (active_nodes[c[k]])
						coal[c[k]] = true;
				}
			}
		}
		MaxMatching<FilterNodes<ListGraph>> coal_m(FilterNodes<ListGraph>(g, coal));
		coal_m.run();
		v[i] = 2 * coal_m.matchingSize();
		if (disp)
			cout << "coal_" << i + 1 << ": " << v[i] << endl;
	}
	FilterNodes<ListGraph> sg(g, active_nodes);
	MaxMatching<FilterNodes<ListGraph>> grand_coal(sg);
	grand_coal.run();
	v[S] = 2 * grand_coal.matchingSize();
	if (unique_imputation) {
		double surplus = v[S];
		for (unsigned short int i = 0; i < N; i++)
			surplus -= v[pow(2, i) - 1];
		if (surplus > prec)
			unique_impu = false;
		else
			unique_impu = true;
	}
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

void norm_banzhaf(vector<double>& banz, vector<unsigned short int>& v, unsigned short int& n, unsigned int& s) {
	double w = 1 / pow(2, n - 1);
	vector<double> expo(n, 1);
	for (unsigned short int j = 1; j < n; j++)
		expo[j] = pow(2, j);
	vector<bool> a(n, false);
	unsigned short int k = 0;
	for (unsigned short int i = 0; i < n; i++)
		banz[i] = w * (double)(v[expo[i] - 1]);
	for (unsigned int i = 0; i < s; i++) {
		de2bi_card(i, a, n, k);
		for (unsigned short int j = 0; j < n; j++) {
			if (!a[j])
				banz[j] += w * (v[i + expo[j]] - v[i]);
		}
	}
	double norm = banz[0];
	for (unsigned short int i = 1; i < n; i++)
		norm += banz[i];
	for (unsigned short int i = 0; i < n; i++)
		banz[i] = banz[i] * v[s] / norm;
	return;
}

void insertion_sort(vector<unsigned short int>& w, vector<double>& deviation, unsigned short int& N) {
	w[0] = 0;
	for (unsigned short int i = 1; i < N; i++) {
		if (deviation[i] <= deviation[w[i - 1]]) {
			w[i] = i;
		}
		else {
			w[i] = w[i - 1];
			if (i == 1) {
				w[0] = i;
			}
			else {
				for (unsigned short int j = i - 2; j >= 0; j--) {
					if (deviation[i] <= deviation[w[j]]) {
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

void de2bi_card(unsigned int& k, vector<bool>& a, unsigned short int& n, unsigned short int& card) {
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

double core_distance(vector<double>& x, vector<unsigned short int>& v, unsigned short int& n, unsigned int& S) {
	double eps = x[0] - v[0];
	if (x[1] - v[1] < eps)
		eps = x[1] - v[1];
	vector<bool> a(n, false);
	double xS = 0;
	for (unsigned int i = 2; i < S; i++) {
		de2bi(i, a, n);
		for (unsigned short int j = 0; j < n; j++)
			if (a[j])
				xS += x[j];
		if (xS - v[i] < eps)
			eps = xS - v[i];
		xS = 0;
	}
	return eps;
}

void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, bool& disp, vector<ListGraph::Node>& c, vector<unsigned short int>& s) {
	for (unsigned short int i = 0; i < N; i++) {
		s[i] = 0;
		for (unsigned short int j = 0; j < Vp; j++) {
			if (leaving[i * Vp + j]) {
				active_nodes[c[i * Vp + j]] = false;
				no_of_active_nodes[i]--;
				leaving[i * Vp + j] = false;
				if (disp)
					cout << "Node " << i * Vp + j << " is leaving at " << Q << " becaused its matched" << endl;
			}
			else {
				if (active_nodes[c[i * Vp + j]] && node_arrives[i * Vp + j] == Q - 4) {
					active_nodes[c[i * Vp + j]] = false;
					no_of_active_nodes[i]--;
					if (disp)
						cout << "Node " << i * Vp + j << " is leaving at " << Q << " because spent 1 year (unmatched) in the program" << endl;
				}
			}
			if (node_arrives[i * Vp + j] == Q) {
				active_nodes[c[i * Vp + j]] = true;
				no_of_active_nodes[i]++;
				if (disp)
					cout << "Node " << i * Vp + j << " is entering at " << Q << endl;
			}
		}
	}
	return;
}

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp) {
	unsigned short int coal = rand() % Vp;
	unsigned short int count = 0;
	if (disp) {
		cout << "Initial pairs:" << endl;
	}
	for (unsigned short int i = 0; i < N; i++) {
		if (disp) {
			cout << "Country " << i << " : ";
		}
		while (count < Vp / 4) {
			if (active_nodes[c[i * Vp + coal]]) {
				coal = rand() % Vp;
			}
			else {
				active_nodes[c[i * Vp + coal]] = true;
				count++;
				if (disp) {
					cout << i * Vp + coal << " (" << count << ") ";
				}
				coal = rand() % Vp;
			}
		}
		if (disp) {
			cout << endl;
		}
		count = 0;
	}
	return;
}

void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods, bool& disp) {
	for (unsigned short int i = 0; i < N; i++) {
		if (disp)
			cout << "Country " << i << " arrivals: ";
		for (unsigned short int j = 0; j < Vp; j++) {
			if (!(active_nodes[c[i * Vp + j]])) {
				node_arrives[i * Vp + j] = rand() % (periods - 1) + 1;
			}
			if (disp)
				cout << node_arrives[i * Vp + j] << " ";
		}
		if (disp)
			cout << endl;
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

void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, bool& disp, unsigned int& m, unsigned short int& no_of_nodes) {
	//disp = true;
	unsigned int l = 6;
	unsigned short int n = 0;
	while (l < line.size() - 7) {
		if (line[l] == '<' && line[l + 1] == 'e') {
			l = l + 17;
			n++;
			if (!is_next_char_digit(line, l)) {
				//node_labels.push_back(char2uint(line[l]));
				node_labels[n - 1] = char2uint(line[l]);
			}
			else {
				if (!is_next_char_digit(line, l + 1)) {
					//node_labels.push_back(10*char2uint(line[l])+char2uint(line[l+1]));
					node_labels[n - 1] = 10 * char2uint(line[l]) + char2uint(line[l + 1]);
					l++;
				}
				else {
					if (!is_next_char_digit(line, l + 2)) {
						//node_labels.push_back(100*char2uint(line[l])+10*char2uint(line[l+1])+char2uint(line[l+2]));
						node_labels[n - 1] = 100 * char2uint(line[l]) + 10 * char2uint(line[l + 1]) + char2uint(line[l + 2]);
						l = l + 2;
					}
					else {
						//if (!is_next_char_digit(line, l + 3)){
						//node_labels.push_back(1000*char2uint(line[l])+100*char2uint(line[l+1])+10*char2uint(line[l+2])+char2uint(line[l+3]));
						node_labels[n - 1] = 1000 * char2uint(line[l]) + 100 * char2uint(line[l + 1]) + 10 * char2uint(line[l + 2]) + char2uint(line[l + 3]);
						l = l + 3;
						//}
					}
				}
			}
			if (n + k - 1 == node_labels[n - 1]) {
				//label_positions.push_back(n - 1);
				label_positions[n + k - 1] = n - 1;
			}
			else {
				while (n + k - 1 < node_labels[n - 1]) {
					//label_positions[n + k -1] = 1410065407;
					label_positions[n + k - 1] = 65535;
					label_positions.push_back(0);
					k++;
					//label_positions.push_back(1410065407);
				}
				//label_positions.push_back(n - 1);
				label_positions[n + k - 1] = n - 1;
			}
			//c.push_back(g.addNode());
			//c[n - 1] = g.addNode();

			c[n - 1] = g.addNode();
			//label[c[14][n - 1]] = node_labels[n - 1];
			if (disp)
				cout << "New node with ID: " << node_labels[n - 1] << endl;
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
			//if (disp)
			//cout << "Donor age: " << donor_ages[n - 1] << endl;
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
			if (disp)
				cout << node_labels[n - 1] << " matches with " << arc_in[m - 1];
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
			//if (disp)
			//cout << " with weight " << arc_weight[m - 1] << endl;
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