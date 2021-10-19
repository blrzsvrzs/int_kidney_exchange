/*
*    int_kidney_exchange
*    shapley_nucl_eqsize.cpp
*    Purpose: computational study for 
*             Benedek et al. (2021) - Computing Balanced Solutions for
*             Large Kidney Exchange Schemes
*             https://arxiv.org/abs/2109.06788
*             Shapley value and nucleolus allocations with equal country sizes
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

double cpuTime();
bool is_next_char_digit(string& line, unsigned int l);
unsigned int char2uint(char& p);
void undi_lemon(unsigned int& m, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, vector<unsigned short int>& label_positions, ListGraph& g, vector<ListGraph::Node>& c, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, ListGraph& g_ideal, vector<unsigned short int>& node_arrives, unsigned short int& no_of_nodes);
void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, bool& dispy, unsigned int& S, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving);
void new_matching(vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<unsigned short int>& w, unsigned short int& p, vector<double>& y, vector<ListGraph::Node>& c, double& opt, vector<unsigned short int>& s, unsigned short int& max_match, double& prec, bool& dispy, vector<double>& target, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& pos, vector<bool>& leaving, vector<double>& credit, unsigned short int& Vp, bool& disp, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2);
void xml_parser(string& line, vector<unsigned short int>& node_labels, vector<unsigned short int>& label_positions, vector<ListGraph::Node>& c, unsigned short int& k, ListGraph& g, vector<unsigned int>& arc_in, vector<unsigned int>& arc_out, bool& disp, unsigned int& m, unsigned short int& no_of_nodes);
void shapley(vector<double>& shapl, vector<unsigned short int>& v, unsigned short int& n, unsigned int& s);
void de2bi_card(unsigned int& k, vector<bool>& a, unsigned short int& n, unsigned short int& card);
void de2bi(unsigned int& k, vector<bool>& a, unsigned short int& n);
void insertion_sort(vector<unsigned short int>& w, vector<double>& y, unsigned short int& N);

void nucl(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& x, vector<unsigned short int>& v, double& prec);
void zeros_mem(vector<bool>& a, unsigned short int& n, unsigned int& s, vector<unsigned short int>& zeros);
void excess_init(vector<double>& exc, vector<bool>& unsettled, vector<double>& x, vector<unsigned short int>& v, unsigned int& s, unsigned short int& n, vector<unsigned short int>& zeros);
void nucl_comp(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& excess, double& prec, vector<bool>& unsettled, unsigned short int& iter, unsigned int& piv, unsigned int& sr, double& t, vector<double>& x, vector<bool>& a, double& t1, vector<double>& singleton_bounds, bool& nlsu, double& min_satisfaction, vector<unsigned short int>& zeros);
void vec_min_uns(double& m, vector<double>& x, vector<bool>& unsettled, unsigned int& s);
void tight_coal2(vector<bool>& T2, vector<double>& x, vector<double>& singleton_bounds, double& prec, unsigned short int& n, vector<unsigned int>& T2_coord, vector<bool>& unsettled_p, unsigned short int& t2_size);
void tight_coal(vector<bool>& T, vector<double>& excess, double& epsi, double& prec, unsigned int& s, vector<unsigned int>& T_coord, vector<bool>& unsettled, unsigned int& t_size);
void pivot(double& epsi, unsigned int& s, vector<double>& excess, double& prec, unsigned short int& n, vector<bool>& a, vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& unsettled, unsigned short int& rank, vector<double>& d, vector<double>& x, bool& disp, vector<vector<bool>>& Asettled, unsigned int& piv, unsigned int& sr_count, unsigned short int& iter, vector<bool>& unsettled_p, vector<double>& singleton_bounds, double& epsi_old, bool& nlsu, vector<unsigned short int>& zeros, vector<bool>& T, vector<unsigned int>& T_coord, vector<bool>& T2, vector<unsigned int>& T2_coord, unsigned int& t_size, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2, double& min_satisfaction);
void subr_upd(vector<vector<double>>& Arref, vector<bool>& J, unsigned int& i, unsigned short int& n, double& prec, vector<bool>& U, vector<bool>& U2, unsigned int& sumt, unsigned int& sumt2, vector<bool>& t, vector<bool>& t2, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, unsigned int& tight_size, unsigned short int& tight2_size, unsigned short int& rank, vector<bool>& unsettled, vector<vector<bool>>& Asettled, bool& disp, unsigned int& s, vector<unsigned int>& T_coord, vector<unsigned int>& T2_coord, double& epsi_old, double& epsi, vector<bool>& unsettled_p, bool& settled, glp_prob*& lp, vector<bool>& ar0pos);
void subroutine(vector<bool>& U, vector<bool>& U2, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<vector<double>>& Arref, vector<bool>& J, double& prec, unsigned short int& n, unsigned int& tight_size, unsigned short int& tight2_size, unsigned short int& rank, bool& disp, vector<vector<bool>>& Asettled, unsigned int& sr_count, bool& u, unsigned int& s, vector<unsigned int>& T_coord, vector<unsigned int>& T2_coord, vector<bool>& unsettled, double& epsi_old, double& epsi, vector<bool>& unsettled_p, bool& settled, bool& nlsu);
void step(vector<bool>& T, vector<bool>& T2, vector<bool>& unsettled, vector<bool>& unsettled_p, unsigned int& s, double& epsi, vector<double>& excess, vector<double>& d, unsigned short int& n, vector<double>& x, vector<double>& singleton_bounds, bool& disp, double& prec, vector<unsigned short int>& zeros, vector<unsigned int>& T_coord, unsigned int& t_size, vector<unsigned int>& T2_coord, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2);
void imprdir(vector<double>& d, unsigned short int& n, unsigned int& t_size, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2, unsigned short int& rank, vector<vector<bool>>& Asettled, bool& disp);
bool binrank2(vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& b, unsigned short int& n, unsigned short int& rank);
bool binrank(vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& b, unsigned short int& n, unsigned short int& rank);
void sc_vec_prod(vector<double>& y, double a, vector<double>& x);
void vec_subtract(vector<double>& z, vector<double>& x, vector<double>& y);
void sum_vecb(unsigned int& s, vector<bool>& x);

void initial_pairs(unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp);
void arrival_times(vector<unsigned short int>& node_arrives, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, unsigned short int& periods, bool& disp);
void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, bool& disp, unsigned short int& Vp, bool& d1);
void changing_nodes(ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, vector<unsigned short int>& node_arrives, unsigned short int& Q, bool& disp, vector<ListGraph::Node>& c, vector<unsigned short int>& s);
void arbitrary_matching(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_nucl, vector<double>& target, vector<unsigned short int>& v, unsigned int& S, double& prec, double& t0, vector<vector<unsigned short int>>& actual_alloc_NC, vector<vector<unsigned short int>>& actual_alloc, vector<vector<double>>& initial_alloc);
void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, bool& target_nucl, vector<double>& target, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, double& t0, bool& d1, unsigned short int& I, vector<vector<double>>& init_alloc, vector<vector<unsigned short int>>& actual_alloc);
void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_nucl, bool& disp, bool& dispy, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, vector<unsigned short int>& s_ideal, double& t0, bool& d1, vector<unsigned short int>& v_ideal, vector<double>& init_alloc_ideal, vector<unsigned short int>& s_ideal_d1, double& ideal_time, double& ideal_d1_time);
void LMC(unsigned short int& Q, unsigned short int& periods, bool& dispy, unsigned short int& N, vector<unsigned short int>& no_of_active_nodes, ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, unsigned int& S, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, unsigned short int& Vp, double& t1, double& game_time, bool& target_nucl, double& t0, vector<double>& target, double& prec, double& target_time, vector<vector<double>>& init_alloc, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, bool& d1, vector<double>& credit, double& matching_time, vector<vector<unsigned short int>>& actual_alloc, vector<unsigned short int>& node_arrives);

int main() {
	cout << "I solemnly swear that I am up to no good." << endl;
	double t0 = cpuTime();
	// input parameters and data
	unsigned short int N = 4;
	unsigned short int inst = 0; // instance number, integer between 0 and 99
	bool target_nucl = false; // true: first to evaluate is the nucleolus
	bool second = true; // false: only the first of Shapley/nucleolus will be evaluated
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

	vector<unsigned short int> v(S + 1, 0);
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
	ideal_matching(S, g_ideal, N, Vp, c, target_nucl, disp, dispy, prec, y, pos, w, p, lb, ub, opt, no_of_nodes, s_ideal, t0, d1, v, init_alloc_ideal, s_ideal_d1, ideal_time, ideal_d1_time);
	ideal_time += init_time + read_time + rand_time + graph_time;
	ideal_d1_time += init_time + read_time + rand_time + graph_time;
	cout << "ideal done... ";
	t0 = cpuTime();

	LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, game_time, target_nucl, t0, target, prec, target_time, init_alloc_LMC, I1, I11, I2, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, matching_time, actual_alloc_LMC, node_arrives);
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
	LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, tmp, target_nucl, tmp, target, prec, tmp, init_alloc_d1C, I, I, I, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, tmp, actual_alloc_d1C, node_arrives);
	t1 = cpuTime();
	double d1C_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	d1 = false;
	cout << "d1+c done... ";

	t0 = cpuTime();
	arbitrary_matching(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, target_nucl, target, v, S, prec, t0, actual_alloc_NC, actual_alloc_rand, init_alloc_rand);
	t1 = cpuTime();
	double arbitrary_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	cout << "arbitrary matching done... ";

	t0 = cpuTime();
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_nucl, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_LM, actual_alloc_LM);
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
	lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_nucl, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_d1, actual_alloc_d1);
	t1 = cpuTime();
	double d1_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
	d1 = false;
	cout << "d1 done!" << endl;

	res.open("results.txt", ofstream::out | ofstream::trunc);
	res << fixed << setprecision(0) << seed << endl << endl;
	res << fixed << setprecision(17) << read_time << endl << graph_time << endl << rand_time << endl << init_time << endl << game_time << endl << target_time << endl << matching_time << endl << endl;
	res << fixed << setprecision(17) << LMC_time << endl << d1C_time << endl << LM_time << endl << d1_time << endl << arbitrary_time << endl << ideal_time << endl << ideal_d1_time << endl << endl;
	res << fixed << setprecision(0) << I1 << endl << I11 << endl << I2 << endl;
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
		if (target_nucl)
			target_nucl = false;
		else
			target_nucl = true;
		I1 = 0;
		I11 = 0;
		I2 = 0;
		Q = 0;
		t0 = cpuTime();
		ideal_matching(S, g_ideal, N, Vp, c, target_nucl, disp, dispy, prec, y, pos, w, p, lb, ub, opt, no_of_nodes, s_ideal, t0, d1, v, init_alloc_ideal, s_ideal_d1, ideal_time, ideal_d1_time);
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

		LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, game_time, target_nucl, t0, target, prec, target_time, init_alloc_LMC, I1, I11, I2, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, matching_time, actual_alloc_LMC, node_arrives);
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
		d1 = true;
		LMC(Q, periods, dispy, N, no_of_active_nodes, g, v, s, c, disp, S, leaving, active_nodes, Vp, t1, tmp, target_nucl, tmp, target, prec, tmp, init_alloc_d1C, I, I, I, edge_card_weight, y, pos, w, p, lb, ub, opt, d1, credit, tmp, actual_alloc_d1C, node_arrives);
		t1 = cpuTime();
		d1C_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		d1 = false;
		cout << "d1+c done... ";

		t0 = cpuTime();
		arbitrary_matching(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, target_nucl, target, v, S, prec, t0, actual_alloc_NC, actual_alloc_rand, init_alloc_rand);
		t1 = cpuTime();
		arbitrary_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		cout << "arbitrary matching done... ";

		t0 = cpuTime();
		lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_nucl, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_LM, actual_alloc_LM);
		t1 = cpuTime();
		LM_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		cout << "lexmin done...";

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
		d1 = true;
		lex_min_without_credits(node_arrives, g, leaving, active_nodes, c, disp, no_of_active_nodes, N, Vp, periods, dispy, s, Q, v, S, target_nucl, target, credit, prec, y, pos, w, p, opt, lb, ub, edge_card_weight, t0, d1, I, init_alloc_d1, actual_alloc_d1);
		t1 = cpuTime();
		d1_time = t1 - t0 + init_time + read_time + rand_time + graph_time;
		d1 = false;
		cout << "d1 done!" << endl;


		res.open("results2.txt", ofstream::out | ofstream::trunc);
		res << fixed << setprecision(0) << seed << endl << endl;
		res << fixed << setprecision(17) << read_time << endl << graph_time << endl << rand_time << endl << init_time << endl << game_time << endl << target_time << endl << matching_time << endl << endl;
		res << fixed << setprecision(17) << LMC_time << endl << d1C_time << endl << LM_time << endl << d1_time << endl << arbitrary_time << endl << ideal_time << endl << ideal_d1_time << endl << endl;
		res << fixed << setprecision(0) << I1 << endl << I11 << endl << I2 << endl;
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

void lex_min_without_credits(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, unsigned short int& Q, vector<unsigned short int>& v, unsigned int& S, bool& target_nucl, vector<double>& target, vector<double>& credit, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, double& opt, vector<unsigned short int>& lb, vector<unsigned short int>& ub, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, double& t0, bool& d1, unsigned short int &I, vector<vector<double>> &init_alloc, vector<vector<unsigned short int>> &actual_alloc) {
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
		coop_game(g, v, s, c, disp, dispy, S, Vp, N, active_nodes, leaving);
		if (target_nucl)
			nucl(disp, N, S, target, v, prec);
		else
			shapley(target, v, N, S);
		init_alloc[Q] = target;
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
		lex_min_matching(N, v[S], S, target, s, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
		actual_alloc[Q] = s;
		Q++;
		changing_nodes(active_nodes, leaving, no_of_active_nodes, N, Vp, node_arrives, Q, disp, c, s);
		if (dispy) {
			cout << endl << "Press enter to continue." << endl;
			cin.get();
		}
	}
	return;
}

void arbitrary_matching(vector<unsigned short int>& node_arrives, ListGraph& g, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, vector<ListGraph::Node>& c, bool& disp, vector<unsigned short int>& no_of_active_nodes, unsigned short int& N, unsigned short int& Vp, unsigned short int& periods, bool& dispy, vector<unsigned short int>& s, bool& target_nucl, vector<double>& target, vector<unsigned short int>& v, unsigned int& S, double& prec, double& t0, vector<vector<unsigned short int>> &actual_alloc_NC, vector<vector<unsigned short int>>& actual_alloc, vector<vector<double>>& initial_alloc) {
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
		coop_game(g, v, s, c, disp, disp, S, Vp, N, active_nodes, leaving);
		if (target_nucl)
			nucl(disp, N, S, target, v, prec);
		else
			shapley(target, v, N, S);
		initial_alloc[Q] = target;
		for (unsigned short int i = 0; i < N; i++)
			actual_alloc_NC[Q][i] = v[pow(2, i) - 1];
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
		actual_alloc[Q] = s;
		if (dispy) {
			cout << "s[" << Q << "]=";
			for (unsigned int i = 0; i < N; i++)
					cout << s[i] << " ";
			cout << endl;
		}
	}
	return;
}

void LMC(unsigned short int& Q, unsigned short int &periods, bool &dispy, unsigned short int &N, vector<unsigned short int> &no_of_active_nodes, ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool &disp, unsigned int &S, vector<bool>& leaving, ListGraph::NodeMap<bool>& active_nodes, unsigned short int& Vp, double& t1, double& game_time, bool &target_nucl, double &t0, vector<double>&target, double &prec, double &target_time, vector<vector<double>> &init_alloc, unsigned short int &I1, unsigned short int& I11, unsigned short int& I2, ListGraph::EdgeMap<unsigned short int> &edge_card_weight, vector<double> &y, vector<bool> &pos, vector<unsigned short int> &w, unsigned short int &p, vector<unsigned short int> &lb, vector<unsigned short int> &ub, double &opt, bool &d1, vector<double> &credit, double &matching_time, vector<vector<unsigned short int>> &actual_alloc, vector<unsigned short int>& node_arrives) {
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
		coop_game(g, v, s, c, disp, dispy, S, Vp, N, active_nodes, leaving);
		t1 = cpuTime();
		game_time += t1 - t0;
		t0 = cpuTime();
		if (target_nucl)
			nucl(disp, N, S, target, v, prec);
		else
			shapley(target, v, N, S);
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
		t1 = cpuTime();
		target_time += t1 - t0;
		init_alloc[Q] = target;
		t0 = cpuTime();
		lex_min_matching(N, v[S], S, target, s, no_of_active_nodes, I1, I11, I2, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
		t1 = cpuTime();
		matching_time += t1 - t0;
		actual_alloc[Q] = s;
		t0 = cpuTime();
		for (unsigned short int i = 0; i < N; i++)
			credit[i] += target[i] - s[i];
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
	return;
}

void ideal_matching(unsigned int& S, ListGraph& g, unsigned short int& N, unsigned short int& Vp, vector<ListGraph::Node>& c, bool& target_nucl, bool& disp, bool& dispy, double& prec, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, unsigned short int& no_of_nodes, vector<unsigned short int>& s_ideal, double& t0, bool& d1, vector<unsigned short int>& v_ideal, vector<double>& target_ideal, vector<unsigned short int>& s_ideal_d1, double& ideal_time, double& ideal_d1_time) {
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
	if (target_nucl)
		nucl(disp, N, S, target_ideal, v_ideal, prec);
	else
		shapley(target_ideal, v_ideal, N, S);
	ideal_time = cpuTime() - t0;
	ideal_d1_time = ideal_time;
	t0 = cpuTime();
	lex_min_matching(N, v_ideal[S], S, target_ideal, s_ideal, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
	ideal_time += cpuTime() - t0;
	t0 = cpuTime();
	d1 = true;
	lex_min_matching(N, v_ideal[S], S, target_ideal, s_ideal_d1, no_of_active_nodes, I, I, I, g, edge_card_weight, c, prec, dispy, y, pos, w, p, lb, ub, opt, active_nodes, leaving, credit, disp, Vp, d1);
	ideal_d1_time += cpuTime() - t0;
	d1 = false;
	return;
}

void lex_min_matching(unsigned short int& N, unsigned short int& grandcoal, unsigned int& S, vector<double>& target, vector<unsigned short int>& s, vector<unsigned short int>& no_of_active_nodes, unsigned short int& I1, unsigned short int& I11, unsigned short int& I2, ListGraph& g, ListGraph::EdgeMap<unsigned short int>& edge_card_weight, vector<ListGraph::Node>& c, double& prec, bool& dispy, vector<double>& y, vector<bool>& pos, vector<unsigned short int>& w, unsigned short int& p, vector<unsigned short int>& lb, vector<unsigned short int>& ub, double& opt, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving, vector<double>& credit, bool& disp, unsigned short int& Vp, bool& d1) {
	p = 0;
	if (disp)
		cout << "w*: " << grandcoal << endl;
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
			if (disp)
				cout << "COUNTRY " << w[i] << " lb: " << lb[w[i]] << " ; ub: " << ub[w[i]] << endl;
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

void coop_game(ListGraph& g, vector<unsigned short int>& v, vector<unsigned short int>& s, vector<ListGraph::Node>& c, bool& disp, bool& dispy, unsigned int& S, unsigned short int& Vp, unsigned short int& N, ListGraph::NodeMap<bool>& active_nodes, vector<bool>& leaving) {
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

void shapley(vector<double>& shapl, vector<unsigned short int>& v, unsigned short int& n, unsigned int& s) {
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

//unsigned int char2uint(char &p){
//	return (int)p-48;
//}

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

double cpuTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}

void nucl(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& x, vector<unsigned short int>& v, double& prec) {
	double min_satisfaction = 0;
	bool nlsu = false;
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	vector<double> excess(s, 0);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	vector<unsigned short int> zeros(s, 0);
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
	vector<bool> a(n, false);
	zeros_mem(a, n, s, zeros);
	excess_init(excess, unsettled, x, v, s, n, zeros);
	nucl_comp(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, a, t1, singleton_bounds, nlsu, min_satisfaction, zeros);
	return;
}

void nucl_comp(bool& disp, unsigned short int& n, unsigned int& s, vector<double>& excess, double& prec, vector<bool>& unsettled, unsigned short int& iter, unsigned int& piv, unsigned int& sr, double& t, vector<double>& x, vector<bool>& a, double& t1, vector<double>& singleton_bounds, bool& nlsu, double& min_satisfaction, vector<unsigned short int>& zeros) {
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
	vec_min_uns(epsi, excess, unsettled, s);
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
	while (rank < n)
		pivot(epsi, s, excess, prec, n, a, Arref, J, unsettled, rank, d, x, disp, Asettled, piv, sr, iter, unsettled_p, singleton_bounds, epsi_old, nlsu, zeros, T, T_coord, T2, T2_coord, t_size, t2_size, Atight, Atight2, U, U2, min_satisfaction);
	//cout << "BNF finished!" << endl;
	if (disp) {
		cout << "The nucleolus solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		//cout << "Time needed: " << t << " seconds" << endl;
		cout << "Iterations needed: " << iter << endl;
		cout << "Pivots needed: " << piv << endl;
		cout << "Subroutine solves needed: " << sr << endl;
	}
	return;
}

void zeros_mem(vector<bool>& a, unsigned short int& n, unsigned int& s, vector<unsigned short int>& zeros) {
	a[0] = true;
	zeros[0] = 1;
	for (unsigned int j = 1; j != s + 1; j++) {
		for (unsigned short int i = 0; i < n; i++) {
			if (!a[i]) {
				zeros[j] = i;
				break;
			}
		}
		a[zeros[j]] = 1;
		for (unsigned short int i = 0; i < zeros[j]; i++)
			a[i] = 0;
	}
	return;
}

void excess_init(vector<double>& exc, vector<bool>& unsettled, vector<double>& x, vector<unsigned short int>& v, unsigned int& s, unsigned short int& n, vector<unsigned short int>& zeros) {
	double ax = x[0];
	vector<double> Ux(n, 0);
	for (unsigned short int i = 0; i < n; i++) {
		Ux[i] = x[i];
		for (unsigned short int j = 0; j < i; j++)
			Ux[i] -= x[j];
	}
	exc[0] = ax - v[0];
	for (unsigned int i = 1; i < s; i++) {
		ax += Ux[zeros[i]];
		if (unsettled[i])
			exc[i] = ax - v[i];
		else
			exc[i] = DBL_MAX;
	}
	return;
}

void vec_min_uns(double& m, vector<double>& x, vector<bool>& unsettled, unsigned int& s) {
	m = DBL_MAX;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i] && x[i] < m)
			m = x[i];
	}
	return;
}

void tight_coal(vector<bool>& T, vector<double>& excess, double& epsi, double& prec, unsigned int& s, vector<unsigned int>& T_coord, vector<bool>& unsettled, unsigned int& t_size) {
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

void tight_coal2(vector<bool>& T2, vector<double>& x, vector<double>& singleton_bounds, double& prec, unsigned short int& n, vector<unsigned int>& T2_coord, vector<bool>& unsettled_p, unsigned short int& t2_size) {
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

void pivot(double& epsi, unsigned int& s, vector<double>& excess, double& prec, unsigned short int& n, vector<bool>& a, vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& unsettled, unsigned short int& rank, vector<double>& d, vector<double>& x, bool& disp, vector<vector<bool>>& Asettled, unsigned int& piv, unsigned int& sr_count, unsigned short int& iter, vector<bool>& unsettled_p, vector<double>& singleton_bounds, double& epsi_old, bool& nlsu, vector<unsigned short int>& zeros, vector<bool>& T, vector<unsigned int>& T_coord, vector<bool>& T2, vector<unsigned int>& T2_coord, unsigned int& t_size, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2, double& min_satisfaction) {
	if (disp)
		cout << "Epsilon: " << epsi << endl;
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
		step(T, T2, unsettled, unsettled_p, s, epsi, excess, d, n, x, singleton_bounds, disp, prec, zeros, T_coord, t_size, T2_coord, t2_size, Atight, Atight2, U, U2);
	}
	else {
		if (disp)
			cout << "Min tight set found! Rank increased to: " << rank << endl;
		if (rank == n)
			return;
		if (!nlsu) {
			a[0] = true;
			for (unsigned short int i = 1; i < n; i++)
				a[i] = false;
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					//de2bi(i, a, n);
					if (!(binrank(Arref, J, a, n, rank))) {
						unsettled[i] = false;
						unsettled[s - 1 - i] = false;
					}
				}
				a[zeros[i + 1]] = 1;
				for (unsigned short int j = 0; j < zeros[i + 1]; j++)
					a[j] = 0;
			}
		}
		for (unsigned short int i = 0; i < n; i++)
			if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
				unsettled_p[i] = false;
		vec_min_uns(epsi, excess, unsettled, s);
		T = vector<bool>(s, false);
		T_coord.clear();
		T2 = vector<bool>(n, false);
		T2_coord.clear();
		t_size = 0;
		tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size);
		t2_size = 0;
		if (epsi > prec || epsi < -prec)
			tight_coal2(T2, x, singleton_bounds, prec, n, T2_coord, unsettled_p, t2_size);
		Atight = vector<vector<bool>>(t_size, vector<bool>(n, false));
		for (unsigned int i = 0; i < t_size; i++)
			de2bi(T_coord[i], Atight[i], n);
		Atight2 = vector<vector<bool>>(t2_size, vector<bool>(n, false));
		for (unsigned int i = 0; i < t2_size; i++)
			de2bi(T2_coord[i], Atight2[i], n);
		U = vector<bool>(t_size, true);
		U2 = vector<bool>(t2_size, true);
	}
	return;
}

void subroutine(vector<bool>& U, vector<bool>& U2, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<vector<double>>& Arref, vector<bool>& J, double& prec, unsigned short int& n, unsigned int& tight_size, unsigned short int& tight2_size, unsigned short int& rank, bool& disp, vector<vector<bool>>& Asettled, unsigned int& sr_count, bool& u, unsigned int& s, vector<unsigned int>& T_coord, vector<unsigned int>& T2_coord, vector<bool>& unsettled, double& epsi_old, double& epsi, vector<bool>& unsettled_p, bool& settled, bool& nlsu) {
	unsigned int sumt = 0;
	vector<bool> t(tight_size, false);
	unsigned int sumt2 = 0;
	vector<bool> t2(tight2_size, false);
	glp_prob* lp;
	lp = glp_create_prob();
	//glp_set_prob_name(lp, "sr");
	//glp_set_obj_name(lp, "obj");
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
	vector<int> ia((n + 1) * (tight_size + tight2_size + rank) + 1, 0);
	vector<int> ja((n + 1) * (tight_size + tight2_size + rank) + 1, 0);
	vector<double> ar((n + 1) * (tight_size + tight2_size + rank) + 1, 0);
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
					if (!(binrank(Arref, J, Atight[i], n, rank))) {
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
					if (!(binrank(Arref, J, Atight2[i], n, rank))) {
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

void subr_upd(vector<vector<double>>& Arref, vector<bool>& J, unsigned int& i, unsigned short int& n, double& prec, vector<bool>& U, vector<bool>& U2, unsigned int& sumt, unsigned int& sumt2, vector<bool>& t, vector<bool>& t2, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, unsigned int& tight_size, unsigned short int& tight2_size, unsigned short int& rank, vector<bool>& unsettled, vector<vector<bool>>& Asettled, bool& disp, unsigned int& s, vector<unsigned int>& T_coord, vector<unsigned int>& T2_coord, double& epsi_old, double& epsi, vector<bool>& unsettled_p, bool& settled, glp_prob*& lp, vector<bool>& ar0pos) {
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
			//ar[count - tight_size - tight2_size - rank + i + 1] = 0;
			glp_set_obj_coef(lp, i + 1, 0);
			sumt++;
			unsettled[T_coord[i]] = false;
			unsettled[s - 1 - T_coord[i]] = false;
			if (binrank2(Arref, J, Atight[i], n, rank)) {
				rank++;
				if (epsi > epsi_old) {
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T_coord[i] << " settled as well)" << endl;
				//rowechform(Arref, J, Atight[i], n, rank);
				Asettled[rank - 1] = Atight[i];
				if (disp) {
					cout << "Arref:" << endl;
					for (unsigned short int j = 0; j < n; j++) {
						for (unsigned short int k = 0; k < n; k++) {
							cout << Arref[j][k] << " ";
						}
						cout << endl;
					}
					cout << "J: ";
					for (unsigned short int j = 0; j < n; j++)
						cout << J[j] << " ";
					cout << endl;
				}
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
			if (binrank2(Arref, J, Atight2[i], n, rank)) {
				rank++;
				if (epsi > epsi_old) {
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
					cout << "lambda_" << T2_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T2_coord[i] << " settled as well)" << endl;
				//rowechform(Arref, J, Atight2[i], n, rank);
				Asettled[rank - 1] = Atight2[i];
				if (disp) {
					cout << "Arref:" << endl;
					for (unsigned short int j = 0; j < n; j++) {
						for (unsigned short int k = 0; k < n; k++) {
							cout << Arref[j][k] << " ";
						}
						cout << endl;
					}
					cout << "J: ";
					for (unsigned short int j = 0; j < n; j++)
						cout << J[j] << " ";
					cout << endl;
				}
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

void step(vector<bool>& T, vector<bool>& T2, vector<bool>& unsettled, vector<bool>& unsettled_p, unsigned int& s, double& epsi, vector<double>& excess, vector<double>& d, unsigned short int& n, vector<double>& x, vector<double>& singleton_bounds, bool& disp, double& prec, vector<unsigned short int>& zeros, vector<unsigned int>& T_coord, unsigned int& t_size, vector<unsigned int>& T2_coord, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2) {
	double ad = d[0];
	vector<double> Ud(n, 0);
	double alpha = DBL_MAX;
	double alpha1 = DBL_MAX;
	double alpha2 = DBL_MAX;
	vector<unsigned int> argmin(0, 0);
	vector<unsigned int> argmin2(0, 0);
	for (unsigned short int i = 0; i < n; i++) {
		Ud[i] = d[i];
		for (unsigned short int j = 0; j < i; j++)
			Ud[i] -= d[j];
	}
	if ((ad < 1 - prec) && unsettled[0] && !T2[0] && !T[0]) {
		alpha1 = (epsi - excess[0]) / (ad - 1);
		argmin.push_back(0);
	}
	for (unsigned int j = 1; j < s; j++) {
		ad += Ud[zeros[j]];
		if ((ad < 1 - prec) && unsettled[j] && !T[j])
			if (alpha1 - prec > (epsi - excess[j]) / (ad - 1)) {
				alpha1 = (epsi - excess[j]) / (ad - 1);
				argmin = vector<unsigned int>(1, 0);
				argmin[0] = j;
			}
			else {
				if (abs(alpha1 - (epsi - excess[j]) / (ad - 1)) < prec)
					argmin.push_back(j);
			}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (!T2[i] && unsettled_p[i] && d[i] < -prec) {
			if ((singleton_bounds[i] - x[i]) / d[i] < alpha2 - prec) {
				alpha2 = (singleton_bounds[i] - x[i]) / d[i];
				argmin2 = vector<unsigned int>(1, 0);
				argmin2[0] = i;
			}
			else {
				if (abs((singleton_bounds[i] - x[i]) / d[i] - alpha2) < prec)
					argmin2.push_back(i);
			}
		}
	}
	if (alpha1 < alpha2 - prec)
		alpha = alpha1;
	else {
		if (alpha2 < alpha1 - prec)
			alpha = alpha2;
		else
			alpha = alpha1;
	}
	if (disp)
		cout << "Step size: " << alpha << endl;
	if (disp)
		cout << endl << "  --==  step size obtained  ==--  " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha * d[i];
	if (disp) {
		cout << "New point: " << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	ad = d[0];
	if (unsettled[0]) {
		excess[0] += alpha * ad;
		if (ad > 1 + prec && T[0]) {
			for (unsigned int i = 0; i < t_size; i++) {
				if (T_coord[i] == 0) {
					T[0] = false;
					T_coord.erase(T_coord.begin() + i);
					Atight.erase(Atight.begin() + i);
					U.erase(U.begin() + i);
					t_size--;
					break;
				}
			}
		}
	}
	else {
		if (T[0]) {
			for (unsigned int i = 0; i < t_size; i++) {
				if (T_coord[i] == 0) {
					T[0] = false;
					T_coord.erase(T_coord.begin() + i);
					Atight.erase(Atight.begin() + i);
					U.erase(U.begin() + i);
					t_size--;
					break;
				}
			}
		}
	}
	for (unsigned int j = 1; j < s; j++) {
		ad += Ud[zeros[j]];
		if (unsettled[j]) {
			excess[j] += alpha * ad;
			if (ad > 1 + prec && T[j]) {
				for (unsigned int i = 0; i < t_size; i++) {
					if (T_coord[i] == j) {
						T[j] = false;
						T_coord.erase(T_coord.begin() + i);
						Atight.erase(Atight.begin() + i);
						U.erase(U.begin() + i);
						t_size--;
						break;
					}
				}
			}
		}
		else {
			if (T[j]) {
				for (unsigned int i = 0; i < t_size; i++) {
					if (T_coord[i] == j) {
						T[j] = false;
						T_coord.erase(T_coord.begin() + i);
						Atight.erase(Atight.begin() + i);
						U.erase(U.begin() + i);
						t_size--;
						break;
					}
				}
			}
		}
	}
	for (unsigned short int j = 0; j < n; j++) {
		if (unsettled_p[j] && T2[j] && d[j] > prec) {
			for (unsigned short int i = 0; i < t2_size; i++) {
				if (T2_coord[i] == pow(2, j) - 1) {
					T2[j] = false;
					T2_coord.erase(T2_coord.begin() + i);
					Atight2.erase(Atight2.begin() + i);
					U2.erase(U2.begin() + i);
					t2_size--;
					break;
				}
			}
		}
		else {
			if (!unsettled_p[j] && T2[j]) {
				for (unsigned short int i = 0; i < t2_size; i++) {
					if (T2_coord[i] == pow(2, j) - 1) {
						T2[j] = false;
						T2_coord.erase(T2_coord.begin() + i);
						Atight2.erase(Atight2.begin() + i);
						U2.erase(U2.begin() + i);
						t2_size--;
						break;
					}
				}
			}
		}
	}
	for (unsigned int i = 0; i < t_size; i++)
		U[i] = true;
	for (unsigned short int i = 0; i < t2_size; i++)
		U2[i] = true;
	epsi += alpha;
	if (alpha1 < alpha2 - prec) {
		for (unsigned int j = 0; j < argmin.size(); j++) {
			T[argmin[j]] = true;
			T_coord.push_back(argmin[j]);
			Atight.push_back(vector<bool>(n, false));
			de2bi(argmin[j], Atight[t_size], n);
			U.push_back(true);
			t_size++;
		}
	}
	else {
		if (alpha2 < alpha1 - prec) {
			for (unsigned int j = 0; j < argmin2.size(); j++) {
				T2[argmin2[j]] = true;
				T2_coord.push_back(pow(2, argmin2[j]) - 1);
				Atight2.push_back(vector<bool>(n, false));
				de2bi(T2_coord[t2_size], Atight2[t2_size], n);
				U2.push_back(true);
				t2_size++;
			}
		}
		else {
			for (unsigned int j = 0; j < argmin.size(); j++) {
				T[argmin[j]] = true;
				T_coord.push_back(argmin[j]);
				Atight.push_back(vector<bool>(n, false));
				de2bi(argmin[j], Atight[t_size], n);
				U.push_back(true);
				t_size++;
			}
			for (unsigned int j = 0; j < argmin2.size(); j++) {
				T2[argmin2[j]] = true;
				T2_coord.push_back(pow(2, argmin2[j]) - 1);
				Atight2.push_back(vector<bool>(n, false));
				de2bi(T2_coord[t2_size], Atight2[t2_size], n);
				U2.push_back(true);
				t2_size++;
			}
		}
	}
	return;
}

void imprdir(vector<double>& d, unsigned short int& n, unsigned int& t_size, unsigned short int& t2_size, vector<vector<bool>>& Atight, vector<vector<bool>>& Atight2, vector<bool>& U, vector<bool>& U2, unsigned short int& rank, vector<vector<bool>>& Asettled, bool& disp) {
	glp_prob* dir_lp;
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
	vector<int> ia(n * (t_size + t2_size + rank) + 1, 0);
	vector<int> ja(n * (t_size + t2_size + rank) + 1, 0);
	vector<double> ar(n * (t_size + t2_size + rank) + 1, 0);
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

bool binrank(vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& b, unsigned short int& n, unsigned short int& rank) {
	double prec = pow(10, -10);
	vector<double> B(n, 0);
	for (unsigned short int i = 0; i < n; i++) {
		if (b[i] == true)
			B[i] = 1;
	}
	// m = rank
	// pivot_col[i] = !J[i]
	if (rank >= n)
		return false;
	else {
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
				if (J[k] == true) {
					return true;
				}
				else {
					while (count < k + 1) {
						if (!J[count])
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

bool binrank2(vector<vector<double>>& Arref, vector<bool>& J, vector<bool>& b, unsigned short int& n, unsigned short int& rank) {
	double prec = pow(10, -10);
	vector<double> B(n, 0);
	for (unsigned short int i = 0; i < n; i++) {
		if (b[i] == true)
			B[i] = 1;
	}
	// m = rank
	// pivot_col[i] = !J[i]
	if (rank >= n)
		return false;
	else {
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
				if (J[k] == true) {
					for (unsigned short i = 0; i < k; i++)
						if (!J[i])
							I++;
					for (unsigned short i = 0; i < rank - I; i++)
						Arref[rank - i] = Arref[rank - i - 1];
					Arref[I] = B;
					J[k] = false;
					return true;
				}
				else {
					while (count < k + 1) {
						if (!J[count])
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

void sum_vecb(unsigned int& s, vector<bool>& x) {
	// sums up the values of boolean x
	s = 0;
	for (unsigned int i = 0; i < x.size(); i++)
		s += x[i];
	return;
}

void vec_subtract(vector<double>& z, vector<double>& x, vector<double>& y) {
	// subtracts vector (double) y from vector (double) x
	for (unsigned int i = 0; i != x.size(); i++)
		z[i] = x[i] - y[i];
	return;
}

void sc_vec_prod(vector<double>& y, double a, vector<double>& x) {
	for (unsigned int i = 0; i < x.size(); i++)
		y[i] = a * x[i];
	return;
}
