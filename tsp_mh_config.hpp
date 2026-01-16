#include "tsp_ts.hpp"
#include <sstream>
#include <map>

//////////////////////////////////////////////////////////////////////////////
// Helper functions for metaheuristic configuration
//////////////////////////////////////////////////////////////////////////////

// Reads a CSV config file into a map (param -> value)
// CSV format: param,value,comment (third column ignored)
//
map<string, string> read_config(const string& config_fname) {
    map<string, string> config;
    ifstream file(config_fname);
    if (!file.is_open()) {
        throw runtime_error("Could not open config file: " + config_fname);
    }

    string line;
    getline(file, line); // skip header
    while (getline(file, line)) {
        istringstream iss(line);
        string param, value, comment;
        if (getline(iss, param, ',') && getline(iss, value, ',')) {
            getline(iss, comment); // ignored
            config[param] = value;
        }
    }
    file.close();
    return config;
}

// Sets random seed from config: 0 = time-based, >0 = fixed seed
//
void set_random_seed(const map<string, string>& config) {
    int i_seed = config.count("seed") ? stoi(config.at("seed")) : 0;
    if (i_seed > 0)
        srand(i_seed);
    else
        srand(static_cast<unsigned int>(time(NULL)));
}

// Creates and configures neighbourhood operator from config
// Options: 2opt, swap, ins, multi (with p_2opt, p_swap, p_ins probabilities)
//
TSPnbrOp* create_neighbourhood(const map<string, string>& config, int ai_n, vector<vector<double>>& ad_dist) {
    string s_nhd = config.count("neighbourhood") ? config.at("neighbourhood") : "2opt";

    if (s_nhd == "2opt")
        return new Op2opt(ai_n, ad_dist);
    else if (s_nhd == "swap")
        return new OpSwap(ai_n, ad_dist);
    else if (s_nhd == "ins")
        return new OpIns(ai_n, ad_dist);
    else if (s_nhd == "multi") {
        OpMulti* p_multi = new OpMulti(ai_n, ad_dist);

        double d_p2opt = config.count("p_2opt") ? stod(config.at("p_2opt")) : 0.34;
        double d_pswap = config.count("p_swap") ? stod(config.at("p_swap")) : 0.33;
        double d_pins = config.count("p_ins") ? stod(config.at("p_ins")) : 0.33;

        Op2opt* p_2opt = new Op2opt(ai_n, ad_dist);
        OpSwap* p_swap = new OpSwap(ai_n, ad_dist);
        OpIns* p_ins = new OpIns(ai_n, ad_dist);

        p_multi->m_add_nhd(*p_2opt, d_p2opt);
        p_multi->m_add_nhd(*p_swap, d_pswap);
        p_multi->m_add_nhd(*p_ins, d_pins);

        return p_multi;
    }
    else
        return new Op2opt(ai_n, ad_dist); // default
}

// Creates and configures initial solution from config
// Options: random, nn (nearest neighbour), nn_2opt (nn + 2-opt local search)
//
TSPsoln* create_initial_solution(const map<string, string>& config, int ai_n, vector<vector<double>>& ad_dist) {
    TSPsoln* x0 = new TSPsoln(ai_n);
    string s_init = config.count("init_soln") ? config.at("init_soln") : "random";

    if (s_init == "random")
        x0->md_z = gd_random_tour(ai_n, ad_dist, x0->mi_x);
    else if (s_init == "nn")
        x0->md_z = gd_nearest_neighbour(ai_n, ad_dist, 0, x0->mi_x);
    else if (s_init == "nn_2opt") {
        x0->md_z = gd_nearest_neighbour(ai_n, ad_dist, 0, x0->mi_x);
        Op2opt ls_nhd(ai_n, ad_dist);
        TSPsoln x_best(ai_n);
        gx_local_search(false, ls_nhd, *x0, x_best);
        *x0 = x_best;
    }
    else
        x0->md_z = gd_random_tour(ai_n, ad_dist, x0->mi_x); // default

    return x0;
}

//////////////////////////////////////////////////////////////////////////////
// Simulated Annealing initialisation
//////////////////////////////////////////////////////////////////////////////

// Initialises simulated annealing components from a CSV configuration file.
//
// Parameters:
//   nhd          - [out] neighbourhood operator (created and configured by this function)
//   x0           - [out] initial solution (created and configured by this function)
//   ai_n         - number of cities in the TSP instance
//   ad_dist      - distance matrix for the TSP instance
//   config_fname - path to CSV config file (defaults to "sa_config_default.csv" if empty)
//
// Returns:
//   SActrl struct containing SA hyperparameters
//
// CSV format (3 columns: param, value, comment/options - third column is ignored):
//   param,value,options
//   seed,0,0=time-based|any positive int for fixed seed
//   neighbourhood,multi,2opt|swap|ins|multi
//   init_soln,nn_2opt,random|nn|nn_2opt
//   init_accept_rate,0.9,
//   accept_limit_coef,10,
//   nbr_limit_coef,100,
//   tfact,0.95,
//   stale_step_lim,10,
//   log_every_it,0,0|1
//   p_2opt,0.34,(multi only)
//   p_swap,0.33,(multi only)
//   p_ins,0.33,(multi only)
//
// Config options:
//   seed             - 0: use time-based seed, >0: use specified seed for reproducibility
//   neighbourhood    - 2opt, swap, ins, or multi (combines all three)
//   init_soln        - random: random tour
//                      nn: nearest neighbour heuristic
//                      nn_2opt: nearest neighbour + 2-opt local search (first improvement)
//   init_accept_rate - initial acceptance rate for calculating T0
//   accept_limit_coef - temperature step ends after accept_limit_coef * n acceptances
//   nbr_limit_coef    - temperature step ends after nbr_limit_coef * n neighbours evaluated
//   tfact            - temperature reduction factor (T_new = T * tfact)
//   stale_step_lim   - algorithm stops after this many temp steps without acceptance
//   log_every_it     - 1: log every iteration, 0: log end of temp step only
//   p_2opt/p_swap/p_ins - probabilities for multi-neighbourhood (should sum to 1.0)
//
SActrl initialise_sa(TSPnbrOp*& nhd, TSPsoln*& x0, int ai_n, vector<vector<double>>& ad_dist, string config_fname="") {
    string s_fname = config_fname.empty() ? "sa_config_default.csv" : config_fname;
    map<string, string> config = read_config(s_fname);

    set_random_seed(config);
    nhd = create_neighbourhood(config, ai_n, ad_dist);
    x0 = create_initial_solution(config, ai_n, ad_dist);

    // Configure SA-specific hyperparameters
    SActrl hyper_params;
    hyper_params.md_init_accept_rate = config.count("init_accept_rate") ? stod(config["init_accept_rate"]) : 0.9;
    hyper_params.mi_accept_limit_coef = config.count("accept_limit_coef") ? stoi(config["accept_limit_coef"]) : 10;
    hyper_params.mi_nbr_limit_coef = config.count("nbr_limit_coef") ? stoi(config["nbr_limit_coef"]) : 100;
    hyper_params.md_tfact = config.count("tfact") ? stod(config["tfact"]) : 0.95;
    hyper_params.mi_stale_step_lim = config.count("stale_step_lim") ? stoi(config["stale_step_lim"]) : 10;
    hyper_params.mb_every_it = config.count("log_every_it") ? (stoi(config["log_every_it"]) != 0) : false;

    return hyper_params;
}

//////////////////////////////////////////////////////////////////////////////
// Tabu Search initialisation
//////////////////////////////////////////////////////////////////////////////

// Initialises tabu search components from a CSV configuration file.
//
// Parameters:
//   nhd          - [out] neighbourhood operator (created and configured by this function)
//   x0           - [out] initial solution (created and configured by this function)
//   tabu_list    - [out] tabu list (created and configured by this function)
//   long_mem     - [out] long-term memory (created if long_mem_wt > 0, otherwise nullptr)
//   ai_n         - number of cities in the TSP instance
//   ad_dist      - distance matrix for the TSP instance
//   config_fname - path to CSV config file (defaults to "ts_config_default.csv" if empty)
//
// Returns:
//   TSctrl struct containing TS hyperparameters
//
// CSV format (3 columns: param, value, comment/options - third column is ignored):
//   param,value,options
//   seed,0,0=time-based|any positive int for fixed seed
//   neighbourhood,2opt,2opt|swap|ins|multi
//   init_soln,nn,random|nn|nn_2opt
//   itmax,1000,max iterations
//   tabu_tenure,10,iterations a solution remains tabu
//   hash_length,1000,length of hash table
//   long_mem_wt,0,0=disabled|>0 weight for frequency penalty
//   p_2opt,0.34,(multi only)
//   p_swap,0.33,(multi only)
//   p_ins,0.33,(multi only)
//
// Config options:
//   seed           - 0: use time-based seed, >0: use specified seed for reproducibility
//   neighbourhood  - 2opt, swap, ins, or multi (combines all three)
//   init_soln      - random: random tour
//                    nn: nearest neighbour heuristic
//                    nn_2opt: nearest neighbour + 2-opt local search (first improvement)
//   itmax          - maximum number of iterations
//   tabu_tenure    - number of iterations a solution remains tabu
//   hash_length    - length of hash table for z-value based tabu list
//   long_mem_wt    - 0: no long-term memory, >0: frequency penalty weight
//   p_2opt/p_swap/p_ins - probabilities for multi-neighbourhood (should sum to 1.0)
//
TSctrl initialise_ts(TSPnbrOp*& nhd, TSPsoln*& x0, TabuListZhash*& tabu_list, TSPnhdMem*& long_mem, int ai_n, vector<vector<double>>& ad_dist, string config_fname="") {
    string s_fname = config_fname.empty() ? "ts_config_default.csv" : config_fname;
    map<string, string> config = read_config(s_fname);

    set_random_seed(config);
    nhd = create_neighbourhood(config, ai_n, ad_dist);
    x0 = create_initial_solution(config, ai_n, ad_dist);

    // Configure tabu list
    int i_hash_len = config.count("hash_length") ? stoi(config["hash_length"]) : 1000;
    int i_tenure = config.count("tabu_tenure") ? stoi(config["tabu_tenure"]) : 10;
    tabu_list = new TabuListZhash(i_hash_len, i_tenure);

    // Configure long-term memory (nullptr if weight is 0)
    double d_long_mem_wt = config.count("long_mem_wt") ? stod(config["long_mem_wt"]) : 0.0;
    if (d_long_mem_wt > 0.0)
        long_mem = new TSPnhdMem(*nhd, d_long_mem_wt);
    else
        long_mem = nullptr;

    // Configure TS-specific hyperparameters
    TSctrl hyper_params;
    hyper_params.mi_itmax = config.count("itmax") ? stoi(config["itmax"]) : 1000;

    return hyper_params;
}
