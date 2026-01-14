#include "BasicTSPheuristics.hpp"
#include "get_tsplib.hpp"

void test_tsplib(string tsp_name) {
	vector<CityPair> v_dat;
	vector<vector<double>> d_dist, d_xy;
    vector<int> i_tour;
    
    int n = get_tsplib(tsp_name, v_dat, d_dist, d_xy);
    i_tour.resize(n);
    double d_z0 = gd_random_tour(n, d_dist, i_tour);

    cout << d_z0 << endl;
}

int main() {
    //g_test_heuristics();
    test_tsplib("berlin52");
    return 0;
}