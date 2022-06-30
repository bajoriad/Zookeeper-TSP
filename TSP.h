// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9


#ifndef TSP_hpp
#define TSP_hpp

#include <stdio.h>
#include <vector>
#include <iomanip>
#include <limits>
#include "xcode_redirect.hpp"
#include <getopt.h>
#include <cmath>
#include <algorithm>
using namespace std;

class TSP
{
   private:
   
    struct coordinate
    {
        int x_coordinate = 0;
        int y_coordinate = 0;
        char category = ' '; // 's'- safe zone, 'w' - wild zone, 'c' - cage zone 
    };
    
    struct prism_table_features
    {
        bool visited_status = false;
        int parent_predecessor = 0;
        double distance = numeric_limits<double>::infinity();
    };
    
    vector<coordinate> input_coordinates;
    vector<prism_table_features> prism_table;
    vector<uint32_t> best_path; // best path
    vector<uint32_t> current_path; // currently working on
    double current_length = 0; // currently working on 
    double best_path_length = 0; // best path length
    double optimal_length = 0;
    double distance_help_Part_c_sqrt(size_t a, size_t b);
    double distance_help_Part_c(size_t a , size_t b);
    vector<uint32_t> path_indices_optimal;
    string run_mode = " ";
    int num_cages = 0;
    int num_border_cages = 0;
    int wild_cages = 0;
    int safe_cages = 0;
    double mst_total_cost = 0;
    
    public:
    
    void get_options(int agrc, char **argv);
    void read_file();
    void MST_distance_helper_A(uint32_t a, uint32_t b);
    void MST_mode_implementation_A();
    void print_and_calculate_MST_A();
    void MST_mode_A();
    void random_insertion_B();
    void print_and_calculate_insertion_B();
    void insertion_B();
    double MST_Part_C_implementation(size_t a);
    void genPerms(size_t a);
    bool promising(size_t b);
    void Mst_germperms_Part_c_implementation();
    void print_and_calculate_Part_c();
    void TSP_salesperson_solution();
    void all_modes_implementation();
};

#endif /* TSP_hpp */
