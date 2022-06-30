// Project Identifier: 3E33912F8BAA7542FC4A1585D2DB6FE0312725B9

#include "TSP.h"
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <limits>
#include "xcode_redirect.hpp"
#include <getopt.h>
#include <cmath>
#include <algorithm>
using namespace std;

void TSP::get_options(int argc, char **argv)
{
    int option_index = 0;
    int choice = 0;
    opterr = false;
    option long_options[] = {
        {"mode", required_argument, nullptr, 'm'},
        {"help", no_argument, nullptr, 'h'},
        { nullptr, 0, nullptr, '\0' }
    };
    while ((choice = getopt_long(argc, argv, "m:h", long_options, &option_index)) != -1)
    {
        switch (choice)
        {
             case 'm':
                run_mode = optarg;
                break;
            
             case 'h':
                std::cout << "This program takes in coordinates and runs different modes.\n"
                          << "The modes are MST, OPTTSP and FASTTSP, which is run based on\n"
                          << "the mode that is inputted\n";
                exit(0);
        }
    }
    if (run_mode == " ")
    {
        cerr << "no mode specified\n";
        exit(1);
    }
}

void TSP::read_file()
{
    int x_coordinate_read = 0;
    int y_coordinate_read = 0;
    cin >> num_cages;
    input_coordinates.resize(num_cages);
    prism_table.resize(num_cages);
    
    for (int i = 0; i < num_cages; ++i)
    {
        cin >> x_coordinate_read >> y_coordinate_read;
        input_coordinates[i].x_coordinate = x_coordinate_read;
        input_coordinates[i].y_coordinate = y_coordinate_read;
        if (x_coordinate_read < 0 && y_coordinate_read < 0)
        {
            input_coordinates[i].category = 'w'; // w - wild
            ++wild_cages;
        }
        else if ((x_coordinate_read == 0 && y_coordinate_read <= 0) || (x_coordinate_read <= 0 && y_coordinate_read == 0))
        {
            input_coordinates[i].category = 'c'; // c - cage border
            ++num_border_cages;
        }
        else
        {
            input_coordinates[i].category = 's'; // s - safe zone
            ++safe_cages;
        }
    }
}

void TSP::MST_distance_helper_A(uint32_t a, uint32_t b)
{
    if (!(input_coordinates[a].category == 's' && input_coordinates[b].category == 'w') && !(input_coordinates[a].category == 'w' && input_coordinates[b].category == 's'))
    {
        double distance_x = ((double)input_coordinates[b].x_coordinate - (double)input_coordinates[a].x_coordinate) * ((double)input_coordinates[b].x_coordinate - (double)input_coordinates[a].x_coordinate);
        double distance_y = ((double)input_coordinates[b].y_coordinate - (double)input_coordinates[a].y_coordinate) * ((double)input_coordinates[b].y_coordinate - (double)input_coordinates[a].y_coordinate);
        double distance_compute = distance_x + distance_y;
        if (distance_compute <= prism_table[b].distance)
        {
            prism_table[b].distance = distance_compute;
            prism_table[b].parent_predecessor = a;
        }
    }
}

void TSP::MST_mode_implementation_A()
{
    int smallest_index = 0;
    double smallest_distance = 0;
    prism_table[0].distance = 0;
    for (int i = 0; i < num_cages; ++i)
    {
        smallest_index = 0;
        smallest_distance = numeric_limits<double>::infinity();
        for (int j = 0; j < num_cages; ++j)
        {
            if(prism_table[j].visited_status == false)
            {
                if(prism_table[j].distance <= smallest_distance)
                {
                    smallest_index = j;
                    smallest_distance = prism_table[j].distance;
                }
            }
        }
        prism_table[smallest_index].visited_status =  true;
        mst_total_cost = mst_total_cost + sqrt(prism_table[smallest_index].distance);
        for (int j1 = 0; j1 < num_cages; ++j1)
        {
            if (prism_table[j1].visited_status == false)
            {
                MST_distance_helper_A(smallest_index, j1);
            }
        }
    }
}

void TSP::print_and_calculate_MST_A()
{
    cout << mst_total_cost << '\n';
    for(int j = 1; j < num_cages; ++j)
    {
        if (j > prism_table[j].parent_predecessor)
        {
            cout << prism_table[j].parent_predecessor << " " << j << '\n';
        }
        else
        {
            cout << j << " " << prism_table[j].parent_predecessor << '\n';
        }
    }
}

void TSP::MST_mode_A()
{
    if (wild_cages > 0 && safe_cages > 0 && num_border_cages == 0)
    {
        cerr << "no border cages" << '\n';
        exit(1);
    }
    MST_mode_implementation_A();
    print_and_calculate_MST_A();
}

void TSP::random_insertion_B()
{
    double distance_compute = 0;
    path_indices_optimal.reserve(num_cages + 1);
    path_indices_optimal.emplace_back(0);
    path_indices_optimal.emplace_back(1);
    path_indices_optimal.emplace_back(2);
    path_indices_optimal.emplace_back(0);
    
    double x_01 = input_coordinates[0].x_coordinate - input_coordinates[1].x_coordinate;
    double y_01 = input_coordinates[0].y_coordinate - input_coordinates[1].y_coordinate;
    double xy_01 = sqrt((x_01 * x_01) + (y_01 * y_01));
    
    double x_12 = input_coordinates[1].x_coordinate - input_coordinates[2].x_coordinate;
    double y_12 = input_coordinates[1].y_coordinate - input_coordinates[2].y_coordinate;
    double xy_12 = sqrt((x_12 * x_12) + (y_12 * y_12));
    
    double x_20 = input_coordinates[0].x_coordinate - input_coordinates[2].x_coordinate;
    double y_20 = input_coordinates[0].y_coordinate - input_coordinates[2].y_coordinate;
    double xy_20 = sqrt((x_20 * x_20) + (y_20 * y_20));
    optimal_length = xy_01 + xy_12 + xy_20;
    
    for (int i = 3; i < num_cages; ++i)
    {
        int smallest_index = 0;
        double smallest_distance = numeric_limits<double>::infinity();
        for (int j = 0; j < i; ++j)
        {
            double x_ji = input_coordinates[path_indices_optimal[j]].x_coordinate  - input_coordinates[i].x_coordinate;
            double y_ji = input_coordinates[path_indices_optimal[j]].y_coordinate - input_coordinates[i].y_coordinate;
            double xy_ji = sqrt((x_ji * x_ji) + (y_ji * y_ji));
            double x_j1i = input_coordinates[path_indices_optimal[j + 1]].x_coordinate - input_coordinates[i].x_coordinate;
            double y_j1i = input_coordinates[path_indices_optimal[j + 1]].y_coordinate - input_coordinates[i].y_coordinate;
            double xy_j1i = sqrt((x_j1i * x_j1i) + (y_j1i * y_j1i));
            double x_jj1 = input_coordinates[path_indices_optimal[j]].x_coordinate - input_coordinates[path_indices_optimal[j + 1]].x_coordinate;
            double y_jj1 = input_coordinates[path_indices_optimal[j]].y_coordinate - input_coordinates[path_indices_optimal[j + 1]].y_coordinate;
            double xy_jj1 = sqrt((x_jj1 * x_jj1) + (y_jj1 * y_jj1));
            
            distance_compute = xy_ji + xy_j1i - xy_jj1;
            if (distance_compute <= smallest_distance)
            {
                smallest_distance = distance_compute;
                smallest_index = j + 1;
            }
        }
        optimal_length = optimal_length + smallest_distance;
        path_indices_optimal.insert((path_indices_optimal.begin() + smallest_index), i);
    }
    path_indices_optimal.pop_back();
}

void TSP::print_and_calculate_insertion_B()
{
    cout << optimal_length << '\n';
    for (size_t i = 0; i < path_indices_optimal.size(); ++i)
    {
        cout << path_indices_optimal[i] << " ";
    }
}

void TSP::insertion_B()
{
    random_insertion_B();
    print_and_calculate_insertion_B();
}

double TSP::distance_help_Part_c_sqrt(size_t a, size_t b)
{
   return  sqrt(((input_coordinates[current_path[a]].x_coordinate - input_coordinates[current_path[b]].x_coordinate) * (input_coordinates[current_path[a]].x_coordinate - input_coordinates[current_path[b]].x_coordinate)) + ((input_coordinates[current_path[a]].y_coordinate - input_coordinates[current_path[b]].y_coordinate) * (input_coordinates[current_path[a]].y_coordinate - input_coordinates[current_path[b]].y_coordinate)));
}

double TSP::distance_help_Part_c(size_t a, size_t b)
{
    return  ((input_coordinates[current_path[a]].x_coordinate - input_coordinates[current_path[b]].x_coordinate) * (input_coordinates[current_path[a]].x_coordinate - input_coordinates[current_path[b]].x_coordinate)) + ((input_coordinates[current_path[a]].y_coordinate - input_coordinates[current_path[b]].y_coordinate) * (input_coordinates[current_path[a]].y_coordinate - input_coordinates[current_path[b]].y_coordinate));

}

double TSP::MST_Part_C_implementation(size_t a)
{
   
    double mst_distance = 0;
    vector<prism_table_features> mst_c;
    int smallest_index = 0;
    double smallest_distance = 0;
    int size = (int)(current_path.size() - a);
    mst_c.resize(size);
    mst_c[0].distance = 0;
    for (size_t i = a; i < current_path.size(); ++i)
    {
        smallest_index = 0;
        smallest_distance = numeric_limits<double>::infinity();
        for (int j = 0; j < size; ++j)
        {
            if(mst_c[j].visited_status == false)
            {
                if(mst_c[j].distance <= smallest_distance)
                {
                    smallest_index = j;
                    smallest_distance = mst_c[j].distance;
                }
            }
        }
    
        mst_c[smallest_index].visited_status =  true;
        mst_distance = mst_distance + sqrt(mst_c[smallest_index].distance);
        for (int j1 = 0; j1 < size; ++j1)
       {
            if (mst_c[j1].visited_status == false)
          {
              double distance_compute(distance_help_Part_c(smallest_index + a, j1 + a));
              if (distance_compute <= mst_c[j1].distance)
              {
                  mst_c[j1].distance = distance_compute;
              }
          }
       }
    }
    return  mst_distance;
}


bool TSP::promising(size_t b)
{
   
        if (current_length >= best_path_length)
        {
            return false;
        }
        double mst_min_distance = 0;
        double smallest_distance_first = numeric_limits<double>::infinity();
        double smallest_distance_second = numeric_limits<double>::infinity();
        for (size_t i = b; i < current_path.size(); ++i)
        {
            double smallest_distance_compute_1(distance_help_Part_c(0, i));
            if (smallest_distance_compute_1 <= smallest_distance_first)
            {
                smallest_distance_first = smallest_distance_compute_1;
            }
            double smallest_distance_compute_2(distance_help_Part_c((b-1), i));
            if (smallest_distance_compute_2 <= smallest_distance_second)
            {
                smallest_distance_second = smallest_distance_compute_2;
            }
        }
        mst_min_distance = MST_Part_C_implementation(b);
        double distance_total = sqrt(smallest_distance_first) + sqrt(smallest_distance_second) + mst_min_distance + current_length;
        
        
        if (distance_total < best_path_length)
        {
            return true;
        }
        else
        {
            return false;
        }
}

void TSP::genPerms(size_t permLength)
{
  if (permLength == current_path.size())
  {
      double distance_border_edge((distance_help_Part_c_sqrt(permLength - 1, 0)));
    current_length = current_length + distance_border_edge;
      if (current_length < best_path_length)
      {
          best_path_length = current_length;
          best_path = current_path;
      }
    current_length = current_length - distance_border_edge;
    return;
  }
  if (!promising(permLength))
  {
      return;
  }
  for (size_t i = permLength; i < current_path.size(); ++i)
    {
       swap(current_path[permLength], current_path[i]);
        double distance_new_edge((distance_help_Part_c_sqrt(permLength, permLength - 1)));
        current_length = current_length + distance_new_edge;
        genPerms(permLength + 1);
        current_length = current_length - distance_new_edge;
        swap(current_path[permLength], current_path[i]);
    }
}

void TSP::Mst_germperms_Part_c_implementation()
{
    random_insertion_B();
    best_path_length = optimal_length;
    best_path = path_indices_optimal;
    current_path = path_indices_optimal;
    size_t perm_length_copy = 1;
    genPerms(perm_length_copy);
}

void TSP::print_and_calculate_Part_c()
{
    
    cout << best_path_length << '\n';
    for (size_t i = 0; i < best_path.size(); ++i)
    {
        cout << best_path[i] << " ";
    }
}

void TSP::TSP_salesperson_solution()
{
    Mst_germperms_Part_c_implementation();
    print_and_calculate_Part_c();
}

void TSP::all_modes_implementation()
{
    transform(run_mode.begin(), run_mode.end(), run_mode.begin(), ::toupper);
    if (run_mode == "MST")
    {
        MST_mode_A();
    }
    else if (run_mode == "FASTTSP")
    {
        insertion_B();
    }
    else if (run_mode == "OPTTSP")
    {
        TSP_salesperson_solution();
    }
}
