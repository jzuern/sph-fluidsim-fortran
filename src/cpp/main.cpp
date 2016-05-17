// include basic stuff
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <time.h>
#include <algorithm>

using namespace std;

// include headers of all outsourced functions
#include "sph_funcs.h"
#include "integrate.h"
#include "linked_list.h"

#include <boost/tuple/tuple.hpp>

// This must be defined before the first time that "gnuplot-iostream.h" is included.
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"



int main( int argc, const char* argv[]){

        // parsing input parameters
        cout << argc << endl;
        const char* paramFile;
        if (argc != 2) { // Check the value of argc. If not enough parameters have been passed, inform user and exit.
                cout << "You need to provide a parameter file.\n"; // Inform the user of how to use the program
                exit(0);
        } else {
                cout << "Reading parameter file\n";
                paramFile = argv[1]; // name of file containing simulation parameters
        }

        // read from param file
        ifstream inputFile(paramFile);
        string line;
        vector<string> paramvector;

        while (getline(inputFile, line)) {
                if (line.at(0) != '#') { // read only lines without preceding '#'
                        // trim at whitespaces
                        size_t found_whitespace = line.find_first_of(" ");
                        line.erase(found_whitespace);
                        cout << line << endl;
                        paramvector.push_back(line);
                }
        }

        sim_param_t params = initialize_sim_param_t(paramvector);
        sim_state_t* state = init_particles(&params);


        int **lc =  init_lc(params, state);// linked cell
        int *ll  =  init_ll(params, state);//  linked list

        setup_neighbour_list(&params, state, ll, lc);

        int nFrames = params.nFrames;
        int nStepsFrame = params.nStepsFrame;
        float dt = params.dt;
        int n = state->n;

        ofstream out;
        out.open(params.filename,ios::ate); //write new data always at the end (ate)

        compute_accel(state, &params,ll,lc);
        leapfrog_start(state,dt);
        check_state(state);

        // print_neighour_list(state,params,ll,lc);

        bool initialize = true;


        for (int frame = 1; frame < nFrames; ++frame) { // main simulation loop
                cout << "Calculating Frame " << frame << " of " << nFrames << "\n";
                for (int i = 0; i < nStepsFrame; ++i) { // frame loop
                        compute_accel(state, &params,ll,lc); // update values for accellerations
                        leapfrog_step(state, dt); // update velocities and positions based on previously calculated accelleration
                        check_state(state);
                }
                plotPoints(frame,out, n, state,params);
                usleep(1000);
        }

        free_state(state);

        // Plotting simulation results
        // for (int frame = 1; frame < nFrames; ++frame) {
        //  plot_points(frame);
        //  usleep(25000);
        // }

        // TEST INTERACTIVE VISUALIZATION

        // Gnuplot gp;
        //
        // // Create field of arrows at random locations.
        // std::vector<boost::tuple<double,double,double,double> > arrows;
        // for(size_t i=0; i<100; i++) {
        //         double x = rand() / double(RAND_MAX);
        //         double y = rand() / double(RAND_MAX);
        //         arrows.push_back(boost::make_tuple(x, y, 0, 0));
        // }
        //
        // double mx=0.5, my=0.5;
        // int mb=1;
        // while(mb != 3 && mb >= 0) {
        //         // Make the arrows point towards the mouse click.
        //         for(size_t i=0; i<arrows.size(); i++) {
        //                 double x = arrows[i].get<0>();
        //                 double y = arrows[i].get<1>();
        //                 double dx = (mx-x) * 0.1;
        //                 double dy = (my-y) * 0.1;
        //                 arrows[i] = boost::make_tuple(x, y, dx, dy);
        //         }
        //
        //         gp << "plot '-' with vectors notitle\n";
        //         gp.send1d(arrows);
        //
        //         gp.getMouse(mx, my, mb, "Left click to aim arrows, right click to exit.");
        //         printf("You pressed mouse button %d at x=%f y=%f\n", mb, mx, my);
        //         if(mb < 0) printf("The gnuplot window was closed.\n");
        // }



        char c;
        cin >> c;

        return 0;
}
