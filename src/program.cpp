#include "program.h"
#include "quad_bcs.h"
#include "global_variables.h"
#include "preprocessor.h"
#include <ctime>
#include "quad_bcs_plus.h"
#include "domain_geometry.h"
#include "Mesh.h"
#include "Boundary_Conditions.h"
#include "external_forces.h"
#include "Solution.h"
#include "Solver.h"
#include <fstream>
#include <iostream>
#include "postprocessor.h"
#include <stdlib.h>
#include <stdio.h>
#include "TECIO.h"
#include <sstream>
#include <tecplot_output.h>
#include <boost/filesystem.hpp>
#include "post_processing.h"

namespace fs = boost::filesystem;

program::program()
{
    //ctor
}

program::~program()
{
    //dtor
}


void program::run(char* xml_input){


    quad_bcs_plus bcs;
    global_variables globals;
    domain_geometry domain;
    initial_conditions initial_conds;
    int mg =0; // first multigrid cycle
    int fmg =0; // first full multigrid cycle
    std::clock_t start = clock();
    double duration;


    preprocessor pre_processor;

    //postprocessor post_processor;

    pre_processor.initialise_program_variables(xml_input, globals, domain,initial_conds,bcs);

    copyfile(xml_input,globals.output_file);
     // create Mesh
    Mesh mesh(domain,globals);

    remove_existing_files(globals);

    // create boundary conditions
    Boundary_Conditions bc(mesh.get_num_x(), mesh.get_num_y());
    bc.assign_boundary_conditions(mesh.get_num_x(), mesh.get_num_y(),bcs, globals.testcase);

    // assign external force terms
    external_forces source_term(mesh.get_total_cells());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force

    //create solution
    Solution soln(mesh.get_total_cells());
    Solution residual(mesh.get_total_cells());
    post_processing post_processor(mesh.get_total_cells());

    if ( globals.fmg_levels > 0 ){

        fmg_cycle(fmg,soln,soln,mesh,bcs,initial_conds,globals,domain, bc);



    }else{

        soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                initial_conds.rho_origin_mag,mesh,globals);
        soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                initial_conds.vel_origin_mag,mesh,globals);
        soln.set_average_rho(initial_conds.average_rho);


    }
    // Solvec

    Solver solve;
    tecplot_output grid(globals,mesh,soln,bc,1,0.0,post_processor);

    solve.Uniform_Mesh_Solver_Clean_MK2(mesh,soln,bc,source_term,globals,domain,initial_conds,bcs,
                              mg,residual,fmg,post_processor);



    soln.post_process(globals.pre_conditioned_gamma,mesh, globals,initial_conds);
    soln.output(globals.output_file, globals, domain);
    soln.output_centrelines(globals.output_file,globals,mesh);
    //post_processor.output_vtk_mesh(globals.output_file,globals,domain);

    std::clock_t end = clock();

    duration = double(end-start)/ CLOCKS_PER_SEC;

    output_globals(globals,duration);

    std::cout << duration << std::endl;
    std::cout << "CPU Cycles:" << double(end-start) << std::endl;


}

void program::remove_existing_files(global_variables &globals){

    std::string output_location;

    output_location = globals.output_file;

    fs::path p(output_location);

    if(fs::exists(p) && fs::is_directory(p)){
        fs::directory_iterator end;

            for (fs::directory_iterator it(p); it !=end; ++it){
                try{
                    if( fs::is_regular_file(it->status()) && (it->path().extension().compare(".plt") == 0))
                    {
                        fs:remove(it->path());
                    }
                }catch(const std::exception &ex){
                    ex;
                }
            }
        }

    }

void program::output_globals (global_variables globals,double duration){

    std::string output_location;
    std::string filename;
    std::ofstream globals_txt ;
    std::string globals_file;
    output_location = globals.output_file;
    filename = globals.simulation_name;
    globals_file = output_location + "/globals.txt";
    double cycles;

    cycles = duration * CLOCKS_PER_SEC;
/// Generic Load Case Input

    globals_txt.open(globals_file.c_str(), ios::out);

    globals_txt << "Runtime:" << duration << "s" << endl;
    globals_txt << "CPU Cycles:" << cycles << endl;
    globals_txt << "File:" << filename.c_str()  << endl;

    globals_txt << "Tau:"  << globals.tau << endl;
    globals_txt << "Mach" << globals.max_velocity *sqrt(3)<< endl;
    globals_txt << "Reynolds" << globals.reynolds_number << endl;

    std::string reynolds_text;
    reynolds_text = globals.reynolds_number;

    globals_txt.close();


}


    // copy in binary mode
void program::copyfile( char* SRC,  std::string  DEST)
{
    DEST.append("/input.xml");

    std::ifstream src(SRC, std::ios::binary);
    std::ofstream dest(DEST, std::ios::binary);
    dest << src.rdbuf();

}

void program::fmg_cycle(int &fmg,Solution &residual , Solution &soln,
                                      Mesh &fine_mesh, quad_bcs_plus &bcs,
                                      initial_conditions &initial_conds,
                                      global_variables globals, domain_geometry &fine_domain,
                                        Boundary_Conditions &fine_bcs){

    fmg = fmg +1;

    int mg = 0;


    // create new coarse Mesh with double up dimensions
    domain_geometry coarse_domain = fine_mesh.create_coarse_mesh_domain();

    Mesh coarse_mesh (coarse_domain,globals);

    globals.update_tau(coarse_domain);
    globals.magnify_time_step();
    Boundary_Conditions bc(coarse_mesh.get_num_x(), coarse_mesh.get_num_y());
    bc.assign_boundary_conditions(coarse_mesh.get_num_x(), coarse_mesh.get_num_y(),bcs,globals.testcase );

    Solution coarse_soln( coarse_mesh.get_total_cells());
    Solution coarse_residual( coarse_mesh.get_total_cells());
    // goto coarsest level
    if( fmg < globals.fmg_levels){
            fmg_cycle(fmg,coarse_residual,coarse_soln,coarse_mesh,bcs,initial_conds,globals,
                      coarse_domain,bc);

    }else{

        //apply initial conditions at coarsest level
        coarse_soln.assign_pressure_gradient(initial_conds.rho_gradient,initial_conds.origin_loc,
                                initial_conds.rho_origin_mag,coarse_mesh,globals);
        coarse_soln.assign_velocity_gradient(initial_conds.vel_gradient,initial_conds.origin_loc,
                                initial_conds.vel_origin_mag,coarse_mesh,globals);
        coarse_soln.set_average_rho(initial_conds.average_rho);

    }

    external_forces source_term(coarse_mesh.get_total_cells());
    source_term.set_uniform_force(initial_conds.pressure_gradient); //pressure force applied through external force


    Solver solve_coarse;

    solve_coarse.Uniform_Mesh_Solver_Clean(coarse_mesh,coarse_soln,bc,source_term,globals,
                                     coarse_domain,initial_conds,bcs,mg,coarse_residual,fmg);

    Solution temp_soln(coarse_mesh.get_total_cells());
    // prolongation occurs here

    coarse_soln.update_bcs(bc,coarse_mesh,coarse_domain);

    soln.prolongation( coarse_soln, temp_soln, soln,coarse_mesh, fine_mesh,fine_bcs,true);
    soln.set_average_rho(initial_conds.average_rho);
    soln.update_bcs(fine_bcs,fine_mesh,fine_domain);
    globals.update_tau(fine_domain);
     globals.reduce_time_step();
    fmg = fmg -1;

}

