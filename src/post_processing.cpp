#include <stdlib.h>
#include <math.h>
#include "post_processing.h"
#include <iostream>

post_processing::post_processing()
{
    //ctor
}

post_processing::~post_processing()
{
    //dtor

     delete [] vorticity;
    vorticity = NULL;
    delete [] streamfunction;
    streamfunction = NULL;

}
post_processing::post_processing(int _total_nodes)
{
    //ctor
    total_nodes = _total_nodes;
     //rho = (double*) malloc (sizeof(double)*(total_nodes));
     vorticity = new double [total_nodes +1];
        if (vorticity==NULL) exit (1);
     streamfunction = new double [total_nodes+1];
        if (streamfunction==NULL) exit (1);
        Initialise();

}

void post_processing::Initialise() {

    std::fill_n(vorticity, total_nodes , 0.00);
    std::fill_n(streamfunction, total_nodes, 0.0);



}

void post_processing::calc_vorticity(Solution &x_grad, Solution & y_grad){

    for(int i =0; i < total_nodes; i++){
        vorticity[i] = x_grad.get_v(i) - y_grad.get_u(i);

    }

}

void post_processing::calc_streamfunction(Mesh &mesh, global_variables &globals,
            Boundary_Conditions &bcs){


    double t,w;

    t= cos(globals.PI/mesh.get_num_x()) + cos(globals.PI/mesh.get_num_y());

    w = (8 - sqrt(64-16*pow(t,2)) )/ pow(t,2);

    double residue,r_min;
    bool loop;
    loop = true;

    while (loop == true){
        r_min = 0;

        for(int i =0; i< total_nodes;i++){
            if( ! bcs.get_bc(i)){
                residue = w*( 0.25* ( streamfunction[mesh.get_w_node(i)] + streamfunction[mesh.get_n_node(i)]+
                            streamfunction[mesh.get_s_node(i)] + streamfunction[mesh.get_e_node(i)]
                            + vorticity[i]) -streamfunction[i]);

                    r_min = r_min + fabs(residue);
                    streamfunction[i] = streamfunction[i] + residue;
            }else{

                streamfunction[i] = 0;
            }
        }
        r_min = r_min/total_nodes;
        std::cout << "R-min:" << r_min << std::endl;

     if (r_min < 10e-9){
        loop = false;
     }


    }


}
