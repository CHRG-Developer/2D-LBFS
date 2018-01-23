#ifndef TECPLOT_OUTPUT_H
#define TECPLOT_OUTPUT_H
#include <global_variables.h>
#include <Mesh.h>
#include <Solution.h>
#include <Boundary_Conditions.h>
#include <post_processing.h>


class tecplot_output
{
    public:
        tecplot_output(global_variables &globals, Mesh &Mesh, Solution &Soln,
                             Boundary_Conditions &bcs, int fileType_e, double timestamp,
                             post_processing &pp);
        virtual ~tecplot_output();

    protected:

    private:
};

#endif // TECPLOT_OUTPUT_H
