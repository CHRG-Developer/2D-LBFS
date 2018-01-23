#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H
#include "global_variables.h"
#include "domain_geometry.h"



class postprocessor
{
    public:
        postprocessor();
        virtual ~postprocessor();
        void output_vtk_mesh(std::string output_location, global_variables &globals,
        domain_geometry &geometry);

    protected:

    private:
};

#endif // POSTPROCESSOR_H
