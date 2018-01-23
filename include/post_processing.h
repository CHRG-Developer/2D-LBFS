#include "Solution.h"
#ifndef POST_PROCESSING_H
#define POST_PROCESSING_H


class post_processing
{
    public:
        post_processing();
        virtual ~post_processing();
        post_processing(int _total_nodes);
        void calc_vorticity(Solution &x_grad, Solution &y_grad);
        void calc_streamfunction(Mesh &mesh,global_variables &globals,
            Boundary_Conditions &bcs);
            void Initialise();
        double *vorticity, *streamfunction;
    protected:

    private:
         private:

        int total_nodes;
};

#endif // POST_PROCESSING_H
