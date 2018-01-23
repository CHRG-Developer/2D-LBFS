#ifndef RESIDUALS_H
#define RESIDUALS_H


class residuals
{
    public:
        residuals();
        virtual ~residuals();
        void reset();
        void add_l2_norm_residuals(double f1, double rho, double f2, double u, double f3, double v);
        void l2_norm_rms();
        double max_error();

        void add_ansys_l2_norm_residuals(double f1, double rho, double f2, double u, double f3, double v,
                                          double delta_t);
        void ansys_5_iter_rms (int timestep);
        double rho_rms, u_rms, v_rms, w_rms;
    protected:
    private:
     double rho_rms_num, u_rms_num, v_rms_num, w_rms_num, error;
    double rho_rms_den, u_rms_den, v_rms_den, w_rms_den;

    double max_rho, max_u, max_v;
};

#endif // RESIDUALS_H
