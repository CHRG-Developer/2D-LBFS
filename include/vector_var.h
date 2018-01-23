#ifndef VECTOR_VAR_H
#define VECTOR_VAR_H


class vector_var
{
    public:
        vector_var();
        vector_var(double _x, double _y, double _z);
        virtual ~vector_var(){};
        double x;
        double y;
        double z;
        double Dot_Product(vector_var b);
        double Magnitude();
        double Angle_Between_Vectors( vector_var b);
        void Get_Gradient(double y1, double y2, vector_var x1, vector_var x2 );
        void add(vector_var b);
        void factor(double a);
        void round(int precision);
        vector_var line_magnitude(vector_var intercept, vector_var slope, vector_var displacement);
        void relative_vectors(vector_var origin, vector_var ref1, vector_var ref2, double const2);
    protected:
    private:
};

#endif // VECTOR_VAR_H
