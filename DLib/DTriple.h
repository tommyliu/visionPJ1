//
//  DLib: A simple image processing library.
//
//  David Crandall, 2003-2005
//  crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
//
//
class DTriple
{
    public:
        DTriple(double _r, double _g, double _b) : r(_r), g(_g), b(_b) {}
        DTriple() {}

        void operator+=(const DTriple &other)
        {
            r+=other.r, g+=other.g, b+=other.b;

            
        }
        void operator/=(const double &d)
        {
            r/=d, g/=d, b/=d;
        }

        double r, g, b;
};

