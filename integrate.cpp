#include "integrate.hpp"

using namespace boost::numeric::odeint;

namespace shared{

template<class T>
struct VecToList
{
    static PyObject* convert(const std::vector<T>& vec)
    {
        boost::python::list* l = new boost::python::list();
        for(size_t i = 0; i < vec.size(); i++)
            (*l).append(vec[i]);

        return l->ptr();
    }
};

Integrate::Integrate(double k) {	
}

void Integrate::do_integration(state_type &x, double &t0, double &te, double &step_size) {	
	integrate_const(runge_kutta4<state_type>() ,
		            std::bind(&Integrate::ode , this , pl::_1 , pl::_2 , pl::_3),
		            x , t0 , te , step_size);
}

void Integrate::setup(std::vector<double> &thetas_star, 
		              std::vector<double> &dot_thetas_star, 
		              std::vector<double> &rhos_star) {
	thetas_star_.clear();
	dot_thetas_star_.clear();
	rhos_star_.clear();
	for (size_t i = 0; i < thetas_star.size(); i++) {
		thetas_star_.push_back(thetas_star[i]);
		dot_thetas_star_.push_back(dot_thetas_star[i]);
		rhos_star_.push_back(rhos_star[i]);
	}
}

void Integrate::ode(const state_type &x , state_type &dxdt , double t) const {std::vector<double> terms({x[2], x[3], rhos_star_[0]*(1.0/(0.41*cos(thetas_star_[1]) + 0.5845) + pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845))) - rhos_star_[1]*(0.205*cos(thetas_star_[1]) + 0.1145)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (-dot_thetas_star_[0] + x[2])*(0.41*dot_thetas_star_[0]*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 0.41*dot_thetas_star_[1]*(-1/(0.41*cos(thetas_star_[1]) + 0.5845) - pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)))*sin(thetas_star_[1])) + (-dot_thetas_star_[1] + x[3])*(-1/(0.41*cos(thetas_star_[1]) + 0.5845) - pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)))*(-0.205*dot_thetas_star_[0]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1])) + (-thetas_star_[0] + x[0])*(-9.81*(-1/(0.41*cos(thetas_star_[1]) + 0.5845) - pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)))*(-0.51*sin(thetas_star_[0]) - 0.205*sin(thetas_star_[0] + thetas_star_[1])) + 2.01105*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[0] + thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)) + (-thetas_star_[1] + x[1])*(rhos_star_[0]*(0.41*sin(thetas_star_[1])/pow(0.41*cos(thetas_star_[1]) + 0.5845, 2) + 0.41*pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)*sin(thetas_star_[1])/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*pow(0.41*cos(thetas_star_[1]) + 0.5845, 2)) - 0.41*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)) - 0.08405*pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)*sin(thetas_star_[1])*cos(thetas_star_[1])/(pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2)*(0.41*cos(thetas_star_[1]) + 0.5845))) + 0.205*rhos_star_[1]*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + 0.08405*rhos_star_[1]*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])*cos(thetas_star_[1])/pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2) - 0.205*(0.205*pow(dot_thetas_star_[0], 2)*sin(thetas_star_[1]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 0.08405*(0.205*pow(dot_thetas_star_[0], 2)*sin(thetas_star_[1]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])*cos(thetas_star_[1])/pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2) + (0.205*pow(dot_thetas_star_[0], 2)*cos(thetas_star_[1]) + 2.01105*sin(thetas_star_[0] + thetas_star_[1]))*(0.205*cos(thetas_star_[1]) + 0.1145)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (-1/(0.41*cos(thetas_star_[1]) + 0.5845) - pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)))*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*cos(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*cos(thetas_star_[1]) + 2.01105*sin(thetas_star_[0] + thetas_star_[1])) + (-0.41*sin(thetas_star_[1])/pow(0.41*cos(thetas_star_[1]) + 0.5845, 2) - 0.41*pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)*sin(thetas_star_[1])/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*pow(0.41*cos(thetas_star_[1]) + 0.5845, 2)) + 0.41*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)) + 0.08405*pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)*sin(thetas_star_[1])*cos(thetas_star_[1])/(pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2)*(0.41*cos(thetas_star_[1]) + 0.5845)))*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1]) - 5.0031*cos(thetas_star_[0]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))) + (0.205*pow(dot_thetas_star_[0], 2)*sin(thetas_star_[1]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*(0.205*cos(thetas_star_[1]) + 0.1145)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (-1/(0.41*cos(thetas_star_[1]) + 0.5845) - pow(0.205*cos(thetas_star_[1]) + 0.1145, 2)/((0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)*(0.41*cos(thetas_star_[1]) + 0.5845)))*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1]) - 5.0031*cos(thetas_star_[0]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1])), -rhos_star_[0]*(0.205*cos(thetas_star_[1]) + 0.1145)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + rhos_star_[1]*(0.41*cos(thetas_star_[1]) + 0.5845)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (-dot_thetas_star_[0] + x[2])*(-0.41*dot_thetas_star_[0]*(0.41*cos(thetas_star_[1]) + 0.5845)*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 0.41*dot_thetas_star_[1]*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)) + (-dot_thetas_star_[1] + x[3])*(0.205*cos(thetas_star_[1]) + 0.1145)*(-0.205*dot_thetas_star_[0]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1]))/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (-thetas_star_[0] + x[0])*(-9.81*(-0.51*sin(thetas_star_[0]) - 0.205*sin(thetas_star_[0] + thetas_star_[1]))*(0.205*cos(thetas_star_[1]) + 0.1145)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 2.01105*(0.41*cos(thetas_star_[1]) + 0.5845)*sin(thetas_star_[0] + thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)) + (-thetas_star_[1] + x[1])*(0.205*rhos_star_[0]*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + 0.08405*rhos_star_[0]*(0.205*cos(thetas_star_[1]) + 0.1145)*sin(thetas_star_[1])*cos(thetas_star_[1])/pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2) - 0.41*rhos_star_[1]*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 0.08405*rhos_star_[1]*(0.41*cos(thetas_star_[1]) + 0.5845)*sin(thetas_star_[1])*cos(thetas_star_[1])/pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2) + 0.41*(0.205*pow(dot_thetas_star_[0], 2)*sin(thetas_star_[1]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + 0.08405*(0.205*pow(dot_thetas_star_[0], 2)*sin(thetas_star_[1]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*(0.41*cos(thetas_star_[1]) + 0.5845)*sin(thetas_star_[1])*cos(thetas_star_[1])/pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2) - (0.205*pow(dot_thetas_star_[0], 2)*cos(thetas_star_[1]) + 2.01105*sin(thetas_star_[0] + thetas_star_[1]))*(0.41*cos(thetas_star_[1]) + 0.5845)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (0.205*cos(thetas_star_[1]) + 0.1145)*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*cos(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*cos(thetas_star_[1]) + 2.01105*sin(thetas_star_[0] + thetas_star_[1]))/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 0.205*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1]) - 5.0031*cos(thetas_star_[0]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*sin(thetas_star_[1])/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) - 0.08405*(0.205*cos(thetas_star_[1]) + 0.1145)*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1]) - 5.0031*cos(thetas_star_[0]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*sin(thetas_star_[1])*cos(thetas_star_[1])/pow(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179, 2)) - (0.205*pow(dot_thetas_star_[0], 2)*sin(thetas_star_[1]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))*(0.41*cos(thetas_star_[1]) + 0.5845)/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179) + (0.205*cos(thetas_star_[1]) + 0.1145)*(-0.205*dot_thetas_star_[0]*dot_thetas_star_[1]*sin(thetas_star_[1]) - 0.205*dot_thetas_star_[1]*(dot_thetas_star_[0] + dot_thetas_star_[1])*sin(thetas_star_[1]) - 5.0031*cos(thetas_star_[0]) - 2.01105*cos(thetas_star_[0] + thetas_star_[1]))/(0.042025*pow(sin(thetas_star_[1]), 2) + 0.01179)});} 

BOOST_PYTHON_MODULE(libintegrate) {
    using namespace boost::python;
    
    class_<std::vector<double> > ("v_double")
             .def(vector_indexing_suite<std::vector<double> >());
    
    class_<Integrate>("Integrate", init<double>())
                        .def("doIntegration", &Integrate::do_integration)                                  
    ;
}

}