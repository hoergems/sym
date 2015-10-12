#include "integrate.cpp"

int main(int argc, char** argv) {
	double k;
	shared::Integrate integrator(k);
	std::vector<double> thetas_star({0.0, 0.0});
	std::vector<double> dot_thetas_star({0.0, 0.0});
	std::vector<double> rho_star({0.0, 0.0});
	
	std::vector<double> current_state({0.0, 0.0, 0.0, 0.0});
	double t0 = 0.0;
	double t1 = 0.1;
	double step_size = 0.00001;
    integrator.setup(thetas_star, dot_thetas_star, rho_star);
    
    integrator.do_integration(current_state, t0, t1, step_size);
}