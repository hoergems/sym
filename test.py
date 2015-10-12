
from sympy import *
import numpy as np
import time
import os
from sympy.printing import print_ccode

from scipy.integrate import ode, odeint

class Test:
    def __init__(self):        
        l1, l2 = symbols("ls[0] ls[1]")
        theta1, theta2 = symbols("x[0] x[1]")
        dot_theta1, dot_theta2 = symbols("x[2] x[3]") 
        theta1_star, theta2_star = symbols("thetas_star_[0] thetas_star_[1]")        
        dot_theta1_star, dot_theta2_star = symbols("dot_thetas_star_[0] dot_thetas_star_[1]") 
        rho1, rho2 = symbols("rhos[0] rhos[1]")
        rho1_star, rho2_star = symbols("rhos_star_[0] rhos_star_[1]")
              
        alpha1, alpha2 = symbols("alpha1 alpha2")
        d1, d2, d3 = symbols("d1 d2 d3")
        m1x, m1y, m1z = symbols("m1x m1y m1z")
        m2x, m2y, m2z = symbols("m2x m2y m2z")
        
        
        m1, m2 = symbols("m1 m2")
        g = symbols("g")
        
        """
        Coordinates of the center of masses relative to the link frames
        """
        ms = [Matrix([[l1 / 2.0], 
                      [0.01], 
                      [0.0]]),
              Matrix([[l2 / 2.0], 
                      [-0.001], 
                      [0.0]])]
        
        
        I1x, I1y, I1z = symbols("I1x I1y I1z")
        I2x, I2y, I2z = symbols("I2x I2y I2z")
        I3x, I3y, I3z = symbols("I3x I3y I3z")       
        
        self.q = [theta1, theta2]
        self.q_star = [theta1_star, theta2_star]
        self.dotq = [dot_theta1, dot_theta2]
        self.dotq_star = [dot_theta1_star, dot_theta2_star]
        self.r = [rho1, rho2]
        self.r_star = [rho1_star, rho2_star]
        
        """
        Get the Jacobians of the links expressed in the robot's base frame
        """
        Jvs, Ocs =  self.get_link_jacobians([l1, l2], 
                                            self.q,
                                            [0.0, 0.0],
                                            [0.0, 0.0],
                                            ms)
        
        
        """
        Inertia parameters of the links
        """
        Is = [[I1x, I1y, I1z],
              [I2x, I2y, I2z],
              [I3x, I3y, I3z]]
        M_is = self.construct_link_inertia_matrices([m1, m2], Is)
        M = simplify(self.calc_inertia_matrix(Jvs, M_is, [l1, l2]))        
        C = simplify(self.calc_coriolis_matrix(self.q, self.dotq, M))        
        N = simplify(self.calc_generalized_forces(self.q,
                                         self.dotq, 
                                         Ocs, 
                                         [m1, m2], 
                                         g))
        print "Get dynamic model"
        f = self.get_dynamic_model(M, C, N, self.q, self.dotq, self.r)        
        #f = simplify(f)      
        fot = self.taylor_approximation(f, 
                                        self.q, 
                                        self.dotq, 
                                        self.q_star, 
                                        self.dotq_star, 
                                        self.r,
                                        self.r_star)
        fot = fot.subs([(m1, 0.2), 
                        (m2, 0.41),
                        (l1, 1.0),
                        (l2, 1.0),
                        (I1x, 0.003),
                        (I1y, 0.000012),
                        (I1z, 0.01),
                        (I2x, 0.0035),
                        (I2y, 0.00012),
                        (I2z, 0.012),
                        (g, -9.81)])
        print "Generating c++ code..."
        self.gen_cpp_code(fot)
        #cmd = "cd build && cmake .. && make -j8"
        #os.system(cmd)
        print "Done"
        """
        Substitude all parameters here
        """
        
        """
        Call test_fot2
        """
        
        """
        Generate c++ code
        """
        
    def gen_cpp_code(self, fot):
        lines = list(open("integrate.cpp", 'r'))
        idx = -1       
        for i in xrange(len(lines)):
            if "void Integrate::ode" in lines[i]:
                idx = i        
        cpp_string = "void Integrate::ode(const state_type &x , state_type &dxdt , double t) const {"
        cpp_string += "std::vector<double> terms({"
        for i in xrange(len(fot)):
            cpp_string += ccode(fot[i])
            if i != len(fot) - 1:
                cpp_string += ", "
        cpp_string += "});"
        
        cpp_string += "dxdt.clear();"
        cpp_string += "for(size_t i = 0; i < x.size(); i++) { dxdt.push_back(terms[i]); }"
        cpp_string += "} \n"
                
        if not idx == -1:
            lines[idx] = cpp_string
            
        os.remove("integrate.cpp")
        with open("integrate.cpp", 'a+') as f:
            for line in lines:
                f.write(line)
        
        
    def get_dynamic_model(self, M, C, N, thetas, dot_thetas, rs):
        M_inv = M.inv()
        Thetas = Matrix([[thetas[i]] for i in xrange(len(thetas))])
        Dotthetas = Matrix([[dot_thetas[i]] for i in xrange(len(dot_thetas))])
        Rs = Matrix([[rs[i]] for i in xrange(len(rs))])
        
        m_upper = Matrix([dot_thetas[i] for i in xrange(len(dot_thetas))])
        m_lower = -M_inv * (C * Dotthetas + N) + M_inv * Rs        
        h = m_upper.col_join(m_lower)        
        return h
        
    def taylor_approximation(self, f, thetas, dot_thetas, thetas_star, dot_thetas_star, rs, rs_star):
        A = f.jacobian(thetas)
        B = f.jacobian(dot_thetas)
        C = f.jacobian(rs)
        for i in xrange(len(thetas)):
            A = A.subs(thetas[i], thetas_star[i])
            A = A.subs(dot_thetas[i], dot_thetas_star[i])
            A = A.subs(rs[i], rs_star[i])
            
            B = B.subs(thetas[i], thetas_star[i])
            B = B.subs(dot_thetas[i], dot_thetas_star[i])
            B = B.subs(rs[i], rs_star[i])
            
            C = C.subs(thetas[i], dot_thetas[i])
            C = C.subs(dot_thetas[i], dot_thetas_star[i])
            C = C.subs(rs[i], rs_star[i])
            
            f = f.subs(thetas[i], thetas_star[i])
            f = f.subs(dot_thetas[i], dot_thetas_star[i])
            f = f.subs(rs[i], rs_star[i])        
        
        q = Matrix([[thetas[i]] for i in xrange(len(thetas))])
        dot_q = Matrix([[dot_thetas[i]] for i in xrange(len(dot_thetas))])
        r = Matrix([[rs[i]] for i in xrange(len(rs))])
        
        q_star = Matrix([[thetas_star[i]] for i in xrange(len(thetas_star))])
        dot_q_star = Matrix([[dot_thetas_star[i]] for i in xrange(len(dot_thetas_star))])
        r_star = Matrix([[rs_star[i]] for i in xrange(len(rs_star))])
        
        print "Build taylor approximation"
        fot = f + A * (q - q_star) + B * (dot_q - dot_q_star) #+ C * (r - r_star)
        return fot 
    
    def test_fot2(self, f, thetas, dot_thetas, thetas_star, dot_thetas_star, rs, rs_star):
        self.fot = self.taylor_approximation(f, thetas, dot_thetas, thetas_star, dot_thetas_star, rs, rs_star)       
         
        
    def test_fot(self, f):
        q1, q2, qdot1, qdot2, qdotdot1, qdotdot2 = symbols("q1 q2 qdot1 qdot2 qdotdot1 qdotdot2")
        r1, r2 = symbols("r1 r2")       
        x1_1, x1_2, x2_1, x2_2, x3_1, x3_2 = symbols("x1_1 x1_2 x2_1 x2_2 x3_1 x3_2")
        
        self.q = [q1, q2]
        self.qdot = [qdot1, qdot2]
        self.r = [r1, r2]
        
        self.q_star = [x1_1, x1_2]
        self.qdot_star = [x2_1, x2_2]
        self.r_star = [x3_1, x3_2]
        self.initial = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
        self.fot = self.taylor_approximation(f, 
                                             q, 
                                             qdot, 
                                             q_star, 
                                             qdot_star, 
                                             r, 
                                             r_star)        
        q_star.extend(qdot_star)
        q_star.extend(r_star)
        print simplify(fot)
        for i in xrange(len(q_star)):
            fot = fot.subs(q_star[i], self.initial[0])
        
        
        #print fot.subs([(q, 0.0), (q_dot, 0.0), (r, 1.0)])
        t = np.linspace(0.0, 0.3, 2)
        
        t0 = time.time()
        eq = odeint(self.f, np.array(self.initial[:4]), t)
        print "calc took " + str(time.time() - t0) + " seconds"
        print "================="
        print eq
        
                
    def f(self, y, t):        
        x1_1 = self.initial[0]
        x1_2 = self.initial[1]
        x2_1 = self.initial[2]
        x2_2 = self.initial[3]
        x3_1 = self.initial[4]
        x3_2 = self.initial[5] 
        q1 = y[0]
        q2 = y[1]
        qdot1 = y[2]
        qdot2 = y[3]
        self.fot.subs(q1, )
        return s   
        
    def calc_generalized_forces(self, 
                                thetas, 
                                dot_thetas, 
                                Ocs, 
                                ms, 
                                g):
        print Ocs              
        V = 0.0
        for i in xrange(len(Ocs)):                                 
            V += ms[i] * g * Ocs[i][2] 
            
        N = Matrix([[diff(V, thetas[i])] for i in xrange(len(thetas))])
        return N
        
        
    def calc_coriolis_matrix(self, thetas, dot_thetas, M):
        C = Matrix([[0.0 for m in xrange(len(thetas))] for n in xrange(len(thetas))])
        for i in xrange(len(thetas)):
            for j in xrange(len(thetas)):
                val = 0.0
                for k in xrange(len(thetas)):                                              
                    val += self.calc_christoffel_symbol(i, j, k, thetas, M) * dot_thetas[k]
                C[i, j] = val
        return C   
    
    def calc_christoffel_symbol(self, i, j, k, thetas, M):
        t_i_j_k = 0.5 * (diff(M[i, j], thetas[k]) + 
                         diff(M[i, k], thetas[j]) -
                         diff(M[k, j], thetas[i]))
        return t_i_j_k
    
    def calc_inertia_matrix(self, Jvs, M_is, links): 
        lc2 = links[1] / 2.0
        res = Matrix([[0.0 for n in xrange(len(Jvs))] for m in xrange(len(Jvs))])        
        for i in xrange(len(Jvs)):
            res += Jvs[i].transpose() * M_is[i] * Jvs[i]        
        return res
    
    def construct_link_inertia_matrices(self, ms, Is):
        M_is = [Matrix([[ms[i], 0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, ms[i], 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, ms[i], 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, Is[i][0], 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, Is[i][1], 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, Is[i][2]]]) for i in xrange(len(ms))]
        return M_is
        
    def get_link_jacobians(self, links, thetas, alphas, ds, ms):
        """
        Center points of the links
        """
        midpoints = [Matrix([[links[i] / 2.0],
                             [0.0],
                             [0.0]]) for i in xrange(len(links))]
        
        """
        Vectors form the center points to the center of mass
        """
        mid_to_m_vectors = [ms[i] - midpoints[i] for i in xrange(len(ms))]
        
        """
        Vectors from the center of mass to the next link
        """
        m_to_link_vectors = [Matrix([[links[i]],
                                      [0.0],
                                      [0.0]]) - ms[i] for i in xrange(len(links))]
        
        trans_matrices = [self.denavit_hartenberg(0.0, 0.0, mid_to_m_vectors[i][0], 0.0) for i in xrange(len(mid_to_m_vectors))]        
        trans_matrices2 = [self.denavit_hartenberg(0.0, 0.0, m_to_link_vectors[i][0], 0.0) for i in xrange(len(m_to_link_vectors))]
        
        dhcs = [self.denavit_hartenberg(thetas[i], alphas[i], midpoints[i][0], ds[i]) for i in xrange(len(midpoints))]
        
        """
        Transformations from the link origins to the center of masses
        """
        dhcs = [dhcs[i] * trans_matrices[i] for i in xrange(len(dhcs))]
        
        Os = [Matrix([[0.0],
                      [0.0],
                      [0.0]])]
        zs = [Matrix([[0.0],
                      [0.0],
                      [1.0]])]
        Ocs = []
        zcs = []
        I = Matrix([[1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0]])
        res = I
        for i in xrange(len(links)):
            res *= dhcs[i]
            col3 = res.col(2)
            col4 = res.col(3)
            
            z = Matrix([col3[j] for j in xrange(3)])
            O = Matrix([col4[j] for j in xrange(3)])
            Ocs.append(O)
            zcs.append(z)
            res = res * trans_matrices2[i]
            col3 = res.col(2)
            col4 = res.col(3)
            
            z = Matrix([col3[j] for j in xrange(3)])
            O = Matrix([col4[j] for j in xrange(3)])
            Os.append(O)
            zs.append(z)
                  
        r1 = Matrix([zcs[0].cross(Ocs[1] - Os[1])])
        Jvs = []
        for i in xrange(len(links)):
            Jv = Matrix([[0.0 for m in xrange(len(links))] for n in xrange(6)])
            for k in xrange(i + 1):
                r1 = Matrix(zcs[i].cross(Ocs[i] - Os[k]))
                for t in xrange(3):
                    Jv[t, k] = r1[t, 0]
                    Jv[t + 3, k] = zcs[i][t, 0]
            Jvs.append(simplify(Jv))
        return Jvs, Ocs
        
    def denavit_hartenberg(self, theta, alpha, a, d):
        matrix = Matrix([[cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), a * cos(theta)],
                         [sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta)],
                         [0, sin(alpha), cos(alpha), d],
                         [0, 0, 0, 1.0]])
        return matrix
        
if __name__ == "__main__":
    Test()