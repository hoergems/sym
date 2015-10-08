from sympy import *
import numpy as np

from scipy.integrate import ode, odeint

class Test:
    def __init__(self):
        #self.test_fot()
        #return
        l1, l2, l3 = symbols("l1 l2 l3")
        theta1, theta2, theta3 = symbols("theta1 theta2 theta3")
        dot_theta1, dot_theta2, dot_theta3 = symbols("dot_theta1 dot_theta2 dot_theta3")
        alpha1, alpha2, alpha3 = symbols("alpha1 alpha2 alpha3")
        d1, d2, d3 = symbols("d1 d2 d3")
        m1x, m1y, m1z = symbols("m1x m1y m1z")
        m2x, m2y, m2z = symbols("m2x m2y m2z")
        m3x, m3y, m3z = symbols("m3x m3y m3z")
        
        rho1, rho2 = symbols("r1 r2")
        
        m1, m2, m3 = symbols("m1 m2 m3")
        
        g = symbols("g")
        '''ms = [Matrix([[m1x], 
                      [m1y], 
                      [m1z]]),
              Matrix([[m2x], 
                      [m2y], 
                      [m2z]]),
              Matrix([[m3x], 
                      [m3y], 
                      [m3z]])]'''
        
        ms = [Matrix([[l1 / 2.0], 
                      [0.0], 
                      [0.0]]),
              Matrix([[l2 / 2.0], 
                      [0.0], 
                      [0.0]])]
        
        
        I1x, I1y, I1z = symbols("I1x I1y I1z")
        I2x, I2y, I2z = symbols("I2x I2y I2z")
        I3x, I3y, I3z = symbols("I3x I3y I3z")       
        
        
        
        lc1, lc2 = symbols("lc1 lc2")
        
        O0 = Matrix([[0.0],
                     [0.0],
                     [0.0]])
        dhc1 = self.denavit_hartenberg(theta1, 0.0, lc1, 0.0)
        
        Oc1 = Matrix([dhc1.col(3)[j] for j in xrange(3)])        
        dhl1 = self.denavit_hartenberg(0.0, 0.0, lc1, 0.0)
        r = dhc1 * dhl1
        O1 = Matrix([r.col(3)[j] for j in xrange(3)])              
        
        dhc2 = self.denavit_hartenberg(theta2, 0.0, lc2, 0.0)
        r = dhc1 * dhl1 * dhc2
        Oc2 = Matrix([r.col(3)[j] for j in xrange(3)])
        dhl2 = self.denavit_hartenberg(0.0, 0.0, lc2, 0.0)
        r = dhc1 * dhl1 * dhc2 * dhl2
        On = Matrix([r.col(3)[j] for j in xrange(3)]) 
        
        z_i1 = Matrix([[0.0],
                       [0.0],
                       [1.0]])
        
        r1 = Matrix([z_i1.cross(Oc2 - O1)])         
        #print simplify(r1)
        
        Jvs, Ocs =  self.get_link_jacobians([l1, l2], 
                                            [theta1, theta2],
                                            [0.0, 0.0],
                                            [0.0, 0.0],
                                            ms)
        Is = [[I1x, I1y, I1z],
              [I2x, I2y, I2z],
              [I3x, I3y, I3z]]
        M_is = self.construct_link_inertia_matrices([m1, m2], Is)
        M = simplify(self.calc_inertia_matrix(Jvs, M_is, [l1, l2]))        
        C = self.calc_coriolis_matrix([theta1, theta2], [dot_theta1, dot_theta2], M)
        
        N = self.calc_generalized_forces([theta1, theta2],
                                         [dot_theta1, dot_theta2], 
                                         Ocs, 
                                         [m1, m2], 
                                         g)
        
        M = simplify(M)
        C = simplify(C)
        N = simplify(N)
        self.dynamic(M, C, N, [theta1, theta2], [dot_theta1, dot_theta2], [rho1, rho2])
        
    def dynamic(self, M, C, N, thetas, dot_thetas, rs):
        M_inv = M.inv()
        Thetas = Matrix([[thetas[i]] for i in xrange(len(thetas))])
        Dotthetas = Matrix([[dot_thetas[i]] for i in xrange(len(dot_thetas))])
        Rs = Matrix([[rs[i]] for i in xrange(len(rs))])
        
        k = -M_inv * (C * Dotthetas + N) + M_inv * Rs       
        m1 = Matrix([dot_thetas[i] for i in xrange(len(dot_thetas))])
        m2 = Matrix([0.0 for i in xrange(len(dot_thetas))])
        h = m1.col_join(k)
        print simplify(h)
        
    def test_fot(self):
        q, qdot, qdotdot = symbols("q q_dot q_dot_dot")
        r = symbols("r")        
        x1, x2, x3 = symbols("x1 x2 x3")        
        
        f = Matrix([[qdot],
                    [sin(q) + cos(qdot)]])        
        
        A = f.jacobian([q])
        B = f.jacobian([qdot])
        
        self.initial = [0.0, 1.0]
        
        fot = f.subs([(q, x1), (qdot, x2)]) + A.subs([(q, x1), (qdot, x2)]) * (q - x1) + B.subs([(q, x1), (qdot, x2)]) * (q - x1)
        fot = simplify(fot)        
        
        #print fot.subs([(q, 0.0), (q_dot, 0.0), (r, 1.0)])
        t = np.linspace(0.0, 0.3, 100)
        
        
        eq = odeint(self.f, np.array([0.0, 0.0]), t)
        print "================="
        print eq
        
        #print fot
                
    def f(self, y, t):
        x1 = self.initial[0]
        x2 = self.initial[1]        
        q = y[0]
        qdot = y[1]
         
        s = np.array([q - x1 + x2,
                      (-q + x1)*np.sin(x2) + (q - x1)*np.cos(x1) + np.sin(x1) + np.cos(x2)])    
        
        #print "s " + str(s)
        return s
        
        
           
        
    def calc_generalized_forces(self, 
                                thetas, 
                                dot_thetas, 
                                Ocs, 
                                ms, 
                                g):              
        V = 0.0
        for i in xrange(len(Ocs)):                                 
            V += ms[i] * g * Ocs[i][1] 
            
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
                            
        
    def get_link_jacobian(self, links, thetas, alphas):
        
        Os = [Matrix([[0.0],
                      [0.0],
                      [0.0]])]
        z0s= [Matrix([[0.0],
                      [0.0],
                      [1.0]])]
        I = Matrix([[1.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0]])
        dhs = []        
        for i in xrange(len(links)):
            dhs.append(self.denavit_hartenberg(thetas[i], alphas[i], links[i], 0.0))
            res = I
            for k in xrange(i + 1):
                res *= dhs[k]                           
            col3 = res.col(2)
            col4 = res.col(3)            
            z = Matrix([col3[j] for j in xrange(3)])                
            O = Matrix([col4[j] for j in xrange(3)])
            
            Os.append(O)
            z0s.append(z)        
        js = [] 
            
        for i in xrange(len(links)):
            Od = Os[-1] - Os[i]            
            r1 = Matrix([z0s[i].cross(Os[-1] - Os[i])])            
            m = Matrix([r1, z0s[i]]).transpose()
            
            #print Matrix([m, m]).transpose()
                  
            js.append(m)
        print js[0]
        print "====="
        print js[1]
        print "====="
        print js[2]
        j = Matrix([js[i] for i in xrange(len(js))])        
        return simplify(j.transpose())
            
        
    def denavit_hartenberg(self, theta, alpha, a, d):
        matrix = Matrix([[cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), a * cos(theta)],
                         [sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta)],
                         [0, sin(alpha), cos(alpha), d],
                         [0, 0, 0, 1.0]])
        return matrix
        
if __name__ == "__main__":
    Test()