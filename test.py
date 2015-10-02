from sympy import *

class Test:
    def __init__(self):
        theta1, theta2, theta3 = symbols("theta1 theta2 theta3")
        m1, m2, m3 = symbols("m1 m2 m3")
        a, d, alpha = symbols("a d alpha") 
        
        I1x, I1y, I1z = symbols("I1x I1y I1z")
        I2x, I2y, I2z = symbols("I2x I2y I2z")
        I3x, I3y, I3z = symbols("I3x I3y I3z")
        
        l1, l2, l3 = symbols("l1 l2 l3")
        dh1 = self.denavit_hartenberg(theta1, 0.0, l1, 0.0)
        dh2 = self.denavit_hartenberg(theta2, 0.0, l2, 0.0)  
        
        z1 = Matrix([[0, 0, 0, 1]]).transpose()
           
        
        res =  dh1.multiply(dh2.multiply(z1))
        
        
        
        matr1 = Matrix([[m1, 0, 0, 0, 0, 0],
                        [0, m1, 0, 0, 0, 0],
                        [0, 0, m1, 0, 0, 0],
                        [0, 0, 0, I1x, 0, 0],
                        [0, 0, 0, 0, I1y, 0],
                        [0, 0, 0, 0, 0, I1z]])
        
        matr2 = Matrix([[m2, 0, 0, 0, 0, 0],
                        [0, m2, 0, 0, 0, 0],
                        [0, 0, m2, 0, 0, 0],
                        [0, 0, 0, I2x, 0, 0],
                        [0, 0, 0, 0, I2y, 0],
                        [0, 0, 0, 0, 0, I2z]])
        
        matr3 = Matrix([[m3, 0, 0, 0, 0, 0],
                        [0, m3, 0, 0, 0, 0],
                        [0, 0, m3, 0, 0, 0],
                        [0, 0, 0, I3x, 0, 0],
                        [0, 0, 0, 0, I3y, 0],
                        [0, 0, 0, 0, 0, I3z]])
        
        matr4 = Matrix([matr1, matr2])
        
        
        links = [l1, l2]
        thetas = [theta1, theta2]
        alphas = [0.0, 0.0]
        self.get_link_jacobian(links, thetas, alphas)
        
        
        
        
    def get_link_jacobian(self, links, thetas, alphas):
        dhs = []
        Os = [Matrix([[0.0],
                      [0.0],
                      [0.0]])]
        z0s= []        
        for i in xrange(len(links)):
            dhs.append(self.denavit_hartenberg(thetas[i], alphas[i], links[i], 0.0))
            col3 = dhs[-1].col(2)
            col4 = dhs[-1].col(3)
            z = Matrix([col3[j] for j in xrange(3)])
                             
            O = Matrix([col4[j] for j in xrange(3)]) 
                             
            Os.append(O)
            z0s.append(z)
        js = [] 
              
        for i in xrange(len(links)):
            Od = Os[-1] - Os[i]
            print Od
            r1 = Matrix([z0s[i].cross(Os[-1] - Os[i])])
            #print r1
            m = Matrix([r1, z0s[i]]).transpose()
            #print Matrix([m, m]).transpose()
                  
            js.append(m)
        j = Matrix([js[i] for i in xrange(len(js))])
        return j.transpose()
            
        
    def denavit_hartenberg(self, theta, alpha, a, d):
        matrix = Matrix([[cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha), a * cos(theta)],
                         [sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha), a * sin(theta)],
                         [0, sin(alpha), cos(alpha), d],
                         [0, 0, 0, 1.0]])
        return matrix
        
if __name__ == "__main__":
    Test()