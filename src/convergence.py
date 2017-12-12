import math

class ConvergenceTest(object):
    def __init__(self):
        self.l2 = 1
        """
        class for testing convergence using the L2 engineering norm
        n is the iteration index, i is the vector content index, I is total number of entries in vector.
        """

    def isConverged(self, vec_n, vec_n1, epsilon):
        sum1 = 0
        for i in range(len(vec_n)):
            error1 = ((vec_n[i] - vec_n1[i])/vec_n[i]) ** 2
            sum1 += error1
        I = len(vec_n)
        l2 = math.sqrt(sum1 / I)

        if l2 < epsilon:
            print "Converged! l2 %g\n" %(l2)
            self.returnL2(l2)
            return True
        else:
            print "Not converged; l2 %g\n" %(l2)
            self.returnL2(l2)
            return False

    def returnL2(self, l2):
        self.l2 = l2
