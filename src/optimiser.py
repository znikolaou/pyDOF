from libs import *
from objectiveFunction import *
from reconstructedTransferFunction import *
#-----------------------------------------------------------------------
class Optimiser(Const, Log):
    """
    Class for the optimiser. 

    Attributes: 
        target (Filter or ReconstructedTransferFunction): the target to 
            optimise against. 
        discreteFilter (DiscreteFilter): the discrete filter object 
            whose coeffs we solve for.
        method (string): the optimiser method to use for constrained 
            optimisation. Suggested candidates are  "SLSQP", and 
            "trust-constr". 
        isZeroDampingAtPi (boolean): True, False. If true, then 
            \hat{G}(pi)=0 always. 
        isPositive (boolean): True, False. If true, then 
            \hat{G}(k) > 0 for all k in [0, pi). 
        isAdaptive (boolean): True, False. If true, the discrete filter 
            stencil size will increase until the target optimisation error 
            Const.TOL_OPTIMISER is reached.
        lamda (float): regularisation/smoothness parameter.

    Author: Z. Nikolaou.
    """
    def __init__(self, target, discreteFilter, method="SLSQP",  \
            isZeroDampingAtPi=True, isPositive=True, isMonotonic=False, isAdaptive=True, \
            lamda=0.e0): 
        Log.__init__(self, type(self).__name__) 
        self._target=target
        self._dfilter=discreteFilter
        self._method=method
        self._isZeroDampingAtPi=isZeroDampingAtPi
        self._isPositive=isPositive
        self._isMonotonic=isMonotonic
        self._isAdaptive=isAdaptive
        self._lamda=lamda
        self._objFunc=ObjectiveFunction(self._target, self._dfilter, lamda=self._lamda)
        self._nk=self._dfilter.getNk()
        self._isInverse=not self._dfilter.isForward()
        self._nVanCittert=self.__getNvanCittert()
        self._solution=None
        self._error=None
        self._uncMethod='BFGS'
        self.__printInfo()
 
    def __printInfo(self):
        if(self._isInverse):
            print('')
            print("Optimiser will run for inverse discrete filter")
            print('')
        else:
            print('')
            print("Optimiser will run for forward discrete filter")
            print('')

    def __getNvanCittert(self):
        nv=None
        if(self._isInverse): nv=self._target.getNvanCittert()
        return nv

    def __getUnConstrainedSolution(self, xo): 
        #solution=minimize(self._objFunc.get, xo, method=self._uncMethod)
        solution=minimize(self._objFunc.get, xo, method='BFGS')
        return solution

    def __getConstrainedSolution(self, xo, constraintsList):
        solution=minimize(self._objFunc.get, xo, method=self._method, \
                constraints=constraintsList)
        return solution

    def run(self): 
        isConverged=False
        maxIter=1
        itr=1
        nmp1=self._dfilter.getNmPlusOne()
        errorOld=None
        solOld=None
        if self._isAdaptive: maxIter=Const.MAX_ITERATIONS
        for i in range(maxIter):
            x=self.__getInitialGuessForThreePointForward(nmp1)
            if(self._isInverse): 
                pass
                #0
                #x=self.__getInitialGuess(nmp1)
                #x=self.__getInitialGuessForThreePointForward(nmp1)
                #solUnc=self.__getUnConstrainedSolution(self.__getInitialGuess(nmp1))
                #x=solUnc.x 
                #1 slightly better
                #solUnc=self.__getUnConstrainedSolution(self.__getInitialGuessForThreePointForward(nmp1))
                #x=solUnc.x
                #2 slightly best
                #solUnc=self.__getUnConstrainedSolution(self.__getInitialGuessForThreePointInverse(nmp1, self._nVanCittert))
                #x=solUnc.x
                #x=self.__getInitialGuessForThreePointInverse(nmp1, self._nVanCittert)
                  
            self._solution=self.__getConstrainedSolution(x, self.__getConstraintsList(nmp1))
            self._error=self._objFunc.getLastError()
            self.__printIterInfo(itr, nmp1, self._error) 
            if self.__isConverged(self._error): 
                isConverged=True
                break
            #TODO: implement more sophisticated error controller
            if self.__isDiverging(itr, self._error, errorOld): 
                isConverged=False
                self._solution=solOld
                self._error=errorOld
                break
            nmp1+=1
            itr+=1
            errorOld=self._error
            solOld=self._solution
            xOld=solOld.x
 
        if isConverged:
            print('\n*Success: converged at iteration ', itr, 'error=', self._error)
        else:
            print('\n*Failure to converge for tolerance=', Const.TOL_OPTIMISER)
            print('\tOptimisation error attained=', self._error)

        self._dfilter.setCoeffs(self.get())
        self._dfilter.setOptError(self._error)
               
    def __isConverged(self, error):
        if error<=Const.TOL_OPTIMISER: return True
        else: return False
    
    def __isDiverging(self, itr, error, errorOld):
        if itr>1:
            if error>Const.DIVERGENCE_FACTOR*errorOld:
                    print('*Increasing error detected, error=', error, 'errorOld=', errorOld)
                    print('*Reverting to last solution ...')
                    return True
            else:
                return False
        else:
            return False

    def get(self): return self._solution.x
 
    def __printIterInfo(self, itr, nmp1, error):
        print("Iteration=", itr, "stencil size (nm+1)=", nmp1, \
                        "error=", error)
    
    def __getInitialGuessForThreePointForward(self, nmp1):
        xo=np.zeros(nmp1)
        xo[0]=Const.HALF
        xo[1]=Const.HALF*Const.HALF
        return xo

    def __getInitialGuessForThreePointInverse(self, nmp1, nVan):
        xo=np.zeros(nmp1)
        xo[0]=nVan*Const.HALF+Const.ONE
        xo[1]=-nVan*Const.HALF*Const.HALF
        return xo

    def __getInitialGuess(self, nmp1):
        xo=np.zeros(nmp1)
        #xo[0]=Const.ONE
        xo[:]=Const.ONE/(2*(nmp1-1)+Const.ONE)
        return xo

    def __getInitialGuessFromOld(self, xo, nsNew):
        xoNew=[]
        for i in xo: xoNew.append(i)
        for j in range(nsNew-len(xo)): xoNew.append(Const.ZERO)
        return xoNew

    def __getConstraintMatrixA(self, nmp1):
        k=self._dfilter.getK()
        A=np.ndarray((self._nk, nmp1))
        A[:, 0]=Const.ONE
        for i in range(self._nk):
            for l in range(1, nmp1):
                A[i,l]=Const.TWO*np.cos(l*k[i])
        return A

    def __getConstraintMatrixD(self, nmp1):
        k=self._dfilter.getK()
        D=np.ndarray((self._nk, nmp1))
        D[:,0]=Const.ZERO
        for i in range(self._nk):
            for l in range(1, nmp1):
                D[i,l]=l*np.sin(l*k[i])
        return D

    def __getConstraintConsistency(self, matrixA):
        return LinearConstraint(matrixA[0, :], Const.ONE, Const.ONE)

    def __getConstraintMonotonicPos(self, matrixD):
        return LinearConstraint(matrixD[:, :], Const.ZERO, np.inf)

    def __getConstraintMonotonicNeg(self, matrixD):
        return LinearConstraint(matrixD[:, :], -np.inf, Const.ZERO)

    def __getConstraintZeroToOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ZERO, Const.ONE)

    def __getConstraintZeroDampingAtPi(self, matrixA):
        return LinearConstraint(matrixA[self._nk-1, :], Const.ZERO, \
                Const.ZERO)

    def __getConstraintFixedValueAtPi(self, matrixA, valueAtPi):
        return LinearConstraint(matrixA[self._nk-1, :], valueAtPi, \
                valueAtPi)

    def __getConstraintLinearOneToNPlusOne(self, matrixA):
        k=self._dfilter.getK()
        gLinear=np.zeros(self._nk)
        for i in range(self._nk):
            gLinear[i]=Const.One+k[i]*(self._nVanCittert+1)/Const.PI
        bl=np.ndarray(self._nk)
        br=np.ndarray(self._nk)
        bl[:]=gLinear[i]
        br[:]=gLinear[i] #TODO: possible divide by zero here
        return LinearConstraint(matrixA[0:self._nk, :], bl[0:self._nk], br[0:self._nk])

    def __getConstraintBetweenOneAndNPlusOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ONE, self._nVanCittert+1)

    def __getConstraintBetweenOneAndInfinity(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk, :], Const.ONE, np.inf)
  
    def __getConstraintBetweenOneAndOneOverHGd(self, matrixA):
        hg=self._target.getOriginalTransferFunction()
        bl=np.ndarray(self._nk)
        br=np.ndarray(self._nk)
        bl[:]=Const.ONE
        br[:]=Const.ONE/hg[:] #TODO: possible divide by zero here
        return LinearConstraint(matrixA[1:self._nk, :], bl[1:self._nk], br[1:self._nk])

    def __getConstraintsList(self, nmp1):
        A=self.__getConstraintMatrixA(nmp1)
        D=self.__getConstraintMatrixD(nmp1)
        constraintsList=[]
        constraintsList.append(self.__getConstraintConsistency(A))
        if(self._isInverse):
            constraintsList.append(self.__getConstraintBetweenOneAndNPlusOne(A))
            if(self._isMonotonic): constraintsList.append(self.__getConstraintMonotonicNeg(D))
        else:
            if(self._isMonotonic): constraintsList.append(self.__getConstraintMonotonicPos(D))
            if(self._isPositive): constraintsList.append(self.__getConstraintZeroToOne(A))
            if(self._isZeroDampingAtPi): constraintsList.append(self.__getConstraintZeroDampingAtPi(A))
        
        return constraintsList

  
