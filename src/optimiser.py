from libs import *
from objectiveFunction import *
from reconstructedTransferFunction import *
from scipy.optimize import Bounds
from scipy.optimize import shgo
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
        isFullDampingAtPi (boolean): True, False. If true, then 
            \hat{G}(pi)=0 always. 
        isZeroToOne (boolean): True, False. If true, then 
            0 < \hat{G}(k) < 1 for all k in (0, pi). 
        isAdaptive (boolean): True, False. If true, the discrete filter 
            stencil size will increase until the target optimisation error 
            self._optError is reached.
        lamda (float): regularisation/smoothness parameter.

    Author: Z. Nikolaou.
    """
    def __init__(self, target, discreteFilter, method="SLSQP", lamda=0.e0, optError=1.0e-6, maxIterations=100, \
            isAdaptive=True, \
            isNoDampingAtZero=False, isFullDampingAtPi=False, isZeroToOne=False, \
            isNoDampingAtPi=False, isFullDampingAtZero=False, \
            isMonotonicNegative=False, isMonotonicPositive=False): 
        Log.__init__(self, type(self).__name__) 
        self._target=target
        self._dfilter=discreteFilter
        self._method=method
        self._isNoDampingAtZero=isNoDampingAtZero
        self._isFullDampingAtPi=isFullDampingAtPi
        self._isNoDampingAtPi=isNoDampingAtPi
        self._isFullDampingAtZero=isFullDampingAtZero
        self._isZeroToOne=isZeroToOne
        self._isMonotonicPositive=isMonotonicPositive
        self._isMonotonicNegative=isMonotonicNegative
        self._isAdaptive=isAdaptive
        self._lamda=lamda
        self._optError=optError
        self._maxIterations=maxIterations
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
            print("Optimiser will run for inverse discrete filter:")
        else:
            print('')
            print("Optimiser will run for forward discrete filter:")

    def __getNvanCittert(self):
        nv=None
        if(self._isInverse): nv=self._target.getNvanCittert()
        return nv

    def __getUnConstrainedSolution(self, xo): 
        return minimize(self._objFunc.get, xo, method=self._uncMethod)

    def __getConstrainedSolution(self, xo, constraintsList):
        return minimize(self._objFunc.get, xo, method=self._method, constraints=constraintsList) 

    def run(self): 
        isConverged=False
        maxIter=1
        itr=1
        nmp1=self._dfilter.getNmPlusOne()
        errorOld=None
        solOld=None
        if self._isAdaptive: maxIter=self._maxIterations
        for i in range(maxIter):
            x=self.__getInitialGuessForThreePointForward(nmp1)
            self._solution=self.__getConstrainedSolution(x, self.__getConstraintsList(nmp1))
            self._error=self._objFunc.getLastError()
            self.__printIterInfo(itr, nmp1, self._error) 
            
            if self.__isConverged(self._error): 
                isConverged=True
                break
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
            print('\n*Failure to converge for error=', self._optError)
            print('\tOptimisation error attained=', self._error)

        self._dfilter.setCoeffs(self.get())
        self._dfilter.setOptError(self._error)
               
    def __isConverged(self, error): return error<=self._optError
          
    def __isDiverging(self, itr, error, errorOld):
        if itr>1:
            if error>Const.DIVERGENCE_FACTOR*errorOld:
                    print('*Increasing error detected, error=', error, 'errorOld=', errorOld)
                    print('\tReverting to last solution ...')
                    return True
            else:
                return False
        else:
            return False

    def get(self): return self._solution.x
 
    def __printIterInfo(self, itr, nmp1, error):
        print("\tIteration=", itr, "stencil size (nm+1)=", nmp1, "error=", error)
   
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

    def __getInitialGuessForConstantAtUnity(self, nmp1):
        xo=np.zeros(nmp1)
        xo[0]=Const.ONE
        return xo

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

    def __getConstraintConsistencyAtPi(self, matrixA):
        return LinearConstraint(matrixA[self._nk-1, :], Const.ONE, Const.ONE)

    def __getConstraintMonotonicNeg(self, matrixD):
        return LinearConstraint(matrixD[:, :], Const.ZERO, np.inf)

    def __getConstraintMonotonicPos(self, matrixD):
        return LinearConstraint(matrixD[:, :], -np.inf, Const.ZERO)

    def __getConstraintZeroToOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ZERO, Const.ONE)

    def __getConstraintFullDampingAtPi(self, matrixA):
        return LinearConstraint(matrixA[self._nk-1, :], Const.ZERO, Const.ZERO)

    def __getConstraintFullDampingAtZero(self, matrixA):
        return LinearConstraint(matrixA[0, :], Const.ZERO, Const.ZERO)

    def __getConstraintBetweenOneAndNPlusOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ONE, self._nVanCittert+1)
 
    def __getConstraintsList(self, nmp1):
        A=self.__getConstraintMatrixA(nmp1)
        D=self.__getConstraintMatrixD(nmp1)
        constraintsList=[]
        if(self._isNoDampingAtZero): constraintsList.append(self.__getConstraintConsistency(A))
        if(self._isFullDampingAtZero): constraintsList.append(self.__getConstraintFullDampingAtZero(A))
        if(self._isFullDampingAtPi): constraintsList.append(self.__getConstraintFullDampingAtPi(A))
        if(self._isNoDampingAtPi): constraintsList.append(self.__getConstraintConsistencyAtPi(A))
        if(self._isMonotonicPositive): constraintsList.append(self.__getConstraintMonotonicPos(D))
        if(self._isMonotonicNegative): constraintsList.append(self.__getConstraintMonotonicNeg(D))
        if(self._isInverse): 
            constraintsList.append(self.__getConstraintBetweenOneAndNPlusOne(A))
        else:
            if(self._isZeroToOne): constraintsList.append(self.__getConstraintZeroToOne(A))

        return constraintsList

  
