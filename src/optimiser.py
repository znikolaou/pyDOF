from libs import *
from utils import *
from objectiveFunction import *
from reconstructedTransferFunction import *
from scipy.optimize import Bounds
from scipy.optimize import shgo
#-----------------------------------------------------------------------
class Optimiser(Const, Log):
    """
    Class for the optimiser. 

    Attributes: 
        target (Filter, ReconstructedTransferFunction): the target to optimise against. 
        discreteFilter (DiscreteFilter): the discrete filter object whose coeffs we solve for.
        method (string): optimisation method. Suggested candidates: "SLSQP", "trust-constr".
        lamda (float): regularisation/smoothness parameter >=0.
        optError (float): adaptive optimiser error.
        maxIterations (int) : max optimiser iterations.
        error (string): error method to use: "msqe", "mabse". 
        adaptive (boolean): True, False. If true, the discrete filter stencil size will 
            increase until the target optimisation error self._optError is reached. 
        zeroToOne (boolean): If true, then 0 <= \hat{G}(k) <= 1 for all k in (0, pi).
        positive (boolean): If true, then 0 <= \hat{G}(k) for all k in (0, pi).
        monotoneNeg (boolean): If true, then 0 <= d\hat{G}(k)/dk <= 0 for all k in [0, pi].
        monotonePos (boolean): If true, then 0 <= d\hat{G}(k)/dk >= 0 for all k in [0, pi].
        fixOne (float []): If not none, then \hat{G}(k)= 1 for all k in list.
        fixZero (float []): If not none, then \hat{G}(k)= 0 for all k in list.
        luBounds ([[, , []]]): List of vlower, vupper, and wavenumbers k where  vlower<=g(k)<=vupper. 

    Author: Z. Nikolaou.
    """
    def __init__(self, target, discreteFilter, method="SLSQP", lamda=0.e0, optError=1.0e-6, maxIterations=100, error="msqe", \
            adaptive=True, zeroToOne=False, positive=False, monotoneNeg=False, monotonePos=False, \
            fixOne=None, fixZero=None, luBounds=None):
        Log.__init__(self, type(self).__name__) 
        self._target=target
        self._dfilter=discreteFilter
        self._method=method
        self._lamda=lamda
        self._optError=optError
        self._maxIterations=maxIterations
        self._error=error
        self._adaptive=adaptive
        self._zeroToOne=zeroToOne
        self._positive=positive
        self._monotonePos=monotonePos
        self._monotoneNeg=monotoneNeg 
        self._fixOne=fixOne
        self._fixZero=fixZero
        self._luBounds=luBounds
        #
        self._objFunc=ObjectiveFunction(self._target, self._dfilter, lamda=self._lamda, error=self._error)
        self._nk=self._dfilter.getNk()
        self._k=self._dfilter.getK()
        self._isInverse=not self._dfilter.isForward()
        self._nVanCittert=self.__getNvanCittert()
        self._solution=None
        self._error=None
        self._uncMethod='BFGS'
        self.__printInfo()
 
    def __printInfo(self):
        if(self._isInverse):
            print("\nOptimiser-inverse filter:")
        else:
            print("\nOptimiser-forward filter:")

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
        errorOld=None
        solOld=None
        maxIter=1
        itr=1
        itrDiv=0
        nmp1=self._dfilter.getNmPlusOne()
        x=self.__getInitialGuess(nmp1)      
        nmp1Old=nmp1
        if self._adaptive: maxIter=self._maxIterations

        for i in range(maxIter): 
            self._solution=self.__getConstrainedSolution(x, self.__getConstraintsList(nmp1))
            self._error=self._objFunc.getLastError()
            if self.__isConverged(self._error): 
                isConverged=True
                break
            if self.__isDiverging(itr, self._error, errorOld, itrDiv, nmp1Old, nmp1):
                itrDiv+=1
                isConverged=False
                self._solution=solOld
                self._error=errorOld
                if itrDiv<=Const.ITR_MX_DIV:
                    x=self.__getPertGuess(xOld, nmp1)     
                    solUnc=self.__getUnConstrainedSolution(x)
                    x=solUnc.x
                else:
                    print('\t*Max internal itr, failed for stencil size ns=', nmp1)
                    print('\t*Reverting to last solution-try initialising from a different stencil size.')
                    print('\t*Terminating ...')
                    break
            else:
                self.__printIterInfo(itr, nmp1, self._error) 
                nmp1Old=nmp1
                itrDiv=0
                nmp1+=1
                itr+=1
                x=self.__getInitialGuess(nmp1)
                errorOld=self._error
                solOld=self._solution
                xOld=solOld.x
        
        if isConverged:
            print('\n*Success: converged: itr', itr, 'error', Form.short(self._error), 'target error', Form.short(self._optError))
        else:
            if self._adaptive:
                print('\n*Failure to converge for error=', self._optError)
                print('\tOptimisation error attained=', self._error)

        self._dfilter.setCoeffs(self.get())
        self._dfilter.setOptError(self._error)
               
    def __isConverged(self, error): return error<=self._optError
          
    def __isDiverging(self, itr, error, errorOld, itrDiv, nmp1Old, nmp1):
        if itr>1:
            if error>errorOld*Const.DIVERGENCE_FACTOR:
                print('\t \titr-div', itrDiv, 'increasing error:', nmp1Old,'-',Form.short(errorOld),'-->',nmp1,'-',Form.short(error),'perturbing ...')
                return True
            else:
                return False
        else:
            return False

    def get(self): return self._solution.x
 
    def __printIterInfo(self, itr, nmp1, error):
        print("\tItr", itr, "(nm+1)", nmp1, "error", Form.short(error))
  
    def __getPertGuess(self, x, nsp):
        ns=np.size(x)
        xp=np.zeros(nsp)
        nc=min(ns, nsp)
        for i in range(nc):
            xp[i]=x[i]+Const.PERT_SIGMA_PARAM*x[i]*np.random.randn()
        return xp

    def __getInitialGuess(self, nmp1):
        if self._isInverse:
            xo=self.__getInitialGuessForThreePointInverse(nmp1)
        else:
            xo=self.__getInitialGuessForThreePointForward(nmp1)
            #xo=self.__getInitialGuessForConstantAtUnity(nmp1) 
        #solUnc=self.__getUnConstrainedSolution(xo)
        return xo
        
    def __getInitialGuessFromOld(self, nmp1, xOld):
        nmp1Old=np.size(xOld)
        xo=np.zeros(nmp1)
        if nmp1>=nmp1Old:
            xo[0:nmp1Old]=xOld
        else:
            xo[0:nmp1]=xOld[0:nmp1]
        return xo

    def __getInitialGuessForThreePointForward(self, nmp1):
        xo=np.zeros(nmp1)
        xo[0]=Const.HALF
        xo[1]=Const.HALF*Const.HALF
        return xo

    def __getInitialGuessForThreePointInverse(self, nmp1):
        xo=np.zeros(nmp1)
        nVan=self._target.getNvanCittert()
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
        for i in range(self._nk): A[i,:]=[Const.TWO*np.cos(l*k[i]) for l in range(nmp1)]
        A[:,0]=Const.ONE        
        return A

    def __getConstraintMatrixD(self, nmp1):
        k=self._dfilter.getK()
        D=np.ndarray((self._nk, nmp1))
        for i in range(self._nk): D[i,:]=[l*np.sin(l*k[i]) for l in range(nmp1)]
        D[:,0]=Const.ZERO
        return D

    def __getConstraintBounded(self, matrixA, kb, vlower, vupper):
        i=self._k.searchsorted(kb)
        return LinearConstraint(matrixA[i, :], vlower, vupper)

    def __getConstraintConsistency(self, matrixA):
        return LinearConstraint(matrixA[0, :], Const.ONE, Const.ONE)

    def __getConstraintMonotonicNeg(self, matrixD):
        return LinearConstraint(matrixD[:, :], Const.ZERO, np.inf)

    def __getConstraintMonotonicPos(self, matrixD):
        return LinearConstraint(matrixD[:, :], -np.inf, Const.ZERO)

    def __getConstraintZeroToOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ZERO, Const.ONE)

    def __getConstraintPositive(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ZERO, np.inf)

    def __getConstraintNoDampingAtK(self, matrixA, k):
        i=self._k.searchsorted(k)
        return LinearConstraint(matrixA[i, :], Const.ONE, Const.ONE)

    def __getConstraintFullDampingAtK(self, matrixA, k):
        i=self._k.searchsorted(k)
        return LinearConstraint(matrixA[i, :], Const.ZERO, Const.ZERO)

    def __getConstraintBetweenOneAndNPlusOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ONE, self._nVanCittert+1)
    
    def __getConstraintsList(self, nmp1):
        A=self.__getConstraintMatrixA(nmp1)
        D=self.__getConstraintMatrixD(nmp1)
        constraintsList=[]
        if(self._monotonePos): constraintsList.append(self.__getConstraintMonotonicPos(D))
        if(self._monotoneNeg): constraintsList.append(self.__getConstraintMonotonicNeg(D))
        if(self._positive): constraintsList.append(self.__getConstraintPositive(A))
        if(self._fixOne is not None):
            for k in self._fixOne: constraintsList.append(self.__getConstraintNoDampingAtK(A, k))
        if(self._fixZero is not None):
            for k in self._fixZero: constraintsList.append(self.__getConstraintFullDampingAtK(A, k))
        if(self._luBounds is not None):
            for itm in self._luBounds:
                krange=itm[2]
                vlower=itm[0]
                vupper=itm[1]
                for k in self._k:
                    if k>=krange[0] and k<=krange[1]:
                        constraintsList.append(self.__getConstraintBounded(A, k, vlower, vupper))
        if(self._isInverse):
            constraintsList.append(self.__getConstraintBetweenOneAndNPlusOne(A))
        else:
            if(self._zeroToOne): constraintsList.append(self.__getConstraintZeroToOne(A))

        return constraintsList

  
