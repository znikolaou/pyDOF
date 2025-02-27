from libs import *
from objectiveFunction import *
from reconstructedTransferFunction import *
#-----------------------------------------------------------------------
class Optimiser(Const, Log):
    """
    Class for the optimiser. 

    Attributes: 
        target: Target, ReconstructedTransferFunction 
            The target to optimise against. 
        discreteFilter: DiscreteFilter 
            The discrete filter object whose coeffs we solve for.
        method: "SLSQP", "trust-constr", ... 
            The optimiser method to use for constrained optimisation. 
            Suggested candidates are  "SLSQP", and "trust-constr". 
        isZeroDampingAtPi: True, False 
            If true then \hat{G}(pi)=0 always. 
        isPositive: True, False 
            If true then \hat{G}(k) > 0 for all k in [0, pi). 
        isAdaptive: True, False 
            If true, the discrete filter stencil size will increase 
            until the target optimisation error Const.TOL_OPTIMISER 
            is reached.

    Author: Z. Nikolaou.
    """
    def __init__(self, target, discreteFilter, method="SLSQP",  \
            isZeroDampingAtPi=True, isPositive=True, isAdaptive=True): 
        Log.__init__(self, type(self).__name__) 
        self._target=target
        self._dfilter=discreteFilter
        self._method=method
        self._isZeroDampingAtPi=isZeroDampingAtPi
        self._isPositive=isPositive
        self._isAdaptive=isAdaptive
        self._objFunc=ObjectiveFunction(self._target, self._dfilter)
        self._nk=self._dfilter.getNk()
        self._isInverse=self.__isInverse(self._target)
        self._nVanCittert=self.__getNvanCittert()
        self._solution=None
        self._error=None
 
    def __isInverse(self, targetOb):
        if(isinstance(targetOb, ReconstructedTransferFunction)):
            print('')
            print("Optimiser will run for inverse discrete filter")
            print('')
            return True
        else:
            print('')
            print("Optimiser will run for forward discrete filter")
            print('')
            return False

    def __getNvanCittert(self):
        nv=None
        if(self._isInverse): nv=self._target.getNvanCittert()
        return nv

    def __getUnConstrainedSolution(self, xo): 
        solution=minimize(self._objFunc.get, xo, method=self._method)
        return solution

    def __getConstrainedSolution(self, xo, constraintsList):
        solution=minimize(self._objFunc.get, xo, method=self._method, \
                constraints=constraintsList)
        return solution

    def run(self):
        itr=1
        nmp1=self._dfilter.getNmPlusOne()
        xo=self.__getInitialGuess(nmp1)

        #get unconstrained solution
        #solUnc=self.__getUnConstrainedSolution(xo)
        #xoUnc=solUnc.x
        self._solution=self.__getConstrainedSolution(xo, \
                self.__getConstraintsList(nmp1))
        error=self._objFunc.getLastError()
        self.__printIterInfo(itr, nmp1, error)
        self._error=error
        if self._isAdaptive: 
            self.__runAdaptive(itr, nmp1, self.get(), error)

        #set dfilter coeffs and optimisation error
        self._dfilter.setCoeffs(self.get())
        self._dfilter.setOptError(self._error)

    def __runAdaptive(self, itro, nmp1o, xo, error):
        errorNew=error
        itr=itro
        nmp1=nmp1o
        while(errorNew > Const.TOL_OPTIMISER):
            itr+=1
            if(itr > Const.MAX_ITERATIONS):
                self.printError("_runAdaptive", "itr>MAX_ITERATIONS")
                exit()
            nmp1+=1
            errorOld=errorNew
            solutionOld=self._solution
            #TODO: some issue here! 
            #xoNew=self.__getInitialGuessFromOld(xo, nmp1)
            xoNew=self.__getInitialGuess(nmp1)
            #
            #
            #solUnc=self.__getUnConstrainedSolution(xoNew)
            #xoNew=solUnc.x
            #run optimiser with new stencil
            self._solution= \
                    self.__getConstrainedSolution(xoNew, \
                    self.__getConstraintsList(nmp1))
     
            errorNew=self._objFunc.getLastError()
            self._error=errorNew
            #xo=xoNew
            self.__printIterInfo(itr, nmp1, errorNew)
        
            if(errorNew>1.05*errorOld): 
                msg="increasing error detected!"+ \
                        "\n try using a different initial stencil"+ \
                        "\n reverting to old solution for now ..."
                self.printWarning("__runAdaptive", msg)
                print("errorNew=", errorNew, "errorOld=", errorOld)
                self._solution=solutionOld
                errorNew=errorOld
                self._error=errorNew
                break
            

    def get(self): return self._solution.x

    def getSolution(self): return self._objFunc.getLastSolution()

    def __printIterInfo(self, itr, nmp1, error):
        print("Iteration=", itr, "stencil size (nm+1)=", nmp1, \
                        "error=", error)
            
    def __getInitialGuess(self, nmp1):
        xo=np.zeros(nmp1)
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

    def __getConstraintConsistency(self, matrixA):
        return LinearConstraint(matrixA[0, :], Const.ONE, Const.ONE)

    def __getConstraintZeroToOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk-1, :], Const.ZERO, Const.ONE)

    def __getConstraintZeroDampingAtPi(self, matrixA):
        return LinearConstraint(matrixA[self._nk-1, :], Const.ZERO, \
                Const.ZERO)

    def __getConstraintBetweenOneAndNPlusOne(self, matrixA):
        return LinearConstraint(matrixA[1:self._nk, :], Const.ONE, \
                    self._nVanCittert+1)
    
    def __getConstraintBetweenOneAndOneOverHGd(self, matrixA):
        hg=self._target.getOriginalTransferFunction()
        bl=np.ndarray(self._nk)
        br=np.ndarray(self._nk)
        bl[:]=Const.ONE
        br[:]=Const.ONE/hg[:]
        return LinearConstraint(matrixA[1:, :], bl[1:], br[1:])

    def __getConstraintsList(self, nmp1):
        A=self.__getConstraintMatrixA(nmp1)
        constraintsList=[]
        constraintsList.append(self.__getConstraintConsistency(A))
        if(self._isInverse):
            constraintsList.append(self.__getConstraintBetweenOneAndOneOverHGd(A))
            return constraintsList
        if(self._isPositive): 
            constraintsList.append(self.__getConstraintZeroToOne(A))
        if(self._isZeroDampingAtPi): 
            constraintsList.append(self.__getConstraintZeroDampingAtPi(A))
        return constraintsList

  
