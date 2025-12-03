from include import *
from utils import *

log=Log("A")

def testError(capsys):
    log.printError("method", "some error")
    stdOut, stdError=capsys.readouterr()
    assert stdOut=="\nError: A.method some error\n"

def testWarning(capsys):
    log.printError("method", "some warning")
    stdOut, stdError=capsys.readouterr()
    assert stdOut=="\nError: A.method some warning\n"

