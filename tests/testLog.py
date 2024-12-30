from include import *
from utils import *

log=Log("A")
log.printError("method", "some error")=="Error: A.method some error"
log.printWarning("method", "warning")=="Warning: A.method warning"
