# ################################################################################
# prints the caseTest results and returns an index for the specification
# input
# GVAR.o   ... a GVAR object
# sigLevel ... the signficance level at which test statistics are based
# output
# testStat ... a vector indicating the specification case (trend, intercept, quadratic trend, etc)
#
#1. row: case IV (H0) vs. case V (H1)
#2. row: case III vs. case IV
#3. row: case II vs. case III
#4. row: case I vs. case II
#


case.testP=function(GVAR.o,sigLevel=0.05){
    X=GVAR.o$caseTest
    pMatrix=sapply(X,function(x) x[,3])
    testStat=apply(pMatrix,2,function(x) which(x<=sigLevel)[1]) # look at first test which is significant
    r=testStat;
    testStat[which(testStat==1)]="V"
    testStat[which(testStat==2)]="IV"
    testStat[which(testStat==3)]="III"
    testStat[which(testStat==4)]="II"
    testStat[which(testStat==4)]="II"
    testStat[which(testStat==4)]="II"
    testStat[which(is.na(testStat))]="I"
    return(testStat)
}




#case.testP(test)

