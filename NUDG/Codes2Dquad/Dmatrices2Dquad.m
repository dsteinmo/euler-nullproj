function [Dr,Ds] = Dmatrices2Dquad(N,r,s,V)
    [DVr,DVs] = GradVandermonde2Dquad(N,r,s);
    Dr = DVr/V;
    Ds = DVs/V;
    
end