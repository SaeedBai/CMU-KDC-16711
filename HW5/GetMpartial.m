function Matrixpartial = Getpartial(M,variable)
%% This function calculates the partial derivate in the coriolis function
Matrixpartial = gradient(M,variable);
end
