function [ Err ] = ComputeEstimatedErrorPatternFromSyndrome( Syndrome, SyndromeMatrix, ErrorPattern )
%COMPUTEESTIMATEDERRORPATTERNFROMSYNDROME Summary of this function goes here
% Error correction depending on the syndrome
Err = zeros(size(Syndrome,1), size(SyndromeMatrix,1));
for ii = 1 : size(Syndrome,1)
    for jj = 1 : size(SyndromeMatrix, 1)
        if Syndrome(ii, :) == SyndromeMatrix(jj, :)
            Err(ii, :) = ErrorPattern(jj, :);
            break;
        end
    end
end
end

