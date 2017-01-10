function [ t ] = hammingbound( n, k )
%HAMMINGBOUND Summary of this function goes here
% Use of the Hamming bound to determine the maximum correcting capability
% given the dimension of the code and number of redundant bits.
t = 0;          % Inital value of error correcting capability
test = 0;       % To test if bound is still satisfied
while test ~= 1
    % Compute the volume of a Hamming sphere of radius t
    S = 0;
    for i=0:t
        S = S + nchoosek(n,i);
    end
    % Compare with the number of syndromes
    if 2^(n-k) >= S
        t = t+1;
    else
        t=t-1; % The maximum correcting capability is t
        test = 1;
    end
end
end

