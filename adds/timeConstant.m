function Td = timeConstant(Tmax, Tmin, Wi)
%Time Constant: Calculates Td
%   Calculates Td based on Tmax, Tmin, Wi

if abs(Wi) > 1
    fprintf('Error: weights must be -1< Wi <1');
    return;
endif

Td = ones(size(Wi)).* Tmax - abs(Wi) .* (Tmax - Tmin);
endfunction

