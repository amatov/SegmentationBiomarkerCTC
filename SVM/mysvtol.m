function tol = svtol(C)
% tolerance for Support Vector Detection
if C==Inf 
    tol = 1e-5;
else
    tol = C*1e-6;
end

end
