% -------------------------------------------------------------------------
% [ sol, res, nite  ] = JP_secant_solver( fun , x , h , toll , nmax)
% -------------------------------------------------------------------------
%
% © Jacopo Prinetto - Aerospace Science and Technology Dept. - PoliMi - 
% Version:      1.0
% Status:       VALIDATED
% Last Update:  11/01/2019
%
% -------------------------------DESCRIPTION-------------------------------
%
% This function solve nonlinear equations using the secant method
%
% ---------------------------------INPUTS----------------------------------  
%
%   fun  = @(x) equation to be solved.             
%   x    = initial value
%   h    = step for forward finite difference computation.
%   toll = absolute tolerance
%   nmax = maximum number of iterations before stop
%
% ---------------------------------OUTPUTS---------------------------------
%
%   sol  = solution
%   nite = numbers of iterations to find sol
%   res  = residual
%
% ----------------------------------NOTES----------------------------------
%
%   Tangent method. derivatives computed using forward finite differences
%   If the initial point is not the solution and shows a derivatives equal 
%   to zero the initial value is shifted: x = x + 10*toll

function [ sol, res, nite  ] = JP_secant_solver( fun , x , h , toll , nmax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nite = 0;

fun_x = fun(x);
dfun_x = (fun(x+h) - fun_x)/h;

% checking if guess is the solution
if fun_x == 0
    sol = x;
    res = 0;
    return;
end

if abs(dfun_x) == 0
    
    x = x + 10*toll;
    
    fun_x = fun(x);
    dfun_x = (fun(x+h) - fun_x)/h;
    
    % checking if guess is the solution
    if fun_x == 0
        sol = x;
        res = 0;
        return;
    end
    
end

%initializing error
error = toll+1;

while error > toll && abs(fun_x) > toll  && nite < nmax
    
    
    x = x - fun_x/dfun_x;
    
    fun_x = fun(x);
    dfun_x = (fun(x+h) - fun_x)/h;
    
    error = abs(fun_x/dfun_x)
    nite = nite + 1;
    
end

sol  = x;
res = fun_x;

% if nite == nmax
%     fprintf('maximum number of iterations reached. error %f ; residual %f \n' , error, res);
% end


end