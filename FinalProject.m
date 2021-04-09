%{
Comparing Polynomial Interpolate Approximations Using the Lagrange
Polynomial and Neville's Algorithm
Numerical Methods Final Project
MATH*2130
Professor M. Demers
Chris Wilson
0784966

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                                                         %
%                          FinalProject.m                                 %
%   This Script compares the result from approximating the value of the   %
%   LaGrange Polynomial at a target value, as well as approximating using % 
%   Neville's Algorithm as well as their runtimes                         %
%                                                                         %
% List of Functions                                                       %
% FinalProject()                      -> Main Function                    %
% initCleanup()                       -> function to clear workspace and  %
%                                        console                          %
% title()                             -> function to print title & points %
% printResults(X,Y)                   -> function to print the various    %
%                                        approximations                   %
% neville(X,Y,val)                    -> function that returns matrix of  %
%                                        Neville approximations           %                          
% laGrange(X,Y)                       -> function to determine the        %
%                                        LaGrange Polynomial              %
% evalLagrange(X, Y,val)              -> function to evaluate LaGrange    %
%                                        Polynomial at the target value   %
% addPointNeville(x,y,nevilleMat,val) -> function that adds a new point   %
%                                         to the matrix returned by the   %
%                                         neville function and            %
%                                         recalculate the approximation   % 
% addPointLaGrange(X,Y,x,y,val)       -> function that computes LaGrange, %
%                                        then adds the new point to the   %
%                                        X and Y vectors,and reevaluates  %
%                                        the LaGrange Polynomial at the   %
%                                        target                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
FinalProject -> Main Function; set points to interpolate here

Calls -> printResults(X,Y)
%}
function FinalProject()
    initCleanup();
    X = [8.1 8.3 8.6 8.7];
    Y = [16.9446 17.56492 18.50515 18.82091];
    
    printResults(X,Y)
end

%{
initCleanup -> Function to clear workspace and console; resets format to
               default
%}
function initCleanup()
    format;
    disp('clearing workspace and console')
    pause(1);
    clear;clc;
end
%{
title -> function to print title and points to interpolate
%}
function title(X,Y)
  disp("Comparing Lagrange Polynomial & Neville's Algorithm")
  disp('---------------------------------------------------')
  disp('points:')
  for i = 1: length(X)
    fprintf("(%0.1f,%f)\n", X(i), Y(i))
  end  
end
%{
printResults -> Main Wrapper Function to call and print results of other
                functions

param X      -> 1 x n matrix of x values for a set of correspoding points
param Y      -> 1 x n matrix of y values for a set of correspoding points
 
Calls        -> title
             -> neville
             -> laGrange
             -> addPointNeville
             -> addPointLaGrange
%}
function printResults(X,Y)
    title(X,Y);
    
    target = 8.4;
    fprintf('\ntargeted x-value: %0.2f\n', target)
    
    fprintf("\nApproximating target using Neville's Algorithm\n\n")
    tic
    approxNeville = neville(X,Y,target);
    toc
    fprintf("Neville Approximation using all points: %0.5f\n\n"...
        , approxNeville(length(X), length(X)+1))
    
    fprintf("\nApproximating target using LaGrange Polynomial\n\n")
    tic
    approxLaGrange = evalLaGrange(X,Y,target);
    toc
    fprintf("Approximation using LaGrange Polynomial of all points: %0.5f\n\n"...
        , approxLaGrange)
    
    fprintf("\nAdding point (8.9,19.12345) and approximating target using Neville's Algorithm\n\n")
    tic
    newNevilleApprox = addPointNeville(8.9, 19.12345, approxNeville, target);
    toc
    fprintf("Neville Approximation of target after adding new point: %0.5f\n\n"...
        , newNevilleApprox)
    
    fprintf("\nAdding point (8.9,19.12345) and approximating target using LaGrange Polynomial\n\n")
    tic
    newLaGrangeApprox = addPointLaGrange(X,Y, 8.9, 19.12345, target);
    toc
    fprintf("Approximation using Lagrange Polynomial after adding new point: %0.5f\n\n"...
        ,newLaGrangeApprox)
end

%{
neville            -> function that calculates approximation of 
                      interpolate at val

param X            -> 1 x n matrix of x values for a set of correspoding points
param Y            -> 1 x n matrix of y values for a set of correspoding points
param val          -> value at which to evaluate interpolate at

returns nevilleMat -> Approximation of Interpolate at x = val
%}
function nevilleMat = neville(X, Y, val)
    n = length(X);
    A = zeros(n,n+1);         % matrix to store approximations
    for i = 1:n
        A(i,1) = X(i);
        A(i,2) = Y(i);
    end
    
    for i = 3:n + 1         % start at 3 since col 1 & 2 are x & y
        for j = i - 1:n
            A(j,i) = (((val - X(j - i + 2)) * A(j, i-1))...
                - ((val - X(j)) * A(j-1,i-1)))...
                / (X(j) - X(j - i + 2));    
        end
    end
    nevilleMat = A;
end

%{
laGrange   -> function that calculates the LaGrange Polynomial; taken from
              lab 6

param X    -> 1 x n matrix of x values for a set of correspoding points
param Y    -> 1 x n matrix of y values for a set of correspoding points

returns    -> LaGrange Polynomial
%}
function laGrangePolynomial = laGrange(X,Y)
    NumPoints = length(X);
    syms x;
    syms L(x); 
    syms currentBasisFn(x);

    L(x) = 0;
    currentBasisFn(x) = 1;

    for i = 1:NumPoints
        for k = 1:NumPoints
            if i~= k
                currentBasisFn(x) = currentBasisFn(x)*(x - X(k))...
                    /(X(i) - X(k));
            end
        end
        L(x) = L(x) + Y(i)*currentBasisFn(x); 
        currentBasisFn(x) = 1;
    end
    laGrangePolynomial = L(x);
end


%{
evalLaGrange  -> function that approximates the LaGrange Polynomial at a
                 given point

param X       -> 1 x n matrix of x values for a set of correspoding points
param Y       -> 1 x n matrix of y values for a set of correspoding points
param val     -> value at which to evaluate interpolate at

returns       -> approximation of LaGrange at val
%}
function evaluatedLaGrangePoly = evalLaGrange(X,Y, val)
    syms x;
    syms L(x);    

    L(x) = laGrange(X,Y);
    evaluatedLaGrangePoly = L(val);
    
end


%{
addPointNeville          -> function to append new row of neville approximations to
                            neville matrix; uses modified neville algo 
                            for even faster results

param x                  -> x-value of new point
param y                  -> y-value of new point
param nevilleMat         -> Matrix of Neville Algo. Approximations
param val                -> value at which to evaluate interpolate at

returns newNevilleApprox -> new Neville Algo approximation of interpolate at val
%}
function newNevilleApprox = addPointNeville(x,y, nevilleMat, val)
   nevilleSize = size(nevilleMat);
   nevilleMat = [nevilleMat zeros(nevilleSize(1),1);zeros(1, nevilleSize(2)+1)];
   nevilleMat(end,1) = x;
   nevilleMat(end,2) = y;
   
   for i = 1:nevilleSize(2) - 1
                           
       nevilleMat(end, i + 2) = nevilleMat(end - 1, i + 1) ...
           + ((((val - nevilleMat(end - i,1))) ...
           * ( nevilleMat( end, i + 1 ) ...
           - nevilleMat( end - 1, i + 1 )))  ...
           / ( nevilleMat(end,1) - nevilleMat(end - i,1) )); 
   end
   
   newNevilleApprox = nevilleMat(end,end);
end

%{
addPointLaGrange          -> function to calculate LaGrange Poly, then add new
                             point and recalculate poly and its approximation

param X                   -> 1 x n matrix of x values for a set of correspoding points
param Y                   -> 1 x n matrix of y values for a set of correspoding points
param x                   -> x-value of new point
param y                   -> y-value of new point
param val                 -> value at which to evaluate interpolate at

returns newLaGrangeApprox -> Approx. of new LaGrange Poly at val  
%}
function newLagrangeApprox = addPointLaGrange(X,Y,x,y, val)
    evalLaGrange(X, Y, val);
    
    X(end + 1) = x;
    Y(end + 1) = y;
    
    newLagrangeApprox = evalLaGrange(X, Y, val);
end

%{
Expected Output

Comparing Lagrange Polynomial & Neville's Algorithm
---------------------------------------------------
points:
(8.1,16.944600)
(8.3,17.564920)
(8.6,18.505150)
(8.7,18.820910)

targeted x-value: 8.40

Approximating target using Neville's Algorithm

Elapsed time is 0.004252 seconds.
Neville Approximation using all points: 17.87709


Approximating target using LaGrange Polynomial

Elapsed time is 1.802849 seconds.
Approximation using LaGrange Polynomial of all points: 17.87709


Adding point (8.9,19.12345) and approximating target using Neville's Algorithm

Elapsed time is 0.001234 seconds.
Neville Approximation of target after adding new point: 17.85633


Adding point (8.9,19.12345) and approximating target using LaGrange Polynomial

Elapsed time is 0.938746 seconds.
Approximation using Lagrange Polynomial after adding new point: 17.85633
%}
