function [gamma,n]=Rodrigues(B,A,b,a)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%given the angle and axes of 2 rotations create their composition
%
%Book: Geometric design linkages
%
%Rotation [C]=[B][A]
%With rotation angles a and b
%and axes A and B

gamma = 2*acos( cos(a/2)*cos(b/2) - sin(a/2)*sin(b/2) *dot(B,A)  );

n = 1/sin(gamma/2)*(...
    sin(b/2)*cos(a/2)*B + sin(a/2)*cos(b/2)*A+...
    sin(b/2)*sin(a/2)*cross(B,A));

end