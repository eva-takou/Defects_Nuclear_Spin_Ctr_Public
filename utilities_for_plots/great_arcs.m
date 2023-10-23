function [V,VV,VVV]=great_arcs(n,l,m)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%given the points n, l , m,  draw arcs
%from l to m and from l to n and from m to n
%


t = linspace(0,2*pi,1e4);



V1 = l;
V2 = m;

D = dot(V1,V2);


b = 1/sqrt(1-D^2);
a=-D*b;

W1 = a*V1+b*V2;

PerpM  = [-m(3); 0 ; m(1)];

for ii=1:length(t)
    

V(:,ii) = cos(t(ii))*V1+sin(t(ii))*W1;

   test2(ii) = abs(PerpM(1)*V(1,ii)+PerpM(2)*V(2,ii)+PerpM(3)*V(3,ii));
   
        if test2(ii)<1e-3

%i need to stop once V belongs in the plane
            break
        end
end

% now we want to draw the line from l to n:


V1 = l;
V2 = n;

D = dot(V1,V2);


b = 1/sqrt(1-D^2);
a=-D*b;

W1 = a*V1+b*V2;

PerpN  = [-n(3); 0 ; n(1)];


for ii=1:length(t)
    

VV(:,ii) = cos(t(ii))*V1+sin(t(ii))*W1;

   test2(ii) = abs(PerpN(1)*VV(1,ii)+PerpN(2)*VV(2,ii)+PerpN(3)*VV(3,ii));
   
        if test2(ii)<1e-3

                %stop once V belongs in the plane
            break
        end
end

%finally draw the line from m to n


V1 = m;
V2 = n;

D = dot(V1,V2);


b = 1/sqrt(1-D^2);
a=-D*b;

W1 = a*V1+b*V2;

PerpN  = [-n(3); 0 ; n(1)];


for ii=1:length(t)
    

VVV(:,ii) = cos(t(ii))*V1+sin(t(ii))*W1;

   test2(ii) = abs(PerpN(1)*VVV(1,ii)+PerpN(2)*VVV(2,ii)+PerpN(3)*VVV(3,ii));
   
        if test2(ii)<1e-3

           %stop once V belongs in the plane
            break
        end
end



end