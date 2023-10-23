function [Param,V]=great_circle(n)
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%Input: rotation axis components
%
%%Output: the parametrization (x,y,z) coordinates
%
%given a rotation axes,
%draw the great circle perpendicular to it.
%call plot(arg1(1,:),arg1(2,:),arg1(3,:)) to plot the great circle.
%
%


%find v1 perp to n
if n(3)~=0 || n(1)~=0
    
    v1 = [-n(3);0;n(1)];
    
elseif  n(2)~=0
    
    v1 = [-n(2);n(1);0];
end

v2 = cross(n,v1);


%if we give 

p = [0;0;0]; %it should also have the center at [0,0,0]

t = linspace(0,2*pi,1e3);
Param = p+1*cos(t).*v1+1*sin(t).*v2;

%In case we want to plot also the plane: 
%This is based on the equation of the plane n*(r-r_0):

%[X,Y]=meshgrid(-1:0.01:1,-1:0.01:1);
%Z = -1/n(3)*(n(1)*X+n(2)*Y);
%hold on
%surf(X,Y,Z,'facecolor','b','edgecolor','none')


%now draw an arc so that it starts from the rotation axis and terminates
%onto the great circle.

theta = acos(n(3));
phi  = atan2(n(2),n(1));


th_Var = linspace(-pi/2,theta,1e4);
for ii=1:length(th_Var)
    
   Point(ii,:) = [sin(th_Var(ii))*cos(phi);sin(th_Var(ii))*sin(phi);cos(th_Var(ii)) ];
    
   
   
   
   test(ii) = abs(n(1)*Point(ii,1)+n(2)*Point(ii,2)+n(3)*Point(ii,3));
   
    
    
end

[~,indx]=min(test,[],'all','linear'); 

PP = Point(indx,:);

%now we want an arc that passes from PP and from n

V1 = n;
V2 = PP;

D = dot(V1,V2);



b = 1/sqrt(1-D^2);
a=-D*b;

W1 = a*V1+b*V2.';

for ii=1:length(t)
    

V(:,ii) = cos(t(ii))*V1+sin(t(ii))*W1;

   test2(ii) = abs(n(1)*V(1,ii)+n(2)*V(2,ii)+n(3)*V(3,ii));
   
        if test2(ii)<1e-3

            %stop once V belongs in the plane
            break
        end










    
    
    
end




