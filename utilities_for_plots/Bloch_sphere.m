function Bloch_sphere
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
%
%Plot a sphere to represent the Bloch-sphere.
%--------------------------------------------------------------------------

%Cartesian axes
X0{1}={[0,1],[0,0],[0,0]};  X0{2}={[0,0],[0,1],[0,0]};  X0{3}={[0,0],[0,0],[0,1]};
for ii=1:3
    hold on
    line(X0{ii}{1},X0{ii}{2},X0{ii}{3},'linewidth',2,'color','k','linestyle','-')    
    hold on
end


[X,Y,Z]=sphere(30); h=surf(X,Y,Z,'Facecolor','k');
set(gcf,'color',[0.1,0.1,0.1]); set(gca,'fontsize',30,'fontname','Microsoft Sans Serif');
xlabel('x');  ylabel('y');  zlabel('z');
h.FaceColor='w';
h.EdgeColor ='k';   h.FaceAlpha = 0.3;   h.EdgeAlpha = 0.4;

axis equal

text(1,0,0,'x','color','w','fontsize',20)
text(0,1,0,'y','color','w','fontsize',20)
text(0,0,1,'z','color','w','fontsize',20)

end
