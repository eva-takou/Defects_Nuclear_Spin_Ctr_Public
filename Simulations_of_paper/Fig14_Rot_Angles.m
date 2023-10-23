function Fig14_Rot_Angles
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------


clc ;clearvars; close all;

wL=314;
A=120;
B=90;
s0=0;
s1=-1;
k=1;


figure(1)
Bloch_sphere
Get_Rot_After_free_Evol(wL,A,B,s0,s1,k,'CPMG')
axis off
zlim([-1,1])

figure(2)
Bloch_sphere
Get_Rot_After_free_Evol(wL,A,B,s0,s1,k,'UDD3')
axis off
zlim([-1,1])

figure(3)
Bloch_sphere
Get_Rot_After_free_Evol(wL,A,B,s0,s1,k,'UDD4')
axis off
zlim([-1,1])

end


function Get_Rot_After_free_Evol(wL,AHF,BHF,s0,s1,k,Sequence)

%======== Colors ========================================================

listcolors={'y','m','c',...
    'r','g','b','Asparagus',...
    'Evergreen','Firebrick','Hot pink',...
    'Indigo','Jade','Nutbrown',...
    'Pear','kumQuat','Sky blue','Tawny',...
    'bUrgundy','Violet','aZure'};

valuecolors={[1.00,1.00,0.00], [1.00,0.00,1.00],[0.00,1.00,1.00],...
    [1.00,0.00,0.00],[0.00,1.00,0.00],[0.00,0.00,1.00],[0.42,0.59,0.24],...
    [0.00,0.50,0.00],[0.70,0.13,0.13],[1.00,0.41,0.71],...
    [0.29,0.00,0.51],[0.00,0.66,0.42],[0.50,0.20,0.00],...
    [0.75,0.75,0.00],[1.00,0.50,0.00],[0.00,0.75,0.75],[0.80,0.34,0.00],...
    [0.50,0.00,0.13],[0.75,0.00,0.75],[0.38,0.74,0.99]};

gray_cmap={[1,1,1],[0.7,0.7,0.7],[0.6,0.6,0.6],[0.5,0.5,0.5],[0.4,0.4,0.4],[0.3,0.3,0.3],[0.1,0.1,0.1]};


%==========================================================================


factor=2*pi*1e-3;

Wj  = @(sj) sqrt((wL+sj*AHF)^2+(sj*BHF)^2)*factor;
mxj = @(sj) sj*BHF*factor/Wj(sj);
mzj = @(sj) (wL+sj*AHF)*factor/Wj(sj);
mj  = @(sj) [mxj(sj);0;mzj(sj)];

t   = 4*pi*(2*k-1)/(Wj(s0)+Wj(s1));

%We will label axes with capital letters (A,B,C,...)
%and angles with small letters (a,b,c,...)                                

switch Sequence
    
    
    case 'CPMG'
        
        q=[1,2,1];
        
        for ii=1:length(q)
           
            times(ii) = q(ii)*t/sum(q);
            
        end
        
        %Total rotation: (A)(B)(A)
        
        a = Wj(s0)*times(1);
        b = Wj(s1)*times(2);
        A  = mj(s0);
        B  = mj(s1);
        
        [c,C]=Rodrigues(B,A,b,a);               %BA compose to C
        [d,D]=Rodrigues(A,C,Wj(s0)*times(3),c); %AC compose to D
        
        Axes_Point={A,B,C,D};
        Axes_Names={'A','B','C','D'};
        
        Arcs={{C,B,A},{D,A,C}};
        
    case 'UDD3'
        
        n=3;

        qend= ( sin( pi*(n+1)/(2*n+2))  )^2 - ( sin(pi*(n)/(2*n+2) ) )^2 ;
        
        for ii=1:n
           
            q(ii)=( sin( pi*ii/(2*n+2) ) )^2 - ( sin( pi*(ii-1)/(2*n+2)) )^2;
            
        end
        
        q=1/2*[q(1),q(2),q(3),(qend+q(1)),q(2),q(3),qend];
        
        for ii=1:length(q)
           
            times(ii)=q(ii)*t/sum(q);
            
        end
        
        %Total rotation: (A)(B)(A)(B)(A)(B)(A)
        
        a  = Wj(s0)*times(1);
        b  = Wj(s1)*times(2);
        A  = mj(s0);
        B  = mj(s1);
        
        [c,C]=Rodrigues(B,A,b,a);               %BA compose to C
        [d,D]=Rodrigues(A,C,Wj(s0)*times(3),c); %AC compose to D
        [e,E]=Rodrigues(B,D,Wj(s1)*times(4),d); %BD compose to E
        [f,F]=Rodrigues(A,E,Wj(s0)*times(5),e); %AE compose to F
        [g,G]=Rodrigues(B,F,Wj(s1)*times(6),f); %BF compose to G
        [h,H]=Rodrigues(A,G,Wj(s0)*times(7),g); %AG compose to H
        
        Axes_Point={A,B,C,D,E,F,G,H};
        Axes_Names={'A','B','C','D','E','F','G','H'};
        Arcs={{C,B,A},{D,A,C},{E,B,D},{F,A,E},{G,B,F},{H,A,G}};
        
    case 'UDD4'
        
        n=4;
        qend= ( sin( pi*(n+1)/(2*n+2))  )^2 - ( sin(pi*(n)/(2*n+2) ) )^2 ;

        for ii=1:n
           
            q(ii)=( sin( pi*ii/(2*n+2) ) )^2 - ( sin( pi*(ii-1)/(2*n+2)) )^2;
            
        end        
        
        q=[q,qend];
        
        for ii=1:length(q)
            
           times(ii)=q(ii)*t/sum(q);  
        end
        
        %Total rotation: (A)(B)(A)(B)(A)
        
        a  = Wj(s0)*times(1);
        b  = Wj(s1)*times(2);
        A  = mj(s0);
        B  = mj(s1);
        
        [c,C]=Rodrigues(B,A,b,a);               %BA compose to C
        [d,D]=Rodrigues(A,C,Wj(s0)*times(3),c); %AC compose to D
        [e,E]=Rodrigues(B,D,Wj(s1)*times(4),d); %BD compose to E
        [f,F]=Rodrigues(A,E,Wj(s0)*times(5),e); %AE compose to F
        
        Axes_Point={A,B,C,D,E,F};
        Axes_Names={'A','B','C','D','E','F'};
        Arcs={{C,B,A},{D,A,C},{E,B,D},{F,A,E}};
end



for ii=1:length(Axes_Names)
    

arrowColors{ii}=listcolors{ii};
colors{ii} = valuecolors{ii};
end    

%PLOT THE AXES WITH ARROWHEADS AND ALSO PLOT THE CIRCLES PERP TO THE AXES
for ii=1:length(Axes_Point)
hold on
line([0,Axes_Point{ii}(1)],[0,Axes_Point{ii}(2)],[0,Axes_Point{ii}(3)],'linewidth',4,'color',colors{ii})
text(Axes_Point{ii}(1),Axes_Point{ii}(2),Axes_Point{ii}(3)+0.1,Axes_Names{ii},'fontsize',30,'color',colors{ii})
arrow3([0;0;0].',Axes_Point{ii}.',arrowColors{ii}) %p1,p2,color
[Great_Perp,~] = great_circle(Axes_Point{ii});
plot3(Great_Perp(1,:),Great_Perp(2,:),Great_Perp(3,:),'linewidth',3,...
    'color',colors{ii},'linestyle','--')
fill3(Great_Perp(1,:),Great_Perp(2,:),Great_Perp(3,:),colors{ii},'FaceAlpha',0.2);
hold on
end


for ii=1:length(Arcs)
    
[V1,V2,V3]=great_arcs(Arcs{ii}{1},Arcs{ii}{2},Arcs{ii}{3});    
V={V1,V2,V3};

for ii=1:3
hold on
plot3(V{ii}(1,:),V{ii}(2,:),V{ii}(3,:),'color',gray_cmap{ii},'linewidth',4)
hold on
end
    
end


end
