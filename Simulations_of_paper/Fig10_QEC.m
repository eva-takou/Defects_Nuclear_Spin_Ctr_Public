function Fig10_QEC
%--------------------------------------------------------------------------
%Created by: Eva Takou
%
%Last modified: Oct 23, 2023
%--------------------------------------------------------------------------
clc; close all;


ket0 = [1;0]; ket1=[0;1]; kety = 1/sqrt(2)*[1; 1i]; 



psi_arb = @(th,phi)  cos(th/2)*ket0  + exp(1i*phi/2)*sin(th/2)*ket1;

CaseNum  = {22,13};
ResNum   = {4,4};
PSI0     = {kety,psi_arb(pi/3,0)};
fignum   = {1,3};

init_State_Label  = {'$|y\rangle$','$|\pi/3\rangle$'};


Error_Type = 'X';

for ii=1:length(CaseNum)
    
    
    
    
    plot_coeffs(PSI0{ii},init_State_Label{ii},CaseNum{ii},ResNum{ii},Error_Type,fignum{ii})
    
    
    
    
end

Plot_Bloch_Sphere_Spin(CaseNum{2},ResNum{2})


end

function plot_coeffs(psi0,init_State_Label,CaseNum,ResNum,Error_Type,fignum)

[PSI,Prob,Prob_arb,g,d]=get_coeffs(psi0,CaseNum,ResNum,Error_Type);

rho_F = PSI{5}*PSI{5}';

rho_red = ptrace(rho_F,[2,3],[2,2,2]);

FntSize=20;

if ~iscell(PSI)
    
    error('Expected cell array for input state.')
    
end

Description_Text={'Initial','Encoded','Error','Decoded','Corrected'};


rows=length(PSI);  cols=1;



if CaseNum==22

    ylims   ={[0,1],[-0.4,0.4],[-0.4,0.4],[-0.7,0.5],[-0.5,0.5]};
    yticklab={{'0','0.5','1'},{'-0.4','0','0.4'},...
              {'-0.4','0','0.4'},{'-0.5','0','0.5'},...
              {'-0.5','0','0.5'}};
    yticks  ={[0,0.5,1],[-0.4,0,0.4],[-0.4,0,0.4],...
              [-0.5,0,0.5],[-0.5,0,0.5]};      
    
elseif CaseNum==13
   
    ylims={[0,1],[-0.5,0.5],[-0.5,0.5],[-0.8,0.6],[0,0.4]};
    yticklab={{'0','0.5','1'},{'0.5','0','0.5'},...
              {'-0.5','0','0.5'},{'-0.6','0','0.6'},...
              {'0','0.2','0.4'}};
    yticks={[0,0.5,1],[-0.5,0,0.5],[-0.5,0,0.5],...
            [-0.6,0,0.6],[0,0.2,0.4]};
end

xx = categorical({'|0\rangle_3','|1\rangle_3','|2\rangle_3','|3\rangle_3',...
                  '|4\rangle_3','|5\rangle_3','|6\rangle_3','|7\rangle_3'});        
              
figure(fignum)


for ii=1:length(PSI)
    
    subplot(rows,cols,ii)
    
    if ii~=5
    bar([real(PSI{ii}),imag(PSI{ii})])
    set(gca,'xticklabel',{})
    else
        bar(xx,[real(PSI{ii}),imag(PSI{ii})])
    end
    set(gca,'ytick',yticks{ii},'yticklabel',yticklab{ii})
    if ii==1
        legend({'$\Re(c_j)$','$\Im(c_j)$'},'interpreter','latex','location','best','color','none',...
            'edgecolor','none','NumColumns',2)
        text(1,0.4,init_State_Label,'interpreter','latex','fontsize',FntSize)
        
    elseif ii==5
        
        text(0.5,0.5,strcat(num2str(round(Prob*100,2)),"%"),'interpreter','tex','fontsize',FntSize)
        
    end
    
    ylabel('$c_j$','interpreter','latex')
    fig_defaults(FntSize)
    ylim(ylims{ii})
  
    %text(0.5,0.5,Description_Text{ii},'fontsize',18,'edgecolor','k')

end

figure(fignum+1)
rows=3;
cols=1;

subplot(rows,cols,1)

imagesc(real(rho_red))
colormap(brewermap([],'PuBu'))
cbh=colorbar;
fig_defaults(FntSize)
title('$\Re(\rho_{el})$','interpreter','latex')
set(gca,'xticklabel',{})
set(gca,'yticklabel',{})

subplot(rows,cols,2)

imagesc(imag(rho_red))
set(gca,'xticklabel',{})
set(gca,'yticklabel',{})
colormap(brewermap([],'PuBu'))
cbh=colorbar;
fig_defaults(FntSize)
title('$\Im(\rho_{el})$','interpreter','latex')

subplot(rows,cols,3)


h=surf(g./(2*pi),d./(2*pi),1-Prob_arb.');
h.EdgeColor='none';
xlim([0,0.5])
view(0,90) 
colormap(brewermap([],'PuBu'))  
cbh=colorbar;
fig_defaults(FntSize)
xlabel('$\gamma/(2\pi)$','interpreter','latex')
ylabel('$\delta/(2\pi)$','interpreter','latex')
title('P_{err}')

end

function [PSI,Prob,Prob_arb,g,d]=get_coeffs(psi0_el,CaseNum,ResNum,Error_Type)
%Output: 

% PSI: State at each part of QEC circuit (cell array 1x5)
% Prob: Probability that we recovered psi0_el electron's initial state
% Prob_arb: Probability that QEC succeeds for arbitrary initial states
% g/d: parameters of the arbitrary state cos(g/2)|0>+exp(i*d)*sin(g/2)|1>

[Three_Partite_Gates,~,~,]=Gates_QEC_circuit(CaseNum,ResNum,Error_Type); %Get the gates of QEC circuit


%======== Defs ===========================================================

I    = eye(2); 
I4   = eye(4); 
ket0 = [1;0]; 
ket1 = [0;1]; 

kron3=@(a,b,c) kron(a,kron(b,c));

psi_arb=@(g,d)  cos(g/2)*ket0 + exp(1i*d)*sin(g/2)*ket1;

Probability = @(psi0_EL,Psi_F) norm(kron(psi0_EL,I4)'*Psi_F)^2;


%======= Get state at each step of QEC and find success Prob =============
psi0 = kron3(psi0_el,ket1,ket1);

tempGate=1;
PSI{1}=psi0;
for ii=1:length(Three_Partite_Gates)
   
    tempGate = Three_Partite_Gates{ii}*tempGate;
    PSI{ii+1}  = tempGate*psi0;
    
end

Prob=Probability(psi0_el,PSI{end});

%========= Get Success Prob for arbitrary initial states ================

All_QEC_Gates = 1;
for ii = 1:length(Three_Partite_Gates)
    All_QEC_Gates=Three_Partite_Gates{ii}*All_QEC_Gates;
end

g=linspace(0,pi,1e3); d=linspace(0,2*pi,1e3);

numd=length(d);
parfor ii=1:length(g)
    for jj=1:numd
        
      Psi =  All_QEC_Gates * kron3(psi_arb(g(ii),d(jj)),ket1,ket1);
      Prob_arb(ii,jj)=Probability(psi_arb(g(ii),d(jj)),Psi);
        
    end
end




end


function [Three_Partite_Gates,Gates_Elec_0,Gates_Elec_1]=Gates_QEC_circuit(CaseNum,ResNum,Error_Type)

I = eye(2); X= [0 1; 1 0]; Y=[0 -1i; 1i 0]; Z=[1 0; 0 -1];
s00 = [1 0 ; 0 0 ]; s11=[ 0 0 ; 0 1];      

Rn=@(a,n) expm(-1i*a/2*(n(1)*X+n(2)*Y+n(3)*Z));
Ry=@(a)   expm(-1i*a/2*Y); Rx=@(a)   expm(-1i*a/2*X);

kron3=@(a,b,c) kron(a,kron(b,c));

B0 = 403; %in Gauss
B0 = B0*1e-4;

gamma_C13 = 6.728284 * 1e7; %rad/T*s
wL        = gamma_C13 * B0 ; 
wL        = wL*1e-6; %2.7 MHz
wL        = wL*1e3; %in KHz
wL        = wL/(2*pi);  %431.5484 kHz

s0 = 0; s1=-1;


load('27SpinsData.mat')

T = Out{CaseNum,ResNum}.t;
N = Out{CaseNum,ResNum}.N;

A = Out{CaseNum,ResNum}.Atarget;
B = Out{CaseNum,ResNum}.Btarget;

Nnuc=1; k=1;

for ii=1:length(A)
    
   Spin     = SubClass_U4Operations(wL,A(ii),B(ii),s0,s1,Nnuc,k,N);
   Spin     = Spin.CPMG(T);
   Spin     = Spin.Rot_Angles_And_Axes;
   phi0(ii) = Spin.angles{1};
   phi1(ii) = Spin.angles{2}; 
   n0(ii,:) = Spin.axes{1};
   n1(ii,:) = Spin.axes{2};
   
end

%============= Construct the multipartite gate ============================
U00 = 1; U11=1;

for ii=1:length(A)
    
   U00 = kron(U00,Rn(phi0(ii),n0(ii,:))) ;
   U11 = kron(U11,Rn(phi1(ii),n1(ii,:))) ;
   
end

CR = kron(s00,U00)+kron(s11,U11);


switch Error_Type
    
    case 'X'
Three_Partite_Gates=...
    {kron3(I,Ry(-pi),Ry(-pi))*CR, kron3(X,I,I),  CR,   Toffoli   };               
         %encoding                  %error     %decoding %correction
Gates_Elec_0.Encoding   = {Ry(-pi)*Rn(phi0(1),n0(1,:)),Ry(-pi)*Rn(phi0(2),n0(2,:))};
Gates_Elec_0.Decoding   = {Rn(phi1(1),n1(1,:)),Rn(phi1(2),n1(2,:))}; %electron was flipped
Gates_Elec_1.Encoding   = {Ry(-pi)*Rn(phi1(1),n1(1,:)),Ry(-pi)*Rn(phi1(2),n1(2,:))};
Gates_Elec_1.Decoding   = {Rn(phi0(1),n0(1,:)),Rn(phi0(2),n0(2,:))}; %electron was flipped
         
    case 'I'
Three_Partite_Gates=...
    {kron3(I,Ry(-pi),Ry(-pi))*CR, kron3(I,I,I),  CR,   Toffoli   };               
         %encoding                  %error     %decoding %correction
Gates_Elec_0.Encoding   = {Ry(-pi)*Rn(phi0(1),n0(1,:)),Ry(-pi)*Rn(phi0(2),n0(2,:))};
Gates_Elec_0.Decoding   = {Rn(phi0(1),n0(1,:)),Rn(phi0(2),n0(2,:))}; %electron was flipped
Gates_Elec_1.Encoding   = {Ry(-pi)*Rn(phi1(1),n1(1,:)),Ry(-pi)*Rn(phi1(2),n1(2,:))};
Gates_Elec_1.Decoding   = {Rn(phi1(1),n1(1,:)),Rn(phi1(2),n1(2,:))}; %electron was flipped
        
end
                     



end



function Ucorr=Toffoli

th1=pi/2;
th2=th1;
th3=th1-3*pi/2;
th4=th3+pi;

I = eye(2); X= [0 1; 1 0]; Y=[0 -1i ; 1i 0]; Z=[1 0; 0 -1]; 
s00=[1 0; 0 0 ]; s11=[0 0;0 1];
 

Ry=@(a) expm(-1i*a/2*Y); Rx=@(a) expm(-1i*a/2*X);

kron3=@(a,b,c) kron(a,kron(b,c));

CRx12 = kron3(s00,Rx(pi/2),I) + kron3(s11,Rx(-pi/2),I);
CRx13 = kron3(s00,I,Rx(pi/2)) + kron3(s11,I,Rx(-pi/2));


Ucorr  = kron3(Rx(th4),I,I ) ...
       * CRx13 * kron3(Ry(th3),I,I) ...
       * CRx12 * kron3(Rx(th2),I,I) ...
       * CRx13 * kron3(Ry(th1),I,I)  ...
       * CRx12 * kron3(I,Ry(-pi/2),Ry(-pi/2));




end


function out=rho_BS(ket,oper)

ket = oper*ket;
rho = ket*ket';

u=2*real(rho(1,2));
v=2*imag(rho(2,1));
w=rho(1,1)-rho(2,2);

out=[u;v;w];



end



function Plot_Bloch_Sphere_Spin(CaseNum,ResNum)

gray_cmap={[1,1,1],[0.7,0.7,0.7],[0.6,0.6,0.6],[0.5,0.5,0.5],[0.4,0.4,0.4],[0.3,0.3,0.3],[0.1,0.1,0.1]};
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

%==========================================================================


load('27SpinsData.mat')
NuclearIndices=Out{CaseNum,ResNum}.TargetNuc;

init_state = [0;1];  O=[0;0;0]; 

Error_Type={'X','I'};

fignum={{5,6},{7,8}};
titles={'Error','NoError'};


for jj=1:length(Error_Type)
    
    [~,Gates_Elec_0,Gates_Elec_1]=Gates_QEC_circuit(CaseNum,ResNum,Error_Type{jj});
    
    for ii=1:2 %2 nuclear spins
      
      init_state_BS      = rho_BS(init_state,eye(2));  
      encoded_state_BS_0 = rho_BS(init_state,Gates_Elec_0.Encoding{ii});
      decoded_state_BS_0 = rho_BS(init_state,Gates_Elec_0.Decoding{ii}*Gates_Elec_0.Encoding{ii});
      encoded_state_BS_1 = rho_BS(init_state,Gates_Elec_1.Encoding{ii});
      decoded_state_BS_1 = rho_BS(init_state,Gates_Elec_1.Decoding{ii}*Gates_Elec_1.Encoding{ii});

      
      figure(fignum{jj}{ii})
      title(titles{jj},'color','w')
      Bloch_sphere;
      
      
      ax = gca; 
      ax.Clipping = 'off';
      axis off
      text(1,1.2,1,strcat('Nucleus C',num2str(NuclearIndices(ii))),'fontsize',30,'color','w')
      %======= Electron in |0> ============================================
      arrow3(O.',init_state_BS.',listcolors{6})
      arrow3(O.',encoded_state_BS_0.',listcolors{5})
      arrow3(O.',decoded_state_BS_0.',listcolors{5})

      line([O(1),init_state_BS(1)],[O(2),init_state_BS(2)],[O(3),init_state_BS(3)],'color',listcolors{6},'linewidth',2)
      line([O(1),encoded_state_BS_0(1)],[O(2),encoded_state_BS_0(2)],[O(3),encoded_state_BS_0(3)],'color',listcolors{5},'linewidth',2)
      line([O(1),decoded_state_BS_0(1)],[O(2),decoded_state_BS_0(2)],[O(3),decoded_state_BS_0(3)],'color',listcolors{5},'linewidth',2)
      Points={encoded_state_BS_0,decoded_state_BS_0};
      
      for ll=1:length(Points)

           text(Points{ll}(1),Points{ll}(2),Points{ll}(3),num2str(ll),'fontsize',20,'color','g') ;
      end

      [V1,V2,V3]=great_arcs([0;0;-1],encoded_state_BS_0,encoded_state_BS_0);
      V={V1,V2,V3};
     
      for ll=1:3
            hold on
            plot3(V{ll}(1,:),V{ll}(2,:),V{ll}(3,:),'color','g','linewidth',3)
            hold on
      end                    
        [V1,V2,V3]=great_arcs(encoded_state_BS_0,decoded_state_BS_0,decoded_state_BS_0);
        V={V1,V2,V3};
      for ll=1:3
            hold on
            plot3(V{ll}(1,:),V{ll}(2,:),V{ll}(3,:),'color','g','linewidth',3)
            hold on
      end            
      
      %======== Electron starts in |1> ====================================
      
      %arrow3(O.',init_state_BS.',listcolors{6})
      arrow3(O.',encoded_state_BS_1.',listcolors{4})
      arrow3(O.',decoded_state_BS_1.',listcolors{4})

      %line([O(1),init_state_BS(1)],[O(2),init_state_BS(2)],[O(3),init_state_BS(3)],'color',listcolors{6},'linewidth',2)
      line([O(1),encoded_state_BS_1(1)],[O(2),encoded_state_BS_1(2)],[O(3),encoded_state_BS_1(3)],'color',listcolors{4},'linewidth',2)
      line([O(1),decoded_state_BS_1(1)],[O(2),decoded_state_BS_1(2)],[O(3),decoded_state_BS_1(3)],'color',listcolors{4},'linewidth',2)
      Points={encoded_state_BS_1,decoded_state_BS_1};
      
      for ll=1:length(Points)

           text(Points{ll}(1),Points{ll}(2),Points{ll}(3),num2str(ll),'fontsize',20,'color',listcolors{4}) ;
      end

      [V1,V2,V3]=great_arcs([0;0;-1],encoded_state_BS_1,encoded_state_BS_1);
      V={V1,V2,V3};
     
      for ll=1:3
            hold on
            plot3(V{ll}(1,:),V{ll}(2,:),V{ll}(3,:),'color',listcolors{4},'linewidth',3)
            hold on
      end                    
        [V1,V2,V3]=great_arcs(encoded_state_BS_1,decoded_state_BS_1,decoded_state_BS_1);
        V={V1,V2,V3};
      for ll=1:3
            hold on
            plot3(V{ll}(1,:),V{ll}(2,:),V{ll}(3,:),'color',listcolors{4},'linewidth',3)
            hold on
      end            
      
      
      
    end
    
    
    
    
    
end








end




