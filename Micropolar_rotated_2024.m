%Manuscript ID: RSFS-2023-0064.R2 - Analysis of micropolar elastic multilaminated composite and its application to bio-ceramic materials for bone reconstruction.
%Journal: Interface Focus 
%Autores: R. Rodríguez-Ramos, Y. Espinosa-Almeyda, D. Guinovart-Sanjuán, H. Camacho-Montes, P. Rodríguez-Bermúdez, H. Brito-Santana, J. A. Otero, and F. J. Sabina.
%Corresponding author email: yoanh.espinosa@uacj.mx

% This program reproduces the numerical results reported in Section 4.2: Micropolar elastic bi-laminated composites.

% For computations, we consider the following data:
% *Tables 1 and 2, section 4.2
% Material data
% c1122_1 = .65;
% c1212_1 = .35;
% c1221_1 = .35; 
% c1111_1=1;
% Entry values:
% V1=0.1:0.1:0.9
% angles=[30]; For Table 1
% angles=[45]; For Table 2
% R1=rotation_matrix([0,0,0]);
% R2=rotation_matrix([0,0,th2]);
% 
% *Figures 1 and 2, section 4.2
% Material data for Figure 1:
% c1122_1 = .65;
% c1212_1 = .35;
% c1221_1 = .30;
% c1111_1=1;
% Material data for Figure 2:
% c1122_1 = -.40;
% c1212_1 = .70;
% c1221_1 = .60;
% c1111_1=1;
% Entry values for both Figures
% V1=0.1:0.1:0.9;
% angles=0:2.5:90;
% R1=rotation_matrix([0,0,th2]);
% R2=rotation_matrix([0,0,0]);
% 
% *Tables 3 and 4, section 4.2
% Material data for Table 3:
% c1122_1 = .65;
% c1212_1 = .35;
% c1221_1 = .30;
% c1111_1=1;
% Material data for Table 4:
% c1122_1 = -.40;
% c1212_1 = .70;
% c1221_1 = .60;
% c1111_1=1;
% Entry values for both Figures
% V1=0.3, 0.5, 0.7;            % fixed V1 volume fractions
% angles=[45, 50, 60, 75];
% R1=rotation_matrix([0,0,35]);
% R2=rotation_matrix([0,0,th2]);

clear
close all
% Entries:
% For stiffness properties we use the data (in GPa): 
c1122_1 = .65;
c1212_1 = .35;
c1221_1 = .30;
% For torques properties we use the data (in N): 
%  c1122_1 = -.40;
%  c1212_1 = .70;
%  c1221_1 = .60;
c1111_1=1;
%c1111_1 = c1122_1 + c1212_1 + c1221_1;
C_1=cubic(c1122_1, c1212_1, c1221_1,c1111_1);

%for  V1=0.1:0.1:0.9            %  <- different values of V1 volume fraction
for  V1=0.5                     %  <- fixed V1 volume fractions
     V=[V1,1-V1];               %  <- volume fractions of both layers
     i=0;
    %angles=0:2.5:90;           % <- rotation angles for Figures 1 and 2
    angles=[45, 50, 60, 75];    % <- rotation angles for  Tables 3 and 4
    for th2=angles
        i=i+1;
        R1=rotation_matrix([0,0,35]);      % <- Rotation first layer:  Put R1=rotation_matrix([0,0,0])   for Tables 1 and 2, whereas R1=rotation_matrix([0,0,35])  for Tables 3 and 4.
        R2=rotation_matrix([0,0,th2]);     % <- Rotation second layer: Put R1=rotation_matrix([0,0,th2]) for Tables 1 and 2, whereas R1=rotation_matrix([0,0,th2]) for Tables 3 and 4.
        C_ro(:,:,:,:,1)=Rotar_C(C_1,R1);
        C_ro_sin_rotar = tensor2matrix(C_ro(:,:,:,:,1));   % <- To display the matrix of constituent 1
        C_sin_rotar_class=classic(C_ro_sin_rotar);
        C_ro(:,:,:,:,2)=Rotar_C(C_1,R2);
        C_ro_rotada=tensor2matrix(C_ro(:,:,:,:,2));        % <- To display the matrix of constituent 2
        C_ro_rotada_class=classic(C_ro_rotada);
        Ejemplo1C=Sumar_tensor(C_ro,V);
        Ejemplo1C_matrix=round(tensor2matrix(Ejemplo1C),8) % <- To display the effective matrix.

        C_eff{1}(i)=Ejemplo1C_matrix(1,1);
        C_eff{2}(i)=Ejemplo1C_matrix(1,2);
        C_eff{3}(i)=Ejemplo1C_matrix(1,3);
        C_eff{4}(i)=Ejemplo1C_matrix(1,4);
        C_eff{5}(i)=Ejemplo1C_matrix(1,5);
        C_eff{6}(i)=Ejemplo1C_matrix(1,6);
        C_eff{7}(i)=Ejemplo1C_matrix(1,7);
        C_eff{8}(i)=Ejemplo1C_matrix(1,8);
        C_eff{9}(i)=Ejemplo1C_matrix(1,9);
        C_eff{10}(i)=Ejemplo1C_matrix(2,2);
        C_eff{11}(i)=Ejemplo1C_matrix(2,3);
        C_eff{12}(i)=Ejemplo1C_matrix(2,4);
        C_eff{13}(i)=Ejemplo1C_matrix(2,5);
        C_eff{14}(i)=Ejemplo1C_matrix(2,6);
        C_eff{15}(i)=Ejemplo1C_matrix(2,7);
        C_eff{16}(i)=Ejemplo1C_matrix(2,8);
        C_eff{17}(i)=Ejemplo1C_matrix(2,9);
        C_eff{18}(i)=Ejemplo1C_matrix(3,3);
        C_eff{19}(i)=Ejemplo1C_matrix(3,4);
        C_eff{20}(i)=Ejemplo1C_matrix(3,5);
        C_eff{21}(i)=Ejemplo1C_matrix(3,6);
        C_eff{22}(i)=Ejemplo1C_matrix(3,7);
        C_eff{23}(i)=Ejemplo1C_matrix(3,8);
        C_eff{24}(i)=Ejemplo1C_matrix(3,9);
        C_eff{25}(i)=Ejemplo1C_matrix(4,4);
        C_eff{26}(i)=Ejemplo1C_matrix(4,5);
        C_eff{27}(i)=Ejemplo1C_matrix(4,6);
        C_eff{28}(i)=Ejemplo1C_matrix(4,7);
        C_eff{29}(i)=Ejemplo1C_matrix(4,8);
        C_eff{30}(i)=Ejemplo1C_matrix(4,9);
        C_eff{31}(i)=Ejemplo1C_matrix(5,5);
        C_eff{32}(i)=Ejemplo1C_matrix(5,6);
        C_eff{33}(i)=Ejemplo1C_matrix(5,7);
        C_eff{34}(i)=Ejemplo1C_matrix(5,8);
        C_eff{35}(i)=Ejemplo1C_matrix(5,9);
        C_eff{36}(i)=Ejemplo1C_matrix(6,6);
        C_eff{37}(i)=Ejemplo1C_matrix(6,7);
        C_eff{38}(i)=Ejemplo1C_matrix(6,8);
        C_eff{39}(i)=Ejemplo1C_matrix(6,9);
        C_eff{40}(i)=Ejemplo1C_matrix(7,7);
        C_eff{41}(i)=Ejemplo1C_matrix(7,8);
        C_eff{42}(i)=Ejemplo1C_matrix(7,9);
        C_eff{43}(i)=Ejemplo1C_matrix(8,8);
        C_eff{44}(i)=Ejemplo1C_matrix(8,9);
        C_eff{45}(i)=Ejemplo1C_matrix(9,9);
    end
    
% Plots of Figures 1 and 2
% !!!!! In each subplot changes the ylabel('---') in relation to effective stiffness and torque properties
    for k=1
       subplot(3,3,1)
        hold on
        plot(angles,C_eff{k}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{1111} (GPa)')   % for effective stiffness properties
        ylabel('D^*_{1111} (N)')      % for effective torques properties
       subplot(3,3,2)
        hold on
        plot(angles,C_eff{k+1}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{1122} (GPa)')
        ylabel('D^*_{1122} (N)')
       subplot(3,3,3)
        hold on
        plot(angles,C_eff{k+2}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12) 
        %ylabel('C^*_{1133} (GPa)')
        ylabel('D^*_{1133} (N)')
       subplot(3,3,4)
        hold on
        plot(angles,C_eff{k+5}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12) 
        %ylabel('C^*_{1112} (GPa)')
        ylabel('D^*_{1112} (N)')
       subplot(3,3,5)
        hold on
        plot(angles,C_eff{k+17}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{3333} (GPa)')
        ylabel('D^*_{3333} (N)')
       subplot(3,3,6)
        hold on
        plot(angles,C_eff{k+24}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{2323} (GPa)')
        ylabel('D^*_{2323} (N)')
       subplot(3,3,7)
        hold on
        plot(angles,C_eff{k+27}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{2332} (GPa)')
        ylabel('D^*_{2332} (N)')
       subplot(3,3,8)
        hold on
        plot(angles,C_eff{k+35}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{1212} (GPa)')
        ylabel('D^*_{1212} (N)')
       subplot(3,3,9)
        hold on
        plot(angles,C_eff{k+38}, 'LineWidth',3, 'Color', [V1,0,0])
        grid on
        xlabel('angle $\theta$', 'Interpreter','latex','FontSize', 12)
        %ylabel('C^*_{1221} (GPa)')
        ylabel('D^*_{1221} (N)')
    end
end

legend('V_{\it{A}}=0.1', 'V_{\itA}=0.2','V_{\itA}=0.3','V_{\itA}=0.4','V_{\itA}=0.5','V_{\itA}=0.6','V_{\itA}=0.7', 'V_{\itA}=0.8', 'V_{\itA}=0.9', 'orientation', 'horizontal','EdgeColor',[1 1 1], 'Fontsize', 12, 'Position',[0.194973167848331 0.00692177597373724 0.678623708676803 0.0491551446475192])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Don't touch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function create a matrix from tensor  
function [C_matrix] = tensor2matrix(C)
C_matrix(1,1)=C(1,1,1,1);
C_matrix(1,2)=C(1,1,2,2);
C_matrix(1,3)=C(1,1,3,3);
C_matrix(1,4)=C(1,1,2,3);
C_matrix(1,5)=C(1,1,3,1);
C_matrix(1,6)=C(1,1,1,2);
C_matrix(1,7)=C(1,1,3,2);
C_matrix(1,8)=C(1,1,1,3);
C_matrix(1,9)=C(1,1,2,1);

C_matrix(2,1)=C(2,2,1,1);
C_matrix(2,2)=C(2,2,2,2);
C_matrix(2,3)=C(2,2,3,3);
C_matrix(2,4)=C(2,2,2,3);
C_matrix(2,5)=C(2,2,3,1);
C_matrix(2,6)=C(2,2,1,2);
C_matrix(2,7)=C(2,2,3,2);
C_matrix(2,8)=C(2,2,1,3);
C_matrix(2,9)=C(2,2,2,1);

C_matrix(3,1)=C(3,3,1,1);
C_matrix(3,2)=C(3,3,2,2);
C_matrix(3,3)=C(3,3,3,3);
C_matrix(3,4)=C(3,3,2,3);
C_matrix(3,5)=C(3,3,3,1);
C_matrix(3,6)=C(3,3,1,2);
C_matrix(3,7)=C(3,3,3,2);
C_matrix(3,8)=C(3,3,1,3);
C_matrix(3,9)=C(3,3,2,1);

C_matrix(4,1)=C(2,3,1,1);
C_matrix(4,2)=C(2,3,2,2);
C_matrix(4,3)=C(2,3,3,3);
C_matrix(4,4)=C(2,3,2,3);
C_matrix(4,5)=C(2,3,3,1);
C_matrix(4,6)=C(2,3,1,2);
C_matrix(4,7)=C(2,3,3,2);
C_matrix(4,8)=C(2,3,1,3);
C_matrix(4,9)=C(2,3,2,1);

C_matrix(5,1)=C(3,1,1,1);
C_matrix(5,2)=C(3,1,2,2);
C_matrix(5,3)=C(3,1,3,3);
C_matrix(5,4)=C(3,1,2,3);
C_matrix(5,5)=C(3,1,3,1);
C_matrix(5,6)=C(3,1,1,2);
C_matrix(5,7)=C(3,1,3,2);
C_matrix(5,8)=C(3,1,1,3);
C_matrix(5,9)=C(3,1,2,1);

C_matrix(6,1)=C(1,2,1,1);
C_matrix(6,2)=C(1,2,2,2);
C_matrix(6,3)=C(1,2,3,3);
C_matrix(6,4)=C(1,2,2,3);
C_matrix(6,5)=C(1,2,3,1);
C_matrix(6,6)=C(1,2,1,2);
C_matrix(6,7)=C(1,2,3,2);
C_matrix(6,8)=C(1,2,1,3);
C_matrix(6,9)=C(1,2,2,1);

C_matrix(7,1)=C(3,2,1,1);
C_matrix(7,2)=C(3,2,2,2);
C_matrix(7,3)=C(3,2,3,3);
C_matrix(7,4)=C(3,2,2,3);
C_matrix(7,5)=C(3,2,3,1);
C_matrix(7,6)=C(3,2,1,2);
C_matrix(7,7)=C(3,2,3,2);
C_matrix(7,8)=C(3,2,1,3);
C_matrix(7,9)=C(3,2,2,1);

C_matrix(8,1)=C(1,3,1,1);
C_matrix(8,2)=C(1,3,2,2);
C_matrix(8,3)=C(1,3,3,3);
C_matrix(8,4)=C(1,3,2,3);
C_matrix(8,5)=C(1,3,3,1);
C_matrix(8,6)=C(1,3,1,2);
C_matrix(8,7)=C(1,3,3,2);
C_matrix(8,8)=C(1,3,1,3);
C_matrix(8,9)=C(1,3,2,1);

C_matrix(9,1)=C(2,1,1,1);
C_matrix(9,2)=C(2,1,2,2);
C_matrix(9,3)=C(2,1,3,3);
C_matrix(9,4)=C(2,1,2,3);
C_matrix(9,5)=C(2,1,3,1);
C_matrix(9,6)=C(2,1,1,2);
C_matrix(9,7)=C(2,1,3,2);
C_matrix(9,8)=C(2,1,1,3);
C_matrix(9,9)=C(2,1,2,1);
end


%% function defined to rotate the layers
function [C_ro]=Rotar_C(C_no_ro,R)
invR=inv(R);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                C_ro(i,j,k,l)=0;
                for m=1:3
                    for s=1:3
                        for p=1:3
                            for q=1:3
                                % C_ro(i,j,k,l)=C_ro(i,j,k,l)+R(i,m)*R(j,s)*...
                                %     C_no_ro(m,s,p,q)*invR(l,q)*invR(k,p);
                                C_ro(i,j,k,l)=C_ro(i,j,k,l)+R(i,m)*R(j,s)*R(k,p)*R(l,q)*...
                                    C_no_ro(m,s,p,q);
                            end
                        end
                    end
                end
            end
        end
    end
end
end



%% Code for R
function [R]= rotation_matrix(th)
%in grades
R(1,1)=cosd(th(2))*cosd(th(3));
R(1,2)=cosd(th(2))*sind(th(3));
R(1,3)=sind(th(2));
R(2,1)=-sind(th(1))*sind(th(2))*cosd(th(3))-cosd(th(1))*sind(th(3));
R(2,2)=cosd(th(1))*cosd(th(3))-sind(th(1))*sind(th(2))*sind(th(3));
R(2,3)=sind(th(1))*cosd(th(2));
R(3,1)=sind(th(1))*sind(th(3))-cosd(th(1))*sind(th(2))*cosd(th(3));
R(3,2)=-cosd(th(1))*sind(th(2))*sind(th(3))-sind(th(1))*cosd(th(3));
R(3,3)=cosd(th(1))*cosd(th(2));
end

%% function to construct the properties of cubic materials
function [Tensor] = cubic(C1122, C1212, C1221, C1111)
delta=eye(3);
for i=1:3
    for j=1:3
        for m=1:3
            for n=1:3
                Tensor(i,j,m,n)=C1122*delta(i,j)*delta(m,n)+C1212*delta(i,m)*delta(j,n)+...
                   C1221*delta(i,n)*delta(j,m)+(C1111-C1122-C1212-C1221)*(i==j&&j==m&&m==n); 
            end
        end
    end
end
end

%% main functions
function [S_prom] = Sumar_tensor(T,V)
%% Creando la matrix de la inversa del promedio de la inversa de C(l,3,k,3) 
prominvMatrix_1=zeros(3,3);
for M=1:length(V)
    for l=1:3
        for k=1:3
            Matrix_1(l,k,M)=T(l,3,k,3,M);
        end
    end
    invMatrix_1(:,:,M)=inv(Matrix_1(:,:,M));
    prominvMatrix_1=prominvMatrix_1+V(M).*invMatrix_1(:,:,M);
end
invprominvMatrix_1=inv(prominvMatrix_1);

%% Creando la matrix del promedio de C^{-1}(k,3,d,3)*C(d,3,p,q) 
promMatrix_2=zeros(3,3,3);
for M=1:length(V)
    for k=1:3
        for p=1:3
            for q=1:3
                Matrix_2(k,p,q,M)=0;
                for d=1:3
                    Matrix_2(k,p,q,M)=Matrix_2(k,p,q,M)+invMatrix_1(k,d,M)*T(d,3,p,q,M);
                end
            end
        end
    end
    promMatrix_2=promMatrix_2+V(M)*Matrix_2(:,:,:,M);
end
%% Sumando todo lo que esta en el parentisis
SUMA=zeros(3,3,3,M);
for M=1:length(V)
    for l=1:3
        for p=1:3
            for q=1:3
                SUMA(l,p,q,M)=SUMA(l,p,q,M)-T(l,3,p,q,M);
                for k=1:3
                    SUMA(l,p,q,M)=SUMA(l,p,q,M)+invprominvMatrix_1(l,k)*promMatrix_2(k,p,q);
                end
            end
        end
    end
end
%% Calculo del coeficiente efectivo
S_prom=zeros(3,3,3,3);
for M=1:length(V)
    for i=1:3
        for j=1:3
            for p=1:3
                for q=1:3
                    ST(i,j,p,q,M)=T(i,j,p,q,M);
                    for m=1:3
                        for l=1:3
                            ST(i,j,p,q,M)=ST(i,j,p,q,M)+T(i,j,m,3,M)*invMatrix_1(m,l,M)*SUMA(l,p,q,M);
                        end
                    end
                end
            end
        end
    end
    S_prom=ST(:,:,:,:,M)*V(M)+S_prom;
end
end
%%%%%%%%%%%%%
% classico
function [C_class] = classic(C)
C_class=C(1:6,1:6);
end
