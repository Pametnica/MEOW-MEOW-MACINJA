close all
clc
dt=0.01;
PopSize=25;
MaxGenerations=150;
s=tf('s');
Km=0.05; %[Nm/A]
J=0.00025; %[kg*m^2]
B=0.0001; %[Nm/rad/s]
L=1.5e-3; %[mH]
R=0.5; % [Om]
G=Km/(J*L*s*s+(J*R+B*L)*s+(B*R+Km*Km));
Desired_R = 0.07;
options=optimoptions(@ga,'PopulationSize',PopSize,'MaxGenerations',MaxGenerations,'OutputFcn',@myfun,'MutationFcn',{@mutationadaptfeasible, 0.8});
%[x,fval]=ga(@(K)pidtest(G,dt,K),3,-eye(3),zeros(3,1),[],[],[],[],[],options);
[x_ga, fval_ga]=ga(@(K)pidtest(G,dt,K,Desired_R),3,-eye(3),zeros(3,1),[],[],[],[],[],options);
disp('===================================================');
disp('  Rezultati od optimizacija so GA');
disp('===================================================');
fprintf('1. Optimalni PID parametri: (Kp, Ki, Kd):\n');
fprintf('   Kp = %.4f\n', x_ga(1));
fprintf('   Ki = %.4f\n', x_ga(2));
fprintf('   Kd = %.4f\n', x_ga(3));
fprintf('\n');
fprintf('2. Minimalna funkcija na cena od GA(J): %.6f\n', fval_ga);

% ???????? ???????????? ?? ????????? ?? ??????????? ?????????
% ???????? ? ???? ???????? ?? ?? ?? ???????? ???? ??????
% ???? Ku ? Tu ????? ?? ?? ???????? ?? ??????? ?? G(s)
% --- ?? ???????? ??????? (mnogu_sakam_matlab.m) ---

% ... (?? ????? ?? GA) ...

% 2. ??????-?????? (ZN) - ????????? ?? ?????????
disp('==================================================');
disp('           Zigler- Nikols ');
disp('==================================================');
%opts = pidtuneOptions('DesignMethod','zieglerNichols');
%C = pidtune(G, 'pid', opts);
    
 %   Kp_zn = C_zn.Kp;
  %  Ki_zn = C_zn.Ki;
   % Kd_zn = C_zn.Kd;
Ku =  750.3666;
Tu =  6.2832e-04;

    % Ziegler-Nichols PID:
    Kp_zn = 0.6*Ku;
    Ti = Tu/2;
    Td = Tu/8;
    Ki_zn=Kp_zn/Ti;
    Kd_zn=Kp_zn*Td;
    %Kp_zn= min(Kp_zn, 5000);
    %Ki_zn = min(Ki_zn, 5000);
    %Kd_zn = min(Kd_zn, 5000);

    C = pid(Kp_zn, Kp_zn/Ti, Kp_zn*Td);
    x_zn = [Kp_zn, Ki_zn, Kd_zn];
    fval_zn = pidtest(G, dt, x_zn, Desired_R);
    disp('   Parametri so  MATLAB ZN .');
    fprintf('   Kp=%.4f, Ki=%.4f, Kd=%.4f\n', x_zn(1), x_zn(2), x_zn(3));

% 3. ????????



fprintf('\n3. Performansi (za w = 0.7 [rad/s]):\n');
detailed_performance(G, dt, Desired_R, x_ga, x_zn, fval_ga, fval_zn); 
