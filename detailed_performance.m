function compare_and_analyze(G, dt, Reference, parms_ga, parms_zn, J_ga, J_zn)
% ?? ????????? ???????? ?? GA ? ZN ????????????, ????????? ??????? ? ????.

    % ??????? ???????? ?? ??????????? ?? ???????
    [y_ga, t, info_ga, MSE_ga, IAE_ga] = get_metrics(G, dt, parms_ga, Reference);
    [y_zn, ~, info_zn, MSE_zn, IAE_zn] = get_metrics(G, dt, parms_zn, Reference);


    % =======================================================
    % 1. ???????? ?? ??????????? ?? ??????
    % =======================================================
    disp('================================================================');
    fprintf('               Sporedba na performansite (R=%.2f)\n', Reference);
    disp('================================================================');
    
    fprintf('%-35s | %-12s | %-12s\n', 'Kriterium', 'GA ', ' ZN ');
    disp('----------------------------------------------------------------');
    
    % ?????????
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Kp', parms_ga(1), parms_zn(1));
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Ki', parms_ga(2), parms_zn(2));
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Kd', parms_ga(3), parms_zn(3));
    disp('----------------------------------------------------------------');

    % ???????? ?? ???? (J)
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Funkcija na cena (J)', J_ga, J_zn);
    disp('----------------------------------------------------------------');

    % ??????
    fprintf('%-35s | %-12.6f | %-12.6f\n', 'MSE', MSE_ga, MSE_zn);
    fprintf('%-35s | %-12.6f | %-12.6f\n', 'IAE', IAE_ga, IAE_zn);
    
    % ?????
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Vreme na porast (Rise Time) [s]', info_ga.RiseTime, info_zn.RiseTime);
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Vreme na smiruvanje (Settling Time) [s]', info_ga.SettlingTime, info_zn.SettlingTime);
    
    % ???????
    fprintf('%-35s | %-12.2f | %-12.2f\n', 'Preskok (Overshoot) [%]', info_ga.Overshoot, info_zn.Overshoot);
    fprintf('%-35s | %-12.4f | %-12.4f\n', 'Pik ', info_ga.Peak, info_zn.Peak);
    disp('================================================================');


    % =======================================================
    % 2. ?????? ?? ?????? (??? ???)
    % =======================================================
    figure;
    
    % ?????? ?? GA (??????)
    plot(t, y_ga, 'linewidth', 2, 'color', [0 0.7 0], 'DisplayName', 'GA');
    hold on;
    grid on;
    % ?????? ?? ZN (??????)
    plot(t, y_zn, 'linewidth', 2, 'color', [0.8 0 0], 'DisplayName', 'ZN');
    
    % ?????????? ?????? (????)
    plot(t, Reference*ones(size(t)), 'b:', 'DisplayName', 'Referentna vrednost');
    
    hold off;
    title('Sporedba na odzivot');
    xlabel('Vreme [s]');
    ylabel('Agolna brzina [rad/s]');
    legend('show', 'Location', 'southeast');
    grid on;
    % ???????? ? ???? ?? ??? ??? (?????????? ????????)


end


% ??????? ????????? ???????? ?? ????????? ?? ???? ???????
function [y, t, info, MSE, IAE] = get_metrics(G, dt, parms, Reference)
    s = tf('s');
    Kp = parms(1); Ki = parms(2); Kd = parms(3);
    tau_d = 0.001;
    K = Kp + Ki/s + Kd*s/(1 + tau_d*s); 

    Loop = series(K, G);
    ClosedLoop = feedback(Loop, 1);

    t_final = 50; 
    t = 0:dt:t_final;

    [y_unscaled, ~] = step(ClosedLoop, t); 
    y = y_unscaled * Reference;
    E = Reference - y; 
    
    % ???????
    MSE = mean(E.^2);
    IAE = dt * sum(abs(E));
    info = stepinfo(y, t, Reference);
end