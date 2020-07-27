function [OS, time, X]=simulateVPwithTimeBreaks(model,l_s,l_r,K,a_rs,a_sr,k,M_death,M_diagnosis,sigma,T_break,T,Cmax,dose,oneDose)
opts_TimetoDeath = odeset('Events', model{2});
% set initial conditions
X_s = (1-sigma)*M_diagnosis;
X_r = sigma*M_diagnosis;
start = 0;

time = [];
X = [];
while true % until death
    % treatment phase
   total_dose = 0;
   while total_dose+oneDose < dose   
        total_dose = total_dose+oneDose;
        C = Cmax;
       
        x = ode45(model{1},[start start+T],[X_s, X_r C],opts_TimetoDeath,l_s,l_r,K,K,a_rs,a_sr,k,M_death);
        x.y(x.y<0) =0;
        X_s = x.y(1,end);
        X_r = x.y(2,end);
        C = x.y(3,end);
        C(C<0) = 0;
        start = start+T;  
        X = [X x.y];
        time = [time x.x];
        if ~isempty(x.xe)
            break;
        end
   end
    if ~isempty( x.xe )
            break
    end
    % drug holiday
    x = ode45(model{1},[start start+T_break],[X_s, X_r C],opts_TimetoDeath,l_s,l_r,K,K,a_rs,a_sr,k,M_death);
    X = [X x.y];
    time = [time x.x];
    X_s = x.y(1,end);
    X_r = x.y(2,end);

        start = start+T_break;  

        if ~isempty( x.xe ) || time(end) > 365*20
            break
        end
end 

OS = time(end)/30;

end