function [OS, response, Time, X] = simulateVP(model,type, l_s,l_r,K,a_rs,a_sr,k,Cmax, ncycles,T, M_diagnosis,M_death,sigma )

% set initial & stopping conditions
X_s = (1-sigma)*M_diagnosis;
X_r = sigma*M_diagnosis;
C = Cmax;
opts_TimetoDeath = odeset('Events', model{2});

% set outputs 
OS = 0;   % overall survival
Response = 0; % Response = tumor volume after treatment/ tumor volume at diagnosis 
Time = [];     % time vector
X = [];        % Amount of Xs, X_r & C over time  

%% == CASE 1: no treatment
if T ==0 && ncycles==0 
    sim = ode45(model{1},[0 20*12*30],[X_s, X_r 0],opts_TimetoDeath,l_s,l_r,K,K,a_rs,a_sr,0,M_death);
    OS  = sim.x(end);
    response = 1;  % no response
    Time = sim.x;
    X = sim.y;
    
    return
end

%% CASE 2: MTD or MT treatment WITHOUT drug holiday
if isequal(type,'MTD') ||isequal(type,'MT')
    
    if T == 0 && ncycles ==1
        x_relapse = ode45(model{1},[0 20*12*30],[X_s, X_r C],opts_TimetoDeath,l_s,l_r,K_r,K_s,a_rs,a_sr,k,M_death);  
        R =  min(sum(x_relapse.y(1:2,:)));
        response = log10((R)/M_diagnosis);
        OS = x_relapse.x(end)/30;
        Time = x_relapse.x;
        X = x_relapse.y;
        return 
    end
    
    % simulate treatment phase
    T_curr = 0;
    for i = 1:ncycles
        C = Cmax;
        x_chemo = ode45(model{1},[T_curr T_curr+T],[X_s, X_r C],opts_TimetoDeath,l_s,l_r,K,K,a_rs,a_sr,k,M_death);
        Time = [Time x_chemo.x];
        X =[X x_chemo.y];
        if ~isempty(x_chemo.xe)
            OS = x_chemo.xe/30;
            response =log10(M_death/M_diagnosis);
            return
        end
        X_s = x_chemo.y(1,end);
        X_r = x_chemo.y(2,end);
        C = x_chemo.y(3,end);    
        T_curr = T_curr+T;  
    end
    % simulate post-treatment phase
    x_relapse = ode45(model{1},[T_curr 20*12*30],[X_s, X_r C],opts_TimetoDeath,l_s,l_r,K,K,a_rs,a_sr,k,M_death);  
    Time = [Time x_relapse.x];
    X =[X x_relapse.y];
    R =  min([sum(x_relapse.y(1:2,:)), 0]);
    if R == 0
       R = .000001; 
    end
    response = log10(R/M_diagnosis);
    OS = x_relapse.x(end)/30;
    return
end

end

