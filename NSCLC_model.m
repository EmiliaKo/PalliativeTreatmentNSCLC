%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model of NSCLC growth                                                   %
% Model definition                                                        %
% Author: Emilia Kozlowska                                                %
% Last update: 03/10/2018                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = NSCLC_model()
out{1}  = @model;
out{2}  = @TimeToDeath;
end

%% the model equations
function dy = model(t,x,lambda_s,lambda_r,K_s,K_r,a_rs,a_sr,k,M_death)
        dy   = zeros(3,1);    
        x_s  = x(1);  % sensitive cells
        x_r  = x(2);  % resistant cells
        C = x(3);     % cisplatin 
        
        dy(1) =  lambda_s*x_s*(1-((x_s+a_rs*x_r)/(K_s)))*(1-C); 
        dy(2) =  lambda_r*x_r*(1-((x_r+a_sr*x_s)/(K_r))); 
        dy(3) = -k*C;
end

%% locate death
function [value,isterminal,direction] = TimeToDeath(t, x, lambda_s,lambda_r,K_s,K_r,a_rs,a_sr,d_chemotherapy,M_death)
    value = ( x(1) + x(2) ) - M_death; 
    isterminal = 1; 
    direction  = 0;  
end




