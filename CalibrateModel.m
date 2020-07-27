function [GMM1, GMM2] = Method2(model,data)
%% bootstrap # of CT cycles(CTC) and time interval between 2 consecutive CT cycles (T)
N = 10000;
cc = zeros(N,1);
tt =  zeros(N,1);
    for idx = 1:N
        patient = randi(numel(data.OS));
        cc(idx) = data.CTC(patient);
        tt(idx) = data.T(patient);
    end
GMM1 = fitgmdist([cc tt],1,'Options', statset('Display','final'));

plotHistogram_T_CT(GMM1,data,N)

%% 2. set constant model parameters and define range for non-constant one
k = 0.2110;
M_death     = (4/3)*pi*10^3*10^9; 
M_diagnosis = (4/3)*pi*2^3*10^9;
K = (4/3)*pi*15^3*10^9; 
a_rs = 0;
a_sr = 0;

for Cmax = 5 : 40
DT_seq          = 1:10:1000;
sigma_seq       = 0:.025:1;
[A,B] = meshgrid(DT_seq, sigma_seq);
A = A(:);
B = B(:);

%% 3. simulate model for each parameter combination (Brute-Force algorithm)
OS = zeros(numel(A),1);
tic
parfor i = 1:numel(A)
    DT_s = A(i);
    DT_r = DT_s/2;
    l_s = log(2)/DT_s;
    l_r  = log(2)/DT_r;

    OS_dash = zeros(100,1);
    n =1;
    while true
        [CTC,T] = generate_CT_T(GMM1);       
        [OS_dash(n)] = simulateVP(model,'MTD',l_s,l_r,K,K,a_rs,a_sr,k,Cmax, CTC,T, M_diagnosis, M_death, B(i));
        n = n+1;
    if n >100
        break
    end

    end

    OS(i) = round(mean(OS_dash));
end
toc
patients = [A B OS];

%% 4. bootstrap the non-constant parameters and fit to GMM model
OS_dash = zeros(10000,1);
for idx = 1:10000
    Patient_OS      = data.OS(randi(numel(data.OS)));
    Virtual_patient = patients(patients(:,3) == Patient_OS,:);    
    Virtual_patient = Virtual_patient(randi(size(Virtual_patient,1)),: );    
    DT(idx) = Virtual_patient(1);
    sigma(idx)     = Virtual_patient(2);
    OS_dash(idx)= Virtual_patient(3);
end
GMM2 = fitgmdist([DT; sigma]',3);

OS_model = zeros(1000,1);
R = zeros(1000,1);
i = 1;
sigma= [];
while true
    disp(i)
    r  = generate_DT_sigma(GMM2); 
    DT_s(i) = r(1);
    DT_r = DT_s(i)/2;
    l_s = log(2)/DT_s(i);
    l_r  = log(2)/DT_r;
    [CC, T]= generate_CT_T(GMM1);

    if(CC==0)
        continue
    end

[OS_model(i), R(i)] = simulateMDTpatient(model,'MTD',l_s,l_r,K,K,a_rs,a_sr,k,Cmax, CC,T, M_diagnosis, M_death,r(2));
    i = i+1;
    if(i > 1000)
       break 
    end
end

%% 6. check fittness to response data    
Resp_VPC = R;
idx_SD = Resp_VPC<log10(1.30) & Resp_VPC>log10(0.90);
idx_CR = Resp_VPC<=log10(.5);
idx_PD = Resp_VPC>log10(1.30);
idx_PR = Resp_VPC > log10(.5) &  Resp_VPC < log10(.90); 
Response_VPC = [sum(idx_PD),sum(idx_SD), sum(idx_PR), sum(idx_CR)];
count = histcounts(data.Response);
Table = [count/numel(data.Response); Response_VPC./1000]';
[h,p]=chi2cont(Table);
end

if p>0.99
    disp('Cmax found:')
    disp(C_max)
end


end

function plotHistogram_T_CT(GMM1,data,N)
%% plot results
CTC_dash = zeros(N,1);
T_dash = zeros(N,1);
i = 1;
while true
    [CTC_dash(i), T_dash(i)] = generate_CT_T(GMM1);
    if CTC_dash(i) ==0 
       continue 
    end
    i = i+1;
    
    if i >N
        break;
    end
    
end

x = [T_dash' data.T_CT];
group = [CTC_dash' data.CTC_CT+.5];
positions = [1 2 3 4 5 6 1.5 2.5 3.5 4.5 5.5 6.5];
boxplot(x,group, 'positions', sort(positions),'MedianStyle','target', 'Symbol','ko');

set(gca,'xtick',[1.25 2.25 3.25 4.25 5.25 6.25 ])
set(gca,'xticklabel',{'1','2','3', '4', '5', '6'})

color = repmat(parula(2),6,1);

h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',0.5);
end

c = get(gca, 'Children');

legend(c(1:2), 'Clinical cohort', 'Virtual cohort' );

xlabel('Number of CT cycles')
ylabel('Time interval between to CT cycles')
set(gca, 'FontSize',14)
end