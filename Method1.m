function  parameters = Method1(model, ranges)

n = 1000;
p = 8;
M_death     = (4/3)*pi*7.5^3*10^9; 
K = (4/3)*pi*15^3*10^9; 
k = 0.2110;

X = lhsdesign(n,p);
DT = X(:,1).*(ranges(2,1)-ranges(1,1))+ranges(1,1);
a_sr = X(:,2).*(ranges(2,2)-ranges(1,2))+ranges(1,2);
a_rs = X(:,3).*(ranges(2,3)-ranges(1,3))+ranges(1,3);
M_diagnosis = X(:,4).*(ranges(2,4)-ranges(1,4))+ranges(1,4);
sigma = X(:,5).*(ranges(2,5)-ranges(1,5))+ranges(1,5);
CC = X(:,6).*(ranges(2,6)-ranges(1,6))+ranges(1,6);
T =  X(:,7).*(ranges(2,7)-ranges(1,7))+ranges(1,7);
Cmax = X(:,8)*(ranges(2,8)-ranges(1,8))+ranges(1,8);


i = 1;
counter = 1;
while 1
    disp(counter)
    DT_r = 2*DT(i);
    l_s = log(2)/DT(i);
    l_r = log(2)/DT_r;
    
    X_s = (1-sigma(i))*M_diagnosis(i);
    X_r = sigma(i)*M_diagnosis(i);
    
    OS_model(i)= simulateVP(model,'MTD',l_s,l_r,K,a_rs(i),a_sr(i),k,Cmax(i), CC(i),T(i), M_diagnosis(i),M_death,sigma(i));

    Y(i,:) = [DT(counter) a_sr(counter) a_rs(counter) M_diagnosis(counter) sigma(counter) CC(counter) T(counter) Cmax(counter)];
    i = i+1;
    counter = counter+1;
    if counter > n
       break 
    end
    
end

p = {'DT','a_{sr}','a_{rs}', 'M_{diagnosis}', '\sigma', 'CT_{cycles}','T','Cmax'};
figure 
for i = 1:8
   subplot(3,3,i)
   plot(Y(:,i),OS_model,'k.')
   title(p{i}) 
end



figure
p = {'DT','a_{sr}','a_{rs}','M_{diagnosis}','\sigma',  'CC', 'T' 'Cmax',};
S = zeros(numel(p),1);
for i = 1:size(X,2)
   S(i) = corr(Y(:,i),OS_model'); 
end

[~, I_sort] = sort(abs(S),'descend');
S_sorted = S(I_sort);
p1=bar(S_sorted);
set(p1,'FaceColor',[100 100 100]./255);
hold on
ax = gca;
ax.XTickLabel =p;
ax.XTickLabel = ax.XTickLabel(I_sort);
ax = gca;
ax.XTickLabelRotation=90;
ylabel('Sensivity coefficient')
set(gca, 'FontSize',14)

parameters = find(S>0.2);


end


