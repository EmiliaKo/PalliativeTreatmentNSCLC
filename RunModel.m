ranges = [
 20 1000   
0 500
0 500
10^9 10^11
0.005 .999
0 6
14 40
0 40
]';

model = NSCLC_model();

 p = Method1(model, ranges);
 
data.OS_NH =[1 2 25 1 5 3 1 20 3 4 2 8 1 27 1 9 2];
data.CTC_NH = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
data.T_NH = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];


data.CTC_CT = [6 4 3 3 3 3 3 2 2 3 5 1 2 2 4 3 3 1 4 2 3 1 1 1 2];
data.T_CT  = [30 30 23 16 17 19 20 15 36 16 30 0 11 29 16 14 14 0 26 26 14 0 0 0 14];
data.OS_CT = [8 14 4 4 8 10 4 4 20 4 10 13 6 81 3 9 19 15 15 14 7 4 6 1 4];
data.Response = [2 0 0 1 2 1 0 0 2 1 0 0 2 1 0 3 0 0 2 0 1 0 0 0 0];

data.OS = [data.OS_NH data.OS_CT];
data.T = [data.T_NH data.T_CT];
data.CTC = [data.CTC_NH data.CTC_CT];

 [GMM1] = Method2(model,data);