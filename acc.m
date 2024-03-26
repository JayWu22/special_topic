clear;
load('inputdata.mat')
load('outputdata.mat')

samples = 2400;
K = 3;
N = 10;
P = 100;
ISR = 0.2;

predRatio = zeros(samples, 3);
sumrate = zeros(samples,1);
for sample = 1:samples
    for i = 0:K-1
        max = 0;
        for j = i*10+1:i*10+10 
            if max < pred(sample, j)
                max = pred(sample, j);
                maxj = j - 10*i - 1;
            end
        end
        predRatio(sample, i+1) = maxj;
    end
end
performance = zeros(samples,1);
A =  MRT(1,:,1);
for sample = 1: samples
    if sample < 400
        P = 1;
    elseif sample < 800
        P = 3.16;
    elseif sample < 1200
        P = 10;
    elseif sample < 1600
        P = 31.6;
    elseif sample < 2000
        P = 100;
    else
        P = 316;
    end
    v1 = (predRatio(sample,1)/(N-1))*MRT(sample,:,1)'+((N-1-predRatio(sample,1))/(N-1))*ZF(sample,:,1)';
    v2 = (predRatio(sample,2)/(N-1))*MRT(sample,:,2)'+((N-1-predRatio(sample,2))/(N-1))*ZF(sample,:,2)';
    v3 = (predRatio(sample,3)/(N-1))*MRT(sample,:,3)'+((N-1-predRatio(sample,3))/(N-1))*ZF(sample,:,3)';
    v1 = v1 / norm(v1);
    v2 = v2 / norm(v2);
    v3 = v3 / norm(v3);
    Rate1 = log2(1+(P*(abs(testh1(sample,:,1)*v1))^2)/(1+ISR*P*((abs(testh2(sample,:,1)*v2))^2+(abs(testh3(sample,:,1)*v3))^2)));
    Rate2 = log2(1+(P*(abs(testh2(sample,:,2)*v2))^2)/(1+ISR*P*((abs(testh1(sample,:,2)*v1))^2+(abs(testh3(sample,:,2)*v3))^2)));
    Rate3 = log2(1+(P*(abs(testh3(sample,:,3)*v3))^2)/(1+ISR*P*((abs(testh1(sample,:,3)*v1))^2+(abs(testh2(sample,:,3)*v2))^2)));
    sumrate(sample,1) = Rate1 + Rate2 + Rate3;
    performance(sample,1) = sumrate(sample,1) / testOptimalRate(sample,1);
end

AVGperformance = mean(performance, "all")

baselineMRT = [mean(RateMRT(1:1999,1),"all"), mean(RateMRT(2000:3999,1),"all"), mean(RateMRT(4000:5999,1),"all"),mean(RateMRT(6000:7999,1),"all"),mean(RateMRT(8000:9999,1),"all"),mean(RateMRT(10000:12000,1),"all")];
baselineZF = [mean(RateZF(1:1999,1),"all"), mean(RateZF(2000:3999,1),"all"), mean(RateZF(4000:5999,1),"all"),mean(RateZF(6000:7999,1),"all"),mean(RateZF(8000:9999,1),"all"),mean(RateZF(10000:12000,1),"all")];
MLrate = [mean(sumrate(1:399,1),"all"), mean(sumrate(400:799,1),"all"), mean(sumrate(800:1199,1),"all"),mean(sumrate(1200:1599,1),"all"),mean(sumrate(1600:1999,1),"all"),mean(sumrate(2000:2400,1),"all")];
optimal = [mean(testOptimalRate(1:399,1),"all"), mean(testOptimalRate(400:799,1),"all"), mean(testOptimalRate(800:1199,1),"all"),mean(testOptimalRate(1200:1599,1),"all"),mean(testOptimalRate(1600:1999,1),"all"),mean(testOptimalRate(2000:2400,1),"all")];

y = [0 5 10 15 20 25];
plot(y,optimal,'-s');
hold on
plot(y,baselineMRT,'-v');
plot(y,baselineZF,'-p');
plot(y,MLrate,'-o');
hold off
legend({'Optimal','baselineMRT','baselineZF', 'ML-based'},'Location','northwest')
title('Sum Rate in 3-user interference channel');
yticks([0 5 10 15 20 25])
xlabel('Transmit Power(dB)') 
ylabel('Rate(bits/s/Hz)')

%normalizedRate = [MLrate(1,1)/optimal(1,1), MLrate(1,2)/optimal(1,2), MLrate(1,3)/optimal(1,3) ,MLrate(1,4)/optimal(1,4) ,MLrate(1,5)/optimal(1,5),MLrate(1,6)/optimal(1,6)];
%plot(y, normalizedRate,'-o');
%title('Normalized rate')
%xlabel('Transmit Power(dB)')
%ylabel('Normalized Rate')
%plot(performance);
