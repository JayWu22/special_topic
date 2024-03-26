clear;
K = 3;
Nt = 4;
N = 10;
ISR = 0.2;

%% for training
samples = 12000;
numFeature = 2*Nt*K^2 + 1;
features = zeros(samples, numFeature);
labels = zeros(samples, K * N);
optimalRate = zeros(samples, 1);
RateMRT = zeros(samples,1);
RateZF = zeros(samples,1);
for sample = 1:samples
    h1 = sqrt(1/2) * (randn(Nt,K) + 1j*randn(Nt,K));
    h2 = sqrt(1/2) * (randn(Nt,K) + 1j*randn(Nt,K));
    h3 = sqrt(1/2) * (randn(Nt,K) + 1j*randn(Nt,K));
    num = 1;
    for i = 1:K
        for j = 1:Nt
            features(sample, num) = real(h1(j,i));
            num = num + 1;
            features(sample, num) = imag(h1(j,i));
            num = num + 1;
            features(sample, num) = real(h2(j,i));
            num = num + 1;
            features(sample, num) = imag(h2(j,i));
            num = num + 1;
            features(sample, num) = real(h3(j,i));
            num = num + 1;
            features(sample, num) = imag(h3(j,i));
            num = num + 1;
        end
    end
    if sample < 2000
        P = 1;
    elseif sample < 4000
        P = 3.16;
    elseif sample < 6000
        P = 10;
    elseif sample < 8000
        P = 31.6;
    elseif sample < 10000
        P = 100;
    else 
        P = 316;
    end
    features(sample, num) = P;
    mrt1 = h1(:,1) / norm(h1(:,1));
    mrt2 = h2(:,2) / norm(h2(:,2));
    mrt3 = h3(:,3) / norm(h3(:,3));
    H1 = [h1(:,2), h1(:,3)];
    H2 = [h2(:,1), h2(:,3)];
    H3 = [h3(:,1), h3(:,2)];
    zf1 = (eye(Nt) - H1*((H1'*H1)\H1'))*h1(:,1);
    zf1 = zf1 / norm(zf1(:,1));
    zf2 = (eye(Nt) - H2*((H2'*H2)\H2'))*h2(:,2);
    zf2 = zf2 / norm(zf2(:,1));
    zf3 = (eye(Nt) - H3*((H3'*H3)\H3'))*h3(:,3);
    zf3 = zf3 / norm(zf3(:,1));
    max = 0;
    for i = 0:N-1
        for j = 0:N-1
            for k = 0:N-1
                v1 = (i/(N-1))*mrt1+((N-1-i)/(N-1))*zf1;
                v2 = (j/(N-1))*mrt2+((N-1-j)/(N-1))*zf2;
                v3 = (k/(N-1))*mrt3+((N-1-k)/(N-1))*zf3;
                v1 = v1 / norm(v1);
                v2 = v2 / norm(v2);
                v3 = v3 / norm(v3);
                Rate1 = log2(1+(P*(abs(h1(:,1)'*v1))^2)/(1+ISR*P*((abs(h2(:,1)'*v2))^2+(abs(h3(:,1)'*v3))^2)));
                Rate2 = log2(1+(P*(abs(h2(:,2)'*v2))^2)/(1+ISR*P*((abs(h1(:,2)'*v1))^2+(abs(h3(:,2)'*v3))^2)));
                Rate3 = log2(1+(P*(abs(h3(:,3)'*v3))^2)/(1+ISR*P*((abs(h1(:,3)'*v1))^2+(abs(h2(:,3)'*v2))^2)));
                sumrate = Rate1 + Rate2 + Rate3;
                if sumrate > max
                    max = sumrate;
                    maxi = i;
                    maxj = j;
                    maxk = k;
                end
            end                                                                                                             
        end
    end
    baselineMRT1 = log2(1+(P*(abs(h1(:,1)'*mrt1))^2)/(1+ISR*P*((abs(h2(:,1)'*mrt2))^2+(abs(h3(:,1)'*mrt3))^2)));
    baselineMET2 = log2(1+(P*(abs(h2(:,2)'*mrt2))^2)/(1+ISR*P*((abs(h1(:,2)'*mrt1))^2+(abs(h3(:,2)'*mrt3))^2)));
    baselineMRT3 = log2(1+(P*(abs(h3(:,3)'*mrt3))^2)/(1+ISR*P*((abs(h1(:,3)'*mrt1))^2+(abs(h2(:,3)'*mrt2))^2)));
    RateMRT(sample,1) = baselineMRT1+baselineMET2+baselineMRT3;
    baselineZF1 = log2(1+(P*(abs(h1(:,1)'*zf1))^2)/(1+ISR*P*((abs(h2(:,1)'*zf2))^2+(abs(h3(:,1)'*zf3))^2)));
    baselineZF2 = log2(1+(P*(abs(h2(:,2)'*zf2))^2)/(1+ISR*P*((abs(h1(:,2)'*zf1))^2+(abs(h3(:,2)'*zf3))^2)));
    baselineZF3 = log2(1+(P*(abs(h3(:,3)'*zf3))^2)/(1+ISR*P*((abs(h1(:,3)'*zf1))^2+(abs(h2(:,3)'*zf2))^2)));
    RateZF(sample,1) = baselineZF1+baselineZF2+baselineZF3;
    labels(sample, maxi + 1) = 1;
    labels(sample, N + maxj+1) = 1;
    labels(sample, 2*N + maxk+1) = 1;
    optimalRate(sample, 1) = max;
end
inputdata_file = matfile('inputdata.mat','Writable',true);
inputdata_file.features = features;
inputdata_file.labels = labels;
inputdata_file.optimalRate = optimalRate;
inputdata_file.RateMRT = RateMRT;
inputdata_file.RateZF = RateZF;
%% for testing
samples = 2400;
features = zeros(samples, numFeature);
labels = zeros(samples, K * N);
optimalRate = zeros(samples, 1);
testMRT = zeros(samples, Nt, K);
testZF = zeros(samples, Nt, K);
testh1 = zeros(samples, Nt, K);
testh2 = zeros(samples, Nt, K);
testh3 = zeros(samples, Nt, K);

for sample = 1:samples
    h1 = sqrt(1/2) * (randn(Nt,K) + 1j*randn(Nt,K));
    testh1(sample,:,:) = h1;
    h2 = sqrt(1/2) * (randn(Nt,K) + 1j*randn(Nt,K));
    testh2(sample,:,:) = h2;
    h3 = sqrt(1/2) * (randn(Nt,K) + 1j*randn(Nt,K));
    testh3(sample,:,:) = h3;
    num = 1;
    for i = 1:K
        for j = 1:Nt
            features(sample, num) = real(h1(j,i));
            num = num + 1;
            features(sample, num) = imag(h1(j,i));
            num = num + 1;
            features(sample, num) = real(h2(j,i));
            num = num + 1;
            features(sample, num) = imag(h2(j,i));
            num = num + 1;
            features(sample, num) = real(h3(j,i));
            num = num + 1;
            features(sample, num) = imag(h3(j,i));
            num = num + 1;
        end
    end
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
    features(sample, num) = P;
    mrt1 = h1(:,1) / norm(h1(:,1));
    mrt2 = h2(:,2) / norm(h2(:,2));
    mrt3 = h3(:,3) / norm(h3(:,3));
    mrt = [mrt1 mrt2 mrt3];
    testMRT(sample,:,:) = mrt;
    H1 = [h1(:,2), h1(:,3)];
    H2 = [h2(:,1), h2(:,3)];
    H3 = [h3(:,1), h3(:,2)];
    zf1 = (eye(Nt) - H1*((H1'*H1)\H1'))*h1(:,1);
    zf1 = zf1 / norm(zf1(:,1));
    zf2 = (eye(Nt) - H2*((H2'*H2)\H2'))*h2(:,2);
    zf2 = zf2 / norm(zf2(:,1));
    zf3 = (eye(Nt) - H3*((H3'*H3)\H3'))*h3(:,3);
    zf3 = zf3 / norm(zf3(:,1));
    zf = [zf1 zf2 zf3];
    testZF(sample,:,:) = zf;
    max = 0;
    for i = 0:N-1
        for j = 0:N-1
            for k = 0:N-1
                v1 = (i/(N-1))*mrt1+((N-1-i)/(N-1))*zf1;
                v2 = (j/(N-1))*mrt2+((N-1-j)/(N-1))*zf2;
                v3 = (k/(N-1))*mrt3+((N-1-k)/(N-1))*zf3;
                v1 = v1 / norm(v1);
                v2 = v2 / norm(v2);
                v3 = v3 / norm(v3);
                Rate1 = log2(1+(P*(abs(h1(:,1)'*v1))^2)/(1+ISR*P*((abs(h2(:,1)'*v2))^2+(abs(h3(:,1)'*v3))^2)));
                Rate2 = log2(1+(P*(abs(h2(:,2)'*v2))^2)/(1+ISR*P*((abs(h1(:,2)'*v1))^2+(abs(h3(:,2)'*v3))^2)));
                Rate3 = log2(1+(P*(abs(h3(:,3)'*v3))^2)/(1+ISR*P*((abs(h1(:,3)'*v1))^2+(abs(h2(:,3)'*v2))^2)));
                sumrate = Rate1 + Rate2 + Rate3;
                if sumrate > max
                    max = sumrate;
                    maxi = i;
                    maxj = j;
                    maxk = k;
                end
            end                                                                                                             
        end
    end
    labels(sample, maxi + 1) = 1;
    labels(sample, N + maxj+1) = 1;
    labels(sample, 2*N + maxk+1) = 1;
    optimalRate(sample, 1) = max;
end

inputdata_file.testFeatures = features;
inputdata_file.testLabels = labels;
inputdata_file.testOptimalRate = optimalRate;
inputdata_file.MRT = testMRT;
inputdata_file.ZF = testZF;
inputdata_file.testh1 = testh1;
inputdata_file.testh2 = testh2;
inputdata_file.testh3 = testh3;