%Yilun Jiang
%Prof. Fontaine
%Problem set 3
clear;
%% Problem 1
%Some global parameter for the simulation
N = 250; %Total number of steps for simulation
variance = 0.01; %variance of the increment dW
stepSize = 0.01; %step size of each time step
interest = 0.1/N/3; %Interest rate
%a)
%simulate 1000 paths of length N/2 and report the average
%The function that generates the simulation is in the end of the script
simulateResult = zeros(1,1000);
for i = 1:1000
    path = simulatePath(N/2,stepSize,variance);
    simulateResult(i) = path(N/2);
end
averageResult = sum(simulateResult) / 1000;

%b)
%Simulate 10 paths and graph them
tenPaths = zeros(10,N/2);
for i = 1:10
    path = simulatePath(N/2,stepSize,variance);
    tenPaths(i,:) = path;
end
figure
n = 1:1:125;
plot(n,tenPaths(:,:))
title('10 sumperimposed simulated paths')
xlabel('step')
ylabel('value')
xlim([1 125])

%c)
%The formula to express St in terms of Xt
%S = Xt * e^(r*(250-t)) where r is the interest rate and S is the simulated
%result when we took 125 steps

%Xt is the discounted price, and for each of the above 10 paths we generated in
%part b, we will simulate it from N/2 to N for 1000 times.
xErrorOverPaths = zeros(1,10);
for i = 1:10
    xtSimulate = zeros(1000,N/2);
    %Start represents the starting point of simulation 
    start = tenPaths(i,1); 
    %function simulatedDiscounted in in the bottom of the script
    xtSimulate = simulateDiscouted(start,N/2,variance,interest,stepSize);
    %Use the property of martingale to pull xt back to N/2 and compare it
    %to the given start. Since we start at N/2, there are only N/2-1 steps left
    xPullBack = xtSimulate(:,125) .* exp(-interest*(N/2-1));
    %Average squared error for simulation of one given path over 1000 times
    xError = sum((xPullBack - start).^2) / 1000;
    xErrorOverPaths(i) = xError;
end
%Each col represents the average squared error for simulation of one given path
%over 1000 times
xErrorOverPaths  
%Complete the sentence
%Actually, this becomes a fairly trivial issue, because dX = 0*dt + sigma*X*dW 
%and it is clearly obvious, therefore, that Xi[N] = Xi[N/2] + Y where Y is a random 
%variable that is independent of the filtration up to time N/2 and (under the risk
%neutral measure) has mean value 0, so obviously E[ X[N] | N/2 ] = X[N/2].

%d)
%This is an estimate of the variance no condition is involved
xErrorEstimate = sum(xErrorOverPaths) / 10

%Through couple times of iteation, all average squared errors I get is pretty small. So I
%guess it is safe to say the martingale property is true
%% Problem 2
%Simulate the interest rate using two different models. One is Hull-White
%and the other is Cox-Ingersoll-Ross

%Some global paramter for Problem 2
T = 2.5; %Total time for the simulation
stepSize = 0.01; %Step size in the time 
variance = 0.01; %variance in the dW
t = 0:stepSize:2.5; %vector representing time, step size is 0.01
Kb = 0.1; %Obtained through trial and error until something reasonable
Ksigma = 0.1; %Obtained through trial and error until something reasonable
Ka = 0.1; %Obtained through trial and error until something reasonable
lowerBound = 0.01; %If the simulated result is below the lowerbound, then we clipped it to the lowerbound

%I have chosen parameters(Ka,Ksimga,Kb) so that the simulated result will
%be around 1
%b)
B_t = Kb * (1.1+sin(pi.*t./T));
Sigma_t = Ksigma * (1.1+cos(4*pi.*t./T));
A_t = 0.5 .* Sigma_t.^2 + Ka * (1.1+cos(pi.*t./T));

figure
plot(t, B_t)
hold on
plot(t, Sigma_t)
plot(t, A_t)
xlabel('time')
ylabel('value')
legend('B_t','Sigma_t','A_t')
title('B_t, Sigma_t, A_t') % make sure all these parameters are strictly positive

R_HW = ones(10,251);
for trial = 1:10
    for pos = 2:251
        R_HW(trial,pos) = R_HW(trial, pos-1) + (A_t(pos) - B_t(pos)*R_HW(trial,pos-1))*stepSize ...
                          + Sigma_t(pos) *  normrnd(0,sqrt(variance));
    end
end
%clipped the simulated result so they stay above the lower bound
indexHW = find(R_HW < lowerBound);
R_HW(indexHW(:)) = lowerBound;
%Graph the simulated result
figure
plot(t, R_HW(:,:))
xlabel('time')
ylabel('value')
title('Simulated interest rate(Hull-White)')

%Simulate Interest rate using CIR model
R_CIR = ones(10,251);
for trial = 1:10
    for pos = 2:251
        R_CIR(trial,pos) = R_CIR(trial, pos-1) + (A_t(pos) - B_t(pos)*R_CIR(trial,pos-1))*stepSize ...
                          + Sigma_t(pos) * sqrt(R_CIR(trial,pos-1)) * normrnd(0,sqrt(variance));
    end
end
%clipped the simulated result so they stay above the lower bound
indexCIR = find(R_CIR < lowerBound);
R_CIR(indexCIR(:)) = lowerBound;
%Graph the simulated interest rate using CIR model
figure
plot(t, R_CIR(:,:))
xlabel('time')
ylabel('value')
title('Simulated interest rate(Cox-Ingersoll-Ross)')

%c)
%Estimate of the mean and variance of both HW and COX model at each time step using monte
%carlo mathod with 500 trials
numberOfSimulation = 500; % A parameter that can be changed

%Estimate mean and variance for the HW model at each time step
R_HW_MeanEstimate = ones(1,251);
R_HW_VarEstimate = ones(1,251);
path = ones(numberOfSimulation ,251);
for trial = 1:numberOfSimulation 
    for pos = 2:251
        path(trial,pos) = path(trial,pos-1) + (A_t(pos) - B_t(pos)*path(trial,pos-1))*stepSize ...
                    + Sigma_t(pos) * normrnd(0,sqrt(variance));
    end
end
for pos = 1:251
    mean = sum(path(:,pos))/numberOfSimulation;
    R_HW_MeanEstimate(pos) = mean;
    var = 1/numberOfSimulation * sum( (path(:,pos)-mean).^2 );
    R_HW_VarEstimate(pos) = var;
end

%Estimate mean and variance at each time step for the CIR model
R_CIR_MeanEstimate = ones(1,251);
R_CIR_VarEstimate = ones(1,251);
path = ones(numberOfSimulation ,251);
for trial = 1:numberOfSimulation 
    for pos = 2:251
        path(trial,pos) = path(trial,pos-1) + (A_t(pos) - B_t(pos)*path(trial,pos-1))*stepSize ...
                    + Sigma_t(pos) * sqrt(path(trial,pos-1)) * normrnd(0,sqrt(variance));
    end
end
for pos = 1:251
    mean = sum(path(:,pos))/numberOfSimulation;
    R_CIR_MeanEstimate(pos) = mean;
    var = 1/numberOfSimulation * sum( (path(:,pos)-mean).^2 );
    R_CIR_VarEstimate(pos) = var;
end

figure
pos = 1:1:251;
plot(pos, R_HW_MeanEstimate)
hold on
plot(pos, R_CIR_MeanEstimate)
hold off
xlabel('trial')
xlim([1 251])
ylabel('value')
legend('HW model', 'CIR model')
title('Estimated mean HW VS CIR')

figure
plot(pos, R_HW_VarEstimate)
hold on
plot(pos,  R_CIR_VarEstimate)
hold off
xlabel('trial')
xlim([1 251])
ylabel('value')
legend('HW model', 'CIR model')
title('Estimated variance HW VS CIR')

%d)
%Check how often does the simulation falls below the lower bound and we
%have to clip it
txt = sprintf('Number of exception occurs at HW model: %d', size(indexHW,1));
txt
txt = sprintf('Number of exception occurs at CIR model: %d', size(indexCIR,1));
txt
%Through couple iterations, I found that the exceptions have never occured
%but that's probably due to my choice of parameter(i.e. Kb, Ksigma and Ka)

%% Problem 3
%Simulate the correlated brownian motion

%Some global parameters for problem 3
m = 2; %Number of correlated brownian motion we will simulate in this problem
covariance = [1,0.75;0.75,0.9]; %Covariance matrix for the random vector 
length = 250; %Total number of time step for the simulation

%a) Function written in the bottom of the script
% Generate a set of M Gaussian random vectors that having identical mean 
% vector and covariance matrix. The core is using Cholesky to decompose the
% covariance matrix and then multiply the decomposed result with the
% normal random vector and add on the mean vector to generate the desired vector

%b)
%Generate 10 pairs of correlated brownian motion vectors with the required
%length and graph each pair of them
meanVector = zeros(m,1); %The mean vector is assumed to be column vector
Wmatrix = zeros(20,250);
for trial = 1:10 
    Wmatrix(trial*2-1:trial*2,:) = generateCorreGaussRV(meanVector, covariance, m, length);
    N = 1:1:250;
    figure
    plot(N, Wmatrix(trial*2-1,:))
    hold on
    plot(N, Wmatrix(trial*2,:))
    xlabel('time')
    ylabel('Simulated W')
    legend('W1', 'W2')
end
%Actually stuck on using Matlab to solve R_HW symbolically
%c)
%Using the generated 10 pairs of correlated brownian motion to drive the
%model and graph the output of the model

%Some global parameters for this problem
N = 250; %Total number of time steps for the simulation 
alphaVec = [0.1/N,0.1/N]; %Use the alpha from the previous problem
sigmaMatrix = [0.1/N,0.2/N;0.3/N,0.4/N]; %Arbitrary chosen sigma matrix
stepSize = 0.01; %Step size for the simulation

S = ones(20,N);
for trial = 1:10
    for pos = 2:250
        S(trial*2-1,pos) = S(trial*2-1,pos-1) * (1 + alphaVec(1)*stepSize + sigmaMatrix(1)...
                   *Wmatrix(trial,pos) + sigmaMatrix(3)*Wmatrix(trial,pos));
        S(trial*2,pos) = S(trial*2,pos-1) * (1 + alphaVec(2)*stepSize + sigmaMatrix(2)...
                   *Wmatrix(trial*2,pos) + sigmaMatrix(4)*Wmatrix(trial*2,pos));
    end
    figure
    steps = 1:1:250;
    plot(steps, S(trial*2-1,:))
    hold on
    plot(steps, S(trial*2,:))
    hold off
    xlabel('time')
    ylabel('value')
    legend('Model driven by W1', 'Model driven by W2')
end

%% Problem 4
%The optimization problem is minimize(w'*sigma*C*sigma'*w) with constraint 
%sum(w) = 1. We discuss here only the case of portfolio with two stocks
%Define Cdot = sigma*C*simga'

cdot = sigmaMatrix*covariance*sigmaMatrix';
syms w1 w2 lenda
volatility = w1*(w1*cdot(1)+w2*cdot(2)) + w2*(w1*cdot(3)+w2*cdot(4)); %This is a volatity function about w1 and w2
% After some algebra(the details are on the paper and attached as a pdf),
% The work is mainly take partial derivates and align the equations nicely
eqn1 = 2*w1*cdot(1) + (cdot(2)+cdot(3))*w2 - lenda == 0;
eqn2 = 2*w2*cdot(4) + (cdot(2)+cdot(3))*w1 -lenda == 0;
eqn3 = w1 + w2 == 1;
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [w1, w2, lenda]);
answerVec = linsolve(A,B);
w1 = answerVec(1);
w2 = answerVec(2);
w_Larange = [w1;w2]
%Compute the volatility for three different weight vectors
volatility_Lagarange = sqrt(w_Larange'*sigmaMatrix*covariance*sigmaMatrix'*w_Larange);
w_boundry1 = [0;1];
volatility_boundry1 = sqrt(w_boundry1'*sigmaMatrix*covariance*sigmaMatrix'*w_boundry1);
w_boundry2 = [1;0];
volatility_boundry2 = sqrt(w_boundry2'*sigmaMatrix*covariance*sigmaMatrix'*w_boundry2);

%Volatility through weight found by Lagarange multiplier w = [1.7426, -0.7426] is 1.26e-4,
%however, one of the weight is negative so we have to reject this point
%Volatility through weight defined as [0;1] is 0.0026
%Volatility through weight defined as [1;0] is 0.0011
%So w=[0;1] will maximize the volatility
%and w=[1;0] will minimize the volatility

w_max = [0;1];
w_min = [1;0];
N = 250;
portfolio = zeros(20,N);
step = 1:1:250;
for trial = 1:10
    %S(trial*2-1:trial*2,:) represents the two simulated stock price using
    %one pair of w
    portfolio(trial*2-1,:) = w_max' * S(trial*2-1:trial*2,:);
    portfolio(trial*2,:) = w_min' * S(trial*2-1:trial*2,:);
    figure
    plot(step, portfolio(trial*2-1,:))
    hold on
    plot(step, portfolio(trial*2,:))
    xlabel('time')
    ylabel('Simulated Portfolio result')
    legend('Portfolio with max volatility', 'Portfolio with min volatility')
    title('Portfolio with two stocks')
end
%Through couple iterations, indeed at most cases, the value of portfolio
%with max volatility fluctuates a lot more often than the portfolio with
%min volatility
%% Functions
%Function for simulating geometric brownian motion used in part a and b
function path = simulatePath(L,stepSize,variance)
    path = ones(1,L);
    alpha = 0.1/L;
    sigma = alpha;
    for pos = 2:size(path,2)
        path(pos) = path(pos-1) + alpha*path(pos-1)*stepSize + sigma*path(pos-1)*normrnd(0,sqrt(variance));
        % basically it is x(t) = x(t-1) + dx(t-1)
    end
    return
end

%Function for simulating discounted price used in part c
function OneKPath = simulateDiscouted(start,L,variance,interest,stepSize)
    OneKPath = zeros(1000,L);
    sigma = 1/L;
    %Make all trials begin with the given starting point
    OneKPath(1:1000) = start; 
    %Start of simulation
    for trial = 1:1000
        for pos = 2:L
            %Simulation of Xt involving the drift(interest rate) and
            %randomness(dWt)
            OneKPath(trial,pos) = OneKPath(trial,pos-1) * ...
                                         (1+sigma*normrnd(0,sqrt(variance))) * exp(interest*stepSize);
        end
    end
    return
end

%Function to generate a set of M gaussian random vectors used in Problem 3
%Output is a matrix of size M*length
%mean vector is assumed to be a column vector
function matrix = generateCorreGaussRV(meanVector, covariance, m, length)
    L = chol(covariance, 'lower');
    matrix = normrnd(0,1,m,length);
    for col = 1:length
        matrix(:,col) = meanVector + L * matrix(:,col);
    end
    return
end