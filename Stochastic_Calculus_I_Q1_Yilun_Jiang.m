% Yilun Jiang
% Problem Set II
% Question 1

% Binomial Asset Pricing Modle

%Some global variables that are subject to play around
S0 = 1; %Initial value of the stock
u = 1.005; %The percentage a stock would go up in one step
r = 0.003; %The continuous interest rate
d = 1.002; %The percentage a stock would go down in one step
N = 100; %number of the time steps
L = 1000; %Length of the simulation
%% main function
%Due to the constraints of Matlab, the implementaion of these function must
%be in the end of the script.
%(a)
%Functions are in the end
%(b)
[riskNeutralP,riskNeutralQ] = riskNeutral(u,r,d);
%(c)
%Functions are in the end
%(d)
%Functions are in the end
%(e)
p1 = 0.4;
p2 = 0.6;
pqtilda = riskNeutral(u,r,d);
ptilda = pqtilda(1);
p = [p1, p2, ptilda];
K = S0 * exp(expectRn(N,p,u,d));

for ii = 1:3 %iterate over three strike prices with p1
    simulationP1(ii,:) = simulation(p(1), d, u, N, L, K(ii)); %First Column is call price, second column is put price
end

for ii = 1:3 %iterate over three strike prices with p2
    simulationP2(ii,:) = simulation(p(2), d, u, N, L, K(ii)); %First Column is call price, second column is put price
end
%(f)
%There is only one risk neutral price, but since we have three different
%strike prices, so we have three estimate for V0. The first column in V0 is
%european call option, the second column is european put option.
[RiskNeutralPrice, V0withStrike1] = RiskNeutralSimulation(u, r, d, N, L, K(1));
[~, V0withStrike2] = RiskNeutralSimulation(u, r, d, N, L, K(2));
[~, V0withStrike3] = RiskNeutralSimulation(u, r, d, N, L, K(3));
%(g)
p = 0.5;
path = binornd(1, 0.5, 3, N);
%The replicating portfolio only cares about risk neutral probability
%We call the wealth function written in part b
delta0 = zeros(3,N);
[xtrial0, deltaTrial0] = wealth(N, K(1), path(1,:), u, r, d);
[xtrial1, deltaTrial1] = wealth(N, K(2), path(2,:), u, r, d);
[xtrial2, deltaTrial2] = wealth(N, K(3), path(3,:), u, r, d);
figure
hold on
n = 1:1:N;
plot(n, deltaTrial0(1,:));
plot(n, deltaTrial0(2,:));
plot(n, deltaTrial0(3,:));
hold off
xlabel('time')
ylabel('delta hedging factor')
legend('max price','call option','put option')
title('trial 1')

figure
hold on
plot(n, deltaTrial1(1,:));
plot(n, deltaTrial1(2,:));
plot(n, deltaTrial1(3,:));
hold off
xlabel('time')
ylabel('delta hedging factor')
legend('max price','call option','put option')
title('trial 2')

figure
hold on
plot(n, deltaTrial2(1,:));
plot(n, deltaTrial2(2,:));
plot(n, deltaTrial2(3,:));
hold off
xlabel('time')
ylabel('delta hedging factor')
legend('max price','call option','put option')
title('trial 3')

figure
hold on
n = 0:1:N;
plot(n, stockPrice(path(1,:), u, d))
plot(n, stockPrice(path(2,:), u, d))
plot(n, stockPrice(path(3,:), u, d))
hold off
xlabel('time')
ylabel('simulated stock price with binomial probabily 1/2')
legend('Simulation 1','Simulation 2','Simulation 3')
title('stock price vs. time')

%(h)
%The first line is with strike price from ptilda, the second line is with
%strike price from p1 and the third line is with strike price from p2
Category = {'StrikePrice';'call option';'put option'};
strikePricePtilda = [K(1);V0withStrike1(1);V0withStrike1(2)];
strikePriceP1 = [K(2);V0withStrike2(1);V0withStrike2(2)];
strikePriceP2 = [K(3);V0withStrike3(1);V0withStrike2(2)];
T = table(Category, strikePricePtilda, strikePriceP1, strikePriceP2)
%% (a)
function expectValue = expectRn(N,p,u,d)
    % n*p is the expected value for binomial(n,p). For the stock price, the
    % np are essentially in the power. We have log(u^np*d^n(1-p))
    expectValue = N.*p*log10(u)+N*(1.-p)*log10(d); 
end

function variance = varianceRn(n,p)
    %np(1-p) is the variance for binomial(n,p). For the stock price, the
    %np(1-p) are essentially in the power. We have log(u^(np(1-p))*d^-(np(1-p)))
    variance = n*p*(1-p)*log10(u/d);
end

%% (b)
function [riskNeutralP,riskNeutralQ] = riskNeutral(u,r,d)
    riskNeutralP = ((1+r)-d)/(u-d);
    riskNeutralQ = 1 - riskNeutralP;
end

%% (c)
% Given the path, return the history of stock price
function stockPriceHistory = stockPrice(path, u, d)
    stockPriceHistory = ones(1,length(path)); %Initialize the stock price history
    for c = 1:length(path)
        if path(c) == 1
            stockPriceHistory(c+1) = stockPriceHistory(c) * u;
        else
            stockPriceHistory(c+1) = stockPriceHistory(c) * d;
        end
    end
end

%calculate the x0 and delta0 with three idfferent portfolios
function [x0, delta0] =wealth(N, K, path, u, r, d)
    stockPriceHistory = stockPrice(path(2:length(path)), u, d);
    [ptilda, qtilda] = riskNeutral(u,r,d); % Get the risk neutral probability
    Vn_Max_head = zeros(1,N);
    Vn_Max_tail = zeros(1,N);
    Vn_Call_head = zeros(1,N);
    Vn_Call_tail = zeros(1,N);
    Vn_Put_head = zeros(1,N);
    Vn_Put_tail = zeros(1,N);
    for c = 1:N
        head_sequence = path(1:c);
        head_sequence(c+1) = 1;%1 represents head
        tail_sequence = path(1:c);
        tail_sequence(c+1) = 0;%0 represents tail
        Vn_Max_head(c) = maxPrice(head_sequence, u, d);
        Vn_Max_tail(c) = maxPrice(tail_sequence, u, d);
        Vn_Call_head(c) = EuropeanCall(head_sequence, u, d, K);
        Vn_Call_tail(c) = EuropeanCall(tail_sequence, u, d, K);
        Vn_Put_head(c) = EuropeanPut(head_sequence, u, d, K);
        Vn_Put_tail(c) = EuropeanPut(tail_sequence, u, d, K);
    end    
        x0_Max= ptilda.*Vn_Max_head ./(1+r) + qtilda.*Vn_Max_tail ./(1+r);
        x0_Call= ptilda.*Vn_Call_head ./(1+r) + qtilda.*Vn_Call_tail ./(1+r);
        x0_Put= ptilda.*Vn_Put_head ./(1+r) + qtilda.*Vn_Put_tail ./(1+r);
        delta0_Max = (Vn_Max_head - Vn_Max_tail) ./ (stockPriceHistory.*u - stockPriceHistory.*d);
        delta0_Call = (Vn_Call_head - Vn_Call_tail) ./ (stockPriceHistory.*u - stockPriceHistory.*d);
        delta0_Put = (Vn_Put_head - Vn_Put_tail) ./ (stockPriceHistory.*u - stockPriceHistory.*d);
        x0 = [x0_Max; x0_Call; x0_Put];
        delta0 = [delta0_Max; delta0_Call; delta0_Put];
end

%% (d)
%caculate the return based on the path, u and d
function price = maxPrice(path, u, d)
    stockPriceHistory = stockPrice(path, u, d);
    price = max(stockPriceHistory);
end

%caculate the return of European call option based on the path, u and d
function callPrice = EuropeanCall(path, u, d, K)
    stockPriceHistory = stockPrice(path, u, d);
    stockPriceAtN = stockPriceHistory(length(stockPriceHistory));%This is the spot price at N
    if  stockPriceAtN > K
        callPrice = stockPriceAtN - K;
    else
        callPrice = 0;
    end
end

%caculate the return of European put option based on the path, u and d
function putPrice = EuropeanPut(path, u, d, K)
    stockPriceHistory = stockPrice(path, u, d);
    stockPriceAtN = stockPriceHistory(length(stockPriceHistory));%This is the spot price at N
    if  K > stockPriceAtN
        putPrice = K - stockPriceAtN;
    else
        putPrice = 0;
    end
end
%% (e)
%receive one value of p and one value of K, output the simulation result
function simulationResult = simulation(p, d, u, N, L, K) 
    path = binornd(1, p, L, N); %Generating L simulations of length N with binomial probabiliy p
    simulationResult = zeros(1,2);
    for c = 1:L
        %The first item in the simulationResult is the sum of return of
        %European call over L trials
        simulationResult(1) = simulationResult(1) + EuropeanCall(path(c,:), u, d, K); 
        %The second item in the simulationResult is the sum of return of 
        %European put over L trials
        simulationResult(2) = simulationResult(2) + EuropeanPut(path(c,:), u, d, K); 
    end
    simulationResult = simulationResult ./ L;
end
%% (f)
%K is a single number and all other regular parameters. It outputs risk
%neutral price and the corresponding portfolio values at time t = 0.
function [RiskNeutralPrice, V0] = RiskNeutralSimulation(u, r, d, N, L, K)
    [riskNeutralP,~] = riskNeutral(u,r,d);
    path = binornd(1, riskNeutralP, L, N);
    RiskNeutralPrice = 0;
    Vn = [0,0];
    for c = 1:L
        stockPriceHistory = stockPrice(path(c,:), u, d);
        RiskNeutralPrice = RiskNeutralPrice + stockPriceHistory(length(stockPriceHistory));
        %Vn is accumulating [return of European call option, return of European put option]
        %K is the strike price
        Vn = Vn + [EuropeanCall(path, u, d, K), EuropeanPut(path, u, d, K)];
    end
    RiskNeutralPrice = RiskNeutralPrice / L;
    Vn = Vn ./ L;
    V0 = Vn / (1+r)^N;
end
