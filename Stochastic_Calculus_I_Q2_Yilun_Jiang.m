% Yilun Jiang
% Problem Set II
% Question 1

% Brownian Motion

%Some global variables that are subject to play around
%When I test run it on my conputer, it keeps busy and I have to stop it. So
%I tuned down the parameter to a smaller scale.(Probably too much for loops
%is slowing down the performance)
N = 10^3; 
L = 1000;
m = [0.5,1,2];
%% main function
%Due to the constraints of Matlab, the implementaion of these function must
%be in the end of the script.
%(a)
%Functions are in the end
%(b)
%Functions are in the end
%(c)
[prob1, condExpect1] = calProb(m(1), N, L); %for M = 0.5
[prob2, condExpect2] = calProb(m(2), N, L); %for M = 1
[prob3, condExpect3] = calProb(m(3), N, L); %for M = 2
%(d)
%The reflected graph may be hidden by the original graph if the prcess
%didn't cross the specified level m
graph(N, 3, m(1));%draw 3 sample paths with m = 0.5
graph(N, 3, m(2));%draw 3 sample paths with m = 1
graph(N, 3, m(3));%draw 3 sample paths with m = 2
%(e)
%(1)
%The reflected graph may be hidden by the original graph if the prcess
%didn't cross either 1 or -1
wallGraph(N, 5, 1, -1);%draw 5 sample paths
%(2)
%It turns out that almost always the path will stay between -1 and 1.
upper = 1;
lower = -1;
[firstRef, avergeNumRef] = calWallProb(N,L,upper,lower);
%(3)
%It turns out that the observation made in (2) is correct. The fact that 
%increasing trials while making the observation more "true" confirms the 
%hypothesis. As you increasing the number of trials in Monte Carlo, due to 
%the law of large number, variance shrinks and estimate gets closer to the 
%true mean.
L = 10^5;
[~, averageNum] = calWallProb(N,L,upper,lower);
averageStay = (L-averageNum)/L; 
%% (a)
%N is the total steps in one simulation.
%L is the total number of simulations you want to generate.
%Since we are using symmetric random walk to approximate Brownian motion,
%the probability of going up and going down are the same, namely 1/2.
%The return of simulation(N,L) is a L by N matrix, each row represents one
%simulation.
function path = simulation(N, L)
    path = binornd(1, 1/2, L, N);
end

%% (b)
%Given a path, a level m and a simulation length N, it determines the first 
%passage time(inf iff it never hits m and the value of the process given the
%simulated path(reflect one time if the process hit the given level m and then 
%leave it alone)
function [passageTime, value, mark] = simulationValue(path, m, N)
    value = zeros(1,length(path));
    step = 1 / sqrt(N);
    mark = 0; %mark is 0 if it does't pass level m and equal to the passage time if passage does occur
    %First generating the process value without reflection
    %The value vector records the process value from n = 0 to n = N
    for c = 2:length(path)
        if path(c) == 1
            value(c) = value(c-1) + step;
            if value(c) > m && mark == 0
                mark = c;
            end
        else
            value(c) = value(c-1) - step;
        end
    end
    
    if mark > 0
        passageTime = mark;
        reflected = reflectedValue(mark, value, m);
        value = reflected;
    else
        passageTime = inf;
    end
end

function reflected = reflectedValue(mark, value, m)
    for c = mark:length(value)
        value(c) = m - (value(c) - m);
    end
    reflected = value;
end
%% (C)
%This function takes in m(the passage level), N(number of steps in one simulation)
%and L(total number of simulations). It returns two things. The first one
%is the probability that the process didn't cross the level and the second
%one is a conditional expectation on the passage time given that the
%process passed the given level
function [prob, condExpect] = calProb(m, N, L)
    path = simulation(N,L);
    totalPassageTime = 0;
    totalFailureTimes = 0; %number of times that processes don't cross the given level
    for c = 1:L
        [passageTime, ~] = simulationValue(path(c,:), m, N);
        if passageTime == inf
            totalFailureTimes = totalFailureTimes + 1;
        else
            totalPassageTime = totalPassageTime + passageTime;
        end
    end
    prob = totalFailureTimes / L;
    condExpect = totalPassageTime / (L - totalFailureTimes);
end

%% (d)
function graph(N, L, m)
    n = (1:1:N)./N;
    path = simulation(N, L);
    figure
    hold on
    for c = 1:L
        [~, value, mark] = simulationValue(path(c,:), m, N);       
        plot(n,value);
        xlabel('Time')
        ylabel('Simulated value')
        if mark > 0
            value = reflectedValue(mark, value, m);
        end
        plot(n,value)      
    end
    legend('original path 1', 'reflected path 1', 'original path 2', 'reflected path 2', 'original path 3', 'reflected path 3')
    hold off
end

%% (e)
%Function doubleBounded takes 3 parameters and return the **exact** number of
%times that a single process hits either the upper or lower bound and the value of
%the process. 
function [value, times] = doubleBounded(value, upper, lower, N)
    mark = 1;
    times = 0;
    while mark < N
        if value(mark) > upper
            value = reflectedValue(mark, value, upper);
            times = times + 1;
        elseif value(mark) < lower
            value = reflectedValue(mark, value, lower);
            times = times + 1;
        end
        mark = mark + 1;
    end
end

%Function wallGraph takes 4 parameters and graph the simulated results
function wallGraph(N, L, upper, lower)
    n = (1:1:N)/N;
    path = simulation(N, L);
    figure
    hold on
    for c = 1:L
        [~, value, mark1] = simulationValue(path(c,:), upper, N);
        [~, ~, mark2] = simulationValue(path(c,:), lower, N);
        %If the process did cross either bound, then we need to fix it
        if mark1  > 0 || mark2 > 0
            [value, ~] = doubleBounded(value, upper, lower, N);
        end
        plot(n,value);
    end
    hold off
    xlabel('time')
    ylabel('Simulated value')
    legend('Trial 1', 'Trial 2', 'Trial 3', 'Trial 4', 'Trial 5')
end

%Function calWallProb takes 4 parameters and return the **average** number of
%times that a process hit either the upper or lower bound and the average 
%time to first reflection.
%It calls function doubleBounded L times and each time when it gets the
%number of crossing for a single process returned from doubleBounded, it
%keeps accumulating the results. When returns get accumulated L times, we
%divide the total times by L to get the average number of crossing for one
%process.
function [firstRef, avergeNumRef] = calWallProb(N,L,upper,lower)
        path = simulation(N, L);
        totalFirstRef = 0;
        totalNumRef = 0;
        totalNonPassTimes = 0; % record the times that process stays within the wall
        for c = 1:L
            [~, value, mark1] = simulationValue(path(c,:), upper, N);
            [~, ~, mark2] = simulationValue(path(c,:), lower, N);
            %test whether the process crosses either side of the wall
            if mark1  > 0 || mark2 > 0 
                [~, times] = doubleBounded(value, upper, lower, N);
                totalNumRef = totalNumRef + times;
                %mark1 represent the first reflection time due to crossing
                %the upper side of the wall, mark2 represents the first
                %reflection time due to crossing the lower side of the wall
                %We pick the smaller one to add to the cumulative result.
                if mark1 > mark2
                    totalFirstRef = totalFirstRef + mark2;
                else
                    totalFirstRef = totalFirstRef + mark1;
                end
            else
                totalNonPassTimes = totalNonPassTimes + 1;
            end
        end
        
        firstRef = totalFirstRef / (L - totalNonPassTimes);
        avergeNumRef = totalNumRef / (L - totalNonPassTimes);
end