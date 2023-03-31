% SIR model with fixed duration quarantine simulation using Euler method
clear;
close all;

% Parameters
b = 1/7;             % recovery rate = 7 days
a = 2.3*b;             % infection rate
c = 1/3;             % recovery rate in quarantine = 3 days

% Initial conditions
N = 10000; % population size
I0 = 0.01; % initial proportion of infected individuals
R0 = 0; % initial proportion of recovered individuals
Q0 = 0.01; % initial proportion of quarantined individuals
S0 = 1-I0-R0-Q0; % initial proportion of susceptible individuals

% Simulation time
tspan = [0 200];        % simulation time span
dt = 1/24;              % time step = one hour

% Quarantine durations to simulate
quarantine_durations = 40:-1:1; %starts from 40, decrease by 1, end with 1, 11 graphs in a movie
n_durations = length(quarantine_durations);
q_array = zeros(1, n_durations);

% Store quarantine rates in q_array to plot Economic Cost vs q
for k = 1:n_durations
    q_array(k) = 1/quarantine_durations(k);
end

% Initialize arrays to store the costs
q_cost_array = zeros(1, n_durations);
i_cost_array = zeros(1, n_durations);
total_cost_array = zeros(1, n_durations);

% Quarantine durations to simulate
for i = 1:n_durations
    quarantine_duration = quarantine_durations(i);
    
    if quarantine_duration == 0
        q = 0;  % no quarantine
    else
        q = 1/quarantine_duration;  % quarantine rate
    end

    % Initialize arrays to store the results
    nsteps = round(diff(tspan)/dt);
    T = linspace(tspan(1), tspan(2), nsteps+1);
    S = zeros(nsteps+1, 1);
    I = zeros(nsteps+1, 1);
    R = zeros(nsteps+1, 1);
    Q = zeros(nsteps+1, 1);

    % Set the initial values
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;
    Q(1) = Q0;

    % Use Euler's method to simulate the system of differential equations
    for n = 1:nsteps
        dSdt = -a*S(n)*I(n);
        dIdt = a*S(n)*I(n) - b*I(n) - q*I(n);
        dRdt = b*I(n) + c*Q(n);
        dQdt = q*I(n) - c*Q(n);

        S(n+1) = S(n) + dt*dSdt;
        I(n+1) = I(n) + dt*dIdt;
        R(n+1) = R(n) + dt*dRdt;
        Q(n+1) = Q(n) + dt*dQdt;
    end
    
    % Create a figure for the time series plots
    figure(1)
    % Create a plot for each quarantine duration
    plot(T, S*N, 'b-', T, I*N, 'r-', T, R*N, 'g-', T, Q*N, 'm-', 'LineWidth', 2.0);
    xlabel('Time (days)');
    ylabel('Number of People');
    legend('Susceptible', 'Infected', 'Recovered', 'Quarantined');
    title(sprintf('%d-day Reaction Time\n', quarantine_duration));

    %pause(0.4)

    % Find the peak of the infected curve
    [peak_infected, peak_index] = max(I);
    % Get the time at which the peak occurs
    peak_time = T(peak_index);
    % Display the peak informationq
    fprintf('%d-day quarantine', quarantine_duration);
    fprintf('Peak infected : %d people at t = %0.1f days\n', round(peak_infected*N), peak_time);

    % Calculate total number of recovered individuals after the pandemic
    total_recovered = sum(R)+sum(Q);
    % Display the total number of recovered individuals after the pandemic
    fprintf('Total number of recovered individuals after the pandemic: %d\n', round(total_recovered));

    % Calculate total impairment to the economy
    integralI = sum(I)*dt;
    integralQ = sum(Q)*dt;
    % Calculate the impairment due to infection
    i_cost = integralI*N;
    % Calculate the impairment due to quarantine
    q_cost = integralQ*N;
    total_cost = i_cost+q_cost;

    % Store the costs in the arrays
    i_cost_array(i) = i_cost;
    q_cost_array(i) = q_cost;
    total_cost_array(i) = total_cost;

    % Display the impairment to the economy
    fprintf('Economic Impairment due to Infection: %d\n', round(i_cost));
    fprintf('Economic Impairment due to Quarantine: %d\n', round(q_cost));
    fprintf('I + Q Total Economic Impairment: %d\n', round(total_cost));
    fprintf('\n');
end

figure(2)
%plot(quarantine_duration,q_cost_array,'b-',quarantine_duration,i_cost_array,'r-',quarantine_duration,total_cost_array,'g-','LineWidth',2.0)
plot(q_array,q_cost_array,'b-',q_array,i_cost_array,'r-',q_array,total_cost_array,'g-','LineWidth',2.0)
xlabel('Quarantine Rate q (1/days)');
ylabel('Impairment to Economy');
legend('Quarantine Cost', 'Infection Cost', 'Total Cost');
title(sprintf('Impairment to Economy'));
