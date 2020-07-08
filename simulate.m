%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   simulate: simulates the tandem queue and derive AoI and system time   %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

% Switches and basic parametersì
running = false;
min_index = 1000; % Parameter to avoid transitory behavior
stepsize = 0.1;
delta = stepsize : stepsize : 25;

% Check whether to rerun the Monte Carlo
if (running)
    % Simulation parameters
    lambda = 0.5;
    mu=[1.6 1];
    num_packets = 10000000;
    
    % Auxiliary parameters
    alpha = mu - lambda;
    rho = lambda ./ mu;
    
    % System time: W1, S1, W2, S2 for each packet
    system_times = zeros(4, num_packets + min_index);
    
    % Compute origin times at the source
    origin_times = cumsum((exprnd(1 / lambda, 1, num_packets + min_index)));
    
    % The first packet has no queue
    system_times(2, :) = exprnd(1 / mu(1), 1, num_packets + min_index);
    system_times(4, :) = exprnd(1 / mu(2), 1, num_packets + min_index);
    
    % Compute waiting time at each relay
    for i = 2 : num_packets + min_index
        system_times(1, i) = max(0, origin_times(i - 1) + sum(system_times(1 : 2, i - 1)) - origin_times(i));
        system_times(3, i) = max(0,  origin_times(i - 1) + sum(system_times(1 : 4, i - 1)) - (origin_times(i) + sum(system_times(1 : 2, i))));
    end
    save('simulated.mat')
else
    load('simulated.mat')
    num_packets = length(origin_times);
end

departure_times = origin_times + sum(system_times, 1);
aoi_sim = departure_times(min_index + 1 : end) - origin_times(min_index : end - 1);

% Case 1A: both systems are full
ind_a = intersect(find(system_times(1, min_index + 1 : end)), find(system_times(3, min_index + 1 : end)));
aoi_a_sim = aoi_sim(ind_a);
t_a_sim = sum(system_times(:, ind_a + min_index));

% Case 1B: the second system is empty
ind_b = intersect(find(system_times(1, min_index + 1 : end)), find(system_times(3, min_index + 1 : end) == 0));
aoi_b_sim = aoi_sim(ind_b);
t_b_sim = sum(system_times(:, ind_b + min_index));

% Case C: the first system is empty
ind_c = intersect(find(system_times(1, min_index + 1 : end) == 0), find(system_times(3, min_index + 1 : end) > 0));
aoi_c_sim = aoi_sim(ind_c);
t_c_sim = sum(system_times(:, ind_c + min_index));


% Case D: both systems are empty
ind_d = intersect(find(system_times(1, min_index + 1 : end) == 0), find(system_times(3, min_index + 1 : end) == 0));
aoi_d_sim = aoi_sim(ind_d);
t_d_sim = sum(system_times(:, ind_d + min_index));

% General case
t_sim = [t_a_sim, t_b_sim, t_c_sim, t_d_sim];
aoi_sim = [aoi_a_sim, aoi_c_sim, aoi_c_sim, aoi_d_sim];
[t_sim_hist, ~] = hist(t_sim, delta);
[aoi_sim_hist, ~] = hist(aoi_sim, delta);

% Compute theoretical results
t_th = system_time(lambda, mu, delta);
aoi_th = peak_aoi(lambda, mu, delta);

% Plot results
f1 = figure(1);
plot(delta, cumsum(t_th) * stepsize, 'b')
hold on
scatter(delta, cumsum(t_sim_hist) / sum(t_sim_hist), 'b')
hold on
xlabel('PAoI \Delta')
ylabel('CDF')
legend('Theoretical system time CDF', 'Empirical system time CDF')

f2 = figure(2);
plot(delta, cumsum(aoi_th) * stepsize, 'b')
hold on
scatter(delta, cumsum(aoi_sim_hist) / sum(aoi_sim_hist), 'b')
hold on
xlabel('PAoI \Delta')
ylabel('CDF')
legend('Theoretical PAoI CDF', 'Empirical PAoI CDF')



