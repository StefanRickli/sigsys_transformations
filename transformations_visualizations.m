%close all;

%% User Parameters
N = 16;                                             % Signal length [samples]
w_noise_var = 5;                                    % Variance of Gaussian white noise
noise_filt_N = 3;                                   % No. of samples over which the Gaussian white noise is smoothed
use_saved_seed = false;                             % Use a saved state for the random number generator
n_plot_periods = 2.2;                               % No. of periods (of length N) over which the x-axis shall be plotted
x_time_limits = ceil(n_plot_periods * N) * [-1 1];  % ^^^ belongs to above definition ^^^
y_time_limits = 15*[-1, 1];                         % Defines Y-Axis limits of time plots


%% Generate and plot finite-energy time signal x[n]
if use_saved_seed && exist('seed.mat', 'file')
    load('seed.mat');
    rng(seed);
end

% Gaussian white noise
white_noise = w_noise_var * randn(1, N);

% Signal that is to be analysed, smoothed by a rectangular filter of
% length noise_filt_N, and containing a small sinoid part.
signal = conv(white_noise, ones(1, noise_filt_N) ./ noise_filt_N, 'same') + ...
         3 * sin(3*2*pi*(0:N-1)/N);

% Plot our signal
t = (1:N) - 1;

f1 = figure('name', 'Finite-energy time signal x(k)');

% Plot main signal
stem(t, signal, 'filled', 'b');

% Add zero padding and vertical lines that separate the periods
hold on;
for ii = 2:ceil(n_plot_periods)
    stem(t + N*(ii-1), zeros(1, length(signal)), 'filled', 'b');
end
for ii = 1:ceil(n_plot_periods)
    stem(t - N*ii, zeros(1, length(signal)), 'filled', 'b');
end

xlim(x_time_limits);
ylim(y_time_limits);
add_vert_lines(gca, n_plot_periods, N);
title('Finite-energy time signal');
xlabel('Time [samples]');
legend('x(k)');


%% DTFT of finite-energy signal
f2 = figure('name', 'DTFT of finite-energy signal');

t = 0:N-1;
omegas = 0:0.001:2*pi;
X_dtft = dtft(signal, t, omegas);

% Plot the DTFT of the finite-energy signal
periodic_plot(omegas, abs(X_dtft), n_plot_periods, 2*pi);

title('DTFT of finite-energy signal');
legend('X(\omega)');

%% Nonperiodic autocorrelation of finite-energy signal
f4 = figure('name', 'Nonperiodic autocorrelation');

taus = -1.5*N:1.5*N;
R_x_nonper = nonperiodic_autocorrelation(signal, taus);
stem(taus, R_x_nonper, 'filled', 'b');

title('Nonperiodic autocorrelation');
xlabel('Time Shift [samples]');
legend('R_x(\tau)');


%% Energy Spectral Density S_X(omega) of finite-energy signal x[n]
f3 = figure('name', 'Energy Spectral Density of finite-energy signal');

subplot(2,1,1);
omegas = 0:0.001:2*pi;
periodic_plot(omegas, abs(X_dtft).^2, n_plot_periods, 2*pi);

title('Energy Spectral Density from DTFT');
legend('S_x(\omega)');

subplot(2,1,2);
taus = -1.5*N:1.5*N;
omegas = 0:0.001:2*pi;
periodic_plot(omegas, abs(dtft(R_x_nonper, taus, omegas)), n_plot_periods, 2*pi);

title('Energy Spectral Density from nonperiodic autocorrelation');
legend('S_x(\omega)');


%% Plot periodic signal x_p[n]
f5 = figure('name', 'Periodic signal');

t = (1:N) - 1;
periodic_stem(t, signal, n_plot_periods, N);

xlim(x_time_limits);
ylim(y_time_limits);
add_vert_lines(gca, n_plot_periods, N);
title('Periodic signal');
xlabel('Time [samples]');
legend('x_p(k)');


%% DFT of periodic signal x_p[n]

% Generate the DFT of x[n].
% Note that, as the DFT assumes the signal to be N-periodic, we can give it
% one (the first) period of the nonperiodic signal.
X_dft = fft(signal);

% Add one period of the DFT to the DTFT plot to show that we align
% correctly in X- and Y-Axis.
figure(f2);
hold on;
stem(linspace(0,2*(N-1)/N*pi, N), abs(X_dft), 'r');
legend('X(\omega)', 'X_N(\omega_n)');
f2_ax = gca;

% Make a new figure with the DFT plot alone, this time with periodicity
f6 = figure('name', 'DFT of periodic signal');
f6_ax = axes();

xlim(xlim(f2_ax));
ylim(ylim(f2_ax));

periodic_stem(linspace(0,2*(N-1)/N*pi, N), abs(X_dft), n_plot_periods, 2*pi);

title('DFT of periodic signal');
legend('X_N(\omega_n)');


%% Periodic Autocorrelation
f8 = figure('name', 'Periodic autocorrelation');

tau = 0:N-1;
periodic_stem(tau, periodic_autocorrelation(signal, tau), n_plot_periods, N);

add_vert_lines(gca, 2, N);
title('Periodic autocorrelation');
xlabel('Time Shift [samples]');
legend('R_x(\tau)');


%% Power Spectral Density PSD
f7 = figure('name', 'Power Spectral Density (PSD)');

subplot(2,1,1);
omega = linspace(0,2*(N-1)/N*pi, N);
tau = 0:N-1;
periodic_stem(omega, abs(dtft(periodic_autocorrelation(signal, tau), tau, omega)), n_plot_periods, 2*pi);

add_vert_lines(gca, 2, 2*pi);
title('Power Spectral Density (PSD) from DFT');
legend('\phi_x(\omega_n)');

subplot(2,1,2);
omega = linspace(0,2*(N-1)/N*pi, N);
periodic_stem(omega, abs(X_dft).^2 / N, n_plot_periods, 2*pi);

add_vert_lines(gca, 2, 2*pi);
title('Power Spectral Density (PSD) from periodic autocorrelation');
legend('\phi_x(\omega_n)');


%% Discrete Time Fourier Transformation Function
function X = dtft(x, ns, omegas)
    % Calculates the Discrete Time Fourier Transformation of the signal x
    % according to
    %
    % X(omega) = sum( x[n] * e^( i*omega*n ) ) for n from -oo to oo
    %
    % Arguments:
    % x:        a vector of discrete time samples of a signal
    % ns:       a vector of time stamps at which the samples have been
    %           taken
    % omegas:   a vector of values for which the function is to be
    %           evaluated
    
    assert(isnumeric(x) && isvector(x));
    assert(isnumeric(ns) && isvector(ns));
    assert(isnumeric(omegas));
    assert(length(x) == length(ns));
    
    X = zeros(1, length(omegas));
    for ii = 1:length(omegas)
        for jj = 1:length(x)
            X(ii) = X(ii) + x(jj) * exp( -1i * ns(jj) * omegas(ii) ); 
        end
    end
end


%% Autocorrelation Function for Nonperiodic Signals
function R = nonperiodic_autocorrelation(x, taus)
    % Calculates the autocorrelation for nonperiodic, finite-energy
    % signals according to
    % 
    % R_x[tau] = sum( x[n] * x[n - tau] ) for n from -oo to oo
    %
    % Arguments:
    % x:    a vector of discrete time samples of a signal
    % taus: a vector of values for tau for which the function is to be
    %     	evaluated
    
    assert(isnumeric(x) && isvector(x));
    assert(isnumeric(taus) && isvector(taus));
    
    R = zeros(1, length(taus));
    for ii = taus(1):taus(end)
        if abs(ii) >= length(x)
            R(ii - taus(1) + 1) = 0;
        else
            x_shift = circshift(x, ii);

            if ii >= 0
                temp = x(1+ii:end) .* x_shift(1+ii:end);
            else
                temp = x(1:end+ii) .* x_shift(1:end+ii);
            end

            R(ii - taus(1) + 1) = sum(temp);
        end
    end
end


%% Autocorrelation Function for Nonperiodic Signals
function R = periodic_autocorrelation(x, taus)
    % Calculates the autocorrelation for nonperiodic, finite-energy
    % signals according to
    % 
    % R_x[tau] = 1/N * sum( x[n] * x[n - tau] ) for n from 0 to N-1
    % where x[n] is interpreted as periodic in N samples
    %
    % Arguments:
    % x:    a vector of discrete time samples of a signal
    % taus: a vector of values for tau for which the function is to be
    %     	evaluated
    
    assert(isnumeric(x) && isvector(x));
    assert(isnumeric(taus) && isvector(taus));
    
    R = zeros(1, length(taus));
    for ii = taus(1):taus(end)
            R(ii - taus(1) + 1) = sum(x .* circshift(x, ii));
    end
    
    R = R ./ length(x);
end

%% Plotting helper functions
function add_vert_lines(ax, n_periods, T)
    % Adds n_periods vertical lines at x = n*T except at zero
    
    hold(ax, 'on');
    
    for ii = 2:ceil(n_periods)
        line(ax,  T * (ii - 1) * [1 1], ylim(ax));
    end
    for ii = 1:ceil(n_periods)
        line(ax, -T * (ii - 1) * [1 1], ylim(ax));
    end
end


function periodic_plot(x, y, n_periods, T)
    
    plot(x, y, 'b');
    hold on;

    dimmed_color = [88, 172, 250] ./ 256;
    for ii = 2:ceil(n_periods)
        plot(x + T*(ii-1), y, 'Color', dimmed_color);
    end
    for ii = 1:ceil(n_periods)
        plot(x - T*ii, y, 'Color', dimmed_color);
    end

    xlim(ceil(n_periods * T) * [-1 1]);
    
    if (T - 2*pi < 2*eps)
        add_radiant_xticks(n_periods);
    end
    
    add_vert_lines(gca, n_periods, T);
end


function periodic_stem(x, y, n_periods, T)
    
    stem(x, y, 'filled', 'b');
    hold on;

    dimmed_color = [88, 172, 250] ./ 256;
    for ii = 2:ceil(n_periods)
        stem(x + T*(ii-1), y, 'filled', 'Color', dimmed_color, 'MarkerEdgeColor', dimmed_color, 'MarkerFaceColor', dimmed_color);
    end
    for ii = 1:ceil(n_periods)
        stem(x - T*ii, y, 'filled', 'Color', dimmed_color, 'MarkerEdgeColor', dimmed_color, 'MarkerFaceColor', dimmed_color);
    end

    xlim(ceil(n_periods * T) * [-1 1]);
    
    if (T - 2*pi < 2*eps)
        add_radiant_xticks(n_periods);
    end
    
    add_vert_lines(gca, n_periods, T);
end


function add_radiant_xticks(n_periods)
    tick_arr = [-2*pi, -pi, 0, pi, 2*pi];
    tick_arr_str = {'-2\pi', '-\pi', '0', '\pi', '2\pi'};

    for ii = 2:floor(n_periods)
        tick_arr = [-ii*2*pi, -(2*(ii-1) + 1) * pi, tick_arr, (2*(ii-1) + 1) * pi, ii*2*pi];
        tick_arr_str = {sprintf('-%d\\pi', ii*2), sprintf('-%d\\pi', (2*(ii-1) + 1)), tick_arr_str{:}, sprintf('%d\\pi', (2*(ii-1) + 1)), sprintf('%d\\pi', ii*2)};
    end
    xticks(tick_arr);
    xticklabels(tick_arr_str);
end

