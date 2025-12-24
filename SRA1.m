clc;
clear;
close all;

[phi_vals, ~, xval_grid] = wavefun('db3', 10);
phi_interp = @(xnum) interp1(xval_grid, phi_vals, xnum, 'pchip', 0);    %db3
W=60;
x_vals = linspace(-4, 4, 500);
t = x_vals;

S = zeros(size(t));
g_1 = @(x)  exp(-x.^2/4) .* sin(2 * pi * x);       %choose any function

for l = -100:100

    arg = W*x_vals - 5*l;

    Theta00 = @(t)-0.0018324229699212648*phi_interp(t) +245.24152109075183*phi_interp(t+4) ...
                    +0.2614259690192384*phi_interp(t+3) -0.4273254019213667*phi_interp(t+2) ...
                    -0.11028510154071818*phi_interp(t+1);

    Theta10 = @(t) -0.012297391131130342*phi_interp(t) -398.2368909789256*phi_interp(t+4) ...
                    +20.0915329951772*phi_interp(t+3) +7.015304108989707*phi_interp(t+2) ...
                    +1.8799924576950435*phi_interp(t+1);

    Theta20 = @(t) 0.022618906374038874*phi_interp(t) +122.92313215171635*phi_interp(t+4) ...
                    +0.867673637473107*phi_interp(t+3) +10.051805230050844*phi_interp(t+2) ...
                    +2.5830241149363533*phi_interp(t+1);

    Theta30 = @(t) 0.8720615623572289*phi_interp(t) +103.93203753708859*phi_interp(t+4) ...
                   -71.36336979810179*phi_interp(t+3) -41.91512204307734*phi_interp(t+2) ...
                    -7.7736416302137945*phi_interp(t+1);

    Theta40 = @(t) 0.11944934536451099*phi_interp(t) -72.8597997936016*phi_interp(t+4) ...
                    +51.14273719699365*phi_interp(t+3) +26.2753381065639*phi_interp(t+2) ...
                    +4.420910159248203*phi_interp(t+1);

    % coefficients
    a = [5, -10, 10, -5, 1] ;
    shifts = [5, 10, 15, 20, 25];

    theta1 = ( a(1)*Theta00(arg-shifts(1)) + a(2)*Theta00(arg-shifts(2)) ...
             + a(3)*Theta00(arg-shifts(3)) + a(4)*Theta00(arg-shifts(4)) ...
             + a(5)*Theta00(arg-shifts(5)) ) * g_1((5*l + 0.02447174)/W);

    theta2 = ( a(1)*Theta10(arg-shifts(1)) + a(2)*Theta10(arg-shifts(2)) ...
             + a(3)*Theta10(arg-shifts(3)) + a(4)*Theta10(arg-shifts(4)) ...
             + a(5)*Theta10(arg-shifts(5)) ) * g_1((5*l + 0.20610737)/W);

    theta3 = ( a(1)*Theta20(arg-shifts(1)) + a(2)*Theta20(arg-shifts(2)) ...
             + a(3)*Theta20(arg-shifts(3)) + a(4)*Theta20(arg-shifts(4)) ...
             + a(5)*Theta20(arg-shifts(5)) ) * g_1((5*l + 0.5)/W);

    theta4 = ( a(1)*Theta30(arg-shifts(1)) + a(2)*Theta30(arg-shifts(2)) ...
             + a(3)*Theta30(arg-shifts(3)) + a(4)*Theta30(arg-shifts(4)) ...
             + a(5)*Theta30(arg-shifts(5)) ) * g_1((5*l + 0.79389263)/W);

    theta5 = ( a(1)*Theta40(arg-shifts(1)) + a(2)*Theta40(arg-shifts(2)) ...
             + a(3)*Theta40(arg-shifts(3)) + a(4)*Theta40(arg-shifts(4)) ...
             + a(5)*Theta40(arg-shifts(5)) ) * g_1((5*l + 0.97552826)/W);

    S = S + theta1 + theta2 + theta3 + theta4 + theta5;
end


% Plot the original function g_1(x)
hold on;
plot(x_vals, g_1(x_vals), 'k', 'LineWidth',1);

% Plot the reconstructed function S(x)
plot(x_vals, S, '--r', 'LineWidth',1);

ylim([-1.5 1.5])
legend('$f$', sprintf('$\\tilde{S}_{%d}f$', W), ...
       'Interpreter', 'latex', ...
       'FontSize', 18);  % Increase font size

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off;

% Compute reconstruction errors
g1_vals = g_1(x_vals);
valid_idx = ~isnan(S) & ~isnan(g1_vals);

S_valid = S(valid_idx);
g1_valid = g1_vals(valid_idx);

max_error = max(abs(S_valid - g1_valid));
mse_error = mean((S_valid - g1_valid).^2);
l2_error = sqrt(sum((S_valid - g1_valid).^2) * (x_vals(2) - x_vals(1)));

% Display error metrics
disp(['Maximum Absolute Error: ', num2str(max_error)]);
disp(['Mean Squared Error: ', num2str(mse_error)]);
disp(['L2 Norm (RMS Error): ', num2str(l2_error)]);
