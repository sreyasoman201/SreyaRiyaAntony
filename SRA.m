clc;
clear;
syms x z t w f(nu)
m = 4; % order of the B-spline

% Define the B-spline Qm(x) of order m
Qm = @(x) (1/factorial(m-1)) * sum(arrayfun(@(j) (-1)^j * nchoosek(m,j) * ...
    ((x-j+abs(x-j))/2)^(m-1), 0:m));

L = 2;
r = 2;
rho = L * r;
A_mat = sym(zeros(rho, rho));
X_values = [0.5,0.75];

% Construct the matrix psi
row_start = 1;

for j = 1:L
    block_rows = r;
    block_cols = rho;
    x_val = X_values(j);
    z_val = exp(2*pi*1i*x_val);
    A_block = sym(zeros(block_rows, block_cols));

    for p = 0:block_rows-1
        for q = 0:block_cols-1
            Qn_p = diff(Qm(x), x, p);
            A_block(p+1, q+1) = subs(Qn_p, x, x_val - q) + ...
                                subs(Qn_p, x, x_val - q + rho) * z;
        end
    end

    A_mat(row_start:row_start + block_rows - 1, 1:rho) = A_block;
    row_start = row_start + r;
end

% Display the symbolic matrix
disp('Symbolic Matrix:');
disp(A_mat);

% Calculate and display the determinant
det_A = simplify(det(A_mat));
disp('Determinant of the matrix:');
disp(det_A);

% Calculate the inverse of the matrix
if det_A ~= 0
    inv_A = simplify(inv(A_mat));
    disp('Inverse of the matrix:');
    disp(inv_A);

    Fourier_inv_A = sym(zeros(size(inv_A)));
    [rA, cA] = size(inv_A);
    for i = 1:rA
        for j = 1:cA
            Fourier_inv_A(i,j) = int(subs(inv_A(i,j), z, exp(2*pi*1i*x)) * ...
                                     exp(-2*pi*1i * w * x), x, 0, 1);
        end
    end

    disp('Fourier transform of Inverse of A_mat:');
    disp(Fourier_inv_A);

    Theta = sym(zeros(L, r));
    for n = 0:L-1
        for p = 0:r-1
            S = 0;
            for q = 0:rho-1
                for nu = -20:20
                    S = S + limit(Fourier_inv_A(q+1, (n*r)+p+1), w, nu) * ...
                            Qm(x - rho*nu - q);
                end
            end
            Theta(n+1, p+1) = S;
        end
    end

    x_vals = linspace(-10, 10, 1000);
    colors = lines(L * r);
    figure;
    hold on;

    legend_entries = {};
    idx = 1;
    for n = 1:L
        for p = 1:r
            y_vals = double(subs(Theta(n, p), x, x_vals));
            plot(x_vals, y_vals, 'LineWidth', 2, 'Color', colors(idx, :));
            legend_entries{end+1} = sprintf('\\Theta_{%d%d}', n-1, p-1); %#ok<SAGROW>
            idx = idx + 1;
        end
    end

    xlabel('x');
    ylabel('\Theta_{np}(x)', 'Interpreter', 'latex');
    title('Plots of \Theta_{np}(x)', 'Interpreter', 'latex');
    legend(legend_entries, 'Interpreter', 'latex');
    grid on;
    hold off;
else
    disp('The matrix is singular and cannot be inverted.');
end

kappa = rho;
X_values = [0.5, 0.75];  % Corresponding to n = 0, 1
l_vals = -10:10;
t_vals = linspace(-5, 5, 500); % Domain of t
j_max = 1;

% Store results: each row for j=0,1; each column for a t-value
Result = zeros(j_max + 1, length(t_vals));

for j = 0:j_max
    for ti = 1:length(t_vals)
        t = t_vals(ti);
        total_sum = 0;

        for i = 0:rho-1
            % Safeguard: binomial(j,i) = 0 if i > j
            if i > j
                bin_coeff = 0;
            else
                bin_coeff = nchoosek(j, i) * factorial(i);
            end

            inner_sum = 0;

            for l = l_vals
                for n = 0:L-1
                    xn = X_values(n+1);
                    shift_term = (xn + kappa*l - t)^(j - i);  % scalar
                    p = mod(i, r);  % column index
                    Theta_val = double(subs(Theta(n+1, p+1), x, t - kappa*l));
                    Theta_val = real(Theta_val);  % remove tiny imaginary parts
                    inner_sum = inner_sum + shift_term * Theta_val;
                end
            end

            total_sum = total_sum + bin_coeff * inner_sum;
        end

        Result(j+1, ti) = total_sum;
    end
end

% Plotting
figure;
plot(t_vals, Result(1,:), 'b', 'LineWidth', 2); hold on;
plot(t_vals, Result(2,:), 'r', 'LineWidth', 2);
yline(1, '--k', 'LineWidth', 1.2);
yline(0, '--k', 'LineWidth', 1.2);

xlabel('t');
ylabel('Summation Value');
title('Numerical Verification of Theta Condition');
legend('\Sigma for j=0', '\Sigma for j=1', '\delta_{j0}', '\delta_{j1}');
grid on;
f_sym = x;  % Change this to any symbolic function, e.g., sin(x), x^2, etc.

max_order = rho - 1;
f_derivs = sym(zeros(1, max_order + 1));
f_handle = cell(1, max_order + 1);

for i = 0:max_order
    f_derivs(i+1) = diff(f_sym, x, i);
    f_handle{i+1} = matlabFunction(f_derivs(i+1), 'Vars', x);
end
f_sym = x;  % Change this to any symbolic function, e.g., sin(x), x^2, etc.

max_order = rho - 1;
f_derivs = sym(zeros(1, max_order + 1));
f_handle = cell(1, max_order + 1);

for i = 0:max_order
    f_derivs(i+1) = diff(f_sym, x, i);
    f_handle{i+1} = matlabFunction(f_derivs(i+1), 'Vars', x);
end

l_vals_SW = -30:30;  % Truncation window (adjust if needed)
S_Wf_vals = zeros(1, length(t_vals));

% Precompute Theta values for speed
Theta_fun = cell(L, r);
for n = 1:L
    for p = 1:r
        Theta_fun{n, p} = matlabFunction(Theta(n, p), 'Vars', x);
    end
end

for ti = 1:length(t_vals)
    t = t_vals(ti);
    total_sum = 0;

epsilon = [4 5 6 7];
a = zeros(1,rho);

for p = 1:rho
    prodv = 1;
    for q = 1:rho
        if q ~= p
            prodv = prodv * epsilon(q)/(epsilon(q)-epsilon(p));
        end
    end
    a(p) = prodv;
end


Theta_tilde = sym(zeros(L,r));

for n = 1:L
    for i = 1:r
        S = 0;
        for p = 1:rho
            S = S + a(p)*subs(Theta(n,i),t,t-epsilon(p));
        end
        Theta_tilde(n,i) = simplify(S);
    end
end
end

f_sym = xs;                % f(x) = x
f_fun = cell(1,r);

for i = 0:r-1
    f_fun{i+1} = matlabFunction(diff(f_sym,xs,i),'Vars',xs);
end


W = 1;                   %scaling vector
Srec = zeros(size(tvals));

for ti = 1:length(tvals)
    tt = tvals(ti);
    S = 0;
    for l = -20:0      
        for i = 0:r-1
            for n = 0:L-1
                fval = f_fun{i+1}((X(n+1)+rho*l)/W);
                thetaval = real(double( ...
                    subs(Theta_tilde(n+1,i+1),t,W*tt-rho*l)));
                S = S + (1/W^i)*fval*thetaval;
            end
        end
    end
    Srec(ti) = S;
end


% Compute original function
f_original_vals = arrayfun(@(t) f_handle{1}(t), t_vals);

% Plot comparison
figure;
plot(t_vals, f_original_vals, 'b--', 'LineWidth', 2); hold on;
plot(t_vals, S_Wf_vals, 'm', 'LineWidth', 2);
xlabel('t'); ylabel('Function Value');
title('Comparison of $f(t)$ and $(S_W f)(t)$', 'Interpreter', 'latex');
legend({'Original $f(t)$', 'Reconstructed $(S_W f)(t)$'}, 'Interpreter', 'latex');
grid on;

% L2 error calculation
squared_diff = (S_Wf_vals - f_original_vals).^2;
L2_error = sqrt(trapz(t_vals, squared_diff));
fprintf('L2 error between reconstructed (S_W f)(t) and original f(t): %.6e\n', L2_error);
