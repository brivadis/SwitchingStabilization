clear all

A = [0, 1; -1, 0];
B = [0; 1];
K = [0, 2];
alpha = 1;

dim = size(A, 1);

f = @(x, u) A*x+B*u;
h = @(x) 1/2*norm(x)^2;

u1 = @(y, xhat) sqrt(y)+1/sqrt(2)*norm(xhat);
u2 = @(xhat) -K*xhat;

Tobs = sqrt(2)/norm(B);
Tstab = 5;

% Number of switches on each figure:
nb_iter = 1000;
nb_iter_2 = 200;

t0 = 0;
x0 = [1, 0];
a = sqrt(2)/2;
xhat0 = [a, sqrt(2*h(x0) - a^2)];
% x0 = [sqrt(2), 0];
% xhat0 = [1, 1];
% x0 = [1, 0];
% xhat0 = [0.9, 1];
eps0 = xhat0-x0;
etat0 = [xhat0, h(x0), eps0, 0]';

Ttot = [];
Zhat = [];
Eps = [];

Change = [];

for k = 1:nb_iter

    [T, Etat] = ode45(@(t, etat) [...
        f(etat(1:dim), u1(h(etat(1:dim)-etat(dim+2:2*dim+1)), etat(1:dim)))...
        - B*u1(h(etat(1:dim)-etat(dim+2:2*dim+1)), etat(1:dim))*etat(end);...
        u1(h(etat(1:dim)-etat(dim+2:2*dim+1)), etat(1:dim))*B'*etat(1:dim) - alpha*etat(end);...
        A*etat(dim+2:2*dim+1) - u1(h(etat(1:dim)-etat(dim+2:2*dim+1)), etat(1:dim))*B*etat(end);
        u1(h(etat(1:dim)-etat(dim+2:2*dim+1)), etat(1:dim))*B'*etat(dim+2:2*dim+1) - alpha*etat(end)...
        ], [t0, t0+Tobs], etat0);

    Ttot = [Ttot; T];
    Change = [Change; size(Ttot, 1), 0];
    Zhat = [Zhat; Etat(:, 1:dim+1)];
    Eps = [Eps; Etat(:, dim+2:end)];
    
    t0 = Ttot(end);
    etat0 = Etat(end, :)';
    
    
%     [T, Etat] = ode45(@(t, etat) [...
%         f(etat(1:dim), u2(etat(1:dim)));...
%         u2(etat(1:dim))*B'*etat(1:dim);...
%         A*etat(dim+2:2*dim+1);...
%         u2(etat(1:dim))*B'*etat(dim+2:2*dim+1)...
%         ], [t0, t0+Tobs], etat0);

    [T, Etat] = ode45(@(t, etat) [...
        f(etat(1:dim), u2(etat(1:dim)))...
        - B*u1(h(etat(1:dim)-etat(dim+2:2*dim+1)), etat(1:dim))*etat(end);...
        u2(etat(1:dim))*B'*etat(1:dim) - alpha*etat(end);...
        A*etat(dim+2:2*dim+1) - u2(etat(1:dim))*B*etat(end);
        u2(etat(1:dim))*B'*etat(dim+2:2*dim+1) - alpha*etat(end)...
        ], [t0, t0+Tobs], etat0);

    Ttot = [Ttot; T];
    Change(end, end) = size(Ttot, 1);
    Zhat = [Zhat; Etat(:, 1:dim+1)];
    Eps = [Eps; Etat(:, dim+2:end)];
    
    t0 = Ttot(end);
    etat0 = Etat(end, :)';
    etat0(dim+1) = h(Etat(end, 1:dim)-Etat(end, dim+2:2*dim+1));
    etat0(end) = 0;
    
end

X = Zhat - Eps;

%% Figs

close all;

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

fig = figure(1);
plot(X(1:floor(end/nb_iter_2), 1), X(1:floor(end/nb_iter_2), 2), 'k-')
hold on
plot(Zhat(1:floor(end/nb_iter_2), 1), Zhat(1:floor(end/nb_iter_2), 2), 'k--')

scatter(X(Change(1:floor(end/nb_iter_2), 1), 1), X(Change(1:floor(end/nb_iter_2), 1), 2), 'k', 'd')
scatter(Zhat(Change(1:floor(end/nb_iter_2), 1), 1), Zhat(Change(1:floor(end/nb_iter_2), 1), 2), 'k', 'd')
scatter(X([1; Change(1:floor(end/nb_iter_2), 2)], 1), X([1; Change(1:floor(end/nb_iter_2), 2)], 2), 'k', 'o')
scatter(Zhat([1; Change(1:floor(end/nb_iter_2), 2)], 1), Zhat([1; Change(1:floor(end/nb_iter_2), 2)], 2), 'k', 'o')

xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
legend({'$(x_1(t), x_2(t))$', '$(\hat z_1(t), \hat z_2(t))$'}, 'Interpreter', 'latex', 'Location', 'northwest')

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

% set text properties
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(1);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))



Err = sqrt(Eps(:, 1).^2 + Eps(:, 2).^2 + + Eps(:, 3).^2);

fig2 = figure(2);

% set(gca,'yscale','log');
xlabel('t')
ylabel('$L^2$-norm of $\varepsilon$', 'Interpreter', 'latex')
hold on

scatter(Ttot([1; Change(1:floor(end/nb_iter_2), 2)]), Err([1; Change(1:floor(end/nb_iter_2), 2)]), 'k', 'o')
scatter(Ttot(Change(1:floor(end/nb_iter_2), 1)), Err(Change(1:floor(end/nb_iter_2), 1)), 'k', 'd')
plot(Ttot(1:floor(end/nb_iter_2)), Err(1:floor(end/nb_iter_2)), 'k-')

legend({'stabilization $\to$ observation', 'observation $\to$ stabilization'}, 'Interpreter', 'latex')

% scaling
fig2.Units               = 'centimeters';
fig2.Position(3)         = opts.width;
fig2.Position(4)         = opts.height;

% set text properties
set(fig2.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(2);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))


% Norm = sqrt(X(:, 1).^2 + X(:, 2).^2);
Norm = (X(:, 1).^2 + X(:, 2).^2 + Zhat(:, 1).^2 + Zhat(:, 2).^2 + Zhat(:, 3).^2);

% figure
% plot(Ttot, Norm)
% title('state norm')
% set(gca,'yscale','log');

[yupper,ylower] = envelope(Norm, 500, 'peak');
fig3 = figure(3);
P1 = polyfit(Ttot(end/2:end), log10(yupper(end/2:end)), 1);
plot(Ttot, 10^(P1(2)) * 10^(P1(1)).^Ttot, 'k--')
hold on
plot(Ttot, yupper, 'k-')
% hold on
% plot(Ttot, Norm)
% plot(Ttot, ylower)
% title('state norm')
set(gca,'yscale','log');
axis([0 Ttot(end) 0.001 10])
xlabel('t', 'Interpreter', 'latex')
ylabel('$|x|^2 + |\hat z|^2$', 'Interpreter', 'latex')
legend({['$y = '  num2str(10^(P1(2))) ' \times ' num2str(10^(P1(1))) '^t$' ]}, 'interpreter', 'latex')

% scaling
fig3.Units               = 'centimeters';
fig3.Position(3)         = opts.width;
fig3.Position(4)         = opts.height;

% set text properties
set(fig3.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
figure(3);
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))



% find(Norm >= 0.2, 1, 'last');
% find(Norm >= 0.15, 1, 'last');
% find(Norm >= 0.1, 1, 'last');
% find(Norm >= 0.05, 1, 'last');

% 
% Err = sqrt(Eps(:, 1).^2 + Eps(:, 2).^2);
% 
% 
% figure
% plot(Ttot, Err)
% title('erreur')
% 
% 
% figure
% plot(Ttot, X(:, 1))
% hold on
% plot(Ttot, X(:, 2))
% plot(Ttot, Zhat(:, 1))
% plot(Ttot, Zhat(:, 2))


