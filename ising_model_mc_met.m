function [E_avg, L_avg] = ising_model_mc_met(t,m)
N = m*m;

ii = 0;
flag = 1;
reinit = 1;
 
while ii < N*1e3 && flag
%% initialization    
if reinit
    % Create mxm lattice q, calculate E0 from H_tot(q)
    q = 2*(rand(m)>0.5) - 1;
    H_tot = -q.*(circshift(q,[1,0]) + circshift(q,[-1,0]) + ...
            circshift(q,[0,-1]) + circshift(q,[0,1]));
    E0 = sum(H_tot(:));

    % precalculate neighbors for each element in q, collapse last two dims for
    % indexing
    nbrs = zeros(4,2,m,m);
    for i = 1:m;
        % control for edge cases
        l = i-1;
        if l<1;
            l = m;
        end
        r = i+1;
        if r>m;
            r = 1;
        end
        for j = 1:m;
            % more control for edge cases
            u = j-1;
            if u<1;
                u = m;
            end
            d = j+1;
            if j+1>m;
                d = 1;
            end     
            nbrs(:,:,i,j) = [l,j; r,j; i,u; i,d];
        end
    end
    E = [];
    L = [];
    x = 1;
end

%% prepare to plot
figure(1);
clf;
subplot(1,2,1)
title(sprintf('Order and Normalized Energy @ t = %.1f',t))
colormap default
subplot(1,2,2)
title(sprintf('Spin Lattice @ t = %.1f',t))
colormap gray

%% Monte Carlo sweep
% use col_ind to index q1(:) and flip spins
index = randperm(N);
for i = 1:N

% make a copy of q and flip a random spin
q1 = q;
q1(index(i)) = -q1(index(i));

% calculate change in energy
nbrs0 = nbrs(:,:,index(i));
dE = 2*(q(index(i))*(q(nbrs0(1,1),nbrs0(1,2)) + q(nbrs0(2,1),nbrs0(2,2)) + ...
                      q(nbrs0(3,1),nbrs0(3,2)) + q(nbrs0(4,1),nbrs0(4,2))));
    
% decide whether to replace q with q1, if so record energy and order
if dE <= 0 || rand() < exp(-dE/t)
    q = q1;
    Nplus = sum(q(:)>0);
    Nminus = sum(q(:)<0);
    E(x) = E0;
    E0 = E0 + dE;
    L(x) = abs(Nplus - Nminus)/N;
    x = x+1;
end

end
subplot(1,2,1)
plot(1:length(E),E/N,1:length(L),L)
subplot(1,2,2)
imagesc(flipud(q))
caxis([-1 1])
drawnow

%% Conditions for continuing, stopping
reinit = 0;

if ii > 5
    %Nplus = sum(q(:)>0);
    %Nminus = sum(q(:)<0);
    %ordered = Nplus == N || Nminus == N;
    ordered = abs(L(end) - 1) < 1.5e-3;
    L_stable = abs(mean(diff(L(end:-m:end-N)))) < 1e-5;
    E_stable = abs(mean(diff(E(end:-m:end-N)/N))) < 1e-4;
    E_avg = mean(E(end:-m:end-N)/N);
    L_avg = mean(L(end:-m:end-N));
    
    if ordered
        flag = 0;
        E_avg = mean(E(end:-1:end-m)/N);
        L_avg = L(end);
    elseif L_stable && E_stable
        flag = 0;
        % reject metastable states
        if t < 2/asinh(1) && E_avg < -1.7 && L_avg < 0.5
            ii = 0;
            flag = 1;
            reinit = 1;
        end
    end
end

ii = ii + 1;
end