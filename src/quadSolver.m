function [ Wout, x] = quadSolver( W, rowSum, columnSum, varargin)
    if((nargin < 3) && (size(W,1) == size(W, 2))); columnSum = rowSum; end
%     if(nargin < 3); columnSum  = []; end
    
    rowSum = double(reshape(rowSum, [], 1));
    columnSum = double(reshape(columnSum, [], 1));
    
    p = inputParser;
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonnan'});
    validVector = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'vector', 'nonnan'});
    addRequired(p, 'W', validNetwork);
    addRequired(p, 'rowSum', validVector);
    addRequired(p, 'columnSum', validVector);
    addParameter(p, 'numRepeat', 1, @isnumeric);
    parse(p, W, rowSum, columnSum, varargin{:});
    param = p.Results;
    
    ignoreRowConstraint = isempty(rowSum);
    ignoreColumnConstraint = isempty(columnSum);
    
    if(~ignoreRowConstraint && (length(rowSum) ~= size(W, 1)))
       error('Length of rowSum must match the number of rows in W.');
    end
    
    if(~ignoreColumnConstraint && (length(columnSum) ~= size(W, 2)))
       error('Length of columnSum must match the number of columns in W.');
    end
    
    useAllEdges = true;
    isNxN = size(W, 1) == size(W, 2);
    for iRepeat = 1:param.numRepeat
        if(iRepeat > 1)
            perm_row = randperm(size(W, 1));
            if(isNxN)
                perm_col = perm_row;
            else
                perm_col = randperm(size(W, 2));
            end
            W = W(perm_row, perm_col);
            rowSum = rowSum(perm_row, :);
            columnSum = columnSum(perm_col, :);
            I_row = zeros(size(perm_row));
            I_row(perm_row) = 1:length(perm_row);
            I_column = zeros(size(perm_col));
            I_column(perm_col) = 1:length(perm_col);
        end
        [m, n] = size(W);
        W = sparse(W);
    %     W = W - diag(diag(W));
%         indices = find(W);
        if(isNxN && ~useAllEdges)
            indices = find(triu(W, 1));
        else
            indices = find(W);
        end
        p = nnz(indices);
        [r, c] = ind2sub(size(W), indices);
%         [r, c] = find(W);
        Aeq1 = sparse(c == sparse(1:n))';
        Aeq2 = sparse(r == sparse(1:m))';
        if(~ignoreColumnConstraint)
            if(~ignoreRowConstraint)
                A = double([Aeq1; Aeq2]);
                Sx = [columnSum; rowSum];
            else
                A = double([Aeq1]);
                Sx = [columnSum];
            end
        else
            if(~ignoreRowConstraint)
                A = double([Aeq2]);
                Sx = [rowSum];
            else
                A = ones(1, p);
                Sx = [p];
            end
        end
    %     A = double([Aeq1; Aeq2]);
    %     Sx = [columnSum; rowSum];
        H = speye(p, p);

        k = size(A, 1);
        Q = [H A'; A sparse(k, k)]; % (p + k) x (p + k)
        Y = [zeros(p, 1); Sx];
%         Q = Q + sparse(1:(p+k), 1:(p+k), 1e-8 + rand(1, p+k) * 1e-5);
%         Y = Y + randn(size(Y)) * 1e-6;
        Q = Q + speye(size(Q))*1e-6;
        
        add_initial_condition = true;
        mldivide_solver = false;
        
        D = sqrt(sum(W, 1));
        x_init = 1 ./ sqrt(D(r) .* D(c))';
        
        %x_init = 1 ./ (D(r) .* D(c))';
%         if(iRepeat > 1)
%            x_init(invalids) = 0;
%         else
%            x_init = 1 ./ (D(r) .* D(c))';
%         end
        x_init(isinf(x_init) | isnan(x_init)) = 1/p;
        gamma = sum(rowSum) / sum(x_init(:));
        x_init = x_init * gamma;
        Winit = sparse(r, c, x_init, m, n);
        
%         for iX = 1:3
%             s1 = sum(Winit, 1);
%             s2 = sum(Winit, 2);
%             sumSqr = sum((s1' - columnSum).^2) + sum((s2 - rowSum).^2)
%             Winit = rowSum .* Winit ./ s1;
%             Winit = columnSum .* Winit ./ s2;
%             Winit = Winit * sum(rowSum) / sum(sum(Winit));
%         end
%         x_init = Winit(indices);
%         fprintf('gamma: %.2f\n', gamma);
        
        
        if(add_initial_condition) % Relaxes the system by adding an estimate for the solution
            Q = [Q; H sparse(p, k)];
            Y = [Y; x_init];
            Q = [Q; ones(1, p) sparse(1, k)];
            Y = [Y; sum(rowSum)];
%             Y = [Y; ones(p, 1)/p];
%             Y = [Y; x_init];
        end
        
%         x_init_q = [x_init; zeros(k, 1)];
        if(iRepeat == 1)
            x_init_q = [x_init; zeros(k, 1)];
            numiter = 1000;
        else
            x_init_q = [x_init; 1*randn(k, 1)];
            numiter = 1000;
        end
        
        override_system = false;
        if(~override_system)
    %         Q = Q + speye
    %         X = Q \ Y;   
            if(mldivide_solver)
                X = Q \ Y;
            else
                X = lsqr(Q, Y, 1e-3, numiter, [], [], x_init_q);
            end
            x = X(1:p);
        else
            Q = A;
            Y = Sx;
            x_init_q = x_init;
            
            Qb = Q' * Y;
            Qm = Q' * Q;
            X = x_init_q;
            
%             tic
%             gradient = Qb - Qm * X;
%             residual = Y - Q * X;
%             alpha = 0.1;
%             X = X - alpha .* residual ./ gradient;
            
%             Q = [Q; H];
%             Y = [Y; x_init_q];
            
%             X = lsqlin(Q,Y,[],[],[],[],sparse(p, 1),[], x_init_q, ...
%                 optimoptions('lsqlin', 'Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIterations', 20));
%             X = lsqnonneg(Q, Y);
%             X = lsqr(Q, Y, 1e-3, 10, [], [], x_init_q);
            disp('abcd');
            toc
            x = X(1:p);

%             x = lsqlin(H, 0, [], [], Aeq, beq, 0, []);
            
%             x = lsqlin(H, d,[] [],Aeq,beq,lb,ub)
        end
        
        invalids = x < 0;
%         x_init_sum = sum(x_init(invalids));
%         invalids = x <= 0;
%         x_init_sum = 0;
%         
%         x(~invalids) = (sum(rowSum(:)) - x_init_sum) ./ sum(x(~invalids));
%         x(invalids) = x_init(invalids);
        x(invalids) = 0;
        
        Wprob = sparse(r, c, x, m, n);
%         Wprob(Wprob < 0) = 0;
%         Wprob(Wprob <= 0) = x_init(Wprob <= 0);
        Wprob = sum(rowSum(:)) * Wprob / sum(sum(Wprob));
        
%         s1 = sum(Wprob, 1);
%         s2 = sum(Wprob, 2);
%         sumSqr = sum((s1' - columnSum).^2) + sum((s2 - rowSum).^2);
%         fprintf('sumSqr: %.1f\n', sumSqr);
        
        if(iRepeat > 1)
            W = W(I_row, I_column);
            rowSum = rowSum(I_row, :);
            columnSum = columnSum(I_column, :);
            Wprob = Wprob(I_row, I_column); 
        end
        if(iRepeat == 1)
           Wout = Wprob ./ param.numRepeat;
        else
           Wout = Wprob ./ param.numRepeat;
%            Wout = Wout + Wprob ./ param.numRepeat;
        end
    end
%     D = sqrt(sum(W, 1));
%     x_init = 1 ./ (D(r) .* D(c))';
%     Dx = sparse(r, c, x_init, size(W, 1), size(W, 2));
%     Dx = (Dx + Dx') * 0.5;
    
    Wout = (Wout + Wout') * 0.5;
    Wout = sum(rowSum) * Wout / sum(sum(Wout));
%     Wout = Dx;
%     Wout = Wout * 0.9 + Dx * 0.1;
    s1 = sum(Wout, 1);
    s2 = sum(Wout, 2);
    sumSqr = sum((s1' - columnSum).^2) + sum((s2 - rowSum).^2);
    fprintf('sumSqr: %.1f\n', sumSqr);
%     Wout = sparse(size(W));
%     Wout(indices) = x;
    
end

