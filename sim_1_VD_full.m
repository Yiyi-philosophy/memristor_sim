
% scale 0.1 time 50     : 65s
% scale 0.1 time 50*2   : 131s
% scale 0.1 time 50*10   : 650s = 11min
% scale 1.0 time 50*10   : 11min* 100

function main  
    clc; clear all;
    parallel = true;
    tic

    if parallel
        % parallel start
        delete(gcp('nocreate'));
        Num_core = feature('numcores');
        disp(['Num_core: ', num2str(Num_core)]);
        parpool(floor(Num_core-2), 'IdleTimeout', 2*480)
    end

    % name = "sim"
    % 三阶网络参数
    tau = 7.5e-4;  
    L = 0.6;  
    C = 68e-9;  
    G0 = 1e-3;  

    % 流控忆阻器参数
    d2 = 100;  
    d0 = 300;  
    a0_prime = -9;  
    a1 = -1e4;  
    b1_prime = 2500;  

    % 环状网络
    scale = 1.0
    time_range = floor(500);
    V_D_sample = floor(40);
    P = floor(30*scale);  
    N = floor(100*scale);  
    M = floor(20*scale); 
    delta = 0.1;  
  
    filename = strcat('neuron_sim_1', ...  
        '-scale_', num2str(scale), ...  
        '-time_range_', num2str(time_range), ...  
        '-V_D_sample_', num2str(V_D_sample), ...  
        '.png');  

    % 初始状态    
    % x1，1个变量
    % iL，N个变量
    % vCi，N个变量
    % xij，N*N个变量
    rng(42); % 固定随机种子
    init_state = -1 + 2*rand(1 + N + N + N*N, 1); % 随机初始化状态  

    % V_D 变量范围  
    V_D_range = linspace(-2.9, -2.1, V_D_sample);  
      
    % 存储SI结果  
    SI_results = zeros(length(V_D_range), 1);  
      
    if parallel
        parfor k = 1:length(V_D_range) 
            display(['proc:', num2str(k), '/', num2str(length(V_D_range))])
            V_D = V_D_range(k);  
            
            % ODE求解  
            [T, Y] = ode45(@(t, y) odefunc(t, y, tau, L, C, G0, V_D, d2, d0, a0_prime, a1, b1_prime, P, N), [0, time_range], init_state);  
            
            % 计算SI  
            SI_results(k) = computeSI(Y, N, M, delta);  
            display(['SI_results: ', num2str(SI_results(k))])
        end  
    else
        for k = 1:length(V_D_range) 
            display(['proc:', num2str(k), '/', num2str(length(V_D_range))])
            V_D = V_D_range(k);  
            
            % ODE求解  
            [T, Y] = ode45(@(t, y) odefunc(t, y, tau, L, C, G0, V_D, d2, d0, a0_prime, a1, b1_prime, P, N), [0, time_range], init_state);  
            
            % 计算SI  
            SI_results(k) = computeSI(Y, N, M, delta);  
            display(['SI_results: ', num2str(SI_results(k))])
        end  
    end
    
      
    % 绘制SI图  
    figure;  
    plot(V_D_range, SI_results, '-o');  
    xlabel('V_D');  
    ylabel('SI');  
    title('Synchronization Index (SI) vs V_D');  

      
    % saveas(gcf, filename, 'png');  
    print(gcf, '-dpng', '-r600',  filename);
    
    % stop parallel
    if parallel
        p = gcp;
        delete(p)
    end

    toc

end  
 
% 这是一个用于解决特定微分方程系统的MATLAB函数。  
% 它定义了系统的状态变量及其随时间变化的动态。  
% 输入参数:  
% t - 当前时间  
% y - 状态变量的向量  
% tau, L, C, G0, V_D, d2, d0, a0_prime, a1, b1_prime, P, N - 模型参数  
% 输出参数:  
% dydt - 状态变量的时间导数向量  
function dydt = odefunc(t, y, tau, L, C, G0, V_D, d2, d0, a0_prime, a1, b1_prime, P, N)  
    % disp("y size:")
    % disp(size(y));  % 这将显示init_state的尺寸 

    % 提取状态变量  
    x1 = y(1);       % 变量x1的当前值  
    iL = y(2:N+1);   % 电流iL的N个组件  
    vCi = y(N+2:2*N+1); % 电容电压vCi的N个组件  
    xij = y(2*N+2:N*N+2*N+1); % 从xij到xNj的N*N个组件  
    % display(size(xij))
    % 初始化dx1dt的微分方程  
    dx1dt = (3 - x1 + mean(vCi)) / tau;  
    % 初始化diLdt的微分方程  
    diLdt = -(V_D - vCi) / L;  
    % 初始化dvCidt和dxijdt的微分方程为零向量  
    dvCidt = zeros(N, 1);  
    dxijdt = zeros(N*N, 1);  
  
    % 遍历所有电容电压vCi  
    for i = 1:N  
        % 求和项初始化为0  
        sumTerm = 0;  
        % 遍历所有邻居节点j  
        for j = 1:N  
            % dist = min(abs(i-j), N-abs(i-j)); % mod
            dist = abs(i-j); % diag
            if i ~= j && dist <= P % 排除自身 % 只考虑P邻域内的节点
                % 确保xij_index在合理范围内  
                xij_index = (i-1)*N + j;  
                % 计算阻抗RB  
                RB = (d2*xij(xij_index)^2 + d0);  
                % 累加邻居节点影响  
                sumTerm = sumTerm + (vCi(j) - vCi(i)) / (2*P*RB);  
                % 计算dxijdt(xij_index)的值，考虑节点i和j间的耦合  
                dxijdt(xij_index) = a1 * (a0_prime + xij(xij_index) + b1_prime * (vCi(j) - vCi(i)) / RB); 
            end  
        end  
        % 计算dvCidt(i)的值，考虑iL、G0、x1、vCi和邻居节点的影响  
        dvCidt(i) = (1/C) * (-iL(i) - G0*x1^2*vCi(i) + sumTerm);  
    end  
    
    % 将所有计算出的微分方程的导数合并为一个向量，作为ode45函数的输出  
    dydt = [dx1dt; diLdt; dvCidt; dxijdt];  
end  
                
  
function SI = computeSI(Y, N, M, delta)  
    % 计算全局平均  
    w1 = Y(:, N+2:2*N+1); % vCi的时间序列  
    w1_avg = mean(w1, 2); % 计算全局平均值  
      
    n = N / M; % 每个箱子的神经元数量  
    sigma = zeros(M, 1); % 初始化局部标准差数组  
      
    for m = 1:M  
        j_start = n*(m-1) + 1;  
        j_end = n*m;  
        % 计算每个箱子的平均值  
        w1_local_avg = mean(w1(:, j_start:j_end), 2);  
        % 计算局部标准差  
        sum_term = sum((w1(:, j_start:j_end) - w1_local_avg).^2, 2);  
        sigma(m) = mean(sqrt(1/n * sum_term));  
    end  
      
    % 计算SI  
    SI = 1 - sum(sigma < delta) / M;  
end  

