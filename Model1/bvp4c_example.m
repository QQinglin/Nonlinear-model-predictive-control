function bvp4c_example
    % 初始猜测解
    solinit = bvpinit(linspace(0, 1, 21), [1 1 1 1 1 1 1 1 1 1]);
    options = bvpset('Stats', 'on', 'RelTol', 1e-1);

    % 解边界值问题
    sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
    t = sol.x;
    y = sol.y;

    % 计算控制变量 u
    u = zeros(5, length(t));  % 初始化 u
    for i = 1:length(t)
        u(:, i) = compute_control(y(:, i));
    end
    
    figure;
    data_x_q = t';  % 创建一个列向量，包含时间点
    hold on;
    plot(t, y(1, :), '-o', 'DisplayName', 'x_1');
    plot(t, y(2, :), '-x', 'DisplayName', 'x_2');
    plot(t, y(3, :), '-x', 'DisplayName', 'x_3');
    plot(t, y(4, :), '-x', 'DisplayName', 'x_4');
    plot(t, y(5, :), '-x', 'DisplayName', 'x_5');
    hold off;

% 将状态变量数据按列堆叠到时间列后面
    for i = 1:5
        data_x_q = [data_x_q, y(i,:)'];  % y(i,:) 是一行向量，需转置为列向量
    end

    writematrix(data_x_q, 'C:\Users\24966\Desktop\project\x_q_BVP.csv');

    xlabel('Time');
    ylabel('State Variables');
    legend('show');
    title('Solution of the BVP: State Variables');
    
    % 绘制协状态变量图
    figure;
    data_lambda_q = t';  % 创建一个列向量，包含时间点
    hold on;
    plot(t, y(6, :), '-o', 'DisplayName', '\lambda_1');
    plot(t, y(7, :), '-x', 'DisplayName', '\lambda_2');
    plot(t, y(8, :), '-x', 'DisplayName', '\lambda_3');
    plot(t, y(9, :), '-x', 'DisplayName', '\lambda_4');
    plot(t, y(10, :), '-x', 'DisplayName', '\lambda_5');
    hold off;

% 将状态变量数据按列堆叠到时间列后面
    for i = 6:10
        data_lambda_q = [data_lambda_q, y(i,:)'];  % y(i,:) 是一行向量，需转置为列向量
    end

    writematrix(data_lambda_q, 'C:\Users\24966\Desktop\project\lambda_q_BVP.csv');

    xlabel('Time');
    ylabel('Adjoint Variables');
    legend('show');
    title('Solution of the BVP: Adjoint Variables');

    % 绘制控制变量 u 的图
    figure;
    data_u_q = t';  % 创建一个列向量，包含时间点
    hold on;
    plot(t, u(1, :), '-o', 'DisplayName', 'u_1');
    plot(t, u(2, :), '-x', 'DisplayName', 'u_2');
    plot(t, u(3, :), '-x', 'DisplayName', 'u_3');
    plot(t, u(4, :), '-x', 'DisplayName', 'u_4');
    plot(t, u(5, :), '-x', 'DisplayName', 'u_5');
    hold off;

% 将状态变量数据按列堆叠到时间列后面
    for i = 1:5
        data_u_q = [data_u_q, u(i,:)'];  % y(i,:) 是一行向量，需转置为列向量
    end

    writematrix(data_u_q, 'C:\Users\24966\Desktop\project\u_q_BVP.csv');
    xlabel('Time');
    ylabel('Input Variables');
    legend('show');
    title('Solution of the BVP: Input Variables');
end

% 计算控制变量 u 的函数
function u = compute_control(y)
    u = [-y(6) * y(1);
         -y(7) * y(2);
         -y(8) * y(3);
         -y(9) * y(4);
         -y(10) * y(5)];

    % 控制变量的限制
    for i = 1:5
        u(i) = max(-5, min(5, u(i)));
    end
end

% ODE for states and costates
function dydt = BVP_ode(t, y)
    % 控制变量 u，基于协状态变量 y(6) 到 y(10)
    u = compute_control(y);

    % 定义微分方程
    dydt = zeros(10,1); % 预先定义 dydt 维度为 10x1 列向量
    
    % 状态变量的微分方程
    dydt(1) = y(1) * u(1) + y(2);
    dydt(2) = y(2) * u(2) + y(1) - y(2) + y(3) + y(5);
    dydt(3) = y(3) * u(3) + y(2) + y(3) + y(4);
    dydt(4) = y(4) * u(4) + y(3) + 2 * y(4) + y(5);
    dydt(5) = y(5) * u(5) + y(2) + y(4) + 3 * y(5);

    % 协状态变量的微分方程
    dydt(6) = -y(1) - u(1) * y(6) - y(7);
    dydt(7) = -y(2) - y(6) - (u(2) - 1) * y(7) - y(8) - y(10);
    dydt(8) = -y(3) - y(7) - (u(3) + 1) * y(8) - y(9);
    dydt(9) = -y(4) - y(8) - (u(4) + 2) * y(9) - y(10);
    dydt(10) = -y(5) - y(7) - y(9) - (3 + u(5)) * y(10);
end

% 边界条件
function res = BVP_bc(ya, yb)
    res = [ya(1)-1;
           ya(2)-1;
           ya(3)-1;
           ya(4)-1;
           ya(5)-1;
           yb(6)-yb(1);
           yb(7)-yb(2);
           yb(8)-yb(3);
           yb(9)-yb(4);
           yb(10)-yb(5)];
end
