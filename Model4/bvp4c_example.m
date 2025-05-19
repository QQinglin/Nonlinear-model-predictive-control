function main
    % 初始猜测解
    solinit = bvpinit(linspace(0, 1, 21), [1 1 1 1 1 1]);
    options = bvpset('Stats', 'on', 'RelTol', 1e-1);

    % 解边界值问题
    sol = bvp4c(@BVP_ode, @BVP_bc, solinit, options);
    t = sol.x;
    y = sol.y;
    
    % 绘制状态变量图
     figure;
    hold on;
    plot(t, y(1, :), '-o', 'DisplayName', 'x_1');
    plot(t, y(2, :), '-x', 'DisplayName', 'x_2');
    plot(t, y(3, :), '-x', 'DisplayName', 'x_3');
    hold off;
    xlabel('Time');
    ylabel('State Variables');
    legend('show');
    title('Solution of the BVP: State Variables');
    
    % 绘制协状态变量图
    % figure;
    % hold on;
    % plot(t, y(4, :), '-s', 'DisplayName', '\lambda_1');
    % plot(t, y(5, :), '-d', 'DisplayName', '\lambda_2');
    % plot(t, y(6, :), '-d', 'DisplayName', '\lambda_2');
    % hold off;
    % xlabel('Time');
    % ylabel('Costate Variables');
    % legend('show');
    % title('Solution of the BVP: Costate Variables');
end

% ODE for states and costates
function dydt = BVP_ode(t, y)
    u = [-y(4) * y(1);
         -y(5) * y(2);
         -y(6) * y(3)];

     u(1) = max(-5,u(1));
     u(1) = min(5,u(1));
     u(2) = max(-5,u(2));
     u(2) = min(5,u(2));
     u(3) = max(-5,u(3));
     u(3) = min(5,u(3));


    dydt = [y(1)*u(1)+y(2);
            y(2)*u(2)+y(1)+y(3);
            +y(3)*u(3)+y(2)+2*y(3);
            -y(1)-y(5)-y(4)*u(1);
            -y(2)-y(4)-y(6)-y(5)*u(2);
            -y(3)-y(5)-2*y(6)-u(3) * y(6)];
end

% 边界条件
function res = BVP_bc(ya, yb)
    res = [ya(1)-1;
           ya(2)-1;
           ya(3)-1;
           yb(4)-yb(1);
           yb(5)-yb(2);
           yb(6)-yb(3)];
end