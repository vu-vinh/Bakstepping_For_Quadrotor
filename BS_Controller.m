for i=1:N-1 
    [ti,X] = ode45(@BS,[t0, t0+Ts], x0); 
    
    psi_d = 0;
    phi_d = asin((ux*sin(psi_d) - uy*cos(psi_d)) / sqrt(ux^2 + uy^2 + (uz + g)^2)); 
    theta_d = atan2(ux*cos(psi_d) + uy*sin(psi_d),uz+g); 

    a_d(:,i+1) = [phi_d; theta_d; psi_d]; 
    a_d_dot(1:2,i+1) = (a_d(1:2,i+1) - a_d(1:2,i))/Ts;
    a_d_2dot(1:2,i+1) = (a_d_dot(1:2,i+1) - a_d_dot(1:2,i))/Ts;

    u(:,i+1) = [T; tau_phi; tau_theta; tau_psi]; 
    w(:,i+1) = sqrt(max(0, M*[T; tau_phi; tau_theta; tau_psi])); 

    xout(:,i+1) = X(end,:)';
    x0 = X(end,:)'; t0 = ti(end);
end

%% plot2

x = xout(1,:);
y = xout(2,:);
z = xout(3,:);

phi = xout(4,:);
theta = xout(5,:);
psi = xout(6,:);

phi_d = a_d(1,:);
theta_d = a_d(2,:);
psi_d = a_d(3,:);

% Tracking errors: desired minus actual
ex = xd - x;
ey = yd - y;
ez = zd - z;
ephi = phi_d - phi;
etheta = theta_d - theta;
epsi = psi_d - psi;

% Save each disturbance channel for plotting
dx = d(1,:);
dy = d(2,:);
dz = d(3,:);
dphi = d(4,:);
dtheta = d(5,:);
dpsi = d(6,:);

save('BS.mat','t','x','y','z','phi','theta','psi', ...
     'xd','yd','zd','phi_d','theta_d','psi_d', ...
     'ex','ey','ez','ephi','etheta','epsi', ...
     'dx','dy','dz','dphi','dtheta','dpsi');

function deq = BS(t, X)
    global m g Ixx Iyy Izz Ir k b l Ax Ay Az
    global M d w i p_d p_d_dot p_d_2dot a_d a_d_dot a_dot a_d_2dot
    global T tau_phi tau_theta tau_psi ux uy uz Ts 
     
    p=X(1:3);
    e_p = p_d(:,i) - p;
    p_dot = X(7:9);
    e_p_dot = p_d_dot(:,i) - p_dot;

    a=X(4:6);
    e_a = a_d(:,i) - a;
    a_dot = X(10:12);
    e_a_dot = a_d_dot(:,i) - a_dot;

    x = p(1);
    y = p(2);
    z = p(3);
    x_dot = p_dot(1);
    y_dot = p_dot(2);
    z_dot = p_dot(3);

    phi = a(1);
    theta = a(2);
    psi = a(3);
    phi_dot = a_dot(1);
    theta_dot = a_dot(2);
    psi_dot = a_dot(3);

    w_Gamma = w(1,i) - w(2,i) + w(3,i) - w(4,i);

    %% Position Controller

    K1 = diag([3 3 3]);
    K2 = diag([4 4 4]);
    Kf = diag([Ax Ay Az]);

    e1 = e_p_dot + K1*e_p;

    u = p_d_2dot + (Kf/m)*p_dot + e_p + K1*e_p_dot + K2*e1;

    ux = u(1);
    uy = u(2);
    uz = u(3);

    T = m*sqrt(ux^2 + uy^2 + (uz + g)^2);

    %% Attitude Controller
    
    G = [1/Ixx    sin(phi)*tan(theta)/Iyy     cos(phi)*tan(theta)/Izz;
         0           cos(phi)/Iyy                  -sin(phi)/Izz;
         0     (sin(phi)/cos(theta))/Iyy   (cos(phi)/cos(theta))/Izz];

    F_eta = [phi_dot*theta_dot*tan(theta) + psi_dot*theta_dot/cos(theta);
             -phi_dot*psi_dot*cos(theta);
             phi_dot*theta_dot/cos(theta) + psi_dot*theta_dot*tan(theta)];

    K1 = diag([16 16 5]);
    K2 = diag([20 20 10]);

    alpha = -K1*e_a;
    e1 = e_a_dot - alpha;

    tau = G\(K1*e_a_dot + a_d_2dot + e_a + K2*e1 - F_eta);

    tau_phi = tau(1);
    tau_theta = tau(2);
    tau_psi = tau(3);

    %% Model
    phi_2dot = theta_dot*psi_dot*(Iyy - Izz)/Ixx - theta_dot*w_Gamma*Ir/Ixx + tau_phi/Ixx + d(4,i);

    theta_2dot = phi_dot*psi_dot*(Izz - Ixx)/Iyy + phi_dot*w_Gamma*Ir/Iyy + tau_theta/Iyy + d(5,i);

    psi_2dot = phi_dot*theta_dot*(Ixx - Iyy)/Izz + tau_psi/Izz + d(6,i);

    x_2dot = (T/m)*( cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)) - (Ax/m)*x_dot + d(1,i);

    y_2dot = (T/m)*( cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)) - (Ay/m)*y_dot + d(2,i);

    z_2dot = (T/m)*cos(phi)*cos(theta) - g - (Az/m)*z_dot + d(3,i);

    deq = zeros(12,1);
    deq(1:3) = p_dot;
    deq(4:6) = a_dot;
    deq(7:9) = [x_2dot; y_2dot; z_2dot];
    deq(10:12) = [phi_2dot; theta_2dot; psi_2dot];
end

function y = sat(x, a, b)
    y = min(max(x, a), b);
end