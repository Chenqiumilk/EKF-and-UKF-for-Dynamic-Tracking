%Copyright---Chen Qiu----2019
%% Compare Extended Kalman Filter and Unscented Kalman Filter
%system setting
clear 
parameters.N=200;
parameters.nofanchor=4;
parameters.nofnode=2;
parameters.time_step=0.2; %s
parameters.v_a=10;%driving noise
parameters.measurement_noise=0.1; %m
parameters.alpha=0.6;
parameters.beta=2;
parameters.kappa=-1;
parameters.nofensemble=50;

topology.anchor=[0,0;50,0;0,50;50,50];
topology.trajectory_start=[2,1;40,25];
topology.trajectory_end=[32,31;15,40];
topology.v_start=ones(2,2);
%% 
%generate trajectory and plot
topology.trajectory=generateTrajectory(parameters,topology);
figure(1)
hold on
for i=1:parameters.nofnode
    h(1)=plot(topology.trajectory(1,:,i),topology.trajectory(2,:,i),'b','LineWidth',1.5);
end
h(2)=scatter(topology.anchor(:,1),topology.anchor(:,2),'r*','LineWidth',1.5);
h(3)=scatter(topology.trajectory_start(:,1),topology.trajectory_start(:,2),'o','LineWidth',1.5);
h(4)=scatter(topology.trajectory_end(:,1),topology.trajectory_end(:,2),'p','LineWidth',1.5);

grid on
xlabel('x(m)');
ylabel('y(m)');
title('topology');

%% 
%generate measurements
measurements=generateMeasurements(parameters,topology);

%% 
error_EKF=zeros(4,parameters.N,2);
error_UKF=zeros(4,parameters.N,2);
runtime=zeros(2,1);
for nofe=1:parameters.nofensemble
%positioning using EKF
tic
[X_EKF,D_EKF]=EKFpositioning(parameters,topology,measurements);
runtime(1)=runtime(1)+toc;
%Copyright---Chen Qiu----2019

%% 
%positioning using UKF
tic
[X_UKF,D_UKF]=UKFpositioning2(parameters, topology, measurements);
runtime(2)=runtime(2)+toc;


%% error
for i=1:parameters.nofnode
    error_EKF(:,:,i)=error_EKF(:,:,i)+X_EKF(:,:,i)-topology.trajectory(:,2:end,i);
    error_UKF(:,:,i)=error_UKF(:,:,i)+X_UKF(:,:,i)-topology.trajectory(:,2:end,i);
end

end
runtime=runtime/parameters.nofensemble;

for i=1:parameters.nofnode
    h(5)=plot(X_EKF(1,:,i),X_EKF(2,:,i),'r-.','LineWidth',1.5)
end
for i=1:parameters.nofnode
    h(6)=plot(X_UKF(1,:,i),X_UKF(2,:,i),'k-','lineWidth',1.5);
end
legend(h(1:6),'true trajectories','Anchors','Start position', 'Destination', 'EKF tracking','UKF tracking')
p_e_ekf=zeros(1,200,2);
v_e_ekf=zeros(1,200,2);
p_e_ukf=zeros(1,200,2);
v_e_ukf=zeros(1,200,2);

for i=1:parameters.nofnode
    %position error
    p_e_ekf(:,:,i)=sqrt((error_EKF(1,:,i)/parameters.nofensemble).^2+(error_EKF(2,:,i)/parameters.nofensemble).^2);
    p_e_ukf(:,:,i)=sqrt((error_UKF(1,:,i)/parameters.nofensemble).^2+(error_UKF(2,:,i)/parameters.nofensemble).^2);
    %velocity error
    v_e_ekf(:,:,i)=sqrt((error_EKF(3,:,i)/parameters.nofensemble).^2+(error_EKF(4,:,i)/parameters.nofensemble).^2);
    v_e_ukf(:,:,i)=sqrt((error_UKF(3,:,i)/parameters.nofensemble).^2+(error_UKF(4,:,i)/parameters.nofensemble).^2);
end
%Copyright---Chen Qiu----2019
for i=1:parameters.nofnode
    figure (2+i)   
    subplot(2,1,1)
    title('Ensemble averaged RMSE of position--node1');
    xlabel('time step/(0.2s)');
    ylabel('position error/(m)');
    grid on
    hold on
    g(1)=plot(1:200,p_e_ekf(:,:,i),'LineWidth',1.5);
    g(2)=plot(1:200,p_e_ukf(:,:,i),'LineWidth',1.5);
    legend(g(1:2),'EKF position error', 'UKF position error');
    subplot(2,1,2)
    title('Ensemble averaged RMSE of velocity--node1');
    hold on
    xlabel('time step/(0.2s)');
    ylabel('velocity error/(m/s)');
    g(3)=plot(1:200,v_e_ekf(:,:,i),'LineWidth',1.5);
    g(4)=plot(1:200,v_e_ukf(:,:,i),'LineWidth',1.5);
    grid on
    legend(g(3:4),'EKF position error', 'UKF position error');
    
end
% %% figure2
% figure (2)
% hold on
% grid on
% scatter(topology.anchor(:,1),topology.anchor(:,2),'r*');
% scatter(topology.trajectory_start(:,1),topology.trajectory_start(:,2),'o');
% line([0,40],[0,25]);
% line([50,40],[0,25]);
% line([0,40],[50,25]);
% line([50,40],[50,25]);
%Copyright---Chen Qiu----2019