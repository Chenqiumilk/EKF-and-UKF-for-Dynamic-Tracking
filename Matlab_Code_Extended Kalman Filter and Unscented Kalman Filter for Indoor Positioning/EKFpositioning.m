%Copyright---Chen Qiu----2019
function [X_result,D_result]=EKFpositioning(parameters,topology,measurements)
N=parameters.N;
%initial state
x0=[15,3;20,39;5,5;5,5];
%initial state covariance matrix
nofnode=parameters.nofnode;
nofa=parameters.nofanchor;
anchor=topology.anchor;
D0=zeros(4,4,nofnode);
for i=1:nofnode
    D0(:,:,i)=diag([4,4,1,1]);
end
%measurement noise
v_m=parameters.measurement_noise*eye(nofa);
%transition matrix
t=parameters.time_step;
phi=[1,0,t,0;0,1,0,t;0,0,1,0;0,0,0,1];
%process noise covariance matrix
tao=[t^2/2,0;0,t^2/2;t,0;0,t];
v_a=parameters.v_a;
Q=tao*v_a*eye(2)*tao';
%save results
X_result=zeros(4,N,nofnode);
D_result=zeros(4,4,nofnode,N);
x_update=zeros(4,nofnode);
D_update=zeros(4,4,nofnode);
%Copyright---Chen Qiu----2019
for i=1:N
    S=zeros(nofa,nofnode);
    H=zeros(nofa,2,nofnode);
    H_EKF=zeros(nofa,4,nofnode);
    v=zeros(nofa,1);
    for j=1:nofnode
        %previous state and covariance
        x_0=x0(:,j);
        D_0=D0(:,:,j);
        %obtain measurements
        z=measurements(:,:,i);
        % H matrix
        for k=1:nofa
            S(k,j)=sqrt((x_0(1)-anchor(k,1))^2+(x_0(2)-anchor(k,2))^2);
            H(k,:,j)=[(x_0(1)-anchor(k,1))/S(k,j),(x_0(2)-anchor(k,2))/S(k,j)];
        end
        H_EKF(:,:,j)=[H(:,:,j),zeros(4,2)];
        z_new=(z(j,:))';
        %one step prediction
        x_onestep=phi*x_0;
        D_onestep=phi*D_0*phi'+Q;
        %kalman gain
        K=D_onestep*(H_EKF(:,:,j))'/(H_EKF(:,:,j)*D_onestep*(H_EKF(:,:,j))'+v_m);
        %residual
        for k=1:nofa
            v(k)=z_new(k)-sqrt((x_onestep(1)-anchor(k,1))^2+(x_onestep(2)-anchor(k,2))^2);
        end
        %update state vector and covariance matrix
        x_update(:,j)=x_onestep+K*v;
        D_update(:,:,j)=(eye(nofa)-K*H_EKF(:,:,j))*D_onestep*(eye(nofa)-K*H_EKF(:,:,j))'+K*v_m*K';
        %save the results
        X_result(:,i,j)=x_update(:,j);
        D_result(:,:,j,i)=D_update(:,:,j);       
    end
    %go to the next time step
     x0=x_update;
     D0=D_update;
end
end
%Copyright---Chen Qiu----2019