%Copyright---Chen Qiu----2019
function [X_UKF,D_UKF]=UKFpositioning2(parameters, topology, measurements)

N=parameters.N;
X_UKF=zeros(4,N,2);
D_UKF=zeros(4,4,N,2);
nofnode=parameters.nofnode;
nofa=parameters.nofanchor;
%anchor
anchor=topology.anchor;
%initial state
x0=[15,39;15,40;5,5;5,5];
%initial state covariance matrix
D0=zeros(4,4,nofnode);
for i=1:nofnode
    D0(:,:,i)=diag([4,4,1,1]);
end
%measurement noise
v_m=parameters.measurement_noise;
%transition matrix
t=parameters.time_step;
phi=[1,0,t,0;0,1,0,t;0,0,1,0;0,0,0,1];
%process noise covariance matrix
tao=[t^2/2,0;0,t^2/2;t,0;0,t];
v_a=parameters.v_a*eye(2);
Q=tao*v_a*tao';

%parameters of ut
alpha=parameters.alpha;
beta=parameters.beta;
kappa=parameters.kappa;
lamda=alpha^2*(4+kappa)-4;
weight_m=zeros(9,1);
weight_m(1)=lamda/(4+lamda);
weight_m(2:9)=lamda/(2*(4+lamda));
weight_c=zeros(9,1);
weight_c(1)=lamda/(4+lamda)+1-alpha^2+beta;
weight_c(2:9)=lamda/(2*(4+lamda));
%Copyright---Chen Qiu----2019
for i=1:N
    x_0=x0(:,1);
    D_0=D0(:,:,1);
    z=measurements(1,:,i);
    %sigma points
    x_0_s=zeros(4,9);
    x_0_s(:,1)=x_0;
    D_0_root=sqrtm(D_0);
    for k=1:4
        x_0_s(:,k+1)=x_0+sqrt(3+lamda)*D_0_root(:,k);
        x_0_s(:,k+5)=x_0-sqrt(3+lamda)*D_0_root(:,k);
    end
    
    %time update
    x_onestep_s=phi*x_0_s;
    x_onestep=x_onestep_s*weight_m;
    D_onestep=weight_c(1)*(x_0_s(:,1)-x_onestep)*(x_0_s(:,1)-x_onestep)';
    for k=1:8
        D_onestep=D_onestep+weight_c(2)*(x_0_s(:,k)-x_onestep)*(x_0_s(:,k)-x_onestep)';
    end
    D_onestep=D_onestep+Q;
    y_s=zeros(4,9);
    for k=1:4
        for k1=1:9
            y_s(k,k1)=sqrt((x_onestep_s(1,k1)-anchor(k,1))^2+(x_onestep_s(2,k1)-anchor(k,2))^2);
        end
    end
    y_mean=y_s*weight_m;
    D_y=weight_c(1)*(y_s(:,1)-y_mean)*(y_s(:,1)-y_mean)';
    for k=1:8
        D_y=D_y+weight_c(2)*(y_s(:,k)-y_mean)*(y_s(:,k)-y_mean)';
    end
    D_y=D_y+v_m*eye(4);%c_yy
    D_xy=weight_c(1)*(x_0_s(:,1)-x_onestep)*(y_s(:,1)-y_mean)';
    for k=1:8
        D_xy=D_xy+weight_c(2)*(x_0_s(:,k)-x_onestep)*(y_s(:,k)-y_mean)';
    end
    %update
    K=D_xy/D_y;
    x_0=x_onestep+K*(z'-y_mean);
    D_0=D_onestep-K*D_y*K';
    
    X_UKF(:,i,1)=x_0;
    D_UKF(:,:,i,1)=D_0;
    
end

alpha=0.5;
beta=2;
kappa=-1;
lamda=alpha^2*(4+kappa)-4;
tao=[t^2/2,0;0,t^2/2;t,0;0,t];
v_a=100*eye(2);
Q=tao*v_a*tao';

weight_m=zeros(9,1);
weight_m(1)=lamda/(4+lamda);
weight_m(2:9)=lamda/(2*(4+lamda));
weight_c=zeros(9,1);
weight_c(1)=lamda/(4+lamda)+1-alpha^2+beta;
weight_c(2:9)=lamda/(2*(4+lamda));
%Copyright---Chen Qiu----2019
for i=1:N
    x_0=x0(:,2);
    D_0=D0(:,:,2);
    z=measurements(2,:,i);
    %sigma points
    x_0_s=zeros(4,9);
    x_0_s(:,1)=x_0;
    D_0_root=sqrtm(D_0);
    for k=1:4
        x_0_s(:,k+1)=x_0+sqrt(3+lamda)*D_0_root(:,k);
        x_0_s(:,k+5)=x_0-sqrt(3+lamda)*D_0_root(:,k);
    end
    
    %time update
    x_onestep_s=phi*x_0_s;
    x_onestep=x_onestep_s*weight_m;
    D_onestep=weight_c(1)*(x_0_s(:,1)-x_onestep)*(x_0_s(:,1)-x_onestep)';
    for k=1:8
        D_onestep=D_onestep+weight_c(2)*(x_0_s(:,k)-x_onestep)*(x_0_s(:,k)-x_onestep)';
    end
    D_onestep=D_onestep+Q;
    y_s=zeros(4,9);
    for k=1:4
        for k1=1:9
            y_s(k,k1)=sqrt((x_onestep_s(1,k1)-anchor(k,1))^2+(x_onestep_s(2,k1)-anchor(k,2))^2);
        end
    end
    y_mean=y_s*weight_m;
    D_y=weight_c(1)*(y_s(:,1)-y_mean)*(y_s(:,1)-y_mean)';
    for k=1:8
        D_y=D_y+weight_c(2)*(y_s(:,k)-y_mean)*(y_s(:,k)-y_mean)';
    end
    D_y=D_y+v_m*eye(1);%c_yy
    D_xy=weight_c(1)*(x_0_s(:,1)-x_onestep)*(y_s(:,1)-y_mean)';
    for k=1:8
        D_xy=D_xy+weight_c(2)*(x_0_s(:,k)-x_onestep)*(y_s(:,k)-y_mean)';
    end
    %update
    K=D_xy/D_y;
    x_0=x_onestep+K*(z'-y_mean);
    D_0=D_onestep-K*D_y*K';
    
    X_UKF(:,i,2)=x_0;
    D_UKF(:,:,i,2)=D_0;
    
end

end
%Copyright---Chen Qiu----2019