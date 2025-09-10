%Copyright---Chen Qiu----2019
function trajectory=generateTrajectory(parameters, topology)
t=parameters.time_step;%time step
T=t*parameters.N;
N=parameters.N+1;
n=parameters.nofnode;
v_start=topology.v_start;
start=topology.trajectory_start;
destination=topology.trajectory_end;
acc=zeros(2,2);
for i=1:n
    for j=1:2
       acc(i,j)=2*(destination(i,j)-start(i,j)-v_start(i,j)*T)/(T^2); 
    end    
end

phi=[1 0 t 0;0,1,0,t;0,0,1,0;0,0,0,1];
tao=[0.5*t^2,0;0,0.5*t^2;t,0;0,t];
trajectory=zeros(4,N,n);
%Copyright---Chen Qiu----2019

for i=1:n
    x=[start(i,1);start(i,2);v_start(i,1);v_start(i,2)];
    for j=1:N
        trajectory(:,j,i)=x;
        x=phi*x+tao*(acc(i,:))';
    end
end

end
%Copyright---Chen Qiu----2019