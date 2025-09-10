%Copyright---Chen Qiu----2019
function measurements=generateMeasurements(parameters, topology)
nofa=parameters.nofanchor;
nofnode=parameters.nofnode;
anchor=topology.anchor;
N=parameters.N;
trajectory=topology.trajectory;
v_m=parameters.measurement_noise;
%Copyright---Chen Qiu----2019
measurements=zeros(nofnode,nofa,N);
for i=1:N
    for j=1:nofnode
        for k=1:nofa
            measurements(j,k,i)=sqrt((trajectory(1,i+1,j)-anchor(k,1))^2+(trajectory(2,i+1,j)-anchor(k,2))^2)+sqrt(v_m)*randn;
        end
    end
end
end
%Copyright---Chen Qiu----2019