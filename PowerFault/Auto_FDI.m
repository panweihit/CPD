% ==============================  Author ============================== 
%   Wei Pan (w.pan11@imperial.ac.uk) 
 
clear all
N=20; %Define number of nodes
delta = .02;  %Define sample interval
tspan = 10;  %Define time span
t1=3;
t2=9;
T_attack1=t1/delta;
T_attack2=t2/delta;
sigma = 1*1e0;  %Define intensity of noise
initial=pi*(1+randn(N,1));    %  initial value: here the initial value is randomly chosen
scale=1;


sparsity=0.24;
W_temp=abs(sprandn(N,N,sparsity));
W_temp=spones(W_temp);

for i=1:N
    if W_temp(i,i)~=0
        W_temp(i,i)=0;
    end
end

% save('W_temp.mat','W_temp')
nnz(W_temp)

% w2 should be small , since it is conductance

% a1=40;b1=60;
% a2=0.015;b2=0.025;
a2=10;b2=100;
a1=1;b1=5;
 weight1 = (a1 + (b1-a1).*rand(N,N)).*W_temp;
 weight2 = (a2 + (b2-a2).*rand(N,N)).*W_temp;
 
weight11=reshape(weight1,N^2,1);
weight12=reshape(weight2,N^2,1);

w_attack11=find(weight11);rnd11=round(1+(size(w_attack11,1)-1)*rand(5,1));
w_attack12=find(weight12);rnd12=rnd11;
weight11(w_attack11(rnd11))=0; 
weight12(w_attack12(rnd12))=0;
attack11=reshape(weight11, N,N);
attack12=reshape(weight12, N,N);
nnz(attack11)

for j=1:1:N
    weight1(j,j)=0;
    weight2(j,j)=0;
    attack11(j,j)=0;
    attack12(j,j)=0;
%     attack21(j,j)=0;
%     attack22(j,j)=0;
    
end

omega=0*rand(N,1);
T=tspan/delta;
T_attack1=t1/delta;
T_attack2=t2/delta;
x=zeros(N,T); % Initialization of state variable
z=zeros(N,T); % Initialization of second derivative
x(:,1)=initial;%  initial value
% y(:,1)=zeros(n,1);%  initial value
%sigma=0.1;
noise=sigma^2*randn(N,T);



%%

X=x;
Z=z;
for t=1:1:T-1
    for i=1:1:N    
        X(i,t+1)=X(i,t)+delta*Z(i,t);
        Z(i,t+1)=Z(i,t)-delta*(Z(i,t)+omega(i)+weight1(i,:)*sin(X(:,t)-ones(N,1)*X(i,t))+weight2(i,:)*cos(X(:,t)-ones(N,1)*X(i,t)));
    end
end


% initial system (no attack)
for t=1:1:T_attack1 -1
    for i=1:1:N     
        x(i,t+1)=x(i,t)+delta*z(i,t);
        z(i,t+1)=z(i,t)-delta*(z(i,t)+omega(i)+weight1(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+weight2(i,:)*cos(x(:,t)-ones(N,1)*x(i,t)))+noise(i,t);
        zz(i,t+1)=weight1(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+weight2(i,:)*cos(x(:,t)-ones(N,1)*x(i,t))+noise(i,t);
        ZZ(i,t+1)=weight1(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+weight2(i,:)*cos(x(:,t)-ones(N,1)*x(i,t));
        output(i,t+1)=zz(i,t+1)-ZZ(i,t+1);
    end
end
% attack 1
for t=T_attack1:1:T
    for i=1:1:N
        x(i,t+1)=x(i,t)+delta*z(i,t);
        z(i,t+1)=z(i,t)-delta*(z(i,t)+omega(i)+attack11(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+attack12(i,:)*cos(x(:,t)-ones(N,1)*x(i,t)))+noise(i,t);
        zz(i,t+1)=attack11(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+attack12(i,:)*cos(x(:,t)-ones(N,1)*x(i,t))+noise(i,t);
        ZZ(i,t+1)=weight1(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+weight2(i,:)*cos(x(:,t)-ones(N,1)*x(i,t));
        output(i,t+1)=zz(i,t+1)-ZZ(i,t+1);
    end
end
% % attack 2
% for t=T_attack2:1:T-1
%     for i=1:1:N
%         x(i,t+1)=x(i,t)+delta*z(i,t);
%         z(i,t+1)=z(i,t)-delta*(z(i,t)+omega(i)+attack21(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+attack22(i,:)*cos(x(:,t)-ones(N,1)*x(i,t)))+noise(i,t);
%         zz(i,t+1)=attack21(i,:)*sin(x(:,t)-ones(N,1)*x(i,t))+attack22(i,:)*cos(x(:,t)-ones(N,1)*x(i,t));
%         ZZ(i,t+1)=weight1(i,:)*sin(x(:,t)-ones(N,1)*X(i,t))+weight2(i,:)*cos(x(:,t)-ones(N,1)*x(i,t));
%         output(i,t+1)=zz(i,t+1)-ZZ(i,t+1);
%     end
% end

%% figure

for i=1:1:N
    XX(i,:)=X(i,:)-floor(X(i,:)./(2*pi))*2*pi;
end
for i=1:1:N
    xx(i,:)=x(i,:)-floor(x(i,:)./(2*pi))*2*pi;
end
figure(1)
plot((0:T_attack1)*delta,xx(mod(rnd11,N),1:T_attack1+1),'g','LineWidth',2); hold on;
plot((T_attack1:T_attack2)*delta,xx(mod(rnd11,N),T_attack1+1:T_attack2+1),'r','LineWidth',2); hold on
figure(2)
plot((0:T_attack1)*delta,XX(mod(rnd11,N),1:T_attack1+1),'g','LineWidth',2); hold on
plot((T_attack1:T_attack2)*delta,XX(mod(rnd11,N),T_attack1+1:T_attack2+1),'r','LineWidth',2); hold on


output_test=output(:,T_attack1:T_attack1+8);
output_test(abs(output_test)<10*sigma)=0;
order=find(output_test(:,end));
weight1-attack11

%%
%  clear all
%   load data1
weight1-attack11
[order,c]=find(weight1-attack11)
order=sort(order)

figure()
% plot((1:T_attack1)*delta,output(order,1:T_attack1),'g','LineWidth',2); hold on
% plot((T_attack1:T)*delta,output(order,T_attack1:T),'LineWidth',2); hold on
plot((1:T)*delta,output(order,1:T),'LineWidth',2); hold on
axis([0 6 -40 40]);
%title('Phase Dynamics of Kuramoto Model'); 
xlabel('time t (s)');
ylabel('$y_i$','Interpreter','LaTex');
legend('bus 2',    'bus 9',    'bus 15', 'bus 18');

figure()
% plot((1:T_attack1)*delta,output(:,1:T_attack1),'g','LineWidth',2); hold on
% plot((T_attack1:T)*delta,output(:,T_attack1:T),'LineWidth',2); hold on
t1=(1:N)';
t2=zeros(N,1);
t2(order)=order;
others=find(t1-t2);


plot((1:T)*delta,output(order,1:T),'LineWidth',2,'LineStyle', '--'); hold on
plot((1:T)*delta,output(others,1:T),'LineWidth',2); hold on
% plot((1:T)*delta,output(order(4),1:T),'LineWidth',2,'LineStyle', '--', 'Color','r'); hold on
% plot((1:T)*delta,output(order(1),1:T),'LineWidth',2,'LineStyle', '--', 'Color','m'); hold on
% plot((1:T)*delta,output(order(3),1:T),'LineWidth',2,'LineStyle', '--', 'Color','y'); hold on
% plot((1:T)*delta,output(order(5),1:T),'LineWidth',2,'LineStyle', '--', 'Color','g'); hold on
% plot((1:T)*delta,output(order(2),1:T),'LineWidth',2,'LineStyle', '--', 'Color','b'); hold on
axis([0 6 -40 40]);
legend('bus 5',    'bus 7',    'bus 11', 'bus 16', 'bus 19');
xlabel('time t (s)');
% ylabel('${e}_i-\hat{e}_i$','Interpreter','LaTex');
ylabel('$y_i$','Interpreter','LaTex');



%% 
% clear all
% load data_0311_1
scale=1e0;

Y_matrix=zeros(T-1,N);
W_matrix=zeros(2*N,N);
PHI_matrix=zeros(T-1,N*(2*N));
for j=1:1:N 
    Y=-scale^2*((Z(j,2:T)-Z(j,1:T-1))/delta+Z(j,1:T-1))';
    PHI=scale*[sin(X(:,1:T-1)-ones(N,1)*X(j,1:T-1)); cos(X(:,1:T-1)-ones(N,1)*X(j,1:T-1))]';
    W=scale*[weight1(j,:), weight2(j,:)]';
    Y_matrix(:,j)=Y;
    PHI_matrix(:,(1+(2*N)*(j-1)):(2*N)*j)=PHI;
    W_matrix(:,j)=W;
end
E=sigma;

T1=T_attack1+1:T_attack1+100;
j=4; norm(Y_matrix(T1,j)-PHI_matrix(T1,(1+(2*N)*(j-1)):(2*N)*j)*W_matrix(:,j))

y_matrix=zeros(T-1,N);
w0_matrix=zeros(2*N,N);
w1_matrix_tmp=zeros(2*N,N);
phi_matrix=zeros(T-1,N*(2*N));
for j=1:1:N 
    y=-scale^2*((z(j,2:T)-z(j,1:T-1))/delta+z(j,1:T-1))';
    phi=scale*[sin(x(:,1:T-1)-ones(N,1)*x(j,1:T-1)); cos(x(:,1:T-1)-ones(N,1)*x(j,1:T-1))]';
    w0=scale*[weight1(j,:), weight2(j,:)]';
    yy=phi*w0;
    w1=scale*[attack11(j,:), attack12(j,:)]';
    y_matrix(:,j)=y;
    yy_matrix(:,j)=yy;
    phi_matrix(:,(1+(2*N)*(j-1)):(2*N)*j)=phi;
    w0_matrix(:,j)=w0;
    w1_matrix_tmp(:,j)=w1;
%     w2_matrix(:,j)=w2;
end
%  y_matrix=(y_matrix-yy_matrix);
y_matrix=zeros(T-1,N);
 y_matrix=output(:,2:end)';
w1_matrix=w1_matrix_tmp-W_matrix;

e=sigma;
T1=T_attack1+1:T_attack1+100;
% T1=1:T_attack1-10;
j=4; norm(y_matrix(T1,j)-phi_matrix(T1,(1+(2*N)*(j-1)):(2*N)*j)*w1_matrix(:,j))



 save data
% (Z(i,2:T)-Z(i,1:T-1))/delta+Z(i,1:T-1)-(omega(i)+weight1(i,:)*sin(X(:,1:T-1)-ones(n,1)*X(i,1:T-1))+weight2(i,:)*cos(X(:,1:T-1)-ones(n,1)*X(i,1:T-1)))
