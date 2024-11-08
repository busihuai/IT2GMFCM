%function[MIhat,AC,t4,tt]=IT2FPCM(n,D,c,v1,data,data0,gnd)
clc
clear all
close all

%%
I = imread('D:\myprogram\FCM\FCM_my\database\storagetanks35.tif');

[m, n ,p] = size(I);
X = reshape(double(I),m*n,p);
load('D:\myprogram\FCM\FCM_my\paper2\testresult\N44.mat')
label_real = label;
label = [];

% data = importdata('Datasets_all30\iris.txt');%K=1.2
% data0 = importdata('Datasets_all30\iris_label.txt');


%%
c = 3;[n,D]=size(X);
iter = 0;
l=10;
nmi = zeros(l,1);
ari = zeros(l,1);
hi = zeros(l,1);
dbi = zeros(l,1);
sc = zeros(l,1);
chi = zeros(l,1);

for l = 1:1
    u1_u=zeros(n,c);
    u2_u=zeros(n,c);
    t1_t=zeros(n,c);
    t2_t=zeros(n,c);
    u_up=zeros(n,c);
    u_down=zeros(n,c);
    t_up=zeros(n,c);
    t_down=zeros(n,c);
    v1_down=zeros(c,D);
    v1_up=zeros(c,D);
    
%         [V0,U0]=fcm(X,c);
        V0 = zeros(c, D);%P为 N行k列
        for i =1:c
            fd = find(label_real == i);
            if isempty(fd)
               break;
            end
            V0(i,:) = X( randi(numel(fd),1) ,:);
        end
    
    %     V0 = abs(rand(k, n));%P为 N行k列
    %     V0 = V0*255./(sum(V0, 2)*ones(1, n)); %sum(P,2)求P每一行的和
    

    v1=V0;
    
    %% 开始循环计时
    mc1=2;
    mc2=4;
    nc1=2;
    nc2=4;
    e=0.0001;
    t4=0;
    err=0.1;
    tt=clock;
    % gnd = [ones(50,1);ones(50,1)*2;ones(50,1)*3];
    while err>e && t4<100
        v=v1;
        for x=1:n
            
            for k=1:c
                d(x,k)=(norm(X(x,:)-v(k,:)))^2+0.0001;
            end
            tep1=0;
            tep2=0;
            for k=1:c
                tep1=tep1+d(x,k)^(-1/(mc1-1));
                tep2=tep2+d(x,k)^(-1/(mc2-1));
            end
            for k=1:c
                u1_u(x,k)=d(x,k)^(-1/(mc1-1))/tep1;
                u2_u(x,k)=d(x,k)^(-1/(mc2-1))/tep2;
            end
            
        end
        
        for k=1:c
            for x=1:n
                
                d(x,k)=(norm(X(x,:)-v(k,:)))^2+0.0001;
                
            end
            tep3=0;
            tep4=0;
            for x=1:n
                
                tep3=tep3+d(x,k)^(-1/(nc1-1));
                tep4=tep4+d(x,k)^(-1/(nc2-1));
                
            end
            for x=1:n
                
                t1_t(x,k)=d(x,k)^(-1/(nc1-1))/tep3;
                t2_t(x,k)=d(x,k)^(-1/(nc2-1))/tep4;
                
            end
            
        end
        
        
        for x=1:n
            
            for k=1:c
                u_up(x,k)=max(u1_u(x,k),u2_u(x,k));
                u_down(x,k)=min(u1_u(x,k),u2_u(x,k));
            end
            
        end
        
        for x=1:n
            
            for k=1:c
                t_up(x,k)=max(t1_t(x,k),t2_t(x,k));
                t_down(x,k)=min(t1_t(x,k),t2_t(x,k));
            end
            
        end
        
        for k=1:c
            tepp5=0;
            tepp6=0;
            tepp7=0;
            tepp8=0;
            for x=1:n
                
                u(x,k)=(u_down(x,k)+u_up(x,k))/2;
                tepp5=tepp5+(u_down(x,k)+t_down(x,k))^mc1*X(x,:);
                tepp6=tepp6+(u_down(x,k)+t_down(x,k))^mc1;
                tepp7=tepp7+(u_up(x,k)+t_up(x,k))^mc1*X(x,:);
                tepp8=tepp8+(u_up(x,k)+t_up(x,k))^mc1;
                
            end
            v1_down(k,:)=tepp5/(tepp6+0.0001);
            v1_up(k,:)=tepp7/(tepp8+0.0001);
            v1=(v1_down+v1_up)/2;
        end
        err=0;
        for k=1:c
            err=err+(v1(k,1)-v(k,1))^2;
        end
        t4=t4+1;
    end
    disp(['运行时间:  ',num2str(etime(clock,tt))]);
%% 计算性能函数

    [~, label] = max(u',[], 1);

    A = label;
    B = label_real';
    [NMI,AMI] = NMI_AMI(A, B);
    [ARI, RI, MI, ~]=RandIndex(label',label_real);
    [HI,~,~] = compute_f(label,label_real');
    DBI = evalclusters(X,label','DaviesBouldin');
    DBI =DBI.CriterionValues;   
%     SC = evalclusters(X,label,'gap');
%     SC =SC.CriterionValues;
    SC = Eva_Entropy(label,X,D);
    CHI = evalclusters(X,label','CalinskiHarabasz');
    CHI =CHI.CriterionValues;
    
    nmi(l) = NMI;
    ari(l) = ARI;  
    hi(l) = HI;
    dbi(l) = DBI;
    sc(l) = SC;
    chi(l) = CHI;
    
end
nmi = getnum(nmi);
ari = getnum(ari);
hi = getnum(hi);
dbi = getnum(dbi);
sc = getnum(sc);
chi = getnum(chi);

%% 写入文件
xlswrite('demo.xlsx',[nmi,ari,hi,dbi,sc,chi]);
%% 画图
[~, label1] = max(u',[], 1);
[m, n ,p] = size(I);
figure
imshow(uint8(reshape(v(label1, :), m, n, p)))
%imwrite(uint8(reshape(v(label, :), m, n, p)),'D:\myprogram\FCM\FCM_my\paper2\testresult\IT2FPCM_storagetanks35.PNG')
%% gray
[m, n ,p] = size(I);
[~, label] = max(u',[], 1);
figure
label1=zeros(m*n,1);
label1(label==1) = 3;
label1(label==2) = 1;
label1(label==3) = 2;
imshow(reshape(uint8((label1-1)*255./(k-1)),m,n));
% imwrite(reshape(uint8((label1-1)*255./(k-1)),m,n),'D:\myprogram\FCM\FCM_my\paper2\testresult\IT2FPCM_N44.PNG')
