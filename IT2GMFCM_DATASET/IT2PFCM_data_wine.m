%function[MIhat,AC,kk4,tt]=IT2PFCM(n,D,c,v1,data,data0,gnd)
clc
clear all
close all

%%


X  = textread('wine.txt');
label_real = textread('wine_label.txt');


% data = importdata('Datasets_all30\iris.txt');%K=1.2
% data0 = importdata('Datasets_all30\iris_label.txt');


%%
c=max(label_real);[n,D]=size(X);

iter = 0;
l=10;
nmi = zeros(l,1);
ari = zeros(l,1);
hi = zeros(l,1);
dbi = zeros(l,1);
sc = zeros(l,1);
chi = zeros(l,1);

for ll = 1:10
    
    m1=2;
    a=0.5;
    b=0.5;
    u1_u=zeros(n,c);
    u2_u=zeros(n,c);
    u_up=zeros(n,c);
    u_down=zeros(n,c);
    t_up=zeros(n,c);
    t_down=zeros(n,c);
    
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
    
     [~,U0]=fcm(X,c);
    
    
    v1=V0;
    t1_t=U0;
    t2_t=U0;
    
    %% 开始循环计时
    mc1=1.2;
    mc2=3;
    nc1=1.2;
    nc2=3;
    e=0.0001;
    kk4=0;
    err=0.1;
    tt=clock;
    
    % gnd = [ones(50,1);ones(50,1)*2;ones(50,1)*3];
    while err>e && kk4<100
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
            tp1=0.0;
            tp2=0.0;
            for x=1:n
                tp1=tp1+t1_t(k,x)^nc1*d(x,k);
                tp2=tp2+t1_t(k,x)'^nc1;
            end
            yita1(k)=0.5*tp1/(tp2+0.0001);
        end
        
        for k=1:c
            tp3=0.0;
            tp4=0.0;
            for x=1:n
                tp3=tp3+t2_t(k,x)^nc2*d(x,k);
                tp4=tp4+t2_t(k,x)^nc2;
            end
            yita2(k)=0.5*tp3/(tp4+0.0001);
        end
        
        for x=1:n
            for k=1:c
                t1_t(x,k)=(1+(d(x,k)/yita1(k))^(1/(nc1-1)))^(-1);
                t2_t(x,k)=(1+(d(x,k)/yita2(k))^(1/(nc2-1)))^(-1);
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
        
        for x=1:n
            
            for k=1:c
                u1(x,k)=a*u_up(x,k)+b*t_up(x,k);
                u2(x,k)=a*u_up(x,k)+b*t_down(x,k);
                u3(x,k)=a*u_down(x,k)+b*t_up(x,k);
                u4(x,k)=a*u_down(x,k)+b*t_down(x,k);
            end
            
        end
        
        for x=1:n
            for k=1:c
                U_up(x,k)= max(max(u1(x,k),u2(x,k)),max(u3(x,k),u4(x,k)));
                U_down(x,k)=min(min(u1(x,k),u2(x,k)),min(u3(x,k),u4(x,k)));
            end
        end
        
        for k=1:c
            tepp1=0;
            tepp2=0;
            for x=1:n
                
                u(x,k)=(U_down(x,k)+U_up(x,k))/2;
                tepp1=tepp1+u(x,k)*X(x,:);
                tepp2=tepp2+u(x,k);
                
            end
            v1(k,:)=tepp1/(tepp2+0.0001);
        end
        
        
        
        %%降型
        
        
        for l=1:D
            sdata(:,l)=sort(X(:,l));
        end
        
        V11=v1;
        V22=v1;
        comparision=0;
        while (comparision==0)
            for k=1:c
                for l=1:D
                    for x=1:n-1
                        if sdata(x+1,l)>=V11(k,l)&&V11(k,l)>=sdata(x,l)
                            N=x;
                        end
                    end
                    for x=1:n
                        if x<=N
                            u_R(x,k,l)=u_down(x,k);
                        else
                            u_R(x,k,l)=u_up(x,k);
                        end
                    end
                    %计算聚类中
                    tep1=0.0;
                    tep2=0.0;
                    for x=1:n
                        tep1=tep1+sdata(x,l)*u_R(x,k,l);
                        tep2=tep2+u_R(x,k,l);
                    end
                    VR(k,l)=tep1/(tep2+0.0001);
                end
            end
            for k=1:c
                if (floor(VR(k,l)) == floor(V11(k,:)))
                    comparision=1;
                else
                    V11(k,:)=VR(k,l);
                end
            end
        end
        
        
        comparision=0;
        while (comparision==0)
            for k=1:c
                for l=1:D
                    for x=1:n-1
                        if sdata(x+1,l)>=V22(k,l)&&V22(k,l)>=sdata(x,l)
                            N=x;
                        end
                    end
                    for x=1:n
                        if x<=N
                            u_L(x,k,l)=u_up(x,k);
                        else
                            u_L(x,k,l)=u_down(x,k);
                        end
                    end
                    %计算聚类中
                    tep3=0.0;
                    tep4=0.0;
                    for x=1:n
                        tep3=tep3+sdata(x,l)*u_L(x,k,l);
                        tep4=tep4+u_L(x,k,l);
                    end
                    VL(k,l)=tep3/(tep4+0.0001);
                end
            end
            for k=1:c
                if (floor(VL(k,l)) == floor(V22(k,:)))
                    comparision=1;
                else
                    V22(k,:)=VL(k,l);
                end
            end
        end
        
        
        
        for k=1:c
            for l=l:D
                v1(k,l)=(VL(k,l)+VR(k,l)/2);
            end
        end
        
        tp0=0.0;
        for l=1:D
            for k=1:c
                tp0=tp0+(v1(k,l)-v(k,l))^2;
            end
        end
        if tp0<0.000001
            err=0.000001;
        end
        kk4=kk4+1;
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
    
    nmi(ll) = NMI;
    ari(ll) = ARI;  
    hi(ll) = HI;
    dbi(ll) = DBI;
    sc(ll) = SC;
    chi(ll) = CHI;
    
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

label1(label==1) = 2;
label1(label==2) = 3;
label1(label==3) = 1;


figure(5)
gscatter(X(:,1),X(:,2),label1,'rgbcmykbrgbcmyk','o+*.xsd^v<>pho+');
box off
l1=legend('Class 1','Class 2','Class 3','Class 4','Class 5','Class 6','Class 7',...
    'Class 8','Class 9','Class 10','Class 11','Class 12','Class 13','Class 14','Class 15');
legend boxoff
set(l1,'Fontname', 'Times New Roman','FontSize',15);
set(gca,'FontSize',20,'FontName','Times New Roman','fontweight','bold');
xlabel('Attribute 1','fontname','times new roman','fontsize',25);
ylabel('Attribute 2','fontname','times new roman','fontsize',25);

saveas(gcf,'D:\myprogram\FCM\FCM_my\paper2\testresult\1_IT2PFCM_wine.png');

% label = bestMap(data0,label1);
% MIhat = MutualInfo(data0,label);
% AC = length(find(gnd == label))/length(gnd);
% disp(['IT2FPCM in the NMF space normalized mutual information:',num2str(MIhat)]);
% disp(['IT2FPCM in the NMF space accuracy:',num2str(AC)]);
% disp(['IT2FPCM运行时间:  ',num2str(etime(clock,tt))]);
% fprintf('IT2FPCM运行次数：%.4f\n',t4);
