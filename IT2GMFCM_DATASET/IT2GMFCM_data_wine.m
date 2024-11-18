clear all;
close all;
clc;

%% Dataset read, using the document dataset_read
X  = textread('wine.txt');
label_real = textread('wine_label.txt');


%% Parameter settings
[m, n] = size(X);
c = max(label_real); max_num = 200;iter = 0;
N = m*n;
b1 = 2; b2 = 4;b3 = (b1+b2)/2;

%% 
l=10;
nmi = zeros(l,1);
ari = zeros(l,1);
hi = zeros(l,1);
f = zeros(l,1);
dbi = zeros(l,1);
sc = zeros(l,1);
chi = zeros(l,1); 

for l = 1:10
    
    U = abs(randn(m, c));
    U = U./(sum(U, 2)*ones(1, c)); 
    V = zeros(c, n);
    for i =1:c
        fd = find(label_real == i);
        if isempty(fd)
           break; 
        end
        V(i,:) = X( randi(numel(fd),1) ,:);
    end
    
    
    J_prev = inf; J = [];
    
    for i = 1:max_num
        iter = iter + 1;

        if i==1
            V1 = inf;
        else
            V1 = V;
        end
        
        

        dist = abs(sum(X.*X, 2)*ones(1, c) + (sum(V.*V, 2)*ones(1, m))'-2*X*V'+eps);
        
        
        t = (dist).^(1/(b1)); 
        t1 = prod(t,2);
        U1 = t1.^(1/c)./t;

        t = (dist).^(1/(b2)); 
        t2 = prod(t,2);
        U2 = t2.^(1/c)./t;
        
        Uup = max(U1,U2);
        Udown = min(U1,U2);
        
        t3 = (Uup.*Udown).^(1/2);
        t3 = t3.^b3; %uij(t+1)的m次方
        V = (X'*t3)'./(sum(t3)');%V(0)
        
        %% 降维
        xx=zeros(c,n);
        yy=zeros(c,n);
        u_new_R=zeros(m,n,c);
        u_new_L=zeros(m,n,c);
        V_R=zeros(c,n);
        V_L=zeros(c,n);
        v11=V;
        v22=V;
        data1=sort(X,1);
%         data2=reshape(data1,m,n);
        
        
        %UR
        comparision=0;
        while (comparision==0)
            for k=1:c
                for y=1:n
                    for x=1:m-1
                        if v11(k,y)>=data1(x,y) && v11(k,y)<=data1(x+1,y)
                            xx(k,y)=x;
                        end
                    end
                end
            end
            
            for k=1:c
                for y=1:n
                    for x=1:m
                        if x<=xx(k,y)
                            u_new_R(x,y,k)=Udown(x,k);
                        else
                            u_new_R(x,y,k)=Uup(x,k);
                        end
                    end
                end
            end
            
            u_new_R1 = sum(u_new_R,2)/n;         
            u_new_R1 = reshape(u_new_R1,m,c);
            
         
            V_R = (X'*u_new_R1.^b3)'./(sum(u_new_R1.^b3)'*ones(1, n));%V(0)
            for k=1:c
                if sum(sum(roundn(V_R(k,:),-1) == roundn(v11(k,:),-1)))/(n) == 1
                    comparision=1;
                else
                    v11=V_R;
                end
            end
        end
        
        %UL
        comparision=0;
        while (comparision==0)
            for k=1:c
                for y=1:n
                    for x=1:m-1
                        if v22(k,y)>=data1(x,y) && v22(k,y)<=data1(x+1,y)
                            yy(k,y)=x;
                        end
                    end
                end
            end
            
            for k=1:c
                for y=1:n
                    for x=1:m
                        if x<=xx(k,y)
                            u_new_L(x,y,k)=Uup(x,k);
                        else
                            u_new_L(x,y,k)=Udown(x,k);
                        end
                    end
                end
            end
            
            u_new_L1 = sum(u_new_L,2)/n;  
            u_new_L1 = reshape(u_new_L1,m,c);
            
           
            V_L = (X'*u_new_L1.^b3)'./(sum(u_new_L1.^b3)'*ones(1, n));%V(0)
            
            for k=1:c
                if sum(sum(roundn(V_L(k,:),-1) == roundn(v22(k,:),-1)))/(n) == 1
                    comparision=1;
                else
                    v22=V_L;
                end
            end
        end
        
        V = (V_L+V_R)./2;
        U = (u_new_L1.*u_new_R1).^(1/2);
        
        
        
        %%
        disti =  norm(V1-V, 2)^2;
        if disti < 1e-4
            break;
        end
        J_cur = sum(sum((U.^b1).*dist))/N;
        J = [J J_cur];
         fprintf('#iteration: %03d, objective function: %f\n, disti: %f\n', iter, J_cur, disti);
        J_prev = J_cur;

    end
    
    %% Calculation of indexs
    [~, label] = min(dist, [], 2);
    
    A = label;
    B = label_real;
    NMI  = compute_nmi (A, B);
    [ARI, RI, MI, HI]=RandIndex(label',label_real);
    [F,~,~] = compute_f(label,label_real);
    DBI = evalclusters(X,label,'DaviesBouldin');
    DBI =DBI.CriterionValues;
    SC = Eva_Entropy(label,X,1);
    CHI = evalclusters(X,label,'CalinskiHarabasz');
    CHI =CHI.CriterionValues;
    
    
    nmi(l) = NMI;
    ari(l) = ARI;
    hi(l) = HI;
    f(l) = F;
    dbi(l) = DBI;
    sc(l) = SC;
    chi(l) = CHI;
end
nmi = getnum(nmi);
ari = getnum(ari);
hi = getnum(hi);
f = getnum(f);
dbi = getnum(dbi);
sc = getnum(sc);
chi = getnum(chi);
%% Write the results of the indicator calculations in an excel sheet
xlswrite('demo.xlsx',[nmi,ari, f,dbi, sc,chi]);
[nmi,ari, f,dbi, sc,chi]

%% 
figure(3)
gscatter(X(:,1),X(:,2),label,'rgbcmykbrgbcmyk','o+*.xsd^v<>pho+');
box off
l1=legend('Class 1','Class 2','Class 3','Class 4','Class 5','Class 6','Class 7',...
    'Class 8','Class 9','Class 10','Class 11','Class 12','Class 13','Class 14','Class 15');
legend boxoff
set(l1,'Fontname', 'Times New Roman','FontSize',15);
set(gca,'FontSize',20,'FontName','Times New Roman','fontweight','bold');
xlabel('Attribute 1','fontname','times new roman','fontsize',25);
ylabel('Attribute 2','fontname','times new roman','fontsize',25);


% label = bestMap(data0,label1);
% MIhat = MutualInfo(data0,label);
% AC = length(find(gnd == label))/length(gnd);
% disp(['IT2FPCM in the NMF space normalized mutual information:',num2str(MIhat)]);
% disp(['IT2FPCM in the NMF space accuracy:',num2str(AC)]);
% disp(['IT2FPCM运行时间:  ',num2str(etime(clock,tt))]);
% fprintf('IT2FPCM运行次数：%.4f\n',t4);
