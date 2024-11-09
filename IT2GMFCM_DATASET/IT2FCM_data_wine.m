clear all;
close all;
clc;

X  = textread('wine.txt');
label_real = textread('wine_label.txt');

[m, n] = size(X);
c = max(label_real); max_num = 500;iter = 0;
N = m*n;
b1 = 2; b2 = 4; b3 = (b1+b2)/2;

%%
l=10;
nmi = zeros(l,1);
ari = zeros(l,1);
hi = zeros(l,1);
dbi = zeros(l,1);
sc = zeros(l,1);
chi = zeros(l,1);

for l = 1:10
    
    U = abs(randn(m, c));%P为 N行k列
    U = U./(sum(U, 2)*ones(1, c)); %sum(P,2)求P每一行的和
%     V = zeros(c, n);%P为 N行k列
%     for i =1:c
%         fd = find(label_real == i);
%         if isempty(fd)
%             break;
%         end
%         V(i,:) = X( randi(numel(fd),1) ,:);
%     end
    
    [V,U] = fcm(X,c);
    
    J_prev = inf; J = [];
    
    for i = 1:max_num
        iter = iter + 1;
        
        V1 = V;
        
        dist = sum(X.*X, 2)*ones(1, c) + (sum(V.*V, 2)*ones(1, m))'-2*X*V'+eps;
        
        t1 = (1./dist).^(1/(b1-1));
        U1 = t1./(sum(t1, 2)*ones(1, c));
        
        t2 = (1./dist).^(1/(b2-1));
        U2 = t2./(sum(t2, 2)*ones(1, c));
        
        Uup = max(U1,U2);
        Udown = min(U1,U2);
        
        t3 = ((Uup+Udown)./2).^b3; %uij(t+1)的m次方
        V = (X'*t3)'./(sum(t3)'*ones(1, n));%V(0)
        
        %     Uup=reshape(Uup,m,n,c);
        %     Udown=reshape(Udown,m,n,c);
        %     [U,V] = compute_v(m,n,c,V,X,reshape(Uup,m,n,c),reshape(Udown,m,n,c));
        %     V = V';
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
                for x=1:m
                    for y=1:n
                        if X(x,y)<=data1(xx(k,y),y)
                            u_new_R(x,y,k)=Udown(x,k);
                        else
                            u_new_R(x,y,k)=Uup(x,k);
                        end
                    end
                end
            end
            
            u_new_R1 = sum(u_new_R,2)/n;
            u_new_R1 = reshape(u_new_R1,m,c);
            
            t = u_new_R1.^b3; %uij(0)的m次方
            V_R = (X'*t)'./(sum(t)'*ones(1, n));%V(0)
            
            if sum(sum(roundn(V_R,0) == roundn(v11,0)))/(c*n) == 1
                comparision=1;
            else
                v11=V_R;
            end
        end
        
        %UL
        comparision=0;
        while (comparision==0)
            for k=1:c
                for y=1:n
                    for x=1:m-1
                        if v22(k,y)>data1(x,y) && v22(k,y)<data1(x+1,y)
                            yy(k,y)=x;
                        end
                    end
                end
            end
            
            for k=1:c
                for x=1:m
                    for y=1:n
                        if X(x,y)<=data1(yy(k,y),y)
                            u_new_L(x,y,k)=Uup(x,k);
                        else
                            u_new_L(x,y,k)=Udown(x,k);
                        end
                    end
                end
            end
            
            u_new_L1 = sum(u_new_L,2)/n;
            u_new_L1 = reshape(u_new_L1,m,c);
            
            t = u_new_L1.^b3; %uij(0)的m次方
            V_L = (X'*t)'./(sum(t)'*ones(1, n));%V(0)
            
            if sum(sum(roundn(V_L,0) == roundn(v22,0)))/(c*n) == 1
                comparision=1;
            else
                v22=V_L;
            end
        end
        
        V = (V_L+V_R)./2;
        U = (u_new_L1+u_new_R1)./2;
        
        %%
        disti =  norm(V1-V, 2)^2;
        if disti < 1e-6
            break;
        end
        J_cur = sum(sum((U.^b1).*dist))/N;
        J = [J J_cur];
        fprintf('#iteration: %03d, objective function: %f\n, disti: %f\n', iter, J_cur, disti);
        J_prev = J_cur;
        
    end
    
    %% 画图sw
    [~, label] = min(dist, [], 2);
    
    A = label;
    B = label_real;
    NMI  = compute_nmi (A, B);
    [ARI, RI, MI, ~]=RandIndex(label',label_real);
    [HI,~,~] = compute_f(label,label_real);
    DBI = evalclusters(X,label,'DaviesBouldin');
    DBI =DBI.CriterionValues;
    SC = Eva_Entropy(label,X,1);
    CHI = evalclusters(X,label,'CalinskiHarabasz');
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
%%
xlswrite('demo.xlsx',[nmi,ari, hi, dbi, sc,chi]);
%%
[~, label1] = max(U,[], 2);

label1(label==1) = 1;
label1(label==2) = 3;
label1(label==3) = 2;
% label1(label==4) = 3;
% label1(label==5) = 6;
% label1(label==6) = 4;
% label1(label==7) = 5;
% label1(label==8) = 7;

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
saveas(gcf,'D:\myprogram\FCM\FCM_my\paper2\testresult\1_IT2FCM_wine.png');




