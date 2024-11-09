clear all;
close all;
clc;

X  = textread('wine.txt');
label_real = textread('wine_label.txt');

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

for l = 1:5
    
%     U = abs(randn(m, c));%P为 N行k列
%     U = U./(sum(U, 2)*ones(1, c)); %sum(P,2)求P每一行的和
    V = zeros(c, n);%P为 N行k列
    for i =1:c
        fd = find(label_real == i);
        if isempty(fd)
           break; 
        end
        V(i,:) = X( randi(numel(fd),1) ,:);
    end
    
      [V,U]=fcm(X,c);
    
    J_prev = inf; J = [];
    
    for i = 1:max_num
        iter = iter + 1;

        if i==1
            V1 = inf;
        else
            V1 = V;
        end
        
%         [U11,~,~]=GMFCM_main(X,label_real,b1);
%         [U22,~,~]=GMFCM_main(X,label_real,b2);
        

         dist = abs(sum(X.*X, 2)*ones(1, c) + (sum(V.*V, 2)*ones(1, m))'-2*X*V'+eps);
        
        
%         t = (dist).^(1/(b1));
%         for j = 1:c
%             t1 = t(:,j);
%             if j == 1
%                 t2 = t1;
%             else
%                 t2 = t2.*t1;%累积
%             end
%         end
%         U1 = (repmat(t2,1,c)).^(1/c)./t;
%         
%         t = (dist).^(1/(b2));
%         for j = 1:c
%             t1 = t(:,j);
%             if j == 1
%                 t2 = t1;
%             else
%                 t2 = t2.*t1;%累积
%             end
%         end
%         U2 = (repmat(t2,1,c)).^(1/c)./t;
        
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
        
        %         Uup=reshape(Uup,m,n,c);
        %         Udown=reshape(Udown,m,n,c);
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
        
        
        %     t3 = (dist./sum(dist,2));
        %
        %     Uup(find(t3<1/3)) = U1(find(t3<1/3)); %#ok<FNDSB>
        %     Uup(find(t3>=1/3)) = U2(find(t3>=1/3)); %#ok<FNDSB>
        %
        %
        %     Udown(find(t3<1/3)) = U2(find(t3<1/3)); %#ok<FNDSB>
        %     Udown(find(t3>=1/3)) = U1(find(t3>=1/3)); %#ok<FNDSB>
        
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
%%
xlswrite('demo.xlsx',[nmi,ari, f,dbi, sc,chi]);
[nmi,ari, f,dbi, sc,chi]
%% 画图
[~, label1] = max(U',[], 1);

label1(label==1) = 1;
label1(label==2) = 3;
label1(label==3) = 2;
% label1(label==4) = 1;
% label1(label==5) = 8;
% label1(label==6) = 4;
% label1(label==7) = 7;
% label1(label==8) = 2;

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

saveas(gcf,'E:\myprogram\FCM\FCM_my\paper2\testresult\IT2GMFCM_3_wine.png');

% label = bestMap(data0,label1);
% MIhat = MutualInfo(data0,label);
% AC = length(find(gnd == label))/length(gnd);
% disp(['IT2FPCM in the NMF space normalized mutual information:',num2str(MIhat)]);
% disp(['IT2FPCM in the NMF space accuracy:',num2str(AC)]);
% disp(['IT2FPCM运行时间:  ',num2str(etime(clock,tt))]);
% fprintf('IT2FPCM运行次数：%.4f\n',t4);
