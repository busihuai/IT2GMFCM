clear all;
close all;
clc;

I = imread('N_44.jpg');
X = rgb2gray(I);
X = double(X);
rng('default');%生成随机数
[m, n] = size(X);

load('D:\myprogram\FCM\FCM_my\paper2\testresult\N44.mat')
label_real = label;
label = [];

c = 3; max_num = 200;iter = 0;
N = m*n;
mc = 2;
%N为样本个数 p为维数

%%
% U = randn(N, k);%P为 N行k列
% U = U./(sum(U, 2)*ones(1, k)); %sum(P,2)求P每一行的和
%%
l=10;
nmi = zeros(l,1);
ari = zeros(l,1);
hi = zeros(l,1);
dbi = zeros(l,1);
sc = zeros(l,1);
chi = zeros(l,1);

for l = 1:1
    
    U = randn(m*n, c);%P为 N行k列
    U = U./(sum(U, 2)*ones(1, c)); %sum(P,2)求P每一行的和
    U = reshape(U,m,n,c);
    
    V = abs(rand(c,1));%P为 N行k列
    V = V*255; %sum(P,2)求P每一行的和
    [V,U] = fcm(reshape(X,m*n,1),c);
    
    J_prev = inf; J = [];
    
    for i = 1:max_num
        iter = iter + 1;
        V1 = V;
        s_up=zeros(m,n,c);
        s_down=zeros(m,n,c);
        u1=zeros(m,n,c);
        u2_u=zeros(m,n,c);
        u2_d=zeros(m,n,c);
        u_up=zeros(m,n,c);
        u_down=zeros(m,n,c);
        
        for x=1:m
            for y=1:n
                for k=1:c
                    dt(x,y,k)=(X(x,y)-V(k))^2+0.0001;
                end
                tp1=0.0;
                for k=1:c
                    tp1=tp1+(dt(x,y,k))^(-1/(mc-1));
                end
                for k=1:c
                    u1(x,y,k)=(dt(x,y,k))^(-1/(mc-1))/tp1;
                end
            end
        end
        
        for k=1:c
            for  x=3:m-2
                for y=3:n-2
                    Dx=zeros(8,1);
                    for x1=-1:1
                        for y1=-1:1
                            if x1~=0 || y1~=0
                                st1=0;
                                Dr=0;
                                for x2=-1:1
                                    for y2=-1:1
                                        if x2~=0 || y2~=0
                                            Dr=Dr+sqrt((x2)^2+(y2)^2);
                                            st1=st1+sqrt((x2)^2+(y2)^2)*abs(X(x+x1+x2,y+y1+y2)-V(k));
                                        end
                                    end
                                end
                                sx=3*(x1+1)+(y1+2);
                                if (sx>=6)
                                    sx=sx-1;
                                end
                                Dx(sx)=st1/Dr;
                            end
                        end
                    end
                    s_up(x,y,k)=max(reshape(Dx,numel(Dx),1));
                    s_down(x,y,k)=min(reshape(Dx,numel(Dx),1));
                end
            end
        end
        
        % 算隶属度
        for x=3:m-2
            for y=3:n-2
                for k=1:c
                    dt(x,y,k)=(X(x,y)-V(k))^2+0.0001;
                end
                tp1=0.0;
                for k=1:c
                    tp1=tp1+(dt(x,y,k))^(-1/(mc-1));
                end
                for k=1:c
                    u2_u(x,y,k)=s_up(x,y,k)^(-2/(mc-1))/tp1;
                    u2_d(x,y,k)=s_down(x,y,k)^(-2/(mc-1))/tp1;
                end
            end
        end
        
        
        for x=1:m
            for y=1:n
                for k=1:c
                    u_up(x,y,k)=max(u1(x,y,k),u2_u(x,y,k));
                    u_down(x,y,k)=min(u1(x,y,k),u2_d(x,y,k));
                end
            end
        end
        
        for x=1:m
            for y=1:n
                sum0=0;
                sum1=0;
                for k=1:c
                    sum0=sum0+u_up(x,y,k);
                    sum1=sum1+u_down(x,y,k);
                end
                for k=1:c
                    u_up(x,y,k)=u_up(x,y,k)/(sum0+0.00001);
                    u_down(x,y,k)=u_down(x,y,k)/(sum1+0.00001);
                end
            end
        end
        
        Uup=reshape(u_up,m,n,c);
        Udown=reshape(u_down,m,n,c);
        %     [U,V] = compute_v(m,n,c,V,X,reshape(Uup,m,n,c),reshape(Udown,m,n,c));
        %     V = V';
        %%
        xx=zeros(c,2);
        yy=zeros(c,2);
        u_new_R=zeros(m,n,c);
        u_new_L=zeros(m,n,c);
        V_R=zeros(c,1);
        V_L=zeros(c,1);
        v11=V;
        v22=V;
        data1=sort(X(:));
        data2=reshape(data1,m,n);
        
        
        %UR
        comparision=0;
        while (comparision==0)
            for k=1:c
                for x=1:m*n-1
                    if v11(k)>data1(x) && v11(k)<data1(x+1)
                        xx(k,1)=x;
                    end
                end
            end
            
            for k=1:c
                for x=1:m
                    for y=1:n
                        if X(x,y)<=data1(xx(k,1))
                            u_new_R(x,y,k)=Udown(x,y,k);
                        else
                            u_new_R(x,y,k)=Uup(x,y,k);
                        end
                    end
                end
            end
            
            for k=1:c
                tpe1=0;
                tpe2=0;
                for x=1:m
                    for y=1:n
                        tpe1=tpe1+u_new_R(x,y,k)*X(x,y);
                        tpe2=tpe2+u_new_R(x,y,k);
                    end
                end
                V_R(k)=tpe1/tpe2;
            end
            if sum(sum(roundn(V_R,-1) == roundn(v11,-1)))/(c) == 1
                comparision=1;
            else
                v11=V_R;
            end
        end
        
        %UL
        comparision=0;
        while (comparision==0)
            for k=1:c
                for x=1:m*n-1
                    if v22(k)>data1(x) && v22(k)<data1(x+1)
                        yy(k,1)=x;
                    end
                end
            end
            
            for k=1:c
                for x=1:m
                    for y=1:n
                        if X(x,y)<=data1(xx(k,1))
                            u_new_L(x,y,k)=Uup(x,y,k);
                        else
                            u_new_L(x,y,k)=Udown(x,y,k);
                        end
                    end
                end
            end
            
            for k=1:c
                tpe1=0;
                tpe2=0;
                for x=1:m
                    for y=1:n
                        tpe1=tpe1+u_new_L(x,y,k)*X(x,y);
                        tpe2=tpe2+u_new_L(x,y,k);
                    end
                end
                V_L(k)=tpe1/tpe2;
            end
            if sum(sum(roundn(V_L,-1) == roundn(v22,-1)))/(c) == 1
                comparision=1;
            else
                v22=V_L;
            end
        end
        
        V = (V_L+V_R)/2;
        U = (u_new_L+u_new_R)./2;
        
        
        %     t3 = (dist./sum(dist,2));
        %
        %     Uup(find(t3<1/3)) = U1(find(t3<1/3)); %#ok<FNDSB>
        %     Uup(find(t3>=1/3)) = U2(find(t3>=1/3)); %#ok<FNDSB>
        %
        %
        %     Udown(find(t3<1/3)) = U2(find(t3<1/3)); %#ok<FNDSB>
        %     Udown(find(t3>=1/3)) = U1(find(t3>=1/3)); %#ok<FNDSB>
        
        %%
        disti = sum((V1-V).^2);
        if disti < 1e-6
            break;
        end
        %         J_cur = sum(sum((reshape(U.^mc.*dt,m*n,c))))/N;
        %         J = [J J_cur];
        %         fprintf('#iteration: %03d, objective function: %f\n, disti: %f\n', iter, J_cur, disti);
        %         J_prev = J_cur;
        
        fprintf('#iteration: %03d, disti: %f\n', iter, disti);
        
    end
    
    %% 画图
    [~, label] = max(reshape(U,m*n,c), [], 2);
    
    A = label;
    B = label_real;
    NMI  = compute_nmi (A, B);
    [ARI, RI, MI, ~]=RandIndex(label',label_real);
    [HI,~,~] = compute_f(label,label_real);
    DBI = evalclusters(reshape(X,m*n,1),label,'DaviesBouldin');
    DBI =DBI.CriterionValues;
    SC = Eva_Entropy(label,reshape(X,m*n,1),1);
    CHI = evalclusters(reshape(X,m*n,1),label,'CalinskiHarabasz');
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
xlswrite('demo.xlsx',[nmi,ari, hi, dbi, sc,chi]);


%%
% U = (Uup+Udown)/2;
% [~, label] = max(reshape(U,m*n,c), [], 2);

% [V,U] = fcm(X,3);
% figure
% imshow(uint8(reshape(V(label, :), m, n, p)))

% [fs,center_p,Num_p,center_lab] = Label_image(I,reshape(label,m,n));
% figure
% imshow(fs);
%imwrite(fs,'D:\myprogram\FCM\FCM_my\paper2\testresult\AIT2FCM_0002.PNG')


%% gray
[~, label] = max(reshape(U,m*n,c), [], 2);
figure
label1=zeros(N,1);
label1(label==1) = 1;
label1(label==2) = 2;
label1(label==3) = 3;
imshow(reshape(uint8((label1-1)*255./(k-1)),m,n));
% imwrite(reshape(uint8((label1-1)*255./(k-1)),m,n),'D:\myprogram\FCM\FCM_my\paper2\testresult\AIT2FCM_N44.PNG')
