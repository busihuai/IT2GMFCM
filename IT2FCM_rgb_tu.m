clear all;
close all;
clc;
%X = imread('6.jpg');
 X = imread('D:\myprogram\FCM\FCM_my\database\baseballdiamond03.tif');
[m, n ,p] = size(X);
X = reshape(double(X),m*n,p);

c = 3; max_num = 500;iter = 0;
N = m*n;
b1 = 2; b2 = 4; b3 = (b1+b2)/2;
%%
for l = 1:1
    
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
        
        dist = sum(X.*X, 2)*ones(1, c) + (sum(V.*V, 2)*ones(1, N))'-2*X*V'+eps;
        
        t1 = (1./dist).^(1/(b1-1));
        U1 = t1./(sum(t1, 2)*ones(1, c));
        
        t2 = (1./dist).^(1/(b2-1));
        U2 = t2./(sum(t2, 2)*ones(1, c));
        
        Uup = max(U1,U2);
        Udown = min(U1,U2);
        
        t3 = ((Uup+Udown)./2).^b3; %uij(t+1)的m次方
        V = (X'*t3)'./(sum(t3)'*ones(1, p));%V(0)
        
        %     Uup=reshape(Uup,m,n,c);
        %     Udown=reshape(Udown,m,n,c);
        %     [U,V] = compute_v(m,n,c,V,X,reshape(Uup,m,n,c),reshape(Udown,m,n,c));
        %     V = V';
        %% 降维
        xx=zeros(c,p);
        yy=zeros(c,p);
        u_new_R=zeros(N,p);
        u_new_L=zeros(N,p);
        V_R=zeros(c,p);
        V_L=zeros(c,p);
        v11=V;
        v22=V;
        data1=sort(X,1);
        
        %UR
        comparision=0;
        while (comparision==0)
            for k=1:c
                for y=1:p
                    for x=1:N-1
                        if v11(k,y)>=data1(x,y) && v11(k,y)<=data1(x+1,y)
                            xx(k,y)=x;
                        end
                    end
                end
            end
            
            for k=1:c
                for x=1:N
                    for y=1:p
                        if X(x,y)<=data1(xx(k,y),y)
                            u_new_R(x,y,k)=Udown(x,k);
                        else
                            u_new_R(x,y,k)=Uup(x,k);
                        end
                    end
                end
            end
            
            u_new_R1 = sum(u_new_R,2)/p;
            u_new_R1 = reshape(u_new_R1,N,c);
            
            t = u_new_R1.^b3; %uij(0)的m次方
            V_R = (X'*t)'./(sum(t)'*ones(1, p));%V(0)
            for ii = 1:p
                if sum(sum(roundn(V_R(ii,:),-1) == roundn(v11(ii,:),-1)))/(c) == 1
                    comparision=1;
                else
                    v11=V_R;
                end
            end
        end
        
        %UL
        %UL
        comparision=0;
        while (comparision==0)
            for k=1:c
                for y=1:p
                    for x=1:N-1
                        if v22(k,y)>=data1(x,y) && v22(k,y)<=data1(x+1,y)
                            yy(k,y)=x;
                        end
                    end
                end
            end
            
            for k=1:c
                for y=1:p
                    for x=1:N
                        if x<=xx(k,y)
                            u_new_L(x,y,k)=Uup(x,k);
                        else
                            u_new_L(x,y,k)=Udown(x,k);
                        end
                    end
                end
            end
            
            u_new_L1 = sum(u_new_L,2)/p;  
            u_new_L1 = reshape(u_new_L1,N,c);
            
           
            V_L = (X'*u_new_L1.^b3)'./(sum(u_new_L1.^b3)'*ones(1, p));%V(0)
            for ii = 1:p
                if sum(sum(roundn(V_L(ii,:),-1) == roundn(v22(ii,:),-1)))/(c) == 1
                    comparision=1;
                else
                    v22=V_L;
                end
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
end

%%
[~, label] = max(U, [], 2);
figure
imshow(uint8(reshape(V(label, :), m, n, p)))
imwrite(uint8(reshape(V(label, :), m, n, p)),'D:\myprogram\FCM\FCM_my\paper2\testresult\IT2FCM_baseballdiamond03.PNG')


