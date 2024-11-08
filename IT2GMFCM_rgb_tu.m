clear all;
close all;
clc;
% X = imread('D:\myprogram\FCM\FCM_my\database\0042.tif');
X = imread('D:\myprogram\FCM\FCM_my\paper3\testresult\airplane35.tif');

[m, n ,p] = size(X);
X = reshape(double(X),m*n,p);


c = 3; max_num = 200;iter = 0;
N = m*n;
b1 = 2; b2 = 5;b3 = (b1+b2)/2;

%%


for l = 1:1
    
%     U = abs(randn(m, c));%P为 N行k列
%     U = U./(sum(U, 2)*ones(1, c)); %sum(P,2)求P每一行的和

%      V = abs(rand(c,p));%P为 N行k列
%      V = V.*255;
     
%      V = zeros(c, p);%P为 N行k列
%     for i =1:c
%         fd = find(label_real == i);
%         if isempty(fd)
%            break; 
%         end
%         V(i,:) = X( randi(numel(fd),1) ,:);
%     end
%     
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
        

        dist = abs(sum(X.*X, 2)*ones(1, c) + (sum(V.*V, 2)*ones(1, N))'-2*X*V'+eps);
        
        
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
        xx=zeros(c,p);
        yy=zeros(c,p);
        u_new_R=zeros(N,p);
        u_new_L=zeros(N,p);
        V_R=zeros(c,p);
        V_L=zeros(c,p);
        v11=V;
        v22=V;
        data1=sort(X,1);
%         data2=reshape(data1,m,n);
        
        
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
                for y=1:p
                    for x=1:N
                        if x<=xx(k,y)
                            u_new_R(x,y,k)=Udown(x,k);
                        else
                            u_new_R(x,y,k)=Uup(x,k);
                        end
                    end
                end
            end
            
            u_new_R1 = sum(u_new_R,2)/p;         
            u_new_R1 = reshape(u_new_R1,N,c);
            
         
            V_R = (X'*u_new_R1.^b3)'./(sum(u_new_R1.^b3)'*ones(1, p));%V(0)
            for ii = 1:p
                if sum(sum(roundn(V_R(ii,:),-1) == roundn(v11(ii,:),-1)))/(c) == 1
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
        if disti < 1e-4
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
% imwrite(uint8(reshape(V(label, :), m, n, p)),'D:\myprogram\FCM\FCM_my\paper2\testresult\IT2GMFCM_0005.PNG')

%% gray
% figure
% label1=zeros(N,1);
% label1(label==1) = 1;
% label1(label==2) = 2;
% label1(label==3) = 3;
% imshow(reshape(uint8((label1-1)*255./(k-1)),m,n));
% imwrite(reshape(uint8((label1-1)*255./(k-1)),m,n),'D:\myprogram\FCM\FCM_my\paper2\testresult\IT2GMFCM_0005.PNG')
%%
% [fs,center_p,Num_p,center_lab] = Label_image( imread('7.jpg'),reshape(label,m,n));
% imshow(fs)
% imwrite(fs,'D:\myprogram\FCM\FCM_my\paper2\testresult\7_5.PNG')