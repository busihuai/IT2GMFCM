clc
clear all
close all

%%
I = imread('D:\myprogram\FCM\FCM_my\database\storagetanks35.tif');

[m, n ,p] = size(I);
X = reshape(double(I),m*n,p);



%%
c = 3;[n,D]=size(X);
iter = 0;

    
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

end

%% 画图
[~, label1] = max(u',[], 1);
[m, n ,p] = size(I);
figure
imshow(uint8(reshape(v(label1, :), m, n, p)))
[fs,center_p,Num_p,center_lab] = Label_image( I,reshape(label1,m,n));
imshow(fs)

