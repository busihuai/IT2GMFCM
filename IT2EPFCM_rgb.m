clc;clear all;

img1 = imread('10.jpg');
img = double(img1);
load('D:\myprogram\FCM\FCM_my\database\10.mat')
label_real = label;
label = [];


c = 3;

m1=2;
mc1=2;
mc2=4;
nc1=2;
nc2=4;

% [v1,U] = fcm(X1,c);


% v1(1)=40;
% v1(2)=150;
% v1(3)=200;
% v1(4)=250;
[m,n,r]=size(img);
% Data=reshape(img,m*n,1);
N=m*n;
Data = reshape(img,N,c);
a=1;%Cf
b=1;%Cp
tt=clock;
[data_points,~]=size(Data);
Tmax=200;
center=zeros(c,r);
MeanK=mean(Data);
for i=1:c
    for j = 1:r
        center(i,j)=( (3*i)/(2*c) )*MeanK(1);
    end
end

mc= 0.5.*(mc1+mc2);
% [u,center_f]=fcm12(m,n,c,r,center,Data,mc);%u-----(c,n*m)  使用fcm初始化
% [t, center_t ]=pcm11(m,n,c,r,v,Data,nc,u);     %使用pcm初始化
% u = u';


     rand('seed',sum(100*clock))
     for ii =1:c
         fd = floor(abs(rand(c, 1))*m*n);
         center_f(ii,:) = Data(fd(ii),:) ;
     end


u = abs(rand(N, c));%P为 N行k列
u = u./(sum(u, 2)*ones(1, c)); 
u = u';

% v=0.5.*(center_f +center_t);
v1 = center_f;
v2 = center_f;
v1_new = v1;
v2_new = v2;

% %循环外的初始化neta
K=1;
neta_old= ones(1,c);
distances_K = abs(sum(Data.*Data, 2)*ones(1, c) + (sum(center.*center, 2)*ones(1, N))'-2*Data*center'+eps);
mf=u.^(mc);
% nf=t.^(nc);
for k=1:c
    %     neta_old(1,k) = ( (mf(k,:).*nf(k,:))*distances_K(:,k) );
    neta_old(1,k) = ( (mf(k,:))*distances_K(:,k) );
end
% neta_old = neta_old./(sum(mf.*nf,2))';
neta_old = neta_old./(sum(mf,2))';
neta = (ones(data_points,1) * neta_old)*K.*2 ; %neta for all clusters

kk5=1;
term_thr = 1e-4;     % Termination threshold
E1 = zeros(Tmax, 1);	% Array for termination measure values
E2 = zeros(Tmax, 1);	% Array for termination measure values

while (kk5<Tmax)    %终止条件
    
    %求解欧式距离
    dist1 = distfcm(v1, Data)+0.001;  %dlist----(c,n*m)
    dist2 = distfcm(v2, Data)+0.001;
    
    %计算可能隶属度---指数型的PFCM
    tep1 = -1.*( b.*(dist1.^2) )./( neta' );
    t1 = nthroot(exp(tep1),nc1);
    tep2 = -1.*( b.*(dist2.^2) )./( neta' );
    t2 = nthroot(exp(tep2),nc2);
    T_up=max(t1,t2);
    T_low=min(t1,t2);
    
    % 计算模糊隶属度
    tmp1 =   ( dist1.^2 ).^(-1/(mc1-1));      % calculate new U, suppose mc != 1
    u1 = tmp1./(ones(c, 1)*sum(tmp1)); %U_new-----(c,n*m)
    tmp2 =   ( dist2.^2 ).^(-1/(mc2-1));      % calculate new U, suppose mc != 1
    u2 = tmp2./(ones(c, 1)*sum(tmp2)); %U_new-----(c,n*m)
    U_up=max(u1,u2);
    U_low=min(u1,u2);
    
    %%混合隶属度
    UT1= a.*U_up + b.*T_up;
    UT2= a.*U_up + b.*T_low;
    UT3= a.*U_low + b.*T_up;
    UT4= a.*U_low + b.*T_low;
    
    UT_up = max(max(UT1,UT2),max(UT3,UT4));
    UT_low= min(min(UT1,UT2),min(UT3,UT4));
    
    
    %降型
        % Karni-Mendel Alg
    for km=1:c
        UT_mean=(UT_up+UT_low)/2;
        maxUT_mean = max(UT_mean);

        index_c = find(UT_mean(km, :) == maxUT_mean);
        Y= Data(index_c,:);

        a_old=UT_up(km,index_c);
        b_old=UT_low(km,index_c);
        a_old=a_old(:);
        b_old=b_old(:);
        F=cat(2,a_old,b_old);
        if size(Y,1)<=1
            a1=Y;
            a2=km;
            a3=index_c;
            a4=UT_mean;
            a5=maxUT_mean;
        end

        [XLeft,XRight,L,R]=KM_Alg(F,Y);

        v1_new(km,:)=XLeft';   
        v2_new(km,:)=XRight';
    end
    

    % check termination condition
    E1(kk5) = norm (v1_new - v1, 1);
    E2(kk5) = norm (v2_new- v2, 1);
     
        if E1(kk5) <= term_thr
            break; 
        else
            v1=v1_new;
%             v2=v2_new;
        end
        if E2(kk5) <= term_thr
            break;
        else
%             v1=v1_new;
            v2=v2_new;  
        end
          kk5=kk5+1;
%     % 终止条件 
%     if   d<0.00001
%         break;
%     else
%         fprintf('IT2EPFCM_x:Iteration count = %d, dv=%f\n', kk,dv);
%         kk=kk+1;
%         v=center_new;
%     end   
end
title('IT2EPFCM');
disp(['IT2EPFCM运行时间：',num2str(etime(clock,tt))]);
timeuse = etime(clock,tt);

%%
u= a.*UT_up + b.*UT_low;
[~,label] = max(u); %找到所属的类
[NMI,ARI, F,DBI, SC,CHI] = index_compute(label',label_real,Data);
xlswrite('demo.xlsx',[NMI,ARI, F,DBI, SC,CHI,timeuse]);

%变化到图像的大小
figure
[fs] = Label_image( img1,reshape(label,m,n));
imshow(fs)
imwrite(fs,'D:\myprogram\FCM\FCM_my\paper4\IT2EPFCM_10.PNG')



