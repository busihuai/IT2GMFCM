clear all;clc;
X  = -1:.05:1;
v1 = -1;
v2 = 1;
V = [v1 v2];
%% 初始化参数
[m, n] = size(X);
k = 2; b = 2;max_num = 300;

name = ["m=2","m=3","m=4","m=5","m=6","m=7","m=8","m=9"];
    



for i = 2:9
    
    b = i;
    dist1 = (X-v1).^2+eps;
    dist2 = (X-v2).^2+eps;
    dist=[dist1;dist2]';
    t = (1./dist).^(1/(b-1));
    
    pij = t./sum(t, 2);
    
   
    
    % a = find(dist == 0);
    % [row,col] = ind2sub(size(X),a);
    % p0 = 2.5;
    % p00 = p0^(-1/(k-1));
    %
    % for j = 1:size(row,1)
    %     pij(row(j),:) = p00;
    %     pij(row(j),col(j)) = p0;
    % end

%     plot(X,log(pij)/log(10),'LineWidth',1.5)
    p(i-1) = plot(X,pij(:,1),'LineWidth',1.5,'DisplayName',name(i-1));
    hold on
end



lgd = legend;
lgd.FontSize = 20;
lgd.Location = 'eastoutside';
grid on;
xlabel('x','fontname','Times New Roman','fontsize',40,'fontweight','bold');
ylabel('u(x)','fontname','Times New Roman','fontsize',40,'fontweight','bold');
box on
set(gca,'FontSize',20,'FontName','Times New Roman');

set(p(1),'LineWidth',3,'color',[0 0.4470 0.7410])
set(p(1),'Marker','o','MarkerSize',3)

set(p(2),'LineWidth',3,'color',[0.8500 0.3250 0.0980])
set(p(2),'Marker','>','MarkerSize',3)

set(p(3),'LineWidth',3,'color',[0.9290 0.6940 0.1250])
set(p(3),'Marker','diamond','MarkerSize',3)

set(p(4),'LineWidth',3,'color',[0.4940 0.1840 0.5560])
set(p(4),'Marker','*','MarkerSize',3)

set(p(5),'LineWidth',3,'color','yellow')
set(p(5),'Marker','^','MarkerSize',3)

set(p(6),'LineWidth',3,'color',[0.3010 0.7450 0.9330])
set(p(6),'Marker','square','MarkerSize',3)

set(p(7),'LineWidth',3,'color',[0.4660 0.6740 0.1880])
set(p(7),'Marker','.','MarkerSize',3)

set(p(8),'LineWidth',3,'color',[0.6350 0.0780 0.1840])
set(p(8),'Marker','pentagram','MarkerSize',3)

