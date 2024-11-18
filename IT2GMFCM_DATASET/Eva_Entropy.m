function entropy_R = Eva_Entropy(Cluster_lable,x,n)
% n is the number of attributes (categorical).

[m,~] = size(x);
t = length( unique( Cluster_lable ) );
R = cell(1,t);
value_lable = unique( Cluster_lable );
number_in_cluster = zeros(1,t);

for i = 1:m
    for j = 1:t
        if Cluster_lable(i) == value_lable(j)
            number_in_cluster(1,j) = number_in_cluster(1,j) + 1;
        end
    end
end

for i = 1:t
    R{i} = zeros(number_in_cluster(1,i),n);
end

%record cluster R
for i = 1:t
    count_x_in_rt = 1;
    for j = 1:m
        if Cluster_lable(j) == value_lable(i)
            for v = 1:n
                R{i}(count_x_in_rt,v) = x(j,v);     
            end
            count_x_in_rt = count_x_in_rt + 1;
        end
    end
end

%number of attr
num_att = zeros(1,n);
for i = 1:n
    num_att(1,i) = length(unique(x(:,i)));
end

%value of attr
value_att = cell(1,n);
for i = 1:n
    value_att{i} = unique(x(:,i));
end

%probabilities,Eq.34
p_ult = cell(1,n);
for i = 1:n
    p_ult{i} = zeros(t,num_att(1,i));
end


for i = 1:n
    for j = 1:num_att(1,i)
        for k = 1:t
            nlt = 0;
            for mmt = 1:number_in_cluster(1,k)
                if value_att{i}(j) == R{k}(mmt,i)
                    nlt = nlt + 1;
                end
            end
            nlt = nlt / number_in_cluster(1,k);
            p_ult{i}(k,j) = nlt / number_in_cluster(1,k);
        end
    end
end


%Hlt,Eq.35
Hlt = zeros(t,n);
for i = 1:t
    for j = 1:n
        for v = 1:num_att(1,j)
            if p_ult{j}(i,v) ~= 0
                Hlt(i,j) = Hlt(i,j) - p_ult{j}(i,v)*log(p_ult{j}(i,v));
            end
        end
    end
end

%Eq.36
Ht = zeros(1,t);
for i = 1:t
    for j = 1:n
        Ht(1,i) = Ht(1,i) + Hlt(i,j);
    end
    Ht(1,i) = Ht(1,i)/n;
end

%Eq.37
entropy_R = 0;
for i = 1:t
    entropy_R = entropy_R + Ht(1,i);
end
entropy_R = entropy_R / t;

