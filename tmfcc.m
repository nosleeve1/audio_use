function [Tmfcc,Road] = tmfcc(a,b,limit,choose)
% DTW法比较两段语音特征距离   Tmfcc为语音mfcc距离值，Road为最短匹配路径

% a，b分别为两段语音的特征矩阵，limit为限路因子，choose为匹配方式选择因子
%% 

a = a';
b = b';
[cc,afn] = size(a); 
bfn = size(b,2);
%% 距离计算
d = ones(afn,bfn)*inf;   %距离
k = bfn/afn;

for m = 1:afn
    for n = max(1,floor(k*m-bfn/limit)):min(ceil(k*m+bfn/limit),bfn)
        dn = (a(:,m)-b(:,n)).^2;    %计算距离
        d(m,n) = sqrt(sum(dn))/cc;   %欧氏距离
    end
end

D = zeros(afn+1,bfn+1);
D(1,:) = inf;        
D(:,1) = inf;         
D(1,1) = 0;
D(2:end, 2:end) = d;
%% 寻找整个过程的最短匹配距离
if choose == 2 %非固定匹配数  
    for m = 1:afn
        for n = max(1,floor(k*m-bfn/limit)):min(ceil(k*m+bfn/limit),bfn)
            dmin = min([D(m,n), D(m,n+1), D(m+1,n)]);
            D(m+1,n+1) = d(m,n)+dmin;
        end
    end    
else %固定匹配数
    for m = 1:afn
        for n = max(1,floor(k*m-bfn/limit)):min(ceil(k*m+bfn/limit),bfn)
            dmin = min([D(m,n)+d(m,n),D(m,n+1),D(m+1,n)]);%斜方向步长取2
            D(m+1,n+1) = d(m,n)+dmin;
        end
    end    
    Count = afn+bfn-2;%匹配数
end
%% 计算匹配路径
i = afn+1;
j = bfn+1;
count = 0;
road = zeros(afn+bfn-2,2);
while i >= 2 && j >= 2
    road(count+1,1) = i-1;
    road(count+1,2) = j-1;
    if(i < 2 && j >= 2)
        loc = 3;
    elseif(i >= 2 && j < 2)
        loc = 2;
    else
        [~,loc] = min([D(i-1,j-1),D(i-1,j),D(i,j-1)]);
    end
    
    if loc == 3
        j = j-1;
    elseif loc == 2
        i = i-1;
    else
        i = i-1;
        j = j-1;
    end
    count = count+1;    
end
%% 输出
Road = road(1:count,:);
% xlabel("a");
% ylabel("b");
if choose == 2    
    Tmfcc = D(afn+1,bfn+1)/count;
else
    Tmfcc = D(afn+1,bfn+1)/Count;
end
end
