function beta = changepoint(y, x)

if size(y,2) == length(x) && size(y,1)>1
    x = repmat(x, [size(y,1) 1]);
end

y = y(:);
x = x(:);
n = length(find(~isnan(y)));

err2min = inf;

grid = [min(x) min(x(x>min(x))) : (max(x)-min(x(x>min(x))))/100 : max(x)];
beta = [0 0 0];

for i = 1:length(grid)
    c = grid(i);
    if i==1
        a = 0;
        b = nanmean(y);
    elseif i==length(grid)
        r = regress(y, [x ones(length(x),1)]);
        a = r(1);
        b = r(2);
    else
        r = regress(y, [min(x,c) ones(length(x),1)]);
        a = r(1);
        b = r(2);
    end
    
    % constrain slope to be positive
    if a < 0
        a = 0;
        b = nanmean(y);
    end

    yhat = a*x + b;
    yhat(x>c) = a*c + b;
    MSS(i) = nansum((yhat-y).^2) / n;
    
    if MSS(i) < err2min 
        beta = [a b c];
        err2min = MSS(i);
    end
end
