function y = sigmoid(xval,param)
    y = param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4))); 
end
