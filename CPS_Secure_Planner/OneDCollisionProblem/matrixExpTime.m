function result = matrixExpTime(t,y)
global A B
result = (expm(A*t(1))*B);  
for i=2:length(t) 
    result = vertcat(result,(expm(A*t(i))*B)) 
end

end