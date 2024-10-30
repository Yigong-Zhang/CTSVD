function lambda = set_lambda_my(Nway,alpha,set) 
N = length(Nway);
lambda = 0;
switch set
    case 'square'
        for i=1:N
            if i ~= N
                [temp1,temp2] = crack_integer(Nway(i),Nway(i+1));
                d = prod(Nway)/(Nway(i)*Nway(i+1));
            else
                [temp1,temp2] = crack_integer(Nway(N),Nway(1));
                d = prod(Nway)/(Nway(N)*Nway(1));
            end
            temp=alpha(i)/sqrt(max(temp1,temp2)*d);
            lambda = lambda + temp;
        end
    case 'normal'
        for i=1:N
            if i ~= N
                temp1 =Nway(i);
                temp2 = Nway(i+1);
                d = prod(Nway)/(Nway(i)*Nway(i+1));
            else
                temp1 =Nway(N);
                temp2 = Nway(1);
                d = prod(Nway)/(Nway(N)*Nway(1));
            end
            temp=alpha(i)/sqrt(max(temp1,temp2)*d);
            lambda = lambda + temp;
        end
end
