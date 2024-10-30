%% generate tensor with CP rank r
function X = gencptensor(r,dim)
N = length(dim);
rng(1)
A1 = normrnd(0,1,int32([dim(1),r]));
rng(2)
A2 = normrnd(0,1,int32([dim(2),r]));
rng(3)
A3 = normrnd(0,1,int32([dim(3),r]));
rng(4)
A4 = normrnd(0,1,int32([dim(4),r]));
% rng(5)
% A5 = normrnd(0,1/dim(5),int32([dim(5),r]));

X = zeros(dim);
for i1 = 1:dim(1)
    for i2 = 1:dim(2)
        for i3 = 1:dim(3)
            for i4 = 1:dim(4)
                X(i1,i2,i3,i4) = sum(A1(i1,:).*A2(i2,:).*A3(i3,:).*A4(i4,:));
            end
        end 
    end
end
%X = (X-min(X(:)))./(max(X(:))-min(X(:)));
% X = X./max(abs(X(:)));


% for i1 = 1:dim(1)
%     for i2 = 1:dim(2)
%         for i3 = 1:dim(3)
%             for i4 = 1:dim(4)
%                 for i5 = 1:dim(5)
%                     X(i1,i2,i3,i4,i5) = sum(A1(i1,:).*A2(i2,:).*A3(i3,:).*A4(i4,:).*A5(i5,:));
%                 end
%             end
%         end 
%     end
% end