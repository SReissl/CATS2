function [permInc] = projectIncome(probs11,probs12,probs21,probs22,pay1,pay2,rate)

mats=zeros(2,2,length(pay1));
mats(1,1,:)=probs11;
mats(1,2,:)=probs12;
mats(2,1,:)=probs21;
mats(2,2,:)=probs22;
pay=zeros(2,1,length(pay1));
pay(1,1,:)=pay1;
pay(2,1,:)=pay2;
perm=pay(:,:,1);
prod=eye(2);
i=2;
while sum(prod*pay(:,:,i)/((1+rate)^(i-1)))>0.0001
    prod=prod*mats(:,:,i-1);
    perm=perm+prod*pay(:,:,i)/((1+rate)^(i-1));
    i=i+1;
end
permInc=perm;

end
