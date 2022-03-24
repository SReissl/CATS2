function [permInc] = projectIncomeMult(probs11,probs12,probs21,probs22,pay1,pay2,rate)

mats=zeros(2,2,length(pay1(:,1)));
mats(1,1,:)=probs11;
mats(1,2,:)=probs12;
mats(2,1,:)=probs21;
mats(2,2,:)=probs22;
prod=zeros(2,2,length(pay1(:,1)));
pay=zeros(2,length(pay1(:,1)),length(pay1(1,:)));
pay(1,:,:)=pay1;
pay(2,:,:)=pay2;
perm=zeros(2,size(pay1,2));
prod(:,:,1)=eye(2);
perms=zeros(2,size(pay1,2),size(pay1,1));
perms(1,:,1)=pay1(1,:);
perms(2,:,1)=0;
for i=2:length(pay1(:,1))
    prod(:,:,i)=prod(:,:,i-1)*mats(:,:,i-1);
    pay(1,i,:)=pay(1,i,:)./((1+rate)^(i-1));
    pays=zeros(2,size(pay1,2));
    pays(:,:)=pay(:,i,:);
    perms(:,:,i)=prod(:,:,i)*pays;
end
p1=zeros(size(pay1,2),size(pay1,1));
p1(:,:)=perms(1,:,:);
perm(1,:)=sum(p1,2);
p2=zeros(size(pay1,2),size(pay1,1));
p2(:,:)=perms(2,:,:);
perm(2,:)=sum(p2,2);
permInc=perm;

end
