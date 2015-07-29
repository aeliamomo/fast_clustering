

%x = [xl(:)'; yl(:)'; zl(:)'];
patches = [1 1 1 1 1;2 2 2 2 2; 4 5 6 7 8; 4 5 6 7 8; 4 5 6 7 8] 
IP = patches' * patches;%inner product
%d = sum((bsxfun(@minus, patches(i,:), patches(j,:)).^2));
d = sqrt(bsxfun(@plus, diag(IP), diag(IP)') - 2 * IP)

noOfPatches = size(patches,1)
n1 = size(patches,1);
n2 = size(patches,1);
D = zeros(n1,n1);

t = cputime;
for i = 1:n1
    for j = 1:n1
        D(i,j) = sum((bsxfun(@minus, patches(i,:), patches(j,:)).^2));
        %D(i,j) = sum((patches(i,:)-patches(j,:)).^2);
    end
end
e = cputime-t;

t2 = cputime;
IP = dot(patches',patches');
D_2 = bsxfun(@plus,IP',IP)-2*(patches*patches');

e2 = cputime-t2;
