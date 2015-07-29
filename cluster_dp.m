function clusterPatch = cluster_dp(xx, percent)
%% This package implement the clustering algorithm described by Alex Rodriguez and Alessandro Laio (2014). 
%  Rodriguez, A., & Laio, A. (2014). Clustering by fast search and find of density peaks. Science, 344(6191), 1492-1496. doi:10.1126/science.1242072

%clear all
close all
disp('The only input needed is a distance matrix file')
disp('The format of this file should be: ')
disp('Column 1: id of element i')
disp('Column 2: id of element j')
disp('Column 3: dist(i,j)')

% Passed the distance matrix file as an input parameter.
% mdist=input('name of the distance matrix file\n','s');
%disp('Reading input distance matrix')
%xx=load(mdist);

ND=max(xx(:,2))
NL=max(xx(:,1));
if (NL>ND)
  ND=NL;
end
%no of elems
N=size(xx,1);
%diagonal zeros for similar elems
for i=1:ND
  for j=1:ND
    dist(i,j)=0;
  end
end
for i=1:N
  ii=xx(i,1);
  jj=xx(i,2);
  dist(ii,jj)=xx(i,3);
  dist(jj,ii)=xx(i,3);
end

% Vary the percent values to alter the cut-off distance (dc).
% Use the input parameter to calculate the cut off distance using the percentage.
%percent=2.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
position=round(N*percent/100);
sda=sort(xx(:,3),'descend');

%sda=sort(xx(:,3));
%max_dist = sda(1); 
%fprintf('maximum distance found in matrix: %12.6f\n', max_dist);
%dc=(max_dist*13/100) 
dc=(sda(position));

fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);

for i=1:ND
  rho(i)=0.;
end
%
% Gaussian kernel
%
for i=1:ND-1
  for j=i+1:ND
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end
%
% "Cut off" kernel
%
%for i=1:ND-1
%  for j=i+1:ND
%    if (dist(i,j)<dc)
%       rho(i)=rho(i)+1.;
%       rho(j)=rho(j)+1.;
%    end
%  end
%end

maxd=max(max(dist));

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

fid = fopen('DECISION_GRAPH', 'w');
for i=1:ND
   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end

disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);

for i=1:ND
  ind(i)=i;
  gamma(i)=rho(i)*delta(i);
end
subplot(2,1,1);
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')


subplot(2,1,1);
rect = getrect(1)
%rect = [135.6989    0.1189   60.2151    0.1653];%01_07
%rect = [214.2473    0.0378  110.2151    0.0383];%01_02
rhomin=rect(1);
deltamin=rect(2);
NCLUST=0;
for i=1:ND
  cl(i)=-1;
end
for i=1:ND
  if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;
     icl(NCLUST)=i;
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')

%assignation
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end
%halo
for i=1:ND
  halo(i)=cl(i);
end
if (NCLUST>1)
  for i=1:NCLUST
    bord_rho(i)=0.;
  end
  for i=1:ND-1
    for j=i+1:ND
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;
        end
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:ND
    if (rho(i)<bord_rho(cl(i)))
      halo(i)=0;
    end
  end
end
for i=1:NCLUST
  nc=0;
  nh=0;
  for j=1:ND
    if (cl(j)==i) 
      nc=nc+1;
    end
    if (halo(j)==i) 
      nh=nh+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
  
end

NCLUST
cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
   
end
dlmwrite('centroids',icl(:));
subplot(2,1,2)

disp('Performing 2D nonclassical multidimensional scaling')

%Y1 = cmdscale(dist);
Y1 = mdscale(dist, 2, 'criterion','metricsstress');

disp('Performing 2D Classical multidimensional scaling')
Y1 = cmdscale(dist);
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Classical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end

% Return the cluster numbers for the patches as an output parameter.
clusterPatch = [];
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (halo(j)==i)
      nn=nn+1;
      clusterPatch = [clusterPatch ;[j,i]];
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  %dlmwrite(sprintf('cluster_%d.txt',i),xx(new_j,:));
  plot(A(1:nn,1),A(1:nn,2),'.','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
savefig(sprintf('/Users/chingyukao/Documents/MATLAB/Multi-Modal-Similarity/tuningParameters/centroid_info_%f.fig',percent));
%save('clusteredPatches.mat','clusterPatch');

%for i=1:ND
%   if (halo(i)>0)
%      ic=int8((halo(i)*64.)/(NCLUST*1.));
%      hold on
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%   end
%end
faa = fopen('CLUSTER_ASSIGNATION', 'w');
%faaa = fopen('faaa', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
    %%icl(cl(i)) represents the centroids
   fprintf(faa, '%i %i %i %i\n',i,cl(i),halo(i),icl(cl(i)));
   %fprintf(faaa, '%i %i %i\n',i,icl(cl(i)),halo(i));
   %icl(i)
end

end